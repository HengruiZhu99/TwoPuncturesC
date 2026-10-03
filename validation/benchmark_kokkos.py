"""Sequential, restartable performance matrix with immutable worker images.

Run inside a Slurm step exposing one GPU and its allocated CPU fraction.
Each variant is a fresh process; no numerical workers overlap. Failed gates
are retained and do not prevent timing the remaining predeclared controls.
"""
import argparse,json,os,platform,re,subprocess,sys,threading,time
from pathlib import Path
import numpy as np
from benchmark_bowen_york import ROOT,digest
from benchmark_krylov_matrix import by_compare,hi_compare,checks


def hardware():
    result=dict(platform=platform.platform(),hostname=platform.node(),python=sys.version,
                affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
                environment={key:os.environ.get(key) for key in ('SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_CPUS_PER_TASK','SLURM_JOB_GPUS','CUDA_VISIBLE_DEVICES','OMP_PROC_BIND','OMP_PLACES','LOADEDMODULES')})
    for key,command in [('cpu',['lscpu']),('gpu',['nvidia-smi','--query-gpu=uuid,name,memory.total,driver_version','--format=csv']),('compiler',['gcc','--version'])]:
        try:result[key]=subprocess.check_output(command,text=True,stderr=subprocess.STDOUT,timeout=15)
        except (OSError,subprocess.SubprocessError) as e:result[key]=str(e)
    return result


def hardware_class(hw):
    """Allow another allocation only with the same measured hardware class."""
    names=('Architecture','Model name','Vendor ID','Thread(s) per core','Core(s) per socket','Socket(s)','L1d cache','L1i cache','L2 cache','L3 cache')
    cpu={name:value.strip() for name,value in re.findall(r'^([^:]+):\s*(.*)$',hw.get('cpu',''),re.M) if name in names}
    gpu=sorted({tuple(v.strip() for v in line.split(',')[1:]) for line in hw.get('gpu','').splitlines()[1:]})
    if not cpu.get('Model name') or not gpu:raise RuntimeError('hardware class cannot be verified')
    return dict(cpu=cpu,gpu=[list(row) for row in gpu],affinity_count=len(hw['affinity']),cpus_per_task=hw['environment'].get('SLURM_CPUS_PER_TASK'))


class GPUProcessMemory:
    """Driver/NVML process bytes at1Hz; reported Kokkos peak is separate."""
    def __init__(self,pid,enabled):self.pid=pid;self.enabled=enabled;self.stop=threading.Event();self.rows=[];self.errors=[];self.thread=None
    def start(self):
        if self.enabled:self.thread=threading.Thread(target=self.poll,daemon=True);self.thread.start()
    def poll(self):
        while not self.stop.is_set():
            try:
                text=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,stderr=subprocess.STDOUT,timeout=5)
                for line in text.splitlines():
                    values=[value.strip() for value in line.split(',')]
                    if len(values)==3 and values[1]==str(self.pid):self.rows.append(dict(seconds=time.time(),uuid=values[0],bytes=int(values[2])*1024**2))
            except (OSError,ValueError,subprocess.SubprocessError) as e:
                if not self.errors:self.errors.append(str(e))
            self.stop.wait(1)
    def finish(self):
        self.stop.set()
        if self.thread:self.thread.join(timeout=6)
        return dict(enabled=self.enabled,period_seconds=1,samples=self.rows,errors=self.errors,
                    sampled_process_peak_bytes=max((row['bytes'] for row in self.rows),default=None),
                    window='whole worker, including untimed verification and snapshot',
                    available=bool(self.rows),
                    note='Discrete process-level driver observations include CUDA context/runtime. Missing samples are unavailable. They are a lower bound on transient peak and differ from exact Kokkos allocation-callback peaks.')


def manifest_images(variant):
    images={}
    for key in ('hispid_library','by_library'):
        path=Path(variant[key]).resolve(strict=True);images[str(path)]=digest(path)
    for path in variant.get('runtime_images',[]):images[str(Path(path).resolve(strict=True))]=digest(path)
    expected=variant.get('resolved_dependency_images')
    if expected is not None:
        actual={}
        for key in ('hispid_library','by_library'):
            listing=subprocess.check_output(['ldd',variant[key]],text=True)
            if 'not found' in listing:raise RuntimeError('unresolved dynamic dependency')
            for path in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)\s+\(',listing,re.M):actual[str(Path(path).resolve(strict=True))]=digest(path)
        if actual!=expected:raise RuntimeError('resolved dependency paths/hashes changed')
        images.update(actual)
    return images


def linear_log(record,text):
    record['linear_history']=[dict(true_l2=float(v),true_relative_l2=float(rel),absolute_target=float(target)) for v,rel,target in re.findall(r'linear_true: ([^\s]+) relative ([^\s]+) target ([^\s]+)',text)]
    record['linear_iterations']=[max(map(int,re.findall(r'^(?:bicgstab|gmres):\s+(\d+)\s',part,re.M)),default=0) for part in re.split(r'(?:bicgstab|gmres):  itmax.*\n',text)[1:]]


def verify_artifacts(result,artifact_root=None):
    root=ROOT if artifact_root is None else Path(artifact_root)
    for label,row in {**result['records'],**result['failures']}.items():
        for key in ('state','log','worker'):
            if key in row and digest(root/row[key])!=row[key+'_sha256']:
                raise RuntimeError(f'retained {key} changed: {label}')
    for path,sha in result.get('input_sha256',{}).items():
        if digest(root/path)!=sha:raise RuntimeError('retained input changed: '+path)


def protocol_checks(record,variant,images,expected_config):
    mode=record['mode'];method=record['krylov'];c=checks(mode,record)
    c['execution']=record['execution']==variant['execution'] and record['cpu_threads']==variant['threads'] and (variant['execution']=='reference' or record['compiled_execution']==variant['space'])
    primary=str(Path(variant['hispid_library' if mode=='hispid' else 'by_library']).resolve())
    dependencies=record['dependency_images']
    c['images']=record['loaded_image_verified'] and record['library_sha256']==images[primary] and all(images.get(path)==sha for path,sha in dependencies.items())
    expected_runtime={str(Path(path).resolve()):images[str(Path(path).resolve())] for path in variant.get('runtime_images',[])}
    c['runtime_images']=record.get('runtime_images')==expected_runtime
    cfg=record['config'];c['configuration']=cfg==expected_config and cfg['tolerance']==1e-12 and cfg['max_newton']==24 and cfg['max_krylov']==2000 and cfg['krylov_restart']==64 and cfg['memory_limit_mib']==32768 and cfg['n']==record['grid'] and record['initial_guess']=='zero' and record['linear_rtol']==.001
    if mode=='by':
        resolved=record.get('resolved_options',{});integers=resolved.get('integers',{});reals=resolved.get('reals',{})
        expected=dict(TP_krylov_solver=0 if method=='gmres' else 1,TP_preconditioner=1,TP_linear_relative=1,TP_krylov_restart=64,TP_krylov_maxit=2000,Newton_maxit=24,give_bare_mass=1,use_external_initial_guess=0,solve_momentum_constraint=0,grid_setup_method=1)
        if variant['execution']=='kokkos':expected.update(TP_execution_backend=1,TP_execution_memory_limit_mib=32768)
        expected.update(zip(('npoints_A','npoints_B','npoints_phi'),cfg['n']))
        expected_reals=dict(par_b=expected_config['hole'][0]['center'][0],par_m_plus=expected_config['hole'][0]['mass'],par_m_minus=expected_config['hole'][1]['mass'],Newton_tol=1e-12,TP_linear_rtol=.001)
        for side,hole in zip(('plus','minus'),expected_config['hole']):
            velocity=np.asarray(hole['velocity']);lorentz=1/np.sqrt(1-velocity@velocity);rest_spin=np.asarray(hole['spin'])
            momentum=hole['mass']*lorentz*velocity;spin=lorentz*rest_spin-lorentz**2/(1+lorentz)*velocity*(velocity@rest_spin)
            for axis in range(3):expected_reals[f'par_P_{side}{axis+1}']=float(momentum[axis]);expected_reals[f'par_S_{side}{axis+1}']=float(spin[axis])
        c['resolved_options']=all(integers.get(k)==v for k,v in expected.items()) and reals==expected_reals
    else:
        actual=record.get('resolved_options',{})
        c['resolved_options']=actual.get('krylov')==method and actual.get('linear_rtol')==.001 and actual.get('preconditioner')=='modal' and (variant['execution']=='reference' or actual.get('native_verified') is True)
    if variant.get('space')=='Cuda':
        device=record.get('device') or {};c['one_visible_gpu']=device.get('visible_count')==1 and device.get('visible_ordinal')==0 and bool(device.get('uuid'))
        observed={row['uuid'] for row in record['driver_memory']['samples']}
        # Device API supplies positive single-device evidence independently of
        # optional NVML sampling. Any observed other UUID is a failure.
        c['telemetry_identity']=not observed or observed=={device.get('uuid')}
    c['passed']=all(v for k,v in c.items() if k!='passed');return c


def finalize(result,expected):
    verify_artifacts(result)
    variants={v['id']:v for v in result['manifest']['variants']}
    result['comparisons']={};missing=[]
    for label,record in result['records'].items():
        record['checks']=protocol_checks(record,variants[record['variant']],result['images'][record['variant']],result['configurations'][':'.join(map(str,record['grid']))])
        grid_id='_'.join(map(str,record['grid']));reference_label=f'{grid_id}_reference_{record["mode"]}_{record["krylov"]}_0'
        if reference_label in result['records']:
            reference=result['records'][reference_label]
            result['comparisons'][label]=(hi_compare if record['mode']=='hispid' else by_compare)(reference,record)
        else:missing.append(label)
    result['comparison_missing_reference']=missing
    result['completed_workers']=len(result['records'])+len(result['failures']);result['expected_workers']=len(expected)
    result['coverage_complete']=set(result['records'])|set(result['failures'])==set(expected)
    binding=result.get('binding',{})
    result['declared_performance_scope']=binding.get('grids')==[[40,80,16],[80,160,16],[128,256,28]] and binding.get('repeats')==3 and set(binding.get('systems',[]))=={'hispid','by'} and set(binding.get('methods',[]))=={'gmres','bicgstab'} and set(variants)=={'reference','serial','openmp1','openmp2','openmp4','openmp8','openmp16','cuda'}
    result['declared_performance_completed']=result['declared_performance_scope'] and result['coverage_complete']
    result['all_stopping_checks_passed']=result['coverage_complete'] and not result['failures'] and bool(result['records']) and all(r['checks']['passed'] for r in result['records'].values())
    result['strict_state_comparisons_passed']=result['coverage_complete'] and bool(expected) and not result['failures'] and not missing and set(result['comparisons'])==set(expected) and all(c['passed'] for c in result['comparisons'].values())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True);p.add_argument('--output',required=True)
    p.add_argument('--grids',default='40:80:16,80:160:16,128:256:28');p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--systems',default='hispid,by');p.add_argument('--methods',default='gmres,bicgstab')
    p.add_argument('--timeout',type=float,default=3600);p.add_argument('--resume',action='store_true');args=p.parse_args()
    manifest=json.loads(Path(args.manifest).read_text());variants=manifest['variants'];grids=[tuple(map(int,s.split(':'))) for s in args.grids.split(',')]
    if args.repeats<1 or not np.isfinite(args.timeout) or args.timeout<=0:p.error('positive repeats and timeout required')
    systems=args.systems.split(',');methods=args.methods.split(',')
    if not systems or len(set(systems))!=len(systems) or not set(systems)<={'hispid','by'}:p.error('unique hispid/by systems required')
    if not methods or len(set(methods))!=len(methods) or not set(methods)<={'gmres','bicgstab'}:p.error('unique gmres/bicgstab methods required')
    protocols=[(system,method) for system in systems for method in methods]
    if not grids or len(set(grids))!=len(grids) or any(len(g)!=3 or min(g)<4 or max(g)>256 or g[2]%2 for g in grids):p.error('unique three-axis grids4..256 with even phi required')
    if len({v['id'] for v in variants})!=len(variants) or sum(v['id']=='reference' for v in variants)!=1:p.error('unique variants with one reference required')
    out=Path(args.output).resolve();raw=ROOT/'validation/raw'/out.stem
    if out.exists() and not args.resume:raise FileExistsError(out)
    if raw.exists() and not args.resume:raise FileExistsError(raw)
    raw.mkdir(parents=True,exist_ok=True)
    resolved_images={v['id']:manifest_images(v) for v in variants}
    sources=[ROOT/'validation'/name for name in ('benchmark_kokkos.py','benchmark_bowen_york.py','benchmark_krylov_matrix.py','benchmark_common_stopping.py')]
    sources+=list((ROOT/'python').glob('*.py'))+list((ROOT/'examples').glob('*.py'))
    binding=dict(grids=[list(g) for g in grids],systems=systems,methods=methods,repeats=args.repeats,timeout=args.timeout,acceptance_sha256=digest(ROOT/'validation/kokkos_acceptance.json'),source_sha256={str(f.relative_to(ROOT)):digest(f) for f in sources})
    expected=[f'{"_".join(map(str,g))}_{v["id"]}_{mode}_{method}_{i}' for g in grids for i in range(args.repeats) for v in variants for mode,method in protocols]
    result=json.loads(out.read_text()) if args.resume and out.exists() else dict(hardware=hardware(),manifest=manifest,manifest_sha256=digest(args.manifest),acceptance_sha256=binding['acceptance_sha256'],binding=binding,images=resolved_images,records={},failures={},comparisons={},input_sha256={})
    if result['images']!=resolved_images or result['manifest_sha256']!=digest(args.manifest) or result.get('binding')!=binding:raise RuntimeError('resume images/protocol/sources changed')
    verify_artifacts(result)
    if (set(result['records'])|set(result['failures']))-set(expected):raise RuntimeError('unexpected retained workers')
    invocation_hardware=hardware();epochs=result.setdefault('allocation_epochs',[])
    compatible=hardware_class(invocation_hardware)
    if 'hardware_class' in result and result['hardware_class']!=compatible:raise RuntimeError('resume hardware/allocated CPU fraction differs; retain separate performance groups')
    result['hardware_class']=compatible
    if not epochs or epochs[-1]!=invocation_hardware:epochs.append(invocation_hardware)
    epoch=len(epochs)-1
    def save():
        temporary=out.with_suffix(out.suffix+'.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(out)
    save()
    # The frozen moderate config is shared with BY, including the declared
    # budget. The Hi worker reconstructs the identical moderate free data.
    from hispid import Backend
    from benchmark_bowen_york import configuration
    backend=Backend(variants[0]['hispid_library'])
    inputs={}
    configurations={}
    for grid in grids:
        c=configuration(backend,'moderate',grid,1e-12);c.memory_limit_mib=32768
        from configs import as_dict
        path=raw/('input_'+'_'.join(map(str,grid))+'.json');content=json.dumps(as_dict(c),indent=2)+'\n'
        if path.exists() and path.read_text()!=content:raise RuntimeError('existing input differs')
        if not path.exists():path.write_text(content)
        inputs[grid]=path;result['input_sha256'][str(path.relative_to(ROOT))]=digest(path)
        configurations[':'.join(map(str,grid))]=as_dict(c)
    del backend
    if 'configurations' in result and result['configurations']!=configurations:raise RuntimeError('resume physical configurations changed')
    result['configurations']=configurations
    save()
    for grid in grids:
        grid_id='_'.join(map(str,grid))
        for repeat in range(args.repeats):
            ordered=variants if repeat%2==0 else list(reversed(variants))
            for variant in ordered:
                for mode,method in protocols if repeat%2==0 else list(reversed(protocols)):
                    label=f'{grid_id}_{variant["id"]}_{mode}_{method}_{repeat}'
                    if label in result['records'] or label in result['failures']:continue
                    # Administrative stop between workers preserves the full
                    # declared protocol and permits another shared allocation.
                    if (raw/'STOP_AFTER_WORKER').exists():
                        finalize(result,expected);save();print('Checkpointed between workers',flush=True);return
                    if manifest_images(variant)!=resolved_images[variant['id']]:raise RuntimeError('worker image changed')
                    worker=raw/(label+'.json');state=raw/(label+'.npz');log=raw/(label+'.log')
                    # Preserve an interrupted attempt which was not yet added
                    # to the atomic checkpoint before rerunning its worker.
                    leftovers=[path for path in (worker,state,log) if path.exists()]
                    if leftovers:
                        orphan=raw/('interrupted_'+str(time.time_ns())+'_'+label);orphan.mkdir()
                        for path in leftovers:path.rename(orphan/path.name)
                    command=[sys.executable,str(ROOT/'validation/benchmark_bowen_york.py'),'--hispid-library',variant['hispid_library'],'--by-library',variant['by_library'],'--case','moderate','--grid',':'.join(map(str,grid)),'--tolerance','1e-12','--output',str(out),'--mode',mode,'--worker-output',str(worker),'--state-output',str(state),'--execution',variant['execution'],'--threads',str(variant['threads']),'--memory-limit-mib','32768','--krylov',method,'--linear-rtol','.001','--physical-check']
                    if mode=='by':command+=['--input',str(inputs[grid]),'--by-verbose','--by-preconditioner','1','--krylov-maxit','2000']
                    env=os.environ.copy();env.update(OMP_NUM_THREADS=str(variant['threads']),OMP_PROC_BIND='close',OMP_PLACES='cores',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');env.update(variant.get('environment',{}))
                    print('Starting',label,flush=True);started=time.time()
                    with log.open('w') as stream:
                        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT);telemetry=GPUProcessMemory(process.pid,variant.get('space')=='Cuda');telemetry.start()
                        try:code=process.wait(timeout=args.timeout)
                        except subprocess.TimeoutExpired:process.kill();code=process.wait();result['failures'][label]=dict(timeout=True)
                        memory=telemetry.finish()
                    if code or not worker.exists():
                        result['failures'].setdefault(label,{}).update(exit_code=code,allocation_epoch=epoch,log=str(log.relative_to(ROOT)),log_sha256=digest(log),command=command,driver_memory=memory,elapsed_seconds=time.time()-started);print('Failed',label,code,flush=True);save();continue
                    record=json.loads(worker.read_text());record.update(variant=variant['id'],mode=mode,repeat=repeat,grid=list(grid),allocation_epoch=epoch,state=str(state.relative_to(ROOT)),state_sha256=digest(state),log=str(log.relative_to(ROOT)),log_sha256=digest(log),worker=str(worker.relative_to(ROOT)),worker_sha256=digest(worker),command=command,driver_memory=memory)
                    if mode=='by':linear_log(record,log.read_text())
                    record['checks']=protocol_checks(record,variant,resolved_images[variant['id']],configurations[':'.join(map(str,grid))])
                    result['records'][label]=record;save();print('Finished',label,record['solve_seconds'],record['max_rss_bytes'],record['checks'],flush=True)
                    reference_label=f'{grid_id}_reference_{mode}_{method}_0'
                    if reference_label in result['records']:
                        reference=result['records'][reference_label]
                        result['comparisons'][label]=(hi_compare if mode=='hispid' else by_compare)(reference,record);save()
    performance={}
    for label,record in result['records'].items():
        key='_'.join(label.split('_')[:-1]);performance.setdefault(key,[]).append(record)
    result['performance']={key:dict(repeats=len(rows),solve_seconds_median=float(np.median([r['solve_seconds'] for r in rows])),solve_seconds_min=min(r['solve_seconds'] for r in rows),solve_seconds_max=max(r['solve_seconds'] for r in rows),ready_seconds_median=float(np.median([r['ready_to_sample_seconds'] for r in rows])),host_rss_bytes_median=float(np.median([r['max_rss_bytes'] for r in rows]))) for key,rows in performance.items()}
    finalize(result,expected)
    save();print('Completed',result['completed_workers'],'of',result['expected_workers'],flush=True)
if __name__=='__main__':main()
