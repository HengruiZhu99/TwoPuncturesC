"""Generate a standalone, data-bound LaTeX performance report (no TeX runtime)."""
import argparse,hashlib,json,re,statistics
from collections import defaultdict
from pathlib import Path
from benchmark_bowen_york import digest
from benchmark_kokkos import verify_artifacts


def tex(value):
    return str(value).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#')

def median(rows,key):return statistics.median(row[key] for row in rows)
def number(value):return '--' if value is None else f'{value:.3g}'
def work(row):
    stats=row['diagnostics'] if row['mode']=='hispid' else row['work_statistics']
    return stats['newton_iterations'],stats['krylov_iterations'],row['work_statistics']['jvp_applications'],row['work_statistics']['preconditioner_applications']

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    payload=Path(args.results).read_bytes();source_sha=hashlib.sha256(payload).hexdigest()
    result=json.loads(payload);verify_artifacts(result);records=result['records'];groups=defaultdict(list)
    for row in records.values():groups[(tuple(row['grid']),row['mode'],row['krylov'],row['variant'])].append(row)
    grids=sorted({key[0] for key in groups});variants=[v['id'] for v in result['manifest']['variants']]
    failed_comparisons=[(label,c) for label,c in result['comparisons'].items() if not c['passed']]
    doc=[r'''\documentclass[11pt]{article}
\usepackage[margin=0.8in]{geometry}
\usepackage{amsmath,booktabs,longtable,array,hyperref}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\setlength{\tabcolsep}{3pt}
\title{Bowen--York and HiSpID: Kokkos performance and numerical qualification}
\author{Isolated TwoPuncturesC implementation}
\date{October 2026}
\begin{document}\maketitle
\begin{abstract}
The performance protocol compares two equation systems, two Krylov methods, CPU execution
and one allocated NVIDIA A100 at matched computational resolution and
stopping norms. Matched grid counts and input seed charges do not equate
physical accuracy or horizon properties. Internal convergence, complete-state
preservation, independent physical checks and performance are separate results.
No new extremal-spin or relativistic binary acceptance is inferred here.
\end{abstract}
\section{Scope and numerical controls}
Bowen--York solves one Hamiltonian equation with analytic momentum data.
HiSpID solves four coupled correction equations: one Hamiltonian and three
momentum constraints for prescribed curved conformal geometry and mean
curvature. The moderate binary has unequal rest masses 0.6/0.4, centers
$\pm3$, generic spin/boost axes, unchanged attenuation and no far filter.
BY input momenta and intrinsic spins match the isolated lab-frame seed
charges; the two seed geometries and their final ADM/horizon quantities differ.
Both use cubic residual row weighting, outer weighted $L^\infty$ tolerance
$10^{-12}$, true RHS-relative inner $L^2$ tolerance $10^{-3}$, restart 64,
maximum 2000 iterations and a zero linear initial guess. GMRES and BiCGStab
share one controller and use the same modal FD approximation within each
system. True inner residuals are recomputed; finite states and failure
counts are checked independently of timing.

The reference is frozen commit \texttt{25ca064}. Its explicit maximum budget
is extended from 8192 to 65536 MiB for this benchmark, with the allocation
formula, equations and stopping arithmetic unchanged. Every variant requests
32768 MiB. The 128$\times$256$\times$28/restart 64 tier was rejected by the
original 8 GiB guard; it does not fit that unchanged limit. The patch and full
build/source/image hashes are bound to the measurement file.

The historical Makefile defaults were checked separately against archived
complete states, work counts and the BY trace on the original platform.
Those controls pass bitwise. Original TwoPunctures fixtures, including
target-mass iteration, pass with both methods. This does not demand bitwise
cross-platform/GPU reductions. The opt-in comparisons retain the declared
$10^{-10}$ scaled unknown/coefficient and physical-field gates. Hi retains
physical metric-gradient and ADM energy/momentum/angular-momentum checks;
BY retains ADM energy and puncture-end masses, without a retained metric-gradient
array or ADM momentum/angular-momentum comparison. PDE-kernel relative-$L^2$
controls use $10^{-10}$; the conditioned
production preconditioner comparison/action control uses a separately named
$10^{-9}$ bound. Independent Cartesian FD constraint differences are bounded
by $\max(10^{-7},5\epsilon_{\rm step})$ with steps 0.004/0.002/0.001.
\section{Hardware and measurement definitions}
''']
    if not result.get('declared_performance_completed'):
        doc[0]=doc[0].replace(r'\begin{abstract}',r'\noindent\textbf{Incomplete measurement snapshot.} The full declared performance protocol is pending.\par'+'\n'+r'\begin{abstract}')
    binding=result.get('binding',{})
    doc.append('Current measurement binding: '+tex(', '.join(binding.get('systems',[])))+', '+str(len(binding.get('grids',[])))+' grid(s), '+str(len(variants))+' variant(s), '+str(binding.get('repeats','unavailable'))+' repeat(s).'+r'\\'+'\n')
    doc.append('Source result: '+r'\texttt{'+tex(Path(args.results).name)+'}.'+r'\\'+'\n')
    hw=result['hardware'];cpu=re.search(r'^Model name:\s*(.+)$',hw.get('cpu',''),re.M)
    doc.append('Initial host: '+tex(hw['hostname'])+'. CPU: '+tex(cpu[1] if cpu else 'unavailable')+'.'+r'\\'+'\n')
    compiler=hw.get('compiler','unavailable').splitlines()[0]
    doc.append('Compiler command version: '+tex(compiler)+'.'+r'\\'+'\n')
    env=hw['environment'];doc.append('Initial Slurm allocation: '+tex(env.get('SLURM_JOB_ID'))+'. Visible GPU setting: '+tex(env.get('CUDA_VISIBLE_DEVICES'))+'.'+r'\\'+'\n')
    for epoch,hardware in enumerate(result.get('allocation_epochs',[])):
        doc.append('Allocation epoch '+str(epoch)+': '+tex(hardware['hostname'])+', Slurm '+tex(hardware['environment'].get('SLURM_JOB_ID'))+'.'+r'\\'+'\n')
    devices={json.dumps({key:row['device'][key] for key in ('name','uuid','total_bytes')},sort_keys=True) for row in records.values() if row.get('device')}
    for description in sorted(devices):
        device=json.loads(description);doc.append('Device: '+tex(device['name'])+', '+number(device['total_bytes']/2**30)+r' GiB; UUID \texttt{'+tex(device['uuid'])+'}.'+r'\\'+'\n')
    doc.append(r'''The one-GPU shared allocation supplies 16 physical CPU cores (32 logical
CPUs). OpenMP uses 1/2/4/8/16 host threads bound to cores; the CUDA variant
uses the same 16-core host fraction. Reference and Serial use one thread.
CPU images share compiler/Release flags; CUDA uses the same host compiler
through nvcc\_wrapper. C remains C99; Kokkos is pinned 4.7.2. No fast-math flag
is requested. Shared Kokkos runtime images, complete CMake caches, actual
concurrency, affinity and source hashes are retained.
Resumed allocations must match the recorded CPU/cache topology, available
core count and GPU model/memory/driver class. Exact node, GPU UUID and
affinity are retained per allocation epoch; incompatible classes are rejected
instead of pooled. Shared-node frequency/load variation appears in the repeats.

Fresh workers run sequentially in interleaved forward/reverse order.
Solve time excludes creation and first sampling. Ready time includes both;
cold time below additionally includes explicit runtime initialization, but
excludes Python import and DSO loading. All device phases are fenced.
Kokkos A/M timings include finite checks; reference BY phase wrappers time
the kernel calls without the controller's finite scans. Linear time includes
A/M. Hi setup
includes native factor construction and device import. BY reports native
modal setup separately from its Kokkos setup. BY verbose iteration
logging and flushing are inside the Newton timer. The F phase covers Newton
evaluations; the final native residual outside Newton is included in ready
time. Setup includes transfers,
which are also reported separately; overlapping phase times must not be added.

Host RSS is the process peak through first sampling, before untimed verification
and snapshot allocation. Reported Kokkos peaks are allocations observed after
owned runtime initialization, including setup overlap; they exclude native
std::vector/GSL storage and initialization/runtime allocations. Pinned memory
is host; unified spaces are classified by their logical device space. Driver
process memory at 1 Hz covers the whole worker, including untimed verification
and snapshot allocation. It includes CUDA context/runtime and is a sampled lower
bound on that window's transient peak. Missing driver samples are unavailable;
the device API independently verifies one visible GPU and its selected UUID.
Logical resident estimates and API copy counters
are diagnostics, not substitutes for measured RAM/VRAM.
The creation/other column below is HiSpID context creation, or for BY the
total data-construction time minus its Newton timer (allocation and native
post-solve work are not individually instrumented). Transfer timings overlap
setup/operator phases; API copy volumes exclude runtime-internal traffic.
\section{Completed measurements}
''')
    completed=len(records)+len(result['failures'])
    expected=len(binding.get('grids',[]))*binding.get('repeats',0)*len(binding.get('systems',[]))*len(binding.get('methods',[]))*len(variants)
    doc.append(f"Completed workers: {completed} of {expected}. Retained worker failures: {len(result['failures'])}. Strict state comparison failures: {len(failed_comparisons)}."+r'\\'+'\n')
    if completed!=expected:doc.append(r'\textbf{This is an incomplete measurement snapshot; pending rows are not estimated.}\par'+'\n')
    if not result.get('declared_performance_completed'):doc.append(r'\textbf{The full declared three-grid, eight-variant, three-repeat protocol is not yet complete. This snapshot cannot substitute for the comprehensive report.}\par'+'\n')
    doc.append(r'\begin{center}\small\begin{tabular}{llrrrr}\toprule Grid & System/method & OMP16 speedup & GPU/OMP16 & GPU cold speedup & OMP RSS saved\\\midrule'+'\n')
    for grid in grids:
        for mode in ('hispid','by'):
            for method in ('gmres','bicgstab'):
                ref=groups.get((grid,mode,method,'reference'),[]);omp=groups.get((grid,mode,method,'openmp16'),[]);gpu=groups.get((grid,mode,method,'cuda'),[])
                if not (ref and omp and gpu) or not all(row['checks']['passed'] for row in ref+omp+gpu):continue
                cold=lambda rows:statistics.median(row['initialization_seconds']+row['ready_to_sample_seconds'] for row in rows)
                values=[r'$\times$'.join(map(str,grid)),('Hi' if mode=='hispid' else 'BY')+'/'+('G' if method=='gmres' else 'B'),number(median(ref,'solve_seconds')/median(omp,'solve_seconds')),number(median(omp,'solve_seconds')/median(gpu,'solve_seconds')),number(cold(ref)/cold(gpu)),number(100*(1-median(omp,'max_rss_bytes')/median(ref,'max_rss_bytes')))+r'\%']
                doc.append(' & '.join(values)+r'\\'+'\n')
    doc.append(r'\bottomrule\end{tabular}\end{center}'+'\n')
    doc.append('Ratios above use internally converged solves; the strict Hi raw-P failure remains a separate failed gate. G/B denotes GMRES/BiCGStab. Negative RAM savings mean increased RSS. Timing tables report medians and observed ranges. A dagger marks a timing group containing a failed stopping/protocol check; those timings are diagnostic.\n')
    for grid in grids:
        doc.append(r'\subsection{Grid '+r'$\times$'.join(map(str,grid))+'}\n')
        doc.append(r'''\footnotesize
\begin{longtable}{llrrrrrr}
\toprule System/method & Variant & $n$ & Solve(s) & Ready(s) & Cold(s) & RSS(GiB) & Device(GiB)\\
\midrule\endhead
''')
        for mode in ('hispid','by'):
            for method in ('gmres','bicgstab'):
                for variant in variants:
                    rows=groups.get((grid,mode,method,variant),[])
                    if not rows:continue
                    device_peaks=[r['execution_statistics']['kokkos_device_peak_bytes']/2**30 for r in rows if r.get('execution_statistics') and r['execution_statistics']['memory_tracking_available']]
                    label=tex(variant)+(r'$\dagger$' if not all(r['checks']['passed'] for r in rows) else '')
                    timing=number(median(rows,'solve_seconds'))+' ['+number(min(r['solve_seconds'] for r in rows))+','+number(max(r['solve_seconds'] for r in rows))+']'
                    cells=[('Hi' if mode=='hispid' else 'BY')+'/'+('G' if method=='gmres' else 'B'),label,str(len(rows)),timing,number(median(rows,'ready_to_sample_seconds')),number(statistics.median(r['initialization_seconds']+r['ready_to_sample_seconds'] for r in rows)),number(median(rows,'max_rss_bytes')/2**30),number(statistics.median(device_peaks)) if device_peaks else '--']
                    doc.append(' & '.join(cells)+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
        doc.append(r'''\scriptsize\begin{longtable}{llrrrrrrr}
\toprule System/method & Variant & Init(s) & Create/other(s) & Sample(s) & Transfer(s) & H2D(GiB) & D2H(GiB) & K host(GiB)\\
\midrule\endhead
''')
        for mode in ('hispid','by'):
            for method in ('gmres','bicgstab'):
                for variant in variants:
                    rows=groups.get((grid,mode,method,variant),[])
                    if not rows:continue
                    stats=[r['execution_statistics'] for r in rows if r.get('execution_statistics') and r['execution']=='kokkos']
                    creation=statistics.median(r['creation_seconds'] if mode=='hispid' else r['setup_and_solve_seconds']-r['solve_seconds'] for r in rows)
                    transfer=[number(statistics.median(s[key] for s in stats)/scale) if stats else '--' for key,scale in (('transfer_seconds',1),('host_to_device_bytes',2**30),('device_to_host_bytes',2**30))]
                    host=[s['kokkos_host_peak_bytes']/2**30 for s in stats if s['memory_tracking_available']]
                    cells=[('Hi' if mode=='hispid' else 'BY')+'/'+('G' if method=='gmres' else 'B'),tex(variant),number(median(rows,'initialization_seconds')),number(creation),number(median(rows,'first_sample_seconds')),*transfer,number(statistics.median(host)) if host else '--']
                    doc.append(' & '.join(cells)+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
        doc.append(r'''\scriptsize\begin{longtable}{llrrrrrrrr}
\toprule System/method & Variant & Newton/K & A/M calls & A(s) & M(s) & Kokkos setup(s) & Native M(s) & F(s) & Driver(GiB)\\
\midrule\endhead
''')
        for mode in ('hispid','by'):
            for method in ('gmres','bicgstab'):
                for variant in variants:
                    rows=groups.get((grid,mode,method,variant),[])
                    if not rows:continue
                    counts=[work(row) for row in rows];count_text=lambda a,b:'/'.join(number(statistics.median(c[i] for c in counts)) for i in (a,b))
                    stats=[r['execution_statistics'] for r in rows if r.get('execution_statistics') and r['execution']=='kokkos']
                    phases=[number(statistics.median(s[key] for s in stats)) if stats else '--' for key in ('operator_seconds','precondition_seconds','setup_seconds')]
                    native=[r['phase_seconds'] for r in rows if r.get('phase_seconds')]
                    if not stats and native:phases[:2]=[number(statistics.median(s[key] for s in native)) for key in ('spectral_jvp','modal_apply')]
                    native_setup=number(statistics.median(s['modal_setup'] for s in native)) if native else '--'
                    residual=number(statistics.median((r['execution_statistics']['nonlinear_seconds'] if r.get('execution_statistics') and r['execution']=='kokkos' else 0)+r['phase_seconds'].get('native_residual',0) for r in rows)) if mode=='by' and native and all('native_residual' in s for s in native) else '--'
                    driver=[r['driver_memory']['sampled_process_peak_bytes']/2**30 for r in rows if r['driver_memory']['sampled_process_peak_bytes'] is not None]
                    doc.append(' & '.join([('Hi' if mode=='hispid' else 'BY')+'/'+('G' if method=='gmres' else 'B'),tex(variant),count_text(0,1),count_text(2,3),*phases,native_setup,residual,number(statistics.median(driver)) if driver else '--'])+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
    doc.append(r'\section{Numerical comparison results}'+'\n')
    doc.append('Worker stopping checks: '+('all passed' if result.get('all_stopping_checks_passed') else 'pending or failed; inspect retained records')+'.\n')
    doc.append(r'''The existing HiSpID weak raw-modal/coefficient preservation failure remains
failed. Agreement of sampled physical metric, curvature, gradients and
charges does not waive that criterion. Kernel qualification and internally
converged timings are therefore reported separately from full port acceptance.
''')
    if failed_comparisons:
        doc.append(r'\small\begin{longtable}{lrrr}\toprule Failed comparison & Raw-P scaled max & Coefficient scaled max & Residual max\\\midrule\endhead'+'\n')
        for label,comparison in failed_comparisons:
            unknown=comparison.get('arrays',{}).get('unknowns',{}).get('max_scaled_difference');coefficient=comparison.get('coefficients_scaled_difference');residual=comparison.get('max_recomputed_weighted_F')
            doc.append(tex(label)+' & '+number(unknown)+' & '+number(coefficient)+' & '+number(residual)+r'\\'+'\n')
        doc.append(r'\bottomrule\end{longtable}\normalsize'+'\n')
    if result['failures']:
        doc.append(r'\paragraph{Failed workers.} '+', '.join(tex(key) for key in result['failures'])+'.\n')
    doc.append(r'''\section{Interpretation and subsequent physical study}
Dense cached differentiation removes repeated transform/trigonometric setup
from BY's linear action and repeated nonlinear residual evaluations. A scoped
workspace reuses immutable seed and differentiation data throughout Newton.
An original-residual confirmation retains the outer stopping criterion;
if needed, subsequent outer residuals use the original routine for polishing.
The F column includes both device evaluations and native confirmation/polishing.
HiSpID avoids replicated modal stencils, duplicate
host factor banks and unused reference derivative storage; CPU shares
geometry, while CUDA discards its staging geometry after upload. Fourier
partners and vector components share factors. Inner Krylov vectors stay in
the selected execution space; CPU seed construction and native block
factorization still limit end-to-end GPU speed. A speed ratio is meaningful
only with the work counts, stopping checks and memory columns alongside it.

The user selected separate equal-rest-mass aligned seed spin $\chi=0.99$ and
nonspinning inward head-on input $\Gamma=10$ investigations. The latter uses
$v=\sqrt{0.99}$ and coordinate separation about 50 times one hole's measured
$M_{\rm irr}=\sqrt{A/(16\pi)}$, rather than total ADM mass. Measured horizon
spins, masses and momenta will be distinguished from seed parameters.
For nonspinning isolated seeds with $m=0.5$, $M_{\rm irr}=0.5$, so the
initial choice is $d\simeq25$. At input $\Gamma=10$ the total isolated seed
energy is about10, giving $d/E_{\rm seed}\simeq2.5$. Binary component masses
and global ADM energy must be measured, and changed-separation calibration
cases retained individually. This differs from the thesis Table4.3
$d/M_{\rm ADM}=100$--400 regime chosen to approximate isolated holes
\cite{thesis}. Failed common-horizon searches do not prove absence.
New physical runs follow the completed performance report. They require
fresh seed controls, at least three resolutions, independent exterior
constraint convergence, charges/covariance, initial-time AthenaK import,
individual/common horizon searches and enclosure of modified regions.
No high-parameter validation is inherited from the performance case.
\section{Reproducibility and sources}
''')
    doc.append('Measurement JSON SHA256:'+r'\par{\footnotesize\ttfamily '+source_sha+r'}\par'+'\n')
    doc.append('Acceptance plan SHA256:'+r'\par{\footnotesize\ttfamily '+tex(result['acceptance_sha256'])+r'}\par'+'\n')
    doc.append(r'''The immutable measurement manifest contains complete build caches,
source/runtime/native-image hashes, commands, resolved options, configurations,
true stopping histories, physical differences, full-state digests and retained
logs. The benchmark is restartable only while those bound images remain
unchanged. Source controls and failed comparisons are retained with the report.

\begin{thebibliography}{9}
\bibitem{hispid} Ruchlin et al., \emph{Puncture Initial Data for Black-Hole
Binaries with High Spins and High Boosts}, \url{https://arxiv.org/abs/1410.8607}.
\bibitem{thesis} I. Ruchlin, \emph{Puncture Initial Data and Evolution of
Black Hole Binaries with High Speed and High Spin}, RIT dissertation,
August2015, Section4.5 and Table4.3 (printed pp.117--119),
\url{https://repository.rit.edu/theses/8797/}.
\bibitem{nersc} NERSC, job policy and interactive resources,
\url{https://docs.nersc.gov/jobs/policy/},
\url{https://docs.nersc.gov/jobs/interactive/}.
\bibitem{arch} NERSC, Perlmutter architecture,
\url{https://docs.nersc.gov/systems/perlmutter/architecture/}.
\bibitem{kokkos} Kokkos, initialization/finalization and reductions,
\url{https://kokkos.org/kokkos-core-wiki/API/core/initialize_finalize/initialize.html},
\url{https://kokkos.org/kokkos-core-wiki/API/core/parallel-dispatch/parallel_reduce.html}.
\end{thebibliography}
\end{document}
''')
    if digest(args.results)!=source_sha:raise RuntimeError('measurement JSON changed during report generation; use an immutable snapshot')
    Path(args.output).write_text(''.join(doc));print(args.output)
if __name__=='__main__':main()
