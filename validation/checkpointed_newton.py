"""Save accepted nonlinear progress without changing native acceptance gates."""
from pathlib import Path


def solve_checkpointed(solution, path, max_calls, krylov, linear_rtol, progress):
    if solution.config.max_newton != 1 or max_calls < 1:
        raise ValueError('checkpointed continuation requires one native Newton attempt per call')
    path = Path(path)
    pending = path.with_name(path.name + '.pending')
    attempts = []
    for _ in range(max_calls):
        diagnostics = solution.solve(krylov=krylov, linear_rtol=linear_rtol)
        attempts.append(dict(diagnostics=diagnostics, linear_rtol=linear_rtol,
                             linear_history=solution.linear_history(), work=solution.work_statistics()))
        receipt = solution.write_checkpoint(pending, 'preliminary' if diagnostics['converged'] else 'diagnostic')
        pending.replace(path)
        receipt['path'] = str(path.resolve())
        progress(attempts, receipt)
        # A limit after one native Newton attempt means its line search
        # accepted the update. Any linear/line-search failure stops here.
        if diagnostics['converged'] or diagnostics['error'] != 'Newton iteration limit':
            break
    return diagnostics, attempts
