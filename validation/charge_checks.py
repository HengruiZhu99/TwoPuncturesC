"""Shared, native-library-free acceptance logic for ADM refinement evidence."""
import numpy as np
from physical import extrapolate, charges


def validate_refinements(radii, quadratures, bound=None):
    if (not radii or not np.isfinite(radii).all() or min(radii) <= 0
            or not np.all(np.diff(radii) > 0) or not quadratures
            or any(nt < 2 or np_ < 4 for nt, np_ in quadratures)):
        raise ValueError('finite increasing radii and positive angular quadratures required')
    if bound is None:
        return
    if not np.isfinite(bound) or not 0 < bound <= 1e-5:
        raise ValueError('qualification bound must be finite, positive and at most 1e-5')
    if len(radii) < 5:
        raise ValueError('qualification requires at least two four-radius fit windows')
    if (len(quadratures) < 3 or len(set(quadratures)) != len(quadratures)
            or any(a > c or b > d or (a, b) == (c, d)
                   for (a, b), (c, d) in zip(quadratures, quadratures[1:]))):
        raise ValueError('qualification requires distinct monotone angular refinements')
    if not (any(a < c and b == d for (a, b), (c, d) in zip(quadratures, quadratures[1:]))
            and any(a == c and b < d for (a, b), (c, d) in zip(quadratures, quadratures[1:]))):
        raise ValueError('qualification requires separate polar-only and azimuthal-only refinements')


def radial_fits(item, radii, values, prefix=''):
    if len(radii) >= 5:
        fits = np.asarray([extrapolate(radii[i:i+4], values[i:i+4])
                           for i in range(len(radii)-3)])
        item[prefix+'radial_fit_windows_EPJ'] = fits.tolist()
        item[prefix+'radial_fit_change_EPJ_linf'] = float(np.max(abs(fits[-1]-fits[-2])))
        item[prefix+'radial_reported_fit_vs_finest_window_linf'] = float(np.max(abs(extrapolate(radii,values)-fits[-1])))


def angular_changes(item, previous, prefixes=('', 'native_'), global_origin=False):
    if previous is None:
        return
    item['refined_polar_only'] = item['ntheta'] > previous['ntheta'] and item['nphi'] == previous['nphi']
    item['refined_azimuthal_only'] = item['ntheta'] == previous['ntheta'] and item['nphi'] > previous['nphi']
    keys = ['native_EPJ', 'independent_EPJ'] + [p+'radial_fit_windows_EPJ' for p in prefixes]
    if global_origin:
        keys += ['global_native_EPJ', 'global_independent_EPJ']
    for key in keys:
        if key in item and key in previous:
            item[key+'_angular_change_linf'] = float(np.max(abs(np.asarray(item[key])-np.asarray(previous[key]))))


def qualify_refinement(records, radii, bound, global_origin=False, sampled_rotation=True):
    """Qualify every independent check and the last separate directional steps.

    This reports charge accuracy only; it never transfers solver or binary
    acceptance. Missing, nonfinite and diagnostic-only evidence all fail.
    """
    validate_refinements(radii, [(r['ntheta'], r['nphi']) for r in records], bound)
    check_keys = ['native_vs_independent_linf',
                  'native_vs_independent_radial_intercepts_linf', 'fixed_origin_extrapolation_linf']
    if sampled_rotation:
        check_keys += ['centered_rotation_linf']
    angular_keys = ['native_EPJ_angular_change_linf', 'independent_EPJ_angular_change_linf',
                    'radial_fit_windows_EPJ_angular_change_linf', 'native_radial_fit_windows_EPJ_angular_change_linf']
    radial_keys = ['radial_fit_change_EPJ_linf', 'native_radial_fit_change_EPJ_linf',
                   'radial_reported_fit_vs_finest_window_linf','native_radial_reported_fit_vs_finest_window_linf']
    if global_origin:
        check_keys += ['global_native_vs_independent_linf', 'global_native_vs_independent_radial_intercepts_linf']
        angular_keys += ['global_radial_fit_windows_EPJ_angular_change_linf',
                         'global_native_radial_fit_windows_EPJ_angular_change_linf',
                         'global_native_EPJ_angular_change_linf', 'global_independent_EPJ_angular_change_linf']
        radial_keys += ['global_radial_fit_change_EPJ_linf', 'global_native_radial_fit_change_EPJ_linf',
                        'global_radial_reported_fit_vs_finest_window_linf','global_native_radial_reported_fit_vs_finest_window_linf']
    def within(row, keys):
        values = [row.get(k, np.inf) for k in keys]
        return bool(np.isfinite(values).all() and min(values) >= 0 and max(values) <= bound)
    independent = bool(records and all(within(r, check_keys) for r in records))
    directional = {};last_steps=[]
    for direction in ('polar', 'azimuthal'):
        candidates = [i for i in range(1,len(records)) if
            (records[i]['ntheta'] > records[i-1]['ntheta'] and records[i]['nphi'] == records[i-1]['nphi']
             if direction == 'polar' else
             records[i]['ntheta'] == records[i-1]['ntheta'] and records[i]['nphi'] > records[i-1]['nphi'])]
        directional[direction] = bool(candidates and within(records[candidates[-1]], angular_keys))
        if candidates:last_steps.append(candidates[-1])
    # A later joint refinement cannot invalidate a directional witness silently.
    subsequent = bool(len(last_steps)==2 and all(within(r,angular_keys) for r in records[min(last_steps):]))
    angular = all(directional.values()) and subsequent
    radial = bool(records and within(records[-1], radial_keys))
    return dict(qualification_bound=bound, independent_checks_verified=independent,
                directional_angular_checks=directional, subsequent_angular_checks_verified=subsequent,angular_verified=angular,
                radial_verified=radial,
                qualified_charge_checks=bool(independent and angular and radial))


def measure_origin_quadrature(solution, radii, ntheta, nphi, center):
    """Independent and native integrals about displaced and global origins.

    No elliptic solve is performed here. The caller owns the freshly solved
    context, binding checks, and incremental failure evidence.
    """
    item = dict(ntheta=ntheta, nphi=nphi)
    center = np.asarray(center)
    values = {}
    for prefix, origin in (('', center), ('global_', np.zeros(3))):
        native = np.asarray([solution.charges(r, center=origin, ntheta=ntheta, nphi=nphi) for r in radii])
        independent = np.asarray([charges(solution.sample, r, center=origin, ntheta=ntheta, nphi=nphi) for r in radii])
        values[prefix] = native, independent
        item[prefix+'native_EPJ'] = native.tolist()
        item[prefix+'independent_EPJ'] = independent.tolist()
        item[prefix+'native_vs_independent_linf'] = float(np.max(abs(native-independent)))
        for data, name in ((independent, ''), (native, 'native_')):
            radial_fits(item, radii, data, prefix+name)
            item[prefix+name+'extrapolated_EPJ'] = extrapolate(radii, data).tolist()
        item[prefix+'native_vs_independent_radial_intercepts_linf'] = float(np.max(abs(
            np.asarray(item[prefix+'radial_fit_windows_EPJ'])-np.asarray(item[prefix+'native_radial_fit_windows_EPJ']))))
    discrepancies = []
    for i in (0, 1):
        expected = extrapolate(radii, values[''][i])
        expected[4:] += np.cross(center, expected[1:4])
        discrepancies.append(np.max(abs(extrapolate(radii, values['global_'][i])-expected)))
    item['fixed_origin_extrapolation_linf'] = float(max(discrepancies))
    return item
