# HiSpID validation report

All numerical jobs use one CPU thread. Configurations, verifier steps, tensor samples, iterations and runtimes are retained in `validation/results.json` and ignored `validation/raw/`. The native BY baseline is unchanged.

The independent verifier differentiates physical gamma/K only; g is metadata for the region bins. Reported H and the physical momentum norm are in total seed rest-mass units. All binary examples here have total seed rest mass 1. RMS means an average over fixed off-grid points, not a volume L2 norm. Residual maxima and normalized ratios remain in JSON; normalization denominator floors are 1e-8. The normalized momentum ratio is uninformative for maximal data, so absolute physical norms define acceptance.

Numerical evidence is tied to each record's native library SHA. Older failed-source cases remain visible. API migration compares fresh, separate native processes with loaded-image checks; checkpoint guards otherwise remain strict. The original same-process comparison was invalid because dyld reused an archived image with an identical install name. It and the affected attempted revalidation are preserved and explicitly withdrawn.

Current Cartesian-regular modal-P binaries remain unaccepted. The leading current source uses mapped Chebyshev coordinates and exact modal FD block elimination; old nodal-V binary results are historical and do not establish axis regularity. All tables retain their own source fingerprints and failed flags.

## Seed controls

| Seed | H RMS | M RMS | Maximum charge error | Passed |
|---|---:|---:|---:|---|
| Schwarzschild | 3.13143e-09 | 0 | 1.37406e-07 | True |
| Kerr | 5.16231e-09 | 6.66172e-11 | 2.36267e-07 | True |
| boosted_Schwarzschild | 1.44699e-08 | 1.42648e-12 | 6.70066e-07 | True |
| boosted_Kerr | 8.43e-09 | 5.61077e-11 | 2.40416e-06 | True |

## Revised isolated seed targets

Seed rest spin chi=.95 and lab speed v=.885, separately and together with generic directions. This does not accept a solved binary. Independent Cartesian constraints use three verifier step sizes; charge fits use increasing extraction radii with separately refined angular quadrature.

| Case | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | Maximum charge error | Passed |
|---|---:|---:|---:|---:|---:|---|
| spin95 | 8.31899e-09 | 5.64663e-13 | 4.35677e-09 | 5.46464e-15 | 2.2496e-09 | True |
| boost885 | 6.89311e-09 | 1.69182e-13 | 6.3375e-09 | 7.33012e-14 | 4.71514e-08 | True |
| spin95_boost885_generic | 1.33124e-08 | 1.66436e-12 | 6.38599e-09 | 8.36343e-14 | 2.22697e-07 | True |

## moderate

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000325103 | 0.000263127 | 4.44205e-06 | 2.51154e-06 | 0.994515 | 7.20466 |
| 40×40×20 | 0.000128889 | 7.1352e-05 | 5.47675e-06 | 1.85797e-06 | 0.994471 | 51.8315 |
| 56×56×28 | 2.55822e-05 | 9.22112e-06 | 2.41992e-06 | 1.52064e-07 | 0.994481 | 205.766 |
| 72×72×32 | 6.60347e-06 | 1.45097e-05 | 1.83275e-07 | 1.22425e-07 | 0.994481 | 572.868 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00136514 | 4.52514e-05 | 0.00132658 | 0.00100694 | 1.19062e-12, 3.9274e-14, 2.43213e-13, 7.4739e-14 |
| 40×40×20 | 0.0013713 | 3.43756e-05 | 0.000567779 | 0.000224773 | 2.49386e-11, 7.63918e-13, 4.34337e-12, 1.56632e-12 |
| 56×56×28 | 0.001373 | 3.32953e-05 | 7.46973e-05 | 4.05882e-05 | 3.2019e-10, 6.76789e-12, 2.95003e-11, 1.43424e-11 |
| 72×72×32 | 0.00137448 | 3.32681e-05 | 3.1667e-05 | 6.05114e-05 | 1.14053e-09, 6.18279e-11, 1.68026e-10, 6.65225e-11 |

Finest source SHA: `ba1694e312e92fb322d81986288784c66cd85e4866809feeb5673575e1450a19`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    72,
    72,
    32
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 40.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944807417572294, 0.010742973411407784, 0.006441496362262738, 0.019544421436106657, 0.03998926725237321, 0.11627048068136969, 0.3377570683800642]`.

## moderate_far0

Preliminary gate: **True**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000325316 | 0.000261737 | 5.16796e-07 | 1.69787e-06 | 0.994498 | 7.22801 |
| 40×40×20 | 0.000119528 | 7.13055e-05 | 8.49785e-08 | 8.71704e-07 | 0.994489 | 48.8798 |
| 56×56×28 | 1.49731e-05 | 9.12169e-06 | 1.13797e-08 | 6.7693e-08 | 0.994489 | 197.904 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00136521 | 4.52545e-05 | 0.00134455 | 0.00100432 | 1.57559e-12, 6.14728e-14, 1.04542e-13, 8.47985e-14 |
| 40×40×20 | 0.00137136 | 3.43781e-05 | 0.000564469 | 0.000225496 | 2.21609e-11, 1.10236e-12, 3.46226e-12, 2.53461e-12 |
| 56×56×28 | 0.00137305 | 3.32975e-05 | 4.77916e-05 | 3.99333e-05 | 2.35495e-10, 8.56549e-12, 3.67172e-11, 2.31978e-11 |

Finest source SHA: `0016203beecd6ff4407d7e6cdbebc2b948dbba72eeaacf3cf87f7e6348533329`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    56,
    56,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944889341908381, 0.010847878260691978, 0.006521486613273054, 0.01971248418084202, 0.03998971736750021, 0.1161460963874341, 0.3377684348634569]`.

## highspin

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 48×48×12 | 0.000527305 | 1.86185e-05 | 4.64286e-07 | 1.21506e-08 | 0.980231 | 40.21 |
| 80×80×16 | 9.6266e-05 | 2.03393e-06 | 9.94255e-09 | 3.05343e-08 | 0.980229 | 82.2513 |
| 112×112×16 | 7.27344e-06 | 3.14436e-06 | 1.44404e-09 | 1.33313e-07 | 0.980229 | 274.196 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 48×48×12 | no support | no support | 0.00199908 | 6.12802e-05 | 4.99824e-09, 2.86305e-13, 3.27015e-12, 6.58186e-13 |
| 80×80×16 | no support | no support | 0.000495192 | 8.24278e-06 | 3.00095e-09, 1.29895e-11, 1.05163e-10, 5.85879e-12 |
| 112×112×16 | no support | no support | 3.40912e-05 | 1.35724e-05 | 1.64302e-09, 2.73657e-12, 8.2647e-13, 2.28e-13 |

Finest source SHA: `0016203beecd6ff4407d7e6cdbebc2b948dbba72eeaacf3cf87f7e6348533329`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.5,
      "center": [
        6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    },
    {
      "mass": 0.5,
      "center": [
        -6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    }
  ],
  "n": [
    112,
    112,
    16
  ],
  "conformal_choice": 0,
  "inner_flatten": 0,
  "omega": [
    0.2,
    0.2
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.0,
    0.0
  ],
  "inner_max": [
    0.0,
    0.0
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9802286763840211, -4.1439550862933305e-09, 3.87855942817875e-08, 4.36365295652967e-08, 7.375625559488769e-08, -7.021529886520288e-06, 0.4948777118417689]`.

## highspin_angular80

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×80×24 | 9.62974e-05 | 1.93369e-06 | 1.01356e-08 | 1.02666e-08 | 0.980229 | 116.329 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×80×24 | no support | no support | 0.000495106 | 8.77315e-06 | 1.12932e-09, 3.18705e-12, 5.3143e-13, 1.38571e-11 |

Finest source SHA: `0016203beecd6ff4407d7e6cdbebc2b948dbba72eeaacf3cf87f7e6348533329`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Reference comparison:

```json
{
  "source": "thesis Table3.1 HS99UU",
  "ADM_energy": 0.980124,
  "absolute_tolerance": 0.0002,
  "energy_error": 0.00010468301373012245,
  "energy_agreement": true,
  "exact_historical_reproduction": false,
  "departure": "modern superposed-metric trace projection; horizon mass/spin unmeasured"
}
```

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.5,
      "center": [
        6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    },
    {
      "mass": 0.5,
      "center": [
        -6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    }
  ],
  "n": [
    80,
    80,
    24
  ],
  "conformal_choice": 0,
  "inner_flatten": 0,
  "omega": [
    0.2,
    0.2
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.0,
    0.0
  ],
  "inner_max": [
    0.0,
    0.0
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9802286830137301, 1.3015722129323795e-09, 1.2698678245914683e-09, 9.848836542664666e-10, 2.128508876574204e-08, -8.380387127556561e-07, 0.49490753418965755]`.

## moderate_far0_stable

Preliminary gate: **True**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000325316 | 0.000261737 | 5.16796e-07 | 1.69787e-06 | 0.994498 | 0.00873458 |
| 40×40×20 | 0.000119528 | 7.13055e-05 | 8.48723e-08 | 8.71704e-07 | 0.994489 | 13.6277 |
| 56×56×28 | 1.49731e-05 | 9.12169e-06 | 1.10712e-08 | 6.7693e-08 | 0.994489 | 53.4101 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00136521 | 4.52545e-05 | 0.00134455 | 0.00100432 | 9.58515e-12, 3.10817e-13, 1.46295e-12, 9.74829e-13 |
| 40×40×20 | 0.00137136 | 3.43781e-05 | 0.000564469 | 0.000225496 | 1.6912e-10, 1.42273e-10, 3.22102e-10, 1.65026e-10 |
| 56×56×28 | 0.00137305 | 3.32975e-05 | 4.77915e-05 | 3.99333e-05 | 4.25335e-10, 1.30587e-10, 1.32595e-10, 5.86782e-11 |

Finest source SHA: `79e96b4c151a15bd178c92a7820a518ff71c15e83e226a29ff7f54ac5542b1e6`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    56,
    56,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 2048
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944889341939372, 0.010847878260680024, 0.006521486613271152, 0.019712484180835252, 0.0399897173671845, 0.11614609638778407, 0.3377684348630647]`.

## highspin_stable

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×80×16 | 9.62592e-05 | 1.82966e-06 | 9.69597e-09 | 2.15961e-10 | 0.980229 | 78.4001 |
| 112×112×16 | 7.2732e-06 | 1.3046e-07 | 1.83243e-09 | 1.29093e-12 | 0.980229 | 184.992 |
| 144×144×16 | 1.77148e-06 | 8.72061e-08 | 1.79708e-09 | 1.04007e-13 | 0.980229 | 206.005 |
| 160×160×16 | 1.5509e-06 | 8.68914e-08 | 1.53984e-09 | 9.83664e-14 | 0.980229 | 358.121 |
| 160×160×24 | 3.17793e-07 | 2.69769e-09 | 1.96005e-09 | 2.7991e-15 | 0.980229 | 604.33 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×80×16 | no support | no support | 0.000495154 | 8.49162e-06 | 3.6725e-08, 3.0397e-12, 2.6675e-11, 2.43838e-12 |
| 112×112×16 | no support | no support | 3.41441e-05 | 5.14659e-07 | 1.9885e-09, 3.91394e-13, 3.77577e-12, 1.38799e-13 |
| 144×144×16 | no support | no support | 4.44377e-06 | 2.42484e-07 | 5.99414e-07, 7.00384e-10, 1.29152e-08, 2.16552e-10 |
| 160×160×16 | no support | no support | 3.55904e-06 | 2.43432e-07 | 6.52199e-08, 2.24582e-11, 8.96724e-11, 1.0085e-11 |
| 160×160×24 | no support | no support | 1.79295e-06 | 9.51411e-09 | 2.50755e-08, 5.10608e-10, 8.80717e-10, 1.99224e-10 |

Finest source SHA: `79e96b4c151a15bd178c92a7820a518ff71c15e83e226a29ff7f54ac5542b1e6`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Reference comparison:

```json
{
  "source": "thesis Table3.1 HS99UU",
  "ADM_energy": 0.980124,
  "absolute_tolerance": 0.0002,
  "energy_error": 0.0001046761166793031,
  "energy_agreement": true,
  "exact_historical_reproduction": false,
  "departure": "modern superposed-metric trace projection; horizon mass/spin unmeasured"
}
```

Combined high-regime gate: **False**. Finest local strict gate: **True**. Original three-grid monotonic gate: **False**.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.5,
      "center": [
        6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    },
    {
      "mass": 0.5,
      "center": [
        -6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.2475
      ],
      "velocity": [
        0.0,
        0.0,
        0.0
      ]
    }
  ],
  "n": [
    160,
    160,
    24
  ],
  "conformal_choice": 0,
  "inner_flatten": 0,
  "omega": [
    0.2,
    0.2
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.0,
    0.0
  ],
  "inner_max": [
    0.0,
    0.0
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 6144
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9802286761166793, -2.2185065529123618e-14, -2.072480955364718e-14, -2.87765220875191e-15, -1.2241035550648884e-13, 1.0794557612498164e-13, 0.4949061137437017]`.

## highboost_far0

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 48×48×8 | 0.00492912 | 0.0181225 | 0.00025473 | 0.00038119 | 2.42555 | 104.212 |
| 80×80×8 | 0.0019164 | 0.00476354 | 0.000168863 | 6.24096e-05 | 2.56572 | 84.9899 |
| 112×112×8 | 0.00151628 | 0.00280584 | 3.05733e-05 | 1.76599e-05 | 2.80211 | 175.713 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 48×48×8 | 0.000624856 | 0.000320971 | 0.0129992 | 0.0430799 | 138.821, 48.1677, 0.605647, 0.650535 |
| 80×80×8 | 0.000328593 | 0.000773704 | 0.00551453 | 0.0167648 | 457.209, 140.655, 1.24362, 1.23176 |
| 112×112×8 | 0.000117745 | 0.000861754 | 0.00562504 | 0.0161933 | 745.113, 227.427, 1.90489, 2.41225 |

Finest source SHA: `79e96b4c151a15bd178c92a7820a518ff71c15e83e226a29ff7f54ac5542b1e6`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Reference comparison:

```json
{
  "source": "specified local Gamma=sqrt5 benchmark",
  "exact_historical_reproduction": false,
  "departure": "thesis Table4.3 lacks complete bare inputs and uses historical step stuffing"
}
```

Combined high-regime gate: **False**. Finest local strict gate: **False**. Original three-grid monotonic gate: **True**.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.5,
      "center": [
        6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.0
      ],
      "velocity": [
        -0.8944271909999159,
        0.0,
        0.0
      ]
    },
    {
      "mass": 0.5,
      "center": [
        -6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.0
      ],
      "velocity": [
        0.8944271909999159,
        0.0,
        0.0
      ]
    }
  ],
  "n": [
    112,
    112,
    8
  ],
  "conformal_choice": 0,
  "inner_flatten": 1,
  "omega": [
    1.0,
    1.0
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.022360679774997897,
    0.022360679774997897
  ],
  "inner_max": [
    0.044721359549995794,
    0.044721359549995794
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 2048
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[2.8021084064757087, -0.04432222445056393, 8.627594038492523e-06, 1.3868003162774863e-05, 0.03225727503977989, 0.02248638513892231, -0.00833118693832135]`.

## highboost_actualop_far0

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 48×48×8 | 0.00322904 | 0.0307012 | 9.54712e-05 | 0.000756258 | 3.0058 | 0.0286617 |
| 80×80×8 | 0.000196529 | 0.00578986 | 4.1288e-06 | 8.97753e-05 | 3.05217 | 167.805 |
| 112×112×8 | 1.24316e-05 | 0.00135046 | 1.42448e-07 | 2.1977e-06 | 3.07724 | 793.382 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 48×48×8 | 0.00018458 | 0.00161039 | 0.0125885 | 0.0866926 | 1.56742e-08, 2.17815e-09, 4.98446e-11, 9.11264e-11 |
| 80×80×8 | 6.15661e-05 | 0.000942688 | 0.000628258 | 0.0236149 | 2.86562e-09, 4.14052e-10, 1.03691e-11, 1.46656e-11 |
| 112×112×8 | 5.8495e-05 | 0.000931577 | 6.67551e-05 | 0.00774977 | 7.19221e-09, 8.98319e-10, 1.04376e-11, 1.72301e-11 |

Finest source SHA: `79e96b4c151a15bd178c92a7820a518ff71c15e83e226a29ff7f54ac5542b1e6`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Reference comparison:

```json
{
  "source": "specified local Gamma=sqrt5 benchmark",
  "exact_historical_reproduction": false
}
```

Combined high-regime gate: **False**. Finest local strict gate: **False**. Original three-grid monotonic gate: **True**.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.5,
      "center": [
        6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.0
      ],
      "velocity": [
        -0.8944271909999159,
        0.0,
        0.0
      ]
    },
    {
      "mass": 0.5,
      "center": [
        -6.0,
        0.0,
        0.0
      ],
      "spin": [
        0.0,
        0.0,
        0.0
      ],
      "velocity": [
        0.8944271909999159,
        0.0,
        0.0
      ]
    }
  ],
  "n": [
    112,
    112,
    8
  ],
  "conformal_choice": 0,
  "inner_flatten": 0,
  "omega": [
    1.0,
    1.0
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.022360679774997897,
    0.022360679774997897
  ],
  "inner_max": [
    0.044721359549995794,
    0.044721359549995794
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 80,
  "memory_limit_mib": 2048
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[3.0772391750112225, 1.6150394823583045e-13, -6.313166045346149e-14, 9.929484896346418e-16, -5.134809203076543e-14, -4.1055353561247923e-13, -3.421495653350521e-12]`.

## moderate_sampler_api

Preliminary gate: **True**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000325316 | 0.000261737 | 5.16796e-07 | 1.69787e-06 | 0.994498 | 0.00668525 |
| 40×40×20 | 0.000119528 | 7.13055e-05 | 8.48705e-08 | 8.71704e-07 | 0.994489 | 10.5622 |
| 56×56×28 | 1.49731e-05 | 9.12169e-06 | 1.10704e-08 | 6.7693e-08 | 0.994489 | 41.6379 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00136521 | 4.52545e-05 | 0.00134455 | 0.00100432 | 1.91088e-11, 5.91809e-13, 2.9139e-12, 1.85532e-12 |
| 40×40×20 | 0.00137136 | 3.43781e-05 | 0.000564469 | 0.000225496 | 1.66068e-10, 1.42039e-10, 3.21978e-10, 1.65849e-10 |
| 56×56×28 | 0.00137305 | 3.32975e-05 | 4.77915e-05 | 3.99333e-05 | 2.46066e-10, 1.30612e-10, 1.48581e-10, 5.4861e-11 |

Finest source SHA: `2dbf450d83d87ad58484fe9b676327dedd93ae8d185759281898322c338b49cb`. Unknown basis: `u=W+(A-1)V, W=sum((1-F)*(psi_seed-1))`; maps: `None`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    56,
    56,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 2048
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944889341928632, 0.010847878260679997, 0.006521486613271172, 0.019712484180835262, 0.03998971736718533, 0.1161460963877938, 0.3377684348630591]`.

## moderate_c2prolate

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.00166547 | 0.0134378 | 0.00150536 | 0.0315229 | 1.00118 | 1.27338 |
| 40×40×20 | 0.00033146 | 0.00922767 | 0.0012517 | 0.0236885 | 0.996097 | 12.1832 |
| 56×56×28 | 0.000329162 | 0.00375352 | 0.000169196 | 0.0275579 | 0.995335 | 62.8489 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.0013598 | 6.17331e-05 | 0.0046602 | 0.077093 | 2.67833e-15, 1.12465e-16, 1.10575e-16, 1.26023e-16 |
| 40×40×20 | 0.00137488 | 4.98361e-05 | 0.00391543 | 0.0550346 | 4.72892e-16, 2.39388e-17, 4.02107e-17, 3.09047e-17 |
| 56×56×28 | 0.00138323 | 4.41192e-05 | 0.00114082 | 0.0679002 | 1.55617e-14, 1.97678e-15, 2.75303e-15, 1.25663e-15 |

Finest source SHA: `5ae1a9ba2076bccf1e577f2d996e14b5c4d179ae3a240025d709d51676c48532`. Unknown basis: `modal_P_C2prolate_v1`; maps: `{'radial_stretch': 1.0, 'angular_stretch': 0.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    56,
    56,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 2048
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9953345312809408, 0.010942415167982366, 0.010540660009669195, 0.02017126897190442, 0.03999620402442297, -0.028043109680303338, -0.9996558097927203]`.

## moderate_c2mapped

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000440002 | 0.00571512 | 0.00026823 | 0.0311885 | 1.00872 | 1.0394 |
| 40×40×20 | 0.000144464 | 0.000922282 | 0.000135276 | 0.0138438 | 0.99481 | 12.972 |
| 56×56×28 | 2.79138e-05 | 0.000531766 | 6.34544e-05 | 0.00162502 | 0.994501 | 62.4589 |
| 80×80×28 | 2.13171e-05 | 0.00024374 | 2.86667e-05 | 0.00245264 | 0.994489 | 177.668 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00138206 | 4.69087e-05 | 0.00186814 | 0.07442 | 6.15522e-14, 1.45902e-15, 3.64652e-15, 1.53006e-15 |
| 40×40×20 | 0.00137465 | 3.29042e-05 | 0.000570397 | 0.0322071 | 2.32137e-15, 4.67484e-17, 1.40052e-16, 6.56172e-17 |
| 56×56×28 | 0.00137494 | 3.37587e-05 | 0.00023602 | 0.00550184 | 7.13837e-14, 7.7142e-16, 2.37396e-15, 9.78866e-16 |
| 80×80×28 | 0.00137476 | 3.33916e-05 | 8.60416e-05 | 0.00655594 | 2.06444e-06, 9.5306e-07, 1.4858e-06, 1.17649e-06 |

Finest source SHA: `d2de7b9fe0de67403a88ff520fa164ab4b3076c84d321c7c36180b892991f652`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    80,
    80,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 12,
  "max_krylov": 600,
  "krylov_restart": 40,
  "memory_limit_mib": 4096
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944888429925259, 0.010779493758228884, 0.006497234743744305, 0.02092104200061416, 0.04184480604301399, 0.19708498930016013, -0.07450431448545698]`.

## moderate_c2cache

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000440002 | 0.00571512 | 0.00026823 | 0.0311885 | 1.00872 | 1.12131 |
| 40×40×20 | 0.000144464 | 0.000922282 | 0.000135276 | 0.0138438 | 0.99481 | 14.5322 |
| 56×56×28 | 2.79139e-05 | 0.000531766 | 6.34542e-05 | 0.00162502 | 0.994501 | 49.0713 |
| 80×80×28 | 2.13271e-05 | 0.000248068 | 2.86933e-05 | 0.00243548 | 0.994492 | 296.947 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00138206 | 4.69087e-05 | 0.00186814 | 0.07442 | 3.3366e-14, 1.08075e-15, 2.2203e-15, 1.09477e-15 |
| 40×40×20 | 0.00137465 | 3.29042e-05 | 0.000570397 | 0.0322071 | 1.2748e-15, 4.81644e-17, 1.38954e-16, 6.65755e-17 |
| 56×56×28 | 0.00137494 | 3.37587e-05 | 0.00023602 | 0.00550184 | 5.35313e-14, 9.22861e-16, 3.62005e-15, 1.7773e-15 |
| 80×80×28 | 0.00137474 | 3.33569e-05 | 8.60485e-05 | 0.00649838 | 5.82536e-15, 2.61247e-16, 7.38453e-16, 3.95576e-16 |

Finest source SHA: `722df13ac2a86cc390b0c7a677d2c2fea92034dd48d5e44f5b48100f08dc9fdd`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    80,
    80,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 80,
  "memory_limit_mib": 4096
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944915447706093, 0.010778773750116393, 0.006494616068893856, 0.020890205419637726, 0.04183562172898036, 0.19854095268331612, -0.06484533114659395]`.

## moderate_c2block

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 24×24×12 | 0.000440002 | 0.00571512 | 0.00026823 | 0.0311885 | 1.00872 | 0.35926 |
| 40×40×20 | 0.000144464 | 0.000922282 | 0.000135276 | 0.0138438 | 0.99481 | 3.01237 |
| 56×56×28 | 2.79139e-05 | 0.000531766 | 6.34545e-05 | 0.00162502 | 0.994501 | 8.45302 |
| 80×80×28 | 2.1327e-05 | 0.000248068 | 2.86931e-05 | 0.00243548 | 0.994492 | 31.3308 |
| 104×104×28 | 5.71749e-06 | 9.49288e-05 | 3.97901e-06 | 0.000231209 | 0.994489 | 42.4317 |
| 128×128×28 | 2.49292e-06 | 7.82026e-05 | 5.5858e-07 | 0.00027654 | 0.994489 | 164.951 |
| 152×152×28 | 9.75803e-07 | 1.66413e-05 | 4.80477e-07 | 0.000174423 | 0.994489 | 144.756 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 24×24×12 | 0.00138206 | 4.69087e-05 | 0.00186814 | 0.07442 | 5.52567e-15, 3.41752e-16, 2.62848e-16, 1.64584e-16 |
| 40×40×20 | 0.00137465 | 3.29042e-05 | 0.000570397 | 0.0322071 | 2.28704e-15, 3.59402e-17, 1.90826e-16, 7.55835e-17 |
| 56×56×28 | 0.00137494 | 3.37587e-05 | 0.000236021 | 0.00550184 | 5.7025e-15, 3.25689e-15, 3.61751e-15, 3.01719e-15 |
| 80×80×28 | 0.00137474 | 3.33569e-05 | 8.60489e-05 | 0.00649838 | 5.2939e-15, 2.3466e-16, 1.0267e-15, 4.12834e-16 |
| 104×104×28 | 0.00137389 | 3.2559e-05 | 2.68589e-05 | 0.000612845 | 1.15421e-14, 4.69233e-16, 1.85525e-15, 7.89531e-16 |
| 128×128×28 | 0.00137351 | 3.23803e-05 | 8.6349e-06 | 0.000789495 | 1.63302e-14, 5.61431e-16, 2.32968e-15, 9.51629e-16 |
| 152×152×28 | 0.00137353 | 3.24807e-05 | 4.93467e-06 | 0.000546329 | 2.46035e-14, 8.3646e-16, 2.57992e-15, 1.63916e-15 |

Finest source SHA: `25ca63784d1d18175bc85b06d20f940577bb6106d1d11792307cfaba1f63ce37`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    152,
    152,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944890545366777, 0.010851866554435543, 0.0065280154165237555, 0.019711651455636074, 0.03995082863045246, 0.10500242710164513, 0.3302949814728183]`.

## moderate_c2block_far40

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×80×28 | 2.15666e-05 | 0.000248099 | 2.81869e-05 | 0.00243607 | 0.994486 | 21.7921 |
| 104×104×28 | 5.75448e-06 | 9.49259e-05 | 3.99992e-06 | 0.000230949 | 0.994482 | 41.4298 |
| 128×128×28 | 3.19933e-06 | 7.82167e-05 | 6.39036e-07 | 0.000276117 | 0.994481 | 172.636 |
| 152×152×28 | 1.29475e-06 | 1.66369e-05 | 4.87401e-07 | 0.00017451 | 0.994481 | 150.07 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×80×28 | 0.00137474 | 3.33547e-05 | 8.38858e-05 | 0.00649946 | 9.98834e-15, 8.42064e-16, 1.13405e-15, 7.6075e-16 |
| 104×104×28 | 0.0013739 | 3.25568e-05 | 2.2213e-05 | 0.000611825 | 1.27005e-14, 4.53393e-16, 1.20792e-15, 7.39528e-16 |
| 128×128×28 | 0.00137351 | 3.23781e-05 | 1.05144e-05 | 0.000787817 | 2.25154e-14, 5.48449e-16, 2.01411e-15, 1.33656e-15 |
| 152×152×28 | 0.00137353 | 3.24785e-05 | 5.07191e-06 | 0.000546629 | 2.08281e-14, 8.75442e-16, 3.87562e-15, 1.38354e-15 |

Finest source SHA: `25ca63784d1d18175bc85b06d20f940577bb6106d1d11792307cfaba1f63ce37`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    152,
    152,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 40.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944813246723054, 0.010748718344948987, 0.006447942549544943, 0.019543546104561743, 0.03995046544162055, 0.10510570111893107, 0.3302411627452412]`.

## moderate_c2diff

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 56×56×28 | 2.79139e-05 | 0.000531766 | 6.34542e-05 | 0.00162502 | 0.994501 | 8.77042 |
| 80×80×28 | 2.13271e-05 | 0.000248068 | 2.86934e-05 | 0.00243548 | 0.994492 | 30.8511 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 56×56×28 | 0.00137494 | 3.37587e-05 | 0.00023602 | 0.00550184 | 3.39595e-15, 2.24242e-15, 1.43121e-15, 1.09673e-15 |
| 80×80×28 | 0.00137474 | 3.33569e-05 | 8.60481e-05 | 0.00649838 | 4.25601e-15, 1.65273e-16, 5.5168e-16, 2.88625e-16 |

Finest source SHA: `c9186bc59fcd7abcd93989b14e380bdb69d40bbff2a2af15a1aec65de3fa8682`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    80,
    80,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 6144
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944915447279237, 0.010778773713855296, 0.006494616145471746, 0.02089020511202079, 0.041835621965717766, 0.1985408285045107, -0.06484530608472855]`.

## moderate_focus05_k3

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×80×28 | 8.55088e-06 | 0.000306247 | 3.03523e-06 | 0.00109514 | 0.994491 | 48.6601 |
| 104×104×28 | 2.00706e-06 | 0.00113485 | 7.75997e-06 | 0.0101254 | 0.994579 | 91.5329 |
| 128×128×28 | 1.46493e-06 | 0.000251054 | 1.34171e-06 | 0.00261068 | 0.994497 | 310.585 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×80×28 | 0.00137411 | 3.27075e-05 | 3.56694e-05 | 0.00362698 | 1.74716e-14, 4.93391e-16, 2.62886e-15, 1.29645e-15 |
| 104×104×28 | 0.00137352 | 3.2466e-05 | 2.57851e-05 | 0.0369336 | 2.60055e-14, 1.17599e-15, 3.21552e-15, 1.75889e-15 |
| 128×128×28 | 0.00137352 | 3.24716e-05 | 7.14575e-06 | 0.00783836 | 4.49742e-14, 2.41655e-15, 7.48883e-15, 4.01978e-15 |

Finest source SHA: `e7a88824a3cc8c231ef58af941d86f77350d36895702d3be6ec6251b9c65a65d`. Unknown basis: `modal_P_C2prolate_map_v3_r0.050000000000000003_k3`; maps: `{'radial_stretch': 0.05, 'angular_stretch': 3.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    128,
    128,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944965211921352, 0.010817840899799492, 0.0069019374067587295, 0.021239124165597176, 0.0395351412886595, -0.014020453451996451, -0.18295199847078986]`.

## moderate_focus_wide_actualop

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 56×56×28 | 0.00637759 | 0.0464328 | 0.033965 | 0.226289 | -0.15861 | 346.391 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 56×56×28 | 0.000766679 | 7.81984e-05 | 0.114521 | 0.794529 | 0.0106518, 1.72769e-06, 1.92891e-06, 1.71675e-06 |

Finest source SHA: `e7a88824a3cc8c231ef58af941d86f77350d36895702d3be6ec6251b9c65a65d`. Unknown basis: `modal_P_C2prolate_map_v3_r0.050000000000000003_k3`; maps: `{'radial_stretch': 0.05, 'angular_stretch': 3.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    56,
    56,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 0,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.05510866900951247,
    0.03698441834070127
  ],
  "inner_max": [
    0.2204346760380499,
    0.14793767336280508
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[-0.1586096344149276, 0.008912973674813342, -0.12910534489879305, -0.5008497130509636, 0.03933514743031652, -736.7741458863702, 735.0936298884321]`.

## moderate_focus_wide_flattened

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 56×56×28 | 6.16208e-05 | 0.000143245 | 3.44861e-07 | 0.000187358 | 0.994496 | 9.3356 |
| 80×80×28 | 1.1197e-05 | 6.80993e-05 | 2.02525e-07 | 0.00035751 | 0.994496 | 49.2063 |
| 104×104×28 | 4.04896e-06 | 5.36288e-05 | 1.82206e-07 | 0.000415962 | 0.994496 | 86.9981 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 56×56×28 | 0.00184652 | 0.000156201 | 0.000282712 | 0.000692254 | 1.13316e-14, 6.63939e-15, 3.63453e-15, 2.3972e-15 |
| 80×80×28 | 0.00188453 | 0.000157476 | 3.94511e-05 | 0.00108669 | 1.80291e-14, 4.40882e-16, 1.92495e-15, 8.72756e-16 |
| 104×104×28 | 0.00186747 | 0.00015696 | 1.88596e-05 | 0.00139586 | 2.25732e-14, 9.01454e-16, 3.67942e-15, 1.46771e-15 |

Finest source SHA: `e7a88824a3cc8c231ef58af941d86f77350d36895702d3be6ec6251b9c65a65d`. Unknown basis: `modal_P_C2prolate_map_v3_r0.050000000000000003_k3`; maps: `{'radial_stretch': 0.05, 'angular_stretch': 3.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    104,
    104,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.05510866900951247,
    0.03698441834070127
  ],
  "inner_max": [
    0.2204346760380499,
    0.14793767336280508
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944961612620147, 0.010463592627560886, 0.005574843725763816, 0.01796740982391691, 0.03990755126426828, 0.11768385314076465, 0.3379262585652986]`.

## moderate_default_square_r128

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×80×28 | 2.13271e-05 | 0.000248068 | 2.86936e-05 | 0.00243548 | 0.994492 | 22.7967 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×80×28 | 0.00137474 | 3.33569e-05 | 8.60489e-05 | 0.00649838 | 9.28758e-15, 1.01057e-15, 1.16792e-15, 8.50917e-16 |

Finest source SHA: `9cbf1108d9060095cade976c27394d8953594f45df97503a90700474f4f822fb`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    80,
    80,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 128,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944915446440036, 0.010778773707088884, 0.006494616453951793, 0.02089020547641751, 0.04183562130693569, 0.1985408896812386, -0.06484543226865393]`.

## moderate_default_polar2_r128

Preliminary gate: **False**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×160×28 | 1.03107e-05 | 4.87957e-05 | 1.79952e-06 | 0.000420865 | 0.994489 | 50.9763 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×160×28 | 0.00137398 | 3.28589e-05 | 5.83639e-05 | 0.00106329 | 1.38446e-14, 1.18692e-15, 2.22098e-15, 9.56894e-16 |

Finest source SHA: `9cbf1108d9060095cade976c27394d8953594f45df97503a90700474f4f822fb`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    80,
    160,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 128,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944894209197925, 0.010798754122580733, 0.006207480166342187, 0.018549083359431305, 0.03922615859313749, -0.0947882812340487, 0.4643616649780857]`.

## moderate_polar_sequence

Preliminary gate: **True**. Stronger gate: **False**. Horizon enclosure: **unverified**.

| Grid | Near H RMS | Near M RMS | Bulk H RMS | Bulk M RMS | ADM E | Solve seconds |
|---|---:|---:|---:|---:|---:|---:|
| 80×160×28 | 1.03108e-05 | 4.87957e-05 | 1.79953e-06 | 0.000420865 | 0.994489 | 1.66055 |
| 104×208×28 | 4.22846e-06 | 1.60529e-05 | 4.77867e-07 | 0.000161743 | 0.994489 | 96.3173 |
| 128×256×28 | 1.94698e-06 | 6.01886e-06 | 1.24645e-07 | 8.01299e-05 | 0.994489 | 209.001 |

| Grid | g<1 H RMS | g<1 M RMS | g=1 H max | g=1 M max | Physical-equivalent collocation maxima H,Mx,My,Mz |
|---|---:|---:|---:|---:|---|
| 80×160×28 | 0.00137398 | 3.28589e-05 | 5.83639e-05 | 0.00106329 | 2.1659e-14, 1.55408e-15, 4.90079e-15, 1.64194e-15 |
| 104×208×28 | 0.00137366 | 3.2564e-05 | 2.4126e-05 | 0.000435066 | 2.58653e-14, 1.01306e-15, 5.58753e-15, 1.87035e-15 |
| 128×256×28 | 0.00137354 | 3.24273e-05 | 1.154e-05 | 0.000207595 | 4.68708e-14, 1.37554e-15, 6.66965e-15, 3.12277e-15 |

Finest source SHA: `126300dc1b2f2a9f8268916f2cff8d3623e13e475121ece9c4b33ad524c37aff`. Unknown basis: `modal_P_C2prolate_mapped_v2`; maps: `{'radial_stretch': 0.2, 'angular_stretch': 2.0}`.

Fully specified finest-grid input:

```json
{
  "hole": [
    {
      "mass": 0.6,
      "center": [
        3.0,
        0.0,
        0.0
      ],
      "spin": [
        0.072,
        0.054,
        0.108
      ],
      "velocity": [
        0.03,
        0.06,
        0.01
      ]
    },
    {
      "mass": 0.4,
      "center": [
        -3.0,
        0.0,
        0.0
      ],
      "spin": [
        -0.032,
        0.048,
        0.016
      ],
      "velocity": [
        -0.02,
        -0.07,
        0.025
      ]
    }
  ],
  "n": [
    128,
    256,
    28
  ],
  "conformal_choice": 1,
  "inner_flatten": 1,
  "omega": [
    0.5,
    0.5
  ],
  "attenuation_power": 4,
  "inner_min": [
    0.045,
    0.03
  ],
  "inner_max": [
    0.09,
    0.06
  ],
  "far_radius": 0.0,
  "tolerance": 1e-14,
  "max_newton": 24,
  "max_krylov": 2000,
  "krylov_restart": 32,
  "memory_limit_mib": 8192
}
```

Charge radii: [100.0, 200.0, 400.0]; quadratic inverse-radius extrapolation. Finest [E,P,J]: `[0.9944890684872375, 0.010847096086811357, 0.006521645371321065, 0.019712897840627703, 0.0399894866430717, 0.11610540248197647, 0.3377515433176252]`.

## Solved coordinate covariance

Gate: **True**. Source case: `moderate_polar_sequence`.

| Grid | gamma relative max | K relative max | Maximum EPJ error |
|---|---:|---:|---:|
| 80×160×28 | 9.37226e-12 | 6.60642e-07 | 0.00415565 |
| 104×208×28 | 3.54093e-11 | 1.23363e-06 | 0.000489784 |
| 128×256×28 | 3.92422e-11 | 2.83925e-06 | 0.000198708 |

Translation errors: `{'gamma': 2.2138817684545064e-15, 'Kij': 3.492347974557706e-13, 'conformal_metric': 6.495658047760385e-16, 'Atilde': 3.2007875803165195e-13, 'psi': 4.813468574275274e-16, 'mean_curvature': 1.241521302417777e-13, 'correction': 3.3969251799675166e-13}`. Rotation errors are compared against five times the measured base-grid truncation difference, with a 1e-9 floor. ADM angular momentum is checked about both translated and fixed origins. The 12×24 surface rule has a rotation-dependent angular quadrature error; see the separate refined-charge evidence below.

## High-spin verifier calibration

The original all-norms monotonic gate remains failed. On the identical 18 bulk points, exact vacuum Brill–Lindquist and both chi=.99 Kerr seeds show the same worsening with smaller Cartesian stencils as the solved binary. This supports a numerical-floor interpretation of the bulk plateau; it does not measure the true bulk residual or change the acceptance gate.

| Exact control or solved grid | h=.016 H RMS | h=.008 | h=.004 | h=.002 | h=.001 |
|---|---:|---:|---:|---:|---:|
| exact analytic Brill-Lindquist | 1.93636e-11 | 6.93961e-11 | 2.96904e-10 | 1.28087e-09 | 4.08757e-09 |
| exact Kerr seed 0 | 2.34083e-11 | 6.75828e-11 | 3.13894e-10 | 1.15766e-09 | 5.0898e-09 |
| exact Kerr seed 1 | 2.29203e-11 | 8.62671e-11 | 3.88131e-10 | 1.2979e-09 | 5.16777e-09 |
| [144, 144, 16] | 2.2004e-11 | 1.0117e-10 | 3.86155e-10 | 1.79708e-09 | 5.1742e-09 |
| [160, 160, 16] | 2.29568e-11 | 6.75056e-11 | 3.33175e-10 | 1.53984e-09 | 5.1736e-09 |
| [160, 160, 24] | 2.83712e-11 | 9.27353e-11 | 4.30533e-10 | 1.96005e-09 | 5.96023e-09 |

Full pointwise/stencil evidence: `validation/verifier_floor_highspin.json`. High-spin angular refinement from160²×16 to160²×24 lowers near H RMS from1.55e-6 to3.18e-7 and M RMS from8.69e-8 to2.70e-9. Near stencil checks leave those results stable. Local strict thresholds pass, but the predeclared aggregate gate does not.

## Current integration and conditioning diagnostics

The matched default-map80×80→80×160×28 polar refinement reduces fixed bulk physical momentum RMS5.79x and identical104×208×28 g1 momentum RMS5.51x. Both binaries remain failed. The new80×160/104×208/128×256×28 sequence is declared in polar_sequence_plan.json; its current records are shown above. Physical thresholds are unchanged.

Analytic metric gradients reduce matched charge integration from13.657s to1.120s (12.20x), with energy differences<2.7e-11 and bit-identical P/J on matching sphere nodes. Global-z polar refinement at fixed phi64 gave misleading stability: phi128 shifted angular momentum by.00935. Aligning the sphere polar axis with the prolate axis yields azimuthal change6.19e-14 and polar change6.35e-10 on the same failed polynomial. Matched independent Cartesian-FD fluxes and centered tensor covariance pass1e-8, with discrepancies<5.6e-11 and1.7e-10. These are integration controls, not binary acceptance.

The optional sin6/(1-t)^6 positive row scaling fails its declared exact-seed far-source floor1e-14: the finest boosted/generic controls reach1.89e-13. No binary has been solved in that norm. The failed artifact is retained. Correcting the stale76-point allocation estimate to the current five-point stencil preserves fields/residual/JVP bit-for-bit and keeps the2048MiB default and8192MiB cap.

The three-grid moderate polar sequence subsequently passes the preliminary physical, refined centered charge and solved coordinate-covariance gates. Stronger accuracy and binary horizon/enclosure remain unverified. Centered extraction uses2Npolar/3Npolar and phi64/128; angular/fit changes are<1.8e-13, and the two finest extrapolated EPJ vectors differ2.33e-7. Ring reuse gives another21.12x speedup on matched64×64 extraction, with differences<5.6e-13; all-mode first-gradient, frame, invalidation and independent centered/off-center FD controls pass.

The covariance fixed-global-origin charges still use coarse12×24 integration: finest error1.99e-4 passes the preliminary.005 bound and is not a stronger quadrature claim. Its process peak is10.32GB over sequential contexts. See docs/HISPID_PERFORMANCE.md for measured serial speedups and the isolated integration handoff.

## Additional reproducible checks

- `validation/asymptotic_inverse_default.json`: control flag=True; binary acceptance is not implied.
- `validation/analytic_charge_method_comparison.json`: control flag=True; binary acceptance is not implied.
- `validation/analytic_charges_operator_equivalence.json`: control flag=True; binary acceptance is not implied.
- `validation/aligned_charges_operator_equivalence.json`: control flag=True; binary acceptance is not implied.
- `validation/aligned_charge_fd_comparison.json`: control flag=True; binary acceptance is not implied.
- `validation/budget_operator_equivalence.json`: control flag=True; binary acceptance is not implied.
- `validation/charge_quadrature_polar2_r128.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/charge_quadrature_polar2_r128_fine.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/charge_quadrature_polar2_r128_analytic.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/charge_quadrature_polar2_r128_phi_refinement.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/charge_quadrature_polar2_aligned.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/far_source_floor_infinity_equilibrated.json`: control flag=False; binary acceptance is not implied.
- `validation/polar_sequence_plan.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/charge_rings_operator_equivalence.json`: control flag=True; binary acceptance is not implied.
- `validation/charge_rings_method_comparison.json`: control flag=True; binary acceptance is not implied.
- `validation/charge_rings_independent_fd.json`: control flag=True; binary acceptance is not implied.
- `validation/charge_rings_private_controls.json`: control flag=True; binary acceptance is not implied.
- `validation/charge_rings_full_tests.json`: control flag=True; binary acceptance is not implied.
- `validation/polar_sequence_refined_charges.json`: control flag=diagnostic only; binary acceptance is not implied.
- `validation/highboost_seed_controls.json`: passed=True.
- `validation/failed_highboost_seed_charge_resolution.json`: passed=False.
- `validation/target_seed_controls.json`: passed=True.
- `validation/target_seed_controls_current.json`: passed=True.
- `validation/target_seed_controls_failed_short_radial_fit.json`: passed=False.
- `validation/charge_quadrature_moderate.json`: independent and refined angular/radial charge evidence.
- `validation/charge_quadrature_highspin.json`: independent and refined angular/radial charge evidence.
- `validation/highboost_continued_far0_continuation.json`: v=0: status=0, Newton=0, GMRES=0; v=0.3: status=0, Newton=4, GMRES=226; v=0.5: status=0, Newton=4, GMRES=246; v=0.65: status=0, Newton=4, GMRES=260; v=0.75: status=0, Newton=4, GMRES=301; v=0.82: status=0, Newton=5, GMRES=409; v=0.86: status=0, Newton=5, GMRES=477; v=0.894427: status=1, Newton=24, GMRES=7134.
- `validation/highboost_continued_fine_far0_continuation.json`: v=0.86: status=0, Newton=0, GMRES=0; v=0.87: status=0, Newton=5, GMRES=592; v=0.88: status=0, Newton=7, GMRES=1174; v=0.885: status=0, Newton=8, GMRES=1003; v=0.89: status=1, Newton=24, GMRES=8520.
- `validation/highboost_actualop_far0_continuation.json`: v=0.885: status=0, Newton=7, GMRES=892; v=0.89: status=0, Newton=7, GMRES=1054; v=0.892: status=0, Newton=5, GMRES=1049; v=0.894427: status=0, Newton=7, GMRES=1120.
- `validation/axis_api_migration.json`: off-axis/equation compatibility=True; solved-axis stencil limitations are recorded separately.
- `validation/failed_axis_basis_trials.json`: passed=False.
- `validation/regular_modes_mapped.json`: passed=True.
- `validation/regular_operators_difference.json`: passed=True.
- `validation/mapped_cache_equivalence.json`: off-axis/equation compatibility=True; solved-axis stencil limitations are recorded separately.
- `validation/modal_block_equivalence.json`: off-axis/equation compatibility=True; solved-axis stencil limitations are recorded separately.
- `validation/difference_derivative_equivalence.json`: off-axis/equation compatibility=True; solved-axis stencil limitations are recorded separately.
- `validation/focused_map_operator_audit.json`: passed=False.
- `validation/private_operator_controls_focus05_k3.json`: passed=True.
- `validation/regular_operators_focus05_k3.json`: passed=True.
- `validation/equation_rows_focus05_k3.json`: independent and refined angular/radial charge evidence.
- `validation/dense_equations_focus05_k3.json`: independent and refined angular/radial charge evidence.
- `validation/dense_equations_phi_focus05_k3.json`: independent and refined angular/radial charge evidence.
- `validation/dense_equations_angular_focus05_k3.json`: independent and refined angular/radial charge evidence.
- `validation/dense_equations_radial_focus05_k3.json`: independent and refined angular/radial charge evidence.
- `validation/collocation_map_api_equivalence.json`: off-axis/equation compatibility=True; solved-axis stencil limitations are recorded separately.
- `validation/anisotropic_inverse_default.json`: passed=False.
- `validation/anisotropic_inverse_restart128.json`: passed=True.
- `validation/dense_equations_default_square_r128.json`: independent and refined angular/radial charge evidence.
- `validation/dense_equations_default_polar2_r128.json`: independent and refined angular/radial charge evidence.
- `validation/angular_refinement_comparison.json`: independent and refined angular/radial charge evidence.

## AthenaK integration

The isolated AthenaK branch starts at PR790 head22baa243970fa1880b2bbc48e88a590069d55e47. Its z4c/hispid pgen loads portable physical gamma/K checkpoints and checks ADM/Z4c round trips. AthenaK FastFlow finds exact isolated Schwarzschild, chi=.95 Kerr, v=.885 boosted Schwarzschild and combined chi=.95/v=.885 Kerr horizons at initial time with zero evolution steps.

The combined seed uses flow alpha=.2,lmax48,ntheta50: expansion RMS9.92e-8, relative area error1.75e-14 and independent sampled relative shape error3.35e-7. Fixed-order ntheta74 quadrature confirms area/RMS stability. The failed alpha1 run and coarse strict flags are retained. The three-level outer mesh constraint checks pass at64³, with combined H/M RMS6.13e-7/2.23e-7. These tests validate isolated dataset import and finder controls, not binary attenuation enclosure.

Current executable/source fingerprints and replay records are in the sibling AthenaK worktree docs/hispid-current-controls.json (native9cbf1108…); historical sampler controls and fixed-order quadrature are in docs/hispid-validation.json. Build and replay instructions are in its docs/hispid.md.

## Failures, scope and limitations

- Initial far40 12²×8,20²×12,28²×16 grids failed badly despite tiny weighted residuals. `validation/failed_moderate_coarse.json` preserves the evidence. The Nyquist derivative defect was independently identified and repaired. Analytic scalar far-filter reparameterization removes the known unresolved shell without changing the equations.
- `validation/failed_seed_verifier_step.json` preserves a verifier-step failure on inner-sheet points; refining the independent stencil restored fourth-order convergence.
- The far40 moderate source remains unsupported by the three-grid physical gate (`failed_moderate_filtered_after_split.json`). The original high-spin source suffered singular seed-divergence cancellation; failed H/M trends are retained. The analytic seed momentum identity repaired that source and restored symmetric charges and decreasing momentum norms.
- Direct Gamma=sqrt5 binary solves failed the nonlinear/Krylov and physical gates; `highboost_far0` retains all three resolutions. Continuation records distinguish solver convergence from independent physical validation. No boost interval is inferred from a few successful Newton stages.
- The declared preliminary moderate threshold is RMS 1e-4 with three converging fine grids and <0.5% charge stability. The stronger near/bulk threshold is RMS 1e-6 and max 1e-4. They are reported separately. Neither certifies production or publication accuracy.
- AthenaK horizon controls pass for exact isolated seeds, but solved-binary horizons, mass/spin measurements and attenuation enclosure remain unverified. Isolated seed radii only screen binary window sizes.
- Historical unconstrained nodal-V fields violate Cartesian axis regularity. Current mapped modal-P fields enforce C2 axis limits and pass independent scalar/vector Cartesian manufactured checks. Fresh moderate binaries still fail physical/charge acceptance; regularity, small internal residuals and future enclosure alone do not establish exterior vacuum accuracy.
- Exact punctures remain excluded. Regular signed-lapse graph formulas support the QI throat and pass independent seed/derivative checks there; high-regime binary accuracy and arbitrary-precision arithmetic remain separate questions.
- Thesis historical step stuffing differs from modern smooth Eq.26, and boosted thesis spin conventions differ from this rest-spin API. Published high-boost head-on descriptions omit bare masses and companion-attenuation widths. No original parameter files were recovered. Fully specified local benchmarks must not be called exact table reproduction; milestone E remains incomplete.

See `docs/HISPID.md` for build/API/conventions, `docs/hispid-review.md` for independent source review, and `docs/HISPID_STATUS.md` for milestone status. Reviewable branch work stays isolated; no main-project merge or shared installation change is performed.
