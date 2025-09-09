## Fenton-Karma Finitewave model

The Fenton-Karma model is a simplified mathematical representation of cardiac action potentials, designed to reproduce the essential excitation–recovery dynamics of cardiac cells. Unlike detailed ionic models, it does not explicitly describe individual ion currents but instead uses a reduced set of variables and parameters to capture the key features of cardiac electrophysiology. This abstraction makes the model computationally efficient while retaining the ability to simulate wave propagation, spiral wave dynamics, and arrhythmia mechanisms in cardiac tissue.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the model dynamics without the need for the entire framework.

### Reference
Fenton, F. H., & Karma, A. (1998). Vortex dynamics in three-dimensional 
continuous myocardium with fiber rotation: Filament instability and fibrillation.
Chaos: An Interdisciplinary Journal of Nonlinear Science, 8(1), 20-47.

DOI: https://doi.org/10.1063/1.166311

### How to use (quickstart)
```bash
python -m example.fenton_karma_example
```

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── fenton_karma/                    # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       # fill with the model equations (pure functions)
├── implementation/                  # 0D Fenton-Karma implementation
│   ├── __init__.py
│   └── fenton_karma_0d.py
├── example/
│   └── fenton_karma_example.py      # minimal script to run a short trace
├── tests/
│   └── fenton_karma_test.py         # smoke test; extend with reproducibility checks
├── .gitignore
├── LICENSE                          # MIT
├── pyproject.toml                   # placeholders to replace
└── README.md                        # this file
```

### Variables
Model state variables: description, units and ranges (optional)
- `u` — Transmembrane potential (mV)
- `v` - Initial recovery variable (dimensionless)
- `w` - Initial activation variable (dimensionless)

### Parameters
Parameters and their default values
- `tau_r   = 130` - Recovery time constant (ms)
- `tau_o   = 12.5` - Activation time constant (ms)
- `tau_d   = 0.172` - Deactivation time constant (ms)
- `tau_si  = 127` - Slow inactivation time constant (ms)
- `tau_v_m = 18.2` - Membrane time constant (ms)
- `tau_v_p = 10` - Potential time constant (ms)
- `tau_w_m = 80` - Activation time constant for `w` (ms)
- `tau_w_p = 1020` - Activation time constant for `w` (ms)
- `k       = 10` - Scaling factor for recovery dynamics
- `u_c     = 0.13` - Threshold for recovery variable (dimensionless)
- `uc_si   = 0.85` - Slow inactivation threshold (dimensionless)


