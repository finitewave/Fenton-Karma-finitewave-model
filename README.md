# Fenton–Karma with ChR2(H134R)

The ionic core uses the porcine ventricular three-variable minimal model in
[Kaboudian et al. (2019)](https://doi.org/10.1126/sciadv.aav6019), pp. 3–4,
with the four-state ChR2 model from
[Williams et al. (2013)](https://doi.org/10.1371/journal.pcbi.1003220),
Eqs. 1–14, Table 1 and Fig. 5F. This is no longer the original 1998 MLR-I
parameter set. The package works as a Finitewave model plugin or standalone.

## Optical stimulation

```python
from implementation.fenton_karma_0d import FentonKarma0D, Stimulation

model = FentonKarma0D(
    dt=0.01,                         # ms
    stimulations=[],                # electrical current pulses
    optical_stimulations=[
        Stimulation(t_start=10.0, duration=20.0, amplitude=1.0),
    ],                              # optical amplitude: mW/mm²
)
model.parameters['temperature'] = 37.0
model.run(t_max=600.0)
voltage = model.history['u']
light = model.light_history
```

Optical pulses add to `parameters['irradiance']`, which defaults to zero.
For Finitewave, `FentonKarma` automatically discovers the new variables and
parameters. Set `model.irradiance` to a scalar or spatial array, updating it
between time steps to apply a light protocol. It is an irradiance parameter,
not an injected electrical current. Negative values are treated as darkness.

## State and units

- `u`, `v`, `w`: original dimensionless voltage and recovery gates.
- `c1`, `o1`, `o2`, `c2`: ChR2 closed/open probabilities; initially `(1,0,0,0)`.
- `p_chr2`: light-activation state, initially zero; distinct from the FK switch `p`.

Time is in ms. ChR2 rates require voltage in **mV**, not normalized `u`:

```text
V = V_rest + V_scale*u
 du/dt = FK_rhs - I_chr2/(C_m*V_scale)
```

The defaults `V_rest=-85 mV`, `V_scale=100 mV`, and `C_m=1 µF/cm²` are explicit
coupling assumptions. Neither paper validates this combined model or supplies
a voltage calibration for this coupling. Adjust them for the intended tissue;
`V_scale` and `C_m` must be positive. `I_chr2` is in µA/cm².

`g_chr2=0.4 mS/cm²` and `gamma_chr2=0.1` follow Table 1. Set `g_chr2=0` to disable
its voltage contribution. With dark initial states and zero irradiance, the
original FK voltage and recovery dynamics are preserved.

The photon and activation constants are exposed through `get_parameters()`:

| Parameter | Default | Units |
|---|---:|---|
| `sigma_ret` | `12e-20` | m² |
| `lambda_nm` | `470` | nm |
| `w_loss` | `0.77` | dimensionless, denominator convention |
| `hc` | `1.986446e-25` | J·m |
| `eps1` | `0.8535` | dimensionless |
| `eps2` | `0.14` | dimensionless |
| `tau_chr2` | `1.3` | ms |

Separate helpers evaluate the Table 1 baseline rates. `ionic_step` applies
Fig. 5F temperature factors and converts irradiance to `theta=100*I`.
`calc_flux` converts `lambda_nm` to metres explicitly; its output is ms⁻¹.
`calc_G_times_V` returns the product `G(V)*V`, not `G(V)` itself, and
`calc_J_chr2` consumes this product without multiplying by voltage again.

`temperature=22°C` reproduces the baseline rates. Setting it to 37°C applies
Fig. 5F's optimized Q10 factors to the ChR2 rates only; the FK kinetics are not
temperature-scaled. The e12/e21 Q10 factors apply to their dark components.
The light activation function uses Table 1's `theta=100*irradiance` and
`tau_chr2=1.3 ms`; photon flux uses SI constants and converts to ms⁻¹.

The paper fixes `E_chr2=0 mV`. We evaluate `G(V)*V` after algebraic cancellation
to remain finite at `V=0`. The empirical fit has a nonzero limiting current
there and crosses zero near +13.65 mV; this implementation preserves that fit.

All states advance by forward Euler, consistent with the FK implementation.
No probability clipping or renormalization is applied. Use a sufficiently small
time step (the tests use 0.01 ms) and check convergence for stronger illumination
or altered parameters. `ionic_step` now returns eight values:
`(rhs, v_new, w_new, c1_new, o1_new, o2_new, c2_new, p_chr2_new)`.

## Validation

```bash
python -m pytest -q
```

Tests cover published rates, temperature scaling, current units, the zero-voltage
limit, dark FK equivalence, probability conservation, optical excitation and
recovery at 22°C/37°C, and time-step convergence. These are equation and numerical
checks, not experimental validation of the combined model.
