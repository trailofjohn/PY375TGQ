# V1 Technical Manual: Relativistic Navigation
**Advanced Implementation**

## 1. Relativistic Dynamics
V1 integrates **General Relativity** into the state dynamics.
$$ \dot{X} = f_{Newton}(X) + f_{GR}(X) $$
Specifically, we track **Proper Time** ($\tau$).
$$ \frac{d\tau}{dt} = 1 - \frac{\mu}{c^2 r} - \frac{v^2}{2 c^2} $$
This differential equation is coupled to the orbital equations. The "Clock Bias" is no longer just a random variable; it includes a deterministic relativistic component.

## 2. Doppler Observables
We introduced a second measurement equation:
$$ z_2 = \dot{\rho} + c \cdot \dot{b} $$
*   $\dot{\rho}$: Range Rate (Geometric).
*   $\dot{b}$: Clock Drift.
This makes velocity **directly observable**.
In V0 (Range only), velocity was unobservable directly; it had to be inferred from the *change* in range over time. This caused "lag" in the velocity estimate. In V1, the velocity estimate is instantaneous.

## 3. Light-Time Solution
We solved the **Retarded Potential** problem.
$$ \vec{r}_{emit} = \vec{r}(t_{recv} - \delta t) $$
$$ \delta t = \frac{|\vec{r}_{emit} - \vec{r}_{recv}|}{c} $$
This is solved via fixed-point iteration in `MeasurementModel.solve_light_time`.
Without this, errors of thousands of kilometers (Earth moves 30km/s * 2400s = 72,000 km in 40 mins) would be introduced in the Jovian scenario.

## 4. Stability Tuning
The inclusion of Relativistic terms ($c \approx 3\times 10^8$) creates a "stiff" system.
*   **Problem**: Doppler noise of $10^{-3}$ m/s vs $c$.
*   **Solution**: We relaxed Doppler noise variance to $0.5$ m/s in the filter tuning ($R$ matrix). This acts as a "damper" to prevent numerical oscillation while still providing high-precision velocity data.
