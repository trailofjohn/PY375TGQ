# V1 Evaluator Questions (Advanced)

## 1. Physics & Relativity
**Q: How much does the simulation improve by adding Proper Time?**
**A:** In the Solar Grazer scenario (Scenario C), the probe moves at 100 km/s. The time dilation is $\gamma \approx 1 + 10^{-7}$. Over a 10-day pass, this accumulates to roughly 0.1 seconds of clock error. Multiplied by $c$, that is **30,000 kilometers** of position error if ignored. V1 captures this.

**Q: Why define the state vector as 10 elements?**
**A:** We needed to store $\tau$ (Proper Time) as a state because it is the result of an integration ($\int \sqrt{g_{\mu\nu} dx^\mu dx^\nu}$). It cannot be calculated analytically from just the instantaneous position/velocity without history.

## 2. Algorithms
**Q: Explain the "Stiffness" issue you encountered.**
**A:** "Stiffness" in numerical ODEs or filters generally refers to systems with vastly different time scales. Here, we had "Orbital Scale" (Year) and "Light Scale" ($\Delta t \times c$). The Jacobian terms for Clock Drift ($\approx c$) are $10^8$ times larger than Position terms. This caused the EKF update step ($K = P H^T S^{-1}$) to be numerically unstable (singular matrix inversion) when $R$ was too small.

**Q: Why is the Light-Time solver iterative?**
**A:** Because the emission time depends on the distance, but the distance depends on the emission time (where the probe was). It's a transcendental equation $t = f(r(t))$. Iteration converges because $v/c \ll 1$ (Contraction Mapping Theorem applies).

## 3. Engineering
**Q: Why did you separate the Config from the Core?**
**A:** Modularity. `simulation_scenarios.py` acts as a "Mission Control" database. We can add a "Neptune Orbiter" (Scenario D) without touching a line of the Physics engine (`simulation_core.py`).
