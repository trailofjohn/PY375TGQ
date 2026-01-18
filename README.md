# Spacetime State Estimator

A relativistic navigation simulator demonstrating **Extended Kalman Filter (EKF)** techniques for deep-space missions. This project visualizes how spacecraft determine their position using radio signals, accounting for the effects of **General Relativity**.

![Martian Standard Dashboard](results/Screenshot%202026-01-18%20160803.png)

## Overview

When a spacecraft ventures into deep space, simple GPS doesn't work. Instead, it must:
1.  Send and receive radio signals from Earth.
2.  Account for the fact that light takes *minutes to hours* to travel.
3.  Correct for Einstein's effects: **gravitational time dilation** and **Shapiro delay**.

This simulator demonstrates these concepts by running three mission scenarios and showing, in real-time, how the EKF estimate (cyan) converges to the true trajectory (green).

---

## Key Features

| Feature | Description |
|:---|:---|
| **Proper Time Integration** | Tracks relativistic time dilation ($d\tau/dt$) as spacecraft velocity and gravitational potential change. |
| **Doppler Observables** | Measures range-rate (velocity along line of sight) for instantaneous velocity estimation, not just position. |
| **Iterative Light-Time Solver** | Corrects for the fact that the probe's signal was emitted minutes ago from a different location. |
| **Shapiro Delay** | Models the extra path length light travels when passing near the Sun's gravitational field. |
| **10-State EKF** | Estimates `[x, y, z, vx, vy, vz, clock_bias, clock_drift, drift_rate, proper_time]`. |

---

## Scenarios

The simulator runs three pre-configured mission scenarios:

### Scenario A: Martian Standard
A classic Hohmann transfer to Mars. Demonstrates nominal EKF convergence over a ~250-day mission.

![Martian Standard](results/Screenshot%202026-01-18%20160803.png)

### Scenario B: Jovian Deep-Dive
An 800-day mission to Jupiter. The key observation is the **Solar Conjunction**—when the probe passes behind the Sun relative to Earth. Watch the **Shapiro Delay** graph spike dramatically.

![Jovian Deep-Dive](results/Screenshot%202026-01-18%20160821.png)

### Scenario C: Solar Grazer
A probe on a highly elliptical orbit that dives to 0.05 AU from the Sun. This scenario demonstrates extreme **relativistic time dilation** and **velocity changes** (over 100 km/s at perihelion).

![Solar Grazer](results/Screenshot%202026-01-18%20160834.png)

---

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/trailofjohn/PY375TGQ.git
    cd PY375TGQ/code
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the simulation:**
    ```bash
    python main.py
    ```
    This launches all three scenarios in parallel, each in its own interactive window.

---

## Usage

```bash
# Run all scenarios (default)
python main.py

# Run a single scenario
python main.py A   # Martian Standard
python main.py B   # Jovian Deep-Dive
python main.py C   # Solar Grazer
```

### Dashboard Controls
- **Speed Slider**: Adjust animation speed from 0.1x to 10x.
- **Pause Button**: Pause/resume the animation.

### Dashboard Interpretation
- **3D Trajectory**: Green = True, Cyan = EKF Estimate, Magenta = Dead Reckoning (no filter).
- **Position Error**: The primary metric. Watch the EKF error drop from ~1000 km to under 1 km.
- **Clock Bias Residual**: Shows the filter estimating the spacecraft's clock offset in meters.
- **Shapiro Delay**: Shows the relativistic "extra path length" of the signal in meters.

---

## Project Structure

```
PY375TGQ/
├── code/
│   ├── main.py                 # Entry point, multiprocessing orchestrator
│   ├── simulation_core.py      # SpaceTimeState, DynamicsModel, MeasurementModel, EKF
│   ├── simulation_scenarios.py # Scenario configurations (A, B, C)
│   ├── visualization.py        # Matplotlib animation dashboard
│   └── requirements.txt        # numpy, matplotlib
├── text/
│   ├── feynman_physics.md      # Intuitive physics explanation
│   ├── textbook.md             # Technical manual
│   ├── walkthrough.md          # User guide
│   ├── code_logic.md           # Code architecture
│   └── IMPROVEMENTS.md         # V2 roadmap
├── results/                    # Dashboard screenshots
└── README.md
```

---

## The Physics

### Relativistic Time Dilation
A clock on a moving spacecraft in a gravitational field ticks slower than a clock on Earth:

$$\frac{d\tau}{dt} = 1 - \frac{GM}{c^2 r} - \frac{v^2}{2c^2}$$

### Shapiro Delay
Light traveling near a massive object takes an extra time to arrive:

$$\Delta t_{Shapiro} = \frac{4GM}{c^3} \ln\left(\frac{r_e + r_p + \rho}{r_e + r_p - \rho}\right)$$

### Light-Time Iteration
The probe's position when the signal was *emitted* is not where it is *now*. We solve:

$$\vec{r}_{emit} = \vec{r}(t_{recv} - \Delta t), \quad \Delta t = \frac{|\vec{r}_{emit} - \vec{r}_{recv}|}{c}$$

---

## Requirements

- Python 3.10+
- NumPy
- Matplotlib

---

## License

MIT License

---

## Acknowledgments

Inspired by the navigational techniques of NASA's Jet Propulsion Laboratory (JPL) for missions like Voyager, Cassini, and the Mars rovers.
