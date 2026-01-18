# Spacetime Estimator V1: Implementation Plan

## Goal Description
Upgrade the Spacetime Estimator to a "Flight-Ready" Relativistic Navigation System.
**Key Upgrades**: Coupling of General Relativistic Physics (Proper Time, Light-Time) with High-Precision Estimation (Doppler Fusion).

## User Review Required
> [!IMPORTANT]
> **Measurement Model Changes**: The system now strictly separates "Emission Time" (Probe frame) from "Reception Time" (Earth frame). This adds computational cost (iteration) but ensures causality.

## Proposed Changes
### 1. Visualization
*   **Log-Scale Error Plot**: To visualize exponential decay of error.
*   **Time Slider**: Interactive scrubbing of the mission timeline.

### 2. Physics Core (`v1/simulation_core.py`)
*   **[NEW] `SpaceTimeState.tau`**: 10th state element for Proper Time.
*   **[MODIFY] `DynamicsModel`**: Add $d\tau/dt$ integration.
*   **[MODIFY] `MeasurementModel`**: Add `solve_light_time` and Doppler partial derivatives.

### 3. Estimator Logic
*   **[MODIFY] `ExtendedKalmanFilter`**: Retuned $Q$ and $R$ matrices to handle the condition number of $c$ (Speed of Light).

## Verification
*   **Scenario C (Solar Grazer)**: Verify $\tau$ accumulation matches GR predictions.
*   **Scenario B (Jupiter)**: Verify Light-Time solver corrects the ~3000km error introduced by neglecting $\sim 40$ min transit time.
