# SABEMMT: Structures And Blade Element Modified Momentum Theory

**SABEMMT** is a MATLAB-based aerodynamic and structural solver designed for the analysis and optimization of airplane propellers, helicopter rotors, and multirotor blades in axial flight.

Unlike standard BEMT codes that fail at low speeds, this project implements an "Engineering Modification" to the Momentum Theory, allowing it to robustly handle **static thrust**, **hover**, and **turbulent wake states**.

## Features

* **Robust Inflow Model**: Uses the Cuerva et al. (2006) modification to handle the "Singularity of Momentum Theory," making it stable for static thrust and vertical descent.
* **Reynolds-Aware Aerodynamics**: interpolation of $C_l$ and $C_d$ polars across varying Reynolds numbers ($Re$) along the blade span.
* **Structural Estimation**: Calculates Centrifugal Force ($F_{cf}$) and Bending Moments ($M_x, M_y$) to estimate stresses in a conservative way.
* **Vectorized Solver**: Optimized MATLAB code for fast execution, suitable for iterative design and optimization loops.
* **Visual Output**: Automatically generates plots for aerodynamic load distribution ($dT/dx, dF_t/dx$), angle of attack ($\alpha$), and structural stress.

## Project Structure

* `main_runner.m`: **Start here.** The driver script that defines the propeller geometry (Chord/Twist), loads airfoil data, and executes the simulation.
* `runSABEMMT.m`: The core BEMT solver function. It handles the iterative induction loops and structural integration.
* `getAirfoilDataE214_MultiRe_cleaned.m`: A sample database containing Eppler E214C-PT airfoil polars (Lift/Drag vs Alpha vs Re).

## Theoretical Basis

This solver is based on **Blade Element Momentum Theory (BEMT)** but replaces the standard momentum equation with a continuous function that bridges the gap between the "Windmill Brake State" and "Propeller State."
For a complete mathematical derivation of the model, including the structural equations, please see the [**Technical Documentation (PDF)**](docs/The_SABEMMT_Aerodynamic_and_Structural_Model_for_Propellers.pdf).

**Reference:**
> Cuerva, A., Sanz-Andrés, A., Meseguer, J., & Espino, J. L. (2006). "An Engineering Modification of the Blade Element Momentum Equation for Vertical Descent". *Journal of the American Helicopter Society*, 51(4), 349-354.

## Usage

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/sabemmt-rotor-analysis.git](https://github.com/miguel-frade/sabemmt.git)
    cd sabemmt-rotor-analysis
    ```

2.  **Open MATLAB** and navigate to the project folder.

3.  **Run the driver script:**
    ```matlab
    main_runner
    ```

4.  **Analyze the Output:**
    The script will output performance metrics to the Command Window and generate a 4-panel figure:
    * **Thrust & Torque Distribution**: Aerodynamic loading along the span.
    * **Angle of Attack**: Checks for stall relative to the airfoil limits.
    * **Structural Stress**: Estimated max stress (useful for material selection).
    * **Geometry**: Visualization of the blade's chord and twist.

## Configuration

To analyze your own propeller, modify `main_runner.m`:

* **Geometry**: Update `R` (radius), `c` (chord distribution), and `theta_deg` (twist distribution).
* **Operating Point**: Change `v_inf` (velocity in m/s) and `rpm`.
* **Airfoil**: Replace `getAirfoilDataE214...` with your own function or data structure. *Ensure your polars cover the necessary Reynolds range.*

## Example Output

For a propeller similar to APC 11x5.5 Sport Propeller used for UAV:

Running BEMT Solver for V=15.0 m/s, RPM=6000...

--- RESULTS ---
Thrust (T):       1.8284 N
Torque (Q):       0.0736 Nm
Power (P):        46.2406 W
Efficiency (eta): 59.31 %
Max Stress:       25.46 MPa
Tip Mach:         0.262

The solver automatically generates performance plots. Below is the output for the APC 11x5.5 propeller at 6000 RPM and 15 m/s:

![SABEMMT Analysis Results](results/SABEMMT_simulation_plot.png)

*Figure: (Top-Left) Aerodynamic load distribution showing smooth gradients due to linear Re-interpolation. (Top-Right) Local angle of attack relative to stall limits. (Bottom-Left) Von Mises stress estimation. (Bottom-Right) Blade geometry.*
