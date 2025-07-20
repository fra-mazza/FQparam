# FQparam: Fluctuating Charge Force Field Parametrization

## Goal

This program aims to optimize the parameters for a Fluctuating Charge (FQ) polarizable force field. The optimization minimizes the difference between interaction energies and multipoles (dipole and quadrupole) of a QM molecule in the presence of a fixed charge, and the same molecule described with FQ charges. The fixed charges used for parametrization are typically chosen from PCM tesserae, but the input can be modified.

## Structure

The project is organized into several Fortran modules, promoting modularity, reusability, and maintainability.

*   **`types_module.f90`**: Defines custom derived types for `FQMol_t` (Fluctuating Charge Molecule), `qm_datapoint_t` (QM reference data), `molecule_t` (molecule definition), `optimization_params_t` (parameters to be optimized), and `optimization_settings_t` (optimization control parameters).
*   **`misc.f90`**: Contains utility functions, including matrix inversion (`inv`) and a basic Gradient Descent (`GD`) optimizer.
*   **`fq_util.f90`**: The core physics engine for FQ calculations. It includes subroutines for:
    *   `InitializeMolecules`: Initializes the FQ molecule structure.
    *   `MakeFQJ`: Constructs and inverts the FQ Jacobian matrix.
    *   `UpdateFQs`: Solves for the fluctuating charges.
    *   `IntEnergy`: Calculates the interaction energy.
    *   `GetPotentialPointCharge`: Calculates the potential from external point charges.
    *   `GetFQsDerivatives`: Computes derivatives of FQ charges with respect to parameters.
    *   `Dipole`, `Quadrupole`: Calculate molecular multipoles from FQ charges.
    *   `UpdateGradients`: Accumulates gradients for the cost function.
*   **`input_module.f90`**: Handles all data input, including molecular geometry, initial parameters, and QM reference data from `FQpar_data.dat`.
*   **`cost_function_module.f90`**: Calculates the total cost function and its analytical gradient based on the difference between FQ and QM properties. This module iterates through each QM data point.
*   **`optimizer_module.f90`**: Manages the overall optimization process, including the main optimization loop and calling the gradient descent algorithm.
*   **`debug_module.f90`**: Contains debugging utilities, specifically `check_gradients` for comparing analytical and numerical gradients.
*   **`main.f90`**: The main program that orchestrates the workflow: reads input, performs gradient checks, runs the optimization, and prints final results.

## Building the Project

To build the executable, navigate to the `src` directory and use the `make` command:

```bash
cd src
make
```

This will compile all Fortran source files and link them into an executable named `optimizer`.

## Running the Program

After successful compilation, you can run the program from the `src` directory:

```bash
./optimizer
```

### Input Data (`FQpar_data.dat`)

The program expects an input file named `FQpar_data.dat` in the `src` directory. This file should contain the QM reference data, including:

*   Number of QM calculations.
*   For each QM calculation:
    *   Coordinates of the fixed external charge.
    *   Intensity of the fixed external charge.
    *   QM interaction energy.
    *   QM dipole moment (3 components).
    *   QM quadrupole moment (6 components).

**Example `FQpar_data.dat` format (conceptual):**

```
1
0.0 0.0 5.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
```
*(Note: The actual format will depend on how `input_module.f90` reads the data. The current implementation reads coordinates, charge, then energy, dipole components, and quadrupole components sequentially.)*

## Gradient Checking

The program includes a built-in gradient check that runs automatically before the main optimization. This feature compares the analytically calculated gradients (from `fq_util.f90` and `cost_function_module.f90`) with numerically computed gradients (using finite differences). This is a crucial debugging tool to ensure the correctness of the gradient implementation.

The output of the gradient check will be printed to the console, showing the analytical, numerical, absolute difference, and relative difference for each parameter.

## Optimization Output

During the optimization, the program will print the iteration number and the current cost function value. Upon convergence or reaching the maximum number of iterations, it will print the final optimized `Chi` and `Eta` parameters.

## Future Improvements

*   **More robust input parsing**: Implement a more flexible input system (e.g., using namelists or a dedicated input file parser) instead of hardcoded values and sequential reads.
*   **Advanced Optimizer**: Integrate more sophisticated optimization algorithms (e.g., L-BFGS, Conjugate Gradient) from external libraries like LAPACK/BLAS or custom implementations for faster and more stable convergence.
*   **Error Handling**: Add more comprehensive error handling and validation for input data and numerical stability issues.
*   **Unit Testing**: Implement dedicated unit tests for each subroutine to ensure correctness and prevent regressions.
*   **Documentation**: Expand inline code comments and provide more detailed mathematical derivations for the FQ model and gradient calculations.