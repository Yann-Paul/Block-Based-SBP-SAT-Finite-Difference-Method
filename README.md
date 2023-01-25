# Block-Based-SBP-SAT-Finite-Difference-Method

The programs shown in this github were used for the master thesis "On Stable Finite-Difference Methods for Convection
Problems on Wavelet Adaptive Block-Based Grids" from Yann-Paul Marsch.

The following modules can be used:

## Calculation of the interpolation Operator:
#### Interpolation_Gleichungssystem_KKT_2te_Ordnung.mlx
Calculates the 2nd order SBP-preserving interpolation operator presented in eq. (2.93) and eq. (2.94)
#### Interpolation_Gleichungssystem_KKT_Matrixform_2te_Ordnung.mlx
Calculates a  2nd order SBP-preserving interpolation operator of adjustable size presented in eq. (2.95)
#### Interpolation_Gleichungssystem_direct_4te_Ordnung_symmetrisch_F2C.mlx
Calculates the 4th order SBP-preserving interpolation operator presented in eq. (B.1)
#### Interpolation_Gleichungssystem_error_calculation_4te_Ordnung.m
Calculates the interpolation error for 4th order SBP-preserving interpolation operator with different stencil size. The results are shown in Table 1 of section 2.5.


## Numerical evaluation of the SBP operator on the block based grid:
#### plot_raw_grids.m
Plots the grid shown in figure 9 and 10 from section 3.1.
#### SBP_SAT_block_based_eigenvalue.m
Calculation of the eigenvalue of the SBP operator with or without the time step method for both grids. The results are shown in figure 11-13 from section 3.1.
#### SBP_SAT_block_based_shell_2d.m
Simulates the advection equation with the initial conditions shown in Figure 14 and plots the H-norm of the numerical solution. The results are shown in figure 15 from section 3.2.
#### convergence_test_block_based_grid.m
Calculates the convergence error shown in Table 2-5 from section 3.2 and plots the convergence error in the 2-norm shown in figure 16 and 17 of secion 3.2.

## Numerical evaluation of the SBP operator on the two blocks of different resolution:
#### convergence_test_two_blocks.m
Calculates the convergence error of the 2-norm for the interface between two blocks of different resolution.
