! some precompiler statements.

! time integration
! 1 Runge-Kutta second order (Heun)
! 2 SSP-Runge-Kutta third order and 4 stages
! 3 Runge-Kutta fourth order
#define TIME_INT 3

! time dependent boundary condition only for time integration 3!
! 1 convential way: inflow condition every substep
! 2 no intermediate boundary condition
! 3 intermediate boundary condition according to carpenter, gottlieb and abarbanel
#define INFLOW_BC 3
