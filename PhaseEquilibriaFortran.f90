!  PhaseEquilibriaFortran.f90 
!
!  FUNCTIONS:
!  PhaseEquilibriaFortran - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: PhaseEquilibriaFortran
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    PROGRAM PhaseEquilibriaFortran
    USE PC_saft_module
    USE root_finding
    IMPLICIT NONE
    !!! para aplicacao do simplex _teste
    REAL(8) :: eta
    !!! simplex
    REAL(8) :: x_r
    !!! PARA USAR O PC_SAFT
    call universal_model_constants()
    
    ncomp_global = 2
    allocate(x_global(ncomp_global), m_global(ncomp_global), sigma_global(ncomp_global), epsilo_global(ncomp_global))
    x_global = [0.4, 0.6]
    m_global = [1.2053, 1.0]
    sigma_global = [3.313, 3.7039]
    epsilo_global = [90.96, 150.03]
    T_global = 200.0d0
    P_global = 30.0d5
    
    x_r = simplex(residuo_pressao, 1.0d-5)
    print*, x_r

    END PROGRAM PhaseEquilibriaFortran

