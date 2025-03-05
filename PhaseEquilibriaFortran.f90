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

    program PhaseEquilibriaFortran

    use PC_saft_module
    implicit none
    integer ncomp
    real(8), allocatable :: mole_fraction(:), numb_seg(:), seg_dia(:), depth_pair(:)
    
    !!! PARA USAR O PC_SAFT
    call universal_model_constants()
    
    ncomp = 2
    allocate(mole_fraction(ncomp), numb_seg(ncomp), seg_dia(ncomp), depth_pair(ncomp))
    mole_fraction = [0.4, 0.6]
    numb_seg = [1.2053, 1.0]
    seg_dia = [3.313, 3.7039]
    depth_pair = [90.96, 150.03]
    saft_T = 200.0
    call initialise_mixture(ncomp, mole_fraction, numb_seg, seg_dia, depth_pair)


    end program PhaseEquilibriaFortran

