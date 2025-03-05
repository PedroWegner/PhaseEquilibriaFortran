module thermo_data
    implicit none
    ! Variáveis do modulo
    integer :: ncomp
    real(8), allocatable :: Tc(:), Pc(:), omega(:), x(:)
    
contains
    subroutine armazenar_dados()
        implicit none
        allocate(Tc(ncomp), Pc(ncomp), omega(ncomp), x(ncomp))
    end subroutine armazenar_dados
    
    subroutine liberar_mem()
        implicit none
        deallocate(Tc, Pc, omega, x)
    end subroutine liberar_mem
    
end module thermo_data