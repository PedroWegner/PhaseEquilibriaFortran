module flash_module
    use thermo_data
    use root_finding
    use EoS_module
    real(8) :: flash_V
    real(8), allocatable :: flash_K(:), flash_z(:), flash_x(:), flash_y(:)
    contains
    
    subroutine init_flash() ! rever aqui!!
        allocate(flash_K(ncomp), flash_z(ncomp), flash_x(ncomp), flash_y(ncomp))
    end subroutine
    
    subroutine flash_solve(ncomp, T, P, tol, max_iter)
        ! aqui entra o do ao max_iter
        real(8), intent(in) :: T, P, tol
        integer, intent(in) :: max_iter, ncomp
        real(8), allocatable :: K_old(:)
        real(8) :: diff_K
        integer :: i, j
        !!! Estimativas iniciais !!!
        call init_flash()
        flash_V = 0.6
        flash_K = [6.0, 3.0, 0.5, 0.3]
        flash_z = [0.3, 0.1, 0.15, 0.45] ! isso nao eh uma estimativa!!!!!
        
        allocate(K_old(ncomp))
        
        do i = 1, max_iter
            flash_V = newton_raphson(f_rachford_rice, df_rachford_rice, flash_V)
            flash_x = set_x()
            flash_y = set_y()
            K_old = flash_K
            flash_K = set_K(T, P, max_iter, tol)
            diff_K = 0.0
            print*, 'diff: ', diff_K
            print*, 'antigo: ', K_old
            print*, 'calculado: ', flash_K
            do j = 1, ncomp
                diff_K = diff_K + (K_old(j) - flash_K(j))**2
            end do
            if (diff_K < tol) exit
        end do
        
        print*, 'O volume final eh: ', flash_V
        print*, 'a composicao liquida eh: ', flash_x
        print*, 'a composicao gas eh: ', flash_y
        
    end subroutine flash_solve
    
    ! funcoes do modulo
    function set_x() result(x_l)
        real(8), allocatable :: x_l(:)
        integer :: i
        if (.NOT. allocated(x_l)) then
            allocate(x_l(ncomp))
        end if
        do i = 1, ncomp
            x_l(i) = flash_z(i) / (1 + flash_V * (flash_K(i) - 1))
        end do
    end function set_x
    
    function set_y() result(y_l)
        real(8), allocatable :: y_l(:)
        integer :: i
        if (.NOT. allocated(y_l)) then
            allocate(y_l(ncomp))
        end if
        do i = 1, ncomp
            y_l(i) = flash_K(i) * flash_x(i)
        end do
    end function set_y
    
    function set_K(T, P, max_iter, tol) result(k_l)
        real(8), intent(in) :: T, P, tol
        integer, intent(in) :: max_iter
        real(8), allocatable :: k_l(:), phi_l(:), phi_g(:)
        integer :: i
        if (.NOT. allocated(k_l) .AND. .NOT. allocated(phi_l) .AND. .NOT. allocated(phi_g)) then
            allocate(k_l(ncomp), phi_l(ncomp), phi_g(ncomp))
        end if
        
        ! Resolve as fases
        x = flash_x
        call solve_Z(T, P, max_iter, tol, .FALSE.)
        call fugacity_solve()
        phi_l = phi
        x = flash_y
        call solve_Z(T, P, max_iter, tol, .TRUE.)
        call fugacity_solve()
        phi_g = phi
        
        do i = 1, ncomp
            k_l(i) = phi_l(i) / phi_g(i)
        end do
    end function set_K
    
    ! Equacoes de Rachford-Rice
    function f_rachford_rice(v) result(f_rr)
        real(8) :: f_rr, f_aux
        integer :: i
        f_aux = 0.0
        do i = 1, ncomp
            f_aux = f_aux + ((flash_z(i) * (flash_K(i) - 1)) / (1 + v * (flash_K(i) - 1)))
        end do
        f_rr = f_aux !! aqui vai dar um problema!!!! provavelmente terei que criar K, z e ncomp LOCAIS no modulo
    end function f_rachford_rice
    
    function df_rachford_rice(v) result(df_rr)
        real(8) :: df_rr, f_aux
        integer :: i
        f_aux = 0.0
        do i = 1, ncomp
            f_aux = f_aux + ((flash_z(i) * (flash_K(i) - 1)**2) / (1 + v * (flash_K(i) - 1))**2)
        end do
        df_rr = - f_aux !! aqui vai dar um problema!!!
    end function df_rachford_rice
    
end module flash_module