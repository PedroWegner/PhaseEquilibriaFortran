module PhaseEquilibria
    use thermo_data
    use root_finding
    use EoS_module
    implicit none
    
    contains
    ! Abaixo eu implemento os modelos
    subroutine bolha_p(ncomp, h, T, tol, max_iter)
        implicit none
        integer, intent(in) :: ncomp, h, max_iter
        real(8), intent(in) :: tol
        real(8), allocatable :: x_space(:), y_space(:), P_space(:), x_l(:,:), y_new(:), phi_liq(:), phi_gas(:)
        real(8) :: T, P_old, P_new, K, cond, y_i ! esse K esta LOCAL!!!! nao deve dar problema no flash.
        real(8), allocatable :: y_old(:)
        integer :: total, i, j, n
        !
        total = total_i(ncomp, h)
        allocate(x_l(total, ncomp), x_space(total), y_space(total), P_space(total), y_new(total))
        allocate(phi_liq(ncomp), phi_gas(ncomp))
        
        x_l = gerar_matriz(ncomp, h, total)
        
        P_old = 3*10**5
        do i=1, total
            ! organizo os vetores de composicao para calcular
            if (.NOT. allocated(y_old)) then
                allocate(y_old(total))
                y_old = x_l(i,:)
                y_old = x_l(i,:)
            end if
            print*, i
            ! A BAIXO PODE SER PENSADO NUMA RECURSAOOOOO
            ! calcula o liquido
            x = x_l(i,:)
            call solve_Z(T, P_old, max_iter, tol, .FALSE.)
            call fugacity_solve()
            phi_liq = phi
            
            ! calcula o gas
            x = y_old
            call solve_Z(T, P_old, max_iter, tol, .TRUE.)
            call fugacity_solve()
            phi_gas = phi
            K = 0.0
            do j = 1,ncomp
                K = K + (phi_liq(j) / phi_gas(j)) * x_l(i,j) 
            end do
            do j = 1, ncomp
                y_old(j) = ((phi_liq(j) / phi_gas(j)) * x_l(i,j)) / K
            end do
            
            do j = 1, 3000
                ! aqui checa a condicao de convergencia de K
                ! essa alteracao de P esta BASTANTE porca.
                if (K > 1) then
                    P_old = P_old*1.0002
                else
                    P_old = P_old*0.9998
                end if
                ! Calcula liquido
                x = x_l(i,:)
                call solve_Z(T, P_old, max_iter, tol, .FALSE.)
                call fugacity_solve()
                phi_liq = phi
                ! calcula gas
                x = y_old
                call solve_Z(T, P_old, max_iter, tol, .TRUE.)
                call fugacity_solve()
                phi_gas = phi
                K = 0.0
                do n = 1,ncomp
                    K = K + (phi_liq(n) / phi_gas(n)) * x_l(i,n) 
                end do
                do n = 1, ncomp
                    y_old(n) = ((phi_liq(n) / phi_gas(n)) * x_l(i,n) ) / K
                end do
                if (ABS(1 - K) < tol) exit
            end do
            x_space(i) = x_l(i,1)
            y_space(i) = y_old(1)
            P_space(i) = (P_old/(10**5))

        end do
        
        ! quero testar salvar um arquivo CSV
        print*, 'comecei a salvar'
        open(unit=10, file='dados2.csv', status='replace', action='write')
            do i = 1, total
                write(10, '(F18.8, ",", F18.8, ",",F10.4)') x_space(i), y_space(i), P_space(i)
            end do
        close(10)
        print*, 'terminei a salvar'

    end subroutine bolha_p

    recursive subroutine gen_comp(ncomp, S, pos, current, solutions, idx)
        implicit none
        integer, intent(in) :: ncomp, S, pos
        integer, intent(inout) :: idx
        integer, intent(inout) :: current(ncomp)
        integer, intent(inout) :: solutions(:,:)
        integer :: x
        if (pos == ncomp) then
            current(pos) = S
            solutions(idx, :) = current(:)
            idx = idx + 1
        else
            do x = 0, S
                current(pos) = x
                call gen_comp(ncomp, S - x, pos + 1, current, solutions, idx)
            end do
        end if
    end subroutine gen_comp
    
    ! Calcula quantas linhas a matriz tera
    function total_i(ncomp, res) result(count_)
        implicit none
        integer, intent(in) :: ncomp, res
        integer :: count_, i
        real(8) :: temp
        temp = 1.0d0
        do i = 1, ncomp - 1
            temp = temp * dble(res + ncomp - i) / dble(i)
        end do
        count_ = int(temp + 0.5d0)
    end function total_i
    
    
    function gerar_matriz(ncomp, h, total) result(mat)
        integer, intent(in) :: ncomp, h, total
        real(8), allocatable :: mat(:,:)
        integer, allocatable :: solutions(:,:)
        integer, allocatable :: current(:)
        integer :: i, j, idx
        real(8) :: soma_t

        ! Determina quantas linhas a matriz de fracoes molares tera
        allocate(solutions(total, ncomp))
        allocate(current(ncomp))
        idx = 1
        call gen_comp(ncomp, h, 1, current, solutions, idx)
        
        ! Aloca a solucao na matriz e a normaliza
        allocate(mat(total, ncomp))
        do i = 1, total
            do j = 1, ncomp
                mat(i, j) = dble(solutions(i, j)) / dble(h)
            end do
        end do
  end function gerar_matriz
    
end module PhaseEquilibria
    