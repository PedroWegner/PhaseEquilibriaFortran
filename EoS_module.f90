module EoS_module
    use thermo_data ! isso aqui eh para conseguirmos pegar os dados salvos no modulo thermo_data
    use root_finding
    implicit none
    real(8) :: R
    integer :: i, j
    real(8), allocatable :: a_ind(:), b_ind(:), a_par(:), b_par(:), q_par(:), phi(:), ln_phi(:)
    real(8), allocatable :: a_ij(:,:), b_ij(:,:), k_ij(:,:)
    real(8) :: a_mix, b_mix, beta, q, e, s, I_, Z

    
contains
    subroutine iniciar_vet_EoS_module()
        implicit none
        allocate(a_ind(ncomp), b_ind(ncomp))
        allocate(a_ij(ncomp,ncomp), b_ij(ncomp,ncomp), k_ij(ncomp,ncomp))
        allocate(a_par(ncomp), b_par(ncomp), q_par(ncomp))
        allocate(ln_phi(ncomp), phi(ncomp))
    
    end subroutine iniciar_vet_EoS_module
    ! a partir daqui, teremos subrotinas para calcular os parametros etc    
    subroutine calcular_a(T)
        implicit none
        real(8), intent(in) :: T
        real(8) :: alpha
        
        do i = 1, ncomp
            alpha = (1 + (0.37464 + 1.54226*omega(i) - 0.26992*omega(i)**2) * (1 - sqrt(T/Tc(i))))**2
            a_ind(i) = (0.45724 * R**2 * Tc(i)**2 * alpha) / Pc(i) 
        end do
    end subroutine calcular_a
    
    subroutine calcular_b()
        implicit none
        do i = 1, ncomp
            b_ind(i) = 0.07780 * R * Tc(i) / Pc(i)            
        end do
    end subroutine calcular_b
    
    subroutine calcular_parametros_cruzados()
        implicit none
        ! calcula os a e b cruzados
        do i = 1, ncomp
            do j = 1, ncomp
                a_ij(i, j) = (a_ind(i)*a_ind(j))**0.5*(1 - k_ij(i, j))
                if (i /= j) then
                    b_ij(i, j) = 0
                else
                    b_ij(i, j) = b_ind(i)
                end if
            end do
        end do
    end subroutine calcular_parametros_cruzados
    
    subroutine calcular_mistura()
        implicit none
        a_mix = 0.0
        b_mix = 0.0
        do i=1, ncomp
            do j=1, ncomp
                a_mix = a_mix + (x(i) * x(j) * a_ij(i,j))
                b_mix = b_mix + x(i) * b_ij(i, j)
            end do
        end do
    end subroutine calcular_mistura
    
    subroutine calcular_param_parcial()
        implicit none
        real(8) :: a_aux
        do i = 1, ncomp
            a_aux = - a_mix
            do j = 1, ncomp
                a_aux = a_aux + 2 * x(j) * (a_ind(i) * a_ind(j))**0.5 * (1 - k_ij(i, j))
            end do
            a_par(i) = a_aux
            b_par(i) = b_ind(i)
            q_par(i) = q*(1 + a_par(i) / a_mix - b_par(i) / b_mix)
        end do

    end subroutine
    
    subroutine solve_Z(T, P, max_iter, tol, vapor)
        implicit none
        real(8), intent(in) :: T, P, tol
        integer, intent(in) :: max_iter
        logical, intent(in) :: vapor
        ! executa as funcoes necessarias para o calculo de Z
        call calcular_a(T)
        call calcular_b()
        call calcular_parametros_cruzados()
        call calcular_mistura()
        ! definicao das variaveis de Peng-Robinson
        beta = b_mix*P / (R*T)
        q = a_mix / (b_mix*R*T)
        e = 1 - sqrt(2.0)
        s = 1 + sqrt(2.0)
        
        ! vou precisar rever aqui
        if (vapor) then
            Z = solve_equation(z_gas, dz_gas, real(1.0,8))
        else
            Z = solve_equation(z_liq, dz_liq, beta)
        end if

    end subroutine
    
    subroutine fugacity_solve()
        implicit none
        call calcular_param_parcial()
        I_ = log((Z + s*beta) / (Z + e*beta)) / (s - e)
        
        do i = 1, ncomp
            ln_phi(i) = (b_par(i) / b_mix) * (Z - 1) - log(Z - beta) - q_par(i) * I_
            phi(i) = exp(ln_phi(i))
        end do

    
        
    end subroutine fugacity_solve
    
    ! Definicao das funcoes do modulo
    function z_gas(Z) result(f_z)
        real(8), intent(in) :: Z
        real(8) :: f_z
        ! Equacao de estado cubica generica para gas
        f_z = Z - (1.0+ beta - q*beta*((Z - beta) / ((Z + e*beta) * (Z + s*beta))))
    end function z_gas
    
    function dz_gas(Z) result(df_z)
        real(8), intent(in) :: Z
        real(8) :: df_z
        ! Derivada da equacao de estado cubica generica para gas
        df_z = 1.0 + q*beta*(((Z + e*beta) * (Z + s*beta)) - (Z - beta)*(2*Z + beta * (e + s))) / ((Z + e*beta) * (Z + s*beta))**2
    end function dz_gas
        
    function z_liq(Z) result(f_z)
        real(8), intent(in) :: Z
        real(8) :: f_z
        ! Equacao de estado cubica generica para liquido
        f_z = Z - (beta + (Z + e*beta) * (Z + s*beta) * (1 + beta - Z) / (q*beta))
    end function z_liq
    
    function dz_liq(Z) result(df_z)
        real(8), intent(in) :: Z
        real(8) :: df_z
        df_z = 1 - ((2*Z + beta * (e + s)) * (1 + beta - Z) - (2*Z + beta * (e + s))) / (q*beta)
    end function dz_liq
        
end module EoS_module