module PC_saft_module
    IMPLICIT NONE
    !!! Variaveis globais do modulo
    !!! Constantes globais
    REAL(8), PARAMETER :: pi = 4.0d0 * atan(1.0d0)
    REAL(8), PARAMETER :: k_boltz = 1.3806e-23 !! m2 kg s-2 K-1
    REAL(8), PARAMETER :: N_avogrado = 6.0221e23 ! mol-1
    REAL(8), PARAMETER :: R_gas = 8.314 ! J mol-1 K -1
    ! Parametros do modelo PC-Saft
    REAL(8), ALLOCATABLE :: saft_a0(:), saft_a1(:), saft_a2(:), saft_b0(:), saft_b1(:), saft_b2(:)
    
    TYPE :: PC_saft_state
        ! Variaveis unicas
        INTEGER :: ncomp
        REAL(8) :: T, P, rho, eta
        REAL(8) :: a_res, a_hs, a_hc, a_disp
        REAL(8) :: Z, Z_hc, Z_disp
        REAL(8) :: I_1, I_2, m2es3, m2e2s3, C_1, C_2
        REAL(8) :: mean_m
        
        ! Array
        REAL(8), ALLOCATABLE :: x(:)        ! Composicao em fracao molar
        REAL(8), ALLOCATABLE :: m(:)        ! Numero de segmentos dos componentes
        REAL(8), ALLOCATABLE :: epsilo(:)   ! Profundidade do potencial
        REAL(8), ALLOCATABLE :: sigma(:)    ! Diametro dos segmentos
        REAL(8), ALLOCATABLE :: d(:)        ! Diametro dependente da temperatura
        REAL(8), ALLOCATABLE :: zeta(:)     ! Empacotamento
        REAL(8), ALLOCATABLE :: am(:)       ! EQ.18
        REAL(8), ALLOCATABLE :: bm(:)       ! EQ.19
        
        ! Matrizes de interacao binaria
        REAL(8), ALLOCATABLE :: g_ij(:, :)  ! Correlacao radial
        REAL(8), ALLOCATABLE :: e_ij(:, :)  ! Energia de interecao (eh o epsilon)
        REAL(8), ALLOCATABLE :: s_ij(:, :)  ! Diametro de interacao (eh o sigma)
    end TYPE PC_saft_state
    
    integer :: saft_ncomp
    REAL(8) :: saft_T, saft_P, saft
    REAL(8), allocatable :: saft_x(:), saft_m(:), saft_epsilon(:), saft_sigma(:), saft_d(:), saft_zeta(:), saft_am(:), saft_bm(:)
    REAL(8), allocatable :: saft_g_ij(:,:), saft_e_ij(:, :), saft_s_ij(:, :)
    REAL(8) :: saft_mean_m, saft_rho, saft_a_hs, saft_a_hc, saft_a_disp, saft_a_res, saft_ehta
    REAL(8) :: saft_I_1, saft_I_2, saft_m2es3, saft_m2e2s3, saft_C1, saft_C2    
    REAL(8) :: saft_Z_hc, saft_Z_disp, saft_Z
    

CONTAINS

    SUBROUTINE universal_model_constants()
        IMPLICIT NONE
        ALLOCATE(saft_a0(7), saft_a1(7), saft_a2(7), saft_b0(7), saft_b1(7), saft_b2(7))
        ALLOCATE(saft_am(7), saft_bm(7))
        !!! Parametros universais do modelo PC-Saft
        saft_a0 = [0.910563145, 0.636128145, 2.686134789, -26.54736249, 97.75920878, -159.5915409, 91.29777408]
        saft_a1 = [-0.308401692, 0.186053116, -2.503004726, 21.41979363, -65.25588533, 83.31868048, -33.74692293]
        saft_a2 = [-0.090614835, 0.452784281, 0.596270073, -1.724182913, -4.130211253, 13.77663187, -8.672847037]
        saft_b0 = [0.724094694, 2.238279186, -4.002584949, -21.00357682, 26.85564136, 206.5513384, -355.6023561]
        saft_b1 = [-0.575549808, 0.699509552, 3.892567339, -17.21547165, 192.6722645, -161.8264617, -165.2076935]
        saft_b2 = [0.097688312, -0.255757498, -9.155856153, 20.64207597, -38.80443005, 93.62677408, -29.66690559]
    END SUBROUTINE universal_model_constants
    
    !!! ABAIXO SAO AS SUBROTINAS PARA MODIFICAR O ESTADO !!!
    SUBROUTINE init_mixture(state, ncomp, x, m, sigma, epsilo)
        implicit none
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER, INTENT(IN) :: ncomp
        REAL(8), INTENT(IN) :: x(:), m(:), sigma(:), epsilo(:)
        
        state%ncomp = ncomp
        ALLOCATE(state%x(ncomp), state%m(ncomp), state%sigma(ncomp), state%epsilo(ncomp))
        state%x = x
        state%m = m
        state%sigma = sigma
        state%epsilo = epsilo
    END SUBROUTINE init_mixture
    
    SUBROUTINE calc_diameter_T(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i
        !!! EQ (A.9)
        DO i = 1, SIZE(state%x)
            state%d(i) = state%sigma(i) * (1.0d0 - 0.12d0 * EXP(-3.0d0 * state%epsilo(i) / state%T))
        END DO
    END SUBROUTINE calc_diameter_T
    
    SUBROUTINE calc_mean_m(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i
        state%mean_m = 0.0d0
        DO i = 1, SIZE(state%x)
            state%mean_m = state%mean_m + state%x(i) * state%m(i)
        END DO
    END SUBROUTINE calc_mean_m
    
    SUBROUTINE calc_combining_rules(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i, j
        !!! EQ (A.14) & (A.15)
        DO i = 1, SIZE(state%x)
            DO j = 1, SIZE(state%x)
                state%s_ij(i, j) = (state%sigma(i) + state%sigma(j)) / 2
                state%e_ij(i, j) = (state%epsilo(i) * state%epsilo(j)) ** 0.5
            END DO
        END DO
    END SUBROUTINE calc_combining_rules
    
    SUBROUTINE calc_zeta(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: soma_aux, zeta_aux
        INTEGER :: n, i
        
        IF (.NOT. ALLOCATED(state%zeta)) THEN
            ALLOCATE(state%zeta(4))
        END IF
        !!! EQ (A.8)
        zeta_aux = pi * state%rho / 6
        DO n = 1, 4
            soma_aux = 0.0d0 
            DO i = 1, SIZE(state%x)
                soma_aux = soma_aux + (state%x(i) * state%m(i) * state%d(i) ** (n - 1))
            END DO
            state%zeta(n) = zeta_aux * soma_aux
        END DO
    END SUBROUTINE calc_zeta
    
    SUBROUTINE calc_hard_sphere(state)
        IMPLICIT NONE    
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3, d_ij
        INTEGER :: i, j
        !!! EQ (A.7)
        aux_1 = 1 - state%zeta(4)
        aux_2 = 3 * state%zeta(3) / aux_1 ** 2
        aux_3 = 2 * state%zeta(3)**2 / aux_1 ** 3
        
        DO i = 1, SIZE(state%x)
            DO j = 1, SIZE(state%x)
                d_ij = (state%d(i) * state%d(j)) / (state%d(i) + state%d(j))
                state%g_ij(i, j) = 1.0d0 / aux_1 + d_ij * aux_2 + d_ij ** 2 * aux_3
            END DO
        END DO
    END SUBROUTINE calc_hard_sphere
    
    
    SUBROUTINE calc_ab_mean_m(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2
        INTEGER :: i
        
        IF (.NOT. ALLOCATED(state%am) .AND. .NOT. ALLOCATED(state%bm)) THEN
            ALLOCATE(state%am(7), state%bm(7))
        END IF
        !!! EQ (A.18) & (A.19)
        aux_1 = (state%mean_m - 1) / state%mean_m
        aux_2 = (state%mean_m -1) * (state%mean_m - 2) / state%mean_m ** 2
        DO i = 1, 7
            state%am(i) = saft_a0(i) * aux_1 * saft_a1(i) + aux_2 * saft_a2(i)
            state%bm(i) = saft_b0(i) * aux_1 * saft_b1(i) + aux_2 * saft_b2(i)
        END DO
    END SUBROUTINE calc_ab_mean_m
    
    SUBROUTINE calc_pertubation_integral(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i
        state%I_1 = 0.0d0
        state%I_2 = 0.0d0
        !!! EQ (A.16) & (A.17)
        DO i = 1, 7
            state%I_1 = state%I_1 + state%am(i) * state%eta ** (i-1)
            state%I_2 = state%I_2 + state%bm(i) * state%eta ** (i-1)
        END DO
    END SUBROUTINE calc_pertubation_integral
    
    SUBROUTINE calc_C_12(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3
        
        !!! EQ (A.11)
        aux_1 = state%mean_m * (8 * state%eta - 2 * state%eta ** 2) / (1 - state%eta) ** 2
        aux_2 = (1 - state%mean_m) * (20 * state%eta -27 * state%eta **2 + 12 * state%eta ** 3 -2 * state%eta ** 40)
        aux_3 = ((1 - state%eta) * (2 - state%eta)) ** 2
        state%C_1 = 1 / (1 + aux_1 + aux_2 / aux_3)

        !!! EQ (A.31)
        aux_1 = state%mean_m * (-4 * state%eta ** 2 + 20 * state%eta + 8) / (1 - state%eta) ** 5 
        aux_2 = (1 - state%mean_m) * (2 * state%eta ** 3 + 12 * state%eta ** 2 - 48 * state%eta + 40)
        aux_3 = ((1 - state%eta) * (2 - state%eta)) ** 3
        state%C_2 = - state%C_1 ** 2 * (aux_1 + aux_2 / aux_3)
    END SUBROUTINE calc_C_12
    
    SUBROUTINE calc_abbr_mes(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i, j
        REAL(8) :: x_ij, m_ij, e_ij
         !!! EQ (A.12) & (A.13)
        
        state%m2es3 = 0
        state%m2e2s3 = 0
        DO i = 1, SIZE(state%x)
            DO j = 1, SIZE(state%x)
                x_ij = state%x(i) * state%x(j)
                m_ij = state%m(i) * state%m(j)
                e_ij = state%e_ij(i, j) / state%T
                state%m2es3 =   state%m2es3 + (x_ij * m_ij * e_ij * state%s_ij(i, j) ** 3)
                state%m2e2s3 =   state%m2es3 + (x_ij * m_ij * e_ij **2 * state%s_ij(i, j) ** 3)
            END DO
        END DO
    END SUBROUTINE calc_abbr_mes
    
    SUBROUTINE calc_rho(state)
        implicit none
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: sum_aux = 0.0d0
        INTEGER :: i
        !!! EQ (A.20)
        DO i = 1, SIZE(state%x)
            sum_aux = sum_aux + state%x(i) * state%m(i) * state%d(i)**3
        end do
        state%rho = (6.0d0 * state%eta / pi) / sum_aux
    END SUBROUTINE calc_rho
    
    
    
    !!!!! ABAIXO NAO FOI REFATORADOOOOOOOOOOO
    subroutine initialise_mixture(ncomp, mole_fraction, number_segment, segment_diameter, depth_potential)
        implicit none
        integer, intent(in) :: ncomp
        real(8), intent(in) :: mole_fraction(:), number_segment(:), depth_potential(:), segment_diameter(:)
        saft_ncomp = ncomp
        
        allocate(saft_x(saft_ncomp), saft_m(saft_ncomp), saft_epsilon(saft_ncomp), saft_sigma(saft_ncomp))
        allocate(saft_d(saft_ncomp))
        allocate(saft_g_ij(saft_ncomp,saft_ncomp), saft_e_ij(saft_ncomp,saft_ncomp), saft_s_ij(saft_ncomp,saft_ncomp))
        saft_x = mole_fraction
        saft_m = number_segment
        saft_epsilon = depth_potential
        saft_sigma = segment_diameter
        !! atualiza valores da array de diametro dependendo da T
        !call calc_diameter_T
        !call calc_mean_m()
        !call calc_combining_rules()
        
    end subroutine initialise_mixture
    
    subroutine up_mixture()
        call calc_diameter_T()
        call calc_rho()
        call calc_mean_m()
        call calc_ab_mean_m()
        call calc_pertubation_integral()
        call calc_zeta()
        call calc_hard_sphere()
        call calc_combining_rules()
    
    end subroutine
    
    subroutine update_mixture(mole_fraction, number_segment, depth_potential, segment_diameter)
        implicit none
        real(8), intent(in) :: mole_fraction, number_segment, depth_potential, segment_diameter
        saft_x = mole_fraction
        saft_m = number_segment
        saft_epsilon = depth_potential
        saft_sigma = segment_diameter
        !! atualiza valores da array de diametro dependendo da T
        call calc_diameter_T
        
    end subroutine update_mixture
    
    subroutine calc_rho()
        implicit none
        real(8) :: sum_aux = 0
        integer :: i
        !!! EQ (A.20)
        do i = 1, saft_ncomp
            sum_aux = sum_aux + saft_x(i) * saft_m(i) * saft_d(i) ** 3
        end do
        saft_rho = (6 * saft_ehta / pi) / sum_aux
    end subroutine calc_rho
    subroutine calc_combining_rules()
        implicit none
        integer :: i, j
        !!! EQ (A.14) & (A.15)
        do i = 1, saft_ncomp
            do j = 1, saft_ncomp
                saft_s_ij(i, j) = (saft_sigma(i) + saft_sigma(j)) / 2
                saft_e_ij(i, j) = (saft_epsilon(i) * saft_epsilon(j))**0.5 !! ignorando o k_ij aqui
            end do
        end do
    end subroutine calc_combining_rules
    
    subroutine calc_mean_m()
        implicit none
        integer :: i
        saft_mean_m = 0
        do i = 1, saft_ncomp
            saft_mean_m = saft_mean_m + saft_x(i) * saft_m(i) 
        end do
    end subroutine calc_mean_m
    
    subroutine calc_diameter_T()
        implicit none
        integer :: i
        !!! EQ (A.9)
        do i = 1, saft_ncomp
            saft_d(i) = saft_sigma(i) * (1 - 0.12 * exp(-3 * saft_epsilon(i) / saft_T))
        end do
    end subroutine calc_diameter_T
    
    subroutine calc_zeta()
        implicit none
        integer :: n, i
        real(8) :: soma_aux, zeta_aux

        if (.NOT. allocated(saft_zeta)) then
            allocate(saft_zeta(4))
        end if
        
        !!! EQ (A.8)
        do n = 1, 4
            zeta_aux = pi * saft_rho / 6
            soma_aux = 0
            do i = 1, saft_ncomp
                soma_aux = soma_aux + (saft_x(i) * saft_m(i) * saft_d(i) ** (n - 1))
            end do
            saft_zeta(n) = zeta_aux * soma_aux
        end do
        call calc_hard_sphere()

    end subroutine calc_zeta
    
    subroutine calc_hard_sphere()
        implicit none
        integer :: i, j
        real(8) :: term_1, term_2, d_comb
        
        !!! EQ (A.7)
        do i = 1, saft_ncomp
            do j = 1, saft_ncomp
                d_comb = (saft_d(i) * saft_d(j)) / (saft_d(i) + saft_d(j))
                term_1 =  d_comb * (3 * saft_zeta(3) / (1 - saft_zeta(4))**2)
                term_2 = d_comb ** 2 * (2 * saft_zeta(3) ** 2 ) / (1 - saft_zeta(4)) ** 3
                saft_g_ij(i, j) = 1 / (1 - saft_zeta(4)) + term_1 + term_2
            end do
        end do
    end subroutine calc_hard_sphere
    
    subroutine calc_ab_mean_m()
        implicit none
        integer :: i
        real(8) :: m_1, m_2
        m_1 = (saft_mean_m - 1) / saft_mean_m
        m_2 = (saft_mean_m - 1) * (saft_mean_m - 2) / saft_mean_m**2
        
        !!! EQ (A.18) & (A.19)
        do i = 1, 7
            saft_am(i) = saft_a0(i) + m_1 * saft_a1(i) + m_2 * saft_a2(i)
            saft_bm(i) = saft_b0(i) + m_1 * saft_b1(i) + m_2 * saft_b2(i)
        end do
    end subroutine calc_ab_mean_m
    
    subroutine calc_pertubation_integral()
        implicit none
        integer :: i
        saft_I_1 = 0
        saft_I_2 = 0
        call calc_ab_mean_m() !!!!!!!!!!!!! REMOVER!!!!!!!!!!!!!!!!!!!!!!!1
        !!! EQ (A.16) & (A.17)
        do i = 1, 7
            saft_I_1 = + saft_I_1 + saft_am(i) * saft_ehta**(i - 1)
            saft_I_2 = + saft_I_2 + saft_bm(i) * saft_ehta**(i - 1)
        end do
    end subroutine calc_pertubation_integral
    
    subroutine helmholtz_hard_chain()
        implicit none
        real(8) :: a_hs_1, a_hs_2, a_hs_3, soma_aux, zeta_aux
        integer :: i
        
        !!! TALVEZ TENHA QUE MUDAR O LUGAR
        call calc_zeta()  !!!!! REMOVER!!!!!
        
        
        !!! EQ (A.6)
        zeta_aux = 1 - saft_zeta(4)
        a_hs_1 = 3 * saft_zeta(2) * saft_zeta(3) / zeta_aux
        a_hs_2 = (saft_zeta(3)**3) / (saft_zeta(4) * zeta_aux**2)
        a_hs_3 = ((saft_zeta(3) ** 3/ saft_zeta(4) ** 2) - saft_zeta(1)) * log(1 - saft_zeta(4))
        saft_a_hs = (1 / saft_zeta(1)) * (a_hs_1 + a_hs_2 + a_hs_3)
        
        soma_aux = 0
        do i = 1, saft_ncomp
            soma_aux = soma_aux + saft_x(i) * (saft_m(i) - 1) * log(saft_g_ij(i, i))
        end do
        !!! EQ (A.4)
        saft_a_hc = saft_mean_m * saft_a_hs - soma_aux
    end subroutine helmholtz_hard_chain
    
    subroutine helmholtz_dispersion()
        implicit none
        integer :: i, j
        real(8) :: C_1, C_2, C_3
        real(8) :: sum_aux_1, sum_aux_2, a_disp_1, a_disp_2
        
        !!! EQ (A.11)
        C_1 = saft_mean_m * (8 * saft_ehta - 2 * saft_ehta ** 2) / (1 - saft_ehta)**4
        C_2 = (1 - saft_mean_m) * (20*saft_ehta - 27*saft_ehta**2 + 12*saft_ehta**3 - 2*saft_ehta**4)
        C_3 = ((1-saft_ehta)*(2-saft_ehta))**2
        saft_C1 = (1 + C_1 + C_2 / C_3)**(-1)
        
        !!! EQ (A.12) & (A.13)
        
        saft_m2es3 = 0
        saft_m2e2s3 = 0
        do i = 1, saft_ncomp
            do j = 1, saft_ncomp
                saft_m2es3 = saft_m2es3 + saft_x(i)*saft_x(j)*saft_m(i)*saft_m(j)*(saft_e_ij(i, j) / saft_T)*saft_s_ij(i, j)**3
                saft_m2e2s3 = saft_m2e2s3 + saft_x(i)*saft_x(j)*saft_m(i)*saft_m(j)*(saft_e_ij(i, j) / saft_T)**2*saft_s_ij(i, j)**3
            end do
        end do
        
        !! aqui tem que colocar para calcular o I1 e I2
        call calc_pertubation_integral()
        a_disp_1 = - 2 * pi * saft_rho * saft_I_1 *  saft_m2es3
        a_disp_2 = - pi * saft_rho * saft_mean_m * saft_C1 * saft_I_2 * saft_m2e2s3
        saft_a_disp = a_disp_1 + a_disp_2
    end subroutine helmholtz_dispersion
    
    subroutine compressibility_hard_chain()
        implicit none
        real(8), allocatable :: prime_g_rho(:)
        real(8) :: prime_term_1, prime_term_2, prime_term_3, comb_d, aux, soma_aux
        real(8) :: z_term_1, z_term_2, z_term_3, z_hs
        integer :: i
        
        if (.NOT. allocated(prime_g_rho)) then
            allocate(prime_g_rho(saft_ncomp))
        end if
        
        !!! EQ (A.27)
        aux = 1 - saft_zeta(4)
        prime_term_1 = saft_zeta(4) / aux**2
        prime_term_2 = 3*saft_zeta(3) / aux ** 2 + 6*saft_zeta(3) * saft_zeta(4) / aux**3
        prime_term_3 = 4*saft_zeta(3)**2 / aux**3 + 6*saft_zeta(3)**2 * saft_zeta(4) / aux**4
        
        soma_aux = 0
        !!! talvez mude aqui depois
        !!! porque no modelo ele calcula o prime_g_ij_rho, eu apenas o  prime_g_ii_rho
        do i = 1, saft_ncomp
            comb_d = saft_d(i) / 2
            prime_g_rho(i) = prime_term_1 + comb_d * prime_term_2 + comb_d**2 * prime_term_3
            soma_aux = soma_aux + (saft_x(i) * (saft_m(i) - 1) * saft_g_ij(i,i) **(-1) * prime_g_rho(i))
        end do
        
        !!! EQ (A.26)
        z_term_1 = saft_zeta(4) / aux
        z_term_2 = 3 * saft_zeta(2) * saft_zeta(3) / (saft_zeta(1) * aux**2)
        z_term_3 = (3*saft_zeta(3)**3 - saft_zeta(4) * saft_zeta(2)**3) / (saft_zeta(1) * aux **3)
        z_hs = z_term_1 + z_term_2 + z_term_3
        !!! EQ (A.25)
        saft_Z_hc = saft_mean_m * z_hs - soma_aux
    end subroutine compressibility_hard_chain
    
    subroutine compressibility_dispersion()
        implicit none
        real(8) :: c2_term_1, c2_term_2
        real(8) :: prime_I1, prime_I2
        real(8) :: z_term_1, z_term_2
        integer :: i

        !!! EQ (A.31)
        c2_term_1 = saft_mean_m * ((-4*saft_ehta**2 +20*saft_ehta + 8) / (1 - saft_ehta)**5)
        c2_term_2 = (1-saft_mean_m) * ((2-saft_ehta**3 + 12*saft_ehta**2 - 48*saft_ehta + 40) / ((1-saft_ehta)*(2-saft_ehta))**3)
        saft_C2 = - saft_C1**2 * (c2_term_1 + c2_term_2)
        
        !!! EQ (A.29) & (A.30)
        prime_I1 = 0
        prime_I2 = 0
        do i = 1, 7
            prime_I1 = prime_I1 + saft_am(i) * i * saft_ehta**(i-1)
            prime_I2 = prime_I2 + saft_bm(i) * i * saft_ehta**(i-1)
        end do
        
        !!! EQ (A.28)
        z_term_1 = -2 * pi * saft_rho * prime_I1 * saft_m2es3
        z_term_2 = -pi * saft_rho * saft_mean_m * (saft_C1 * prime_I2 + saft_C2 * saft_ehta * saft_I_2) * saft_m2e2s3
        saft_Z_disp = z_term_1 + z_term_2
    end subroutine compressibility_dispersion
    
    
    subroutine residual_helmholtz()
        implicit none
        !!! Preenche a_hc e a_disp
        call helmholtz_hard_chain
        call helmholtz_dispersion
        !!! EQ (A.3)
        saft_a_res = saft_a_hc + saft_a_disp
    end subroutine residual_helmholtz
    
    subroutine compressibility()
        implicit none
        !!! Preenche Z_hc e Z_disp
        call compressibility_hard_chain
        call compressibility_dispersion
        !!! EQ (A.24)
        saft_Z = 1 + saft_Z_hc + saft_Z_disp
    end subroutine compressibility
    
    !!! SUBROTINAS PARA MINIMIZACAO
    !!! acredito que para 
    subroutine residuo_pression(m, n, x, fvec, iflag)
        implicit none
        integer, intent(in) :: m, n, iflag
        real(8), intent(in) :: x(n)
        real(8), intent(out) :: fvec(m)
        real(8) :: ehta, P_calc, P_sis
        !!! Abaixo eh definido o parametro que vai ser lido na minimizacao
        ehta = x(1)
        !!! A pressao calculada
        P_calc = calc_pressao(ehta)
        P_sis = 30.0d5       
        !!! Abaixo eh definido o residuo a ser minimizado
        fvec(1) = (P_calc - P_sis)/P_sis
    end subroutine residuo_pression
    
    !!! FUNCOES 
    function calc_pressao(ehta) result(pressao)
        implicit none
        real(8), intent(in) :: ehta
        real(8) :: pressao
        !!! Atribui valor ao ehta da rotina e atualiza parametros dependentes
        saft_ehta = ehta
        call up_mixture()
        !!! Calcula helmholtz residual e a compressibilidade Z
        call residual_helmholtz()
        call compressibility
        !!! EQ (A.23)
        pressao = (saft_Z * k_boltz * saft_rho * saft_T)*1e30
        print*, pressao
    end function calc_pressao
    
end module PC_saft_module