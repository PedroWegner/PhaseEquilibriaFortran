module PC_saft_module
    IMPLICIT NONE
    !!! PARA O TESTE
    INTEGER :: ncomp_global
    REAL(8) :: T_global, P_global
    REAL(8), ALLOCATABLE :: x_global(:)
    REAL(8), ALLOCATABLE :: m_global(:)
    REAL(8), ALLOCATABLE :: sigma_global(:)
    REAL(8), ALLOCATABLE :: epsilo_global(:)
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
        REAL(8) :: Z, Z_hs, Z_hc, Z_disp
        REAL(8) :: I_1, I_2, m2es3, m2e2s3, C_1, C_2
        REAL(8) :: prime_I_1, prime_I_2             ! EQ (A.29) & (A.30)
        REAL(8) :: mean_m
        
        ! Array
        REAL(8), ALLOCATABLE :: x(:)                ! Composicao em fracao molar
        REAL(8), ALLOCATABLE :: m(:)                ! Numero de segmentos dos componentes
        REAL(8), ALLOCATABLE :: epsilo(:)           ! Profundidade do potencial
        REAL(8), ALLOCATABLE :: sigma(:)            ! Diametro dos segmentos
        REAL(8), ALLOCATABLE :: d(:)                ! Diametro dependente da temperatura
        REAL(8), ALLOCATABLE :: zeta(:)             ! Empacotamento
        REAL(8), ALLOCATABLE :: am(:)               ! EQ.18
        REAL(8), ALLOCATABLE :: bm(:)               ! EQ.19
        
        ! Matrizes de interacao binaria
        REAL(8), ALLOCATABLE :: g_ij(:, :)          ! Correlacao radial
        REAL(8), ALLOCATABLE :: e_ij(:, :)          ! Energia de interecao (eh o epsilon)
        REAL(8), ALLOCATABLE :: s_ij(:, :)          ! Diametro de interacao (eh o sigma)
        REAL(8), ALLOCATABLE :: grad_g_ij_rho(:, :) ! Definido pela EQ (A.27)
    end TYPE PC_saft_state
    

CONTAINS

    SUBROUTINE universal_model_constants()
        IMPLICIT NONE
        ALLOCATE(saft_a0(7), saft_a1(7), saft_a2(7), saft_b0(7), saft_b1(7), saft_b2(7))
        !!! Parametros universais do modelo PC-Saft
        saft_a0 = [0.910563145, 0.636128145, 2.686134789, -26.54736249, 97.75920878, -159.5915409, 91.29777408]
        saft_a1 = [-0.308401692, 0.186053116, -2.503004726, 21.41979363, -65.25588533, 83.31868048, -33.74692293]
        saft_a2 = [-0.090614835, 0.452784281, 0.596270073, -1.724182913, -4.130211253, 13.77663187, -8.672847037]
        saft_b0 = [0.724094694, 2.238279186, -4.002584949, -21.00357682, 26.85564136, 206.5513384, -355.6023561]
        saft_b1 = [-0.575549808, 0.699509552, 3.892567339, -17.21547165, 192.6722645, -161.8264617, -165.2076935]
        saft_b2 = [0.097688312, -0.255757498, -9.155856153, 20.64207597, -38.80443005, 93.62677408, -29.66690559]
    END SUBROUTINE universal_model_constants
    
    !!! ABAIXO SAO AS SUBROTINAS PARA MODIFICAR O ESTADO !!!
    SUBROUTINE allocate_state_arrays(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: n
        
        n = state%ncomp
        !!! Aloca os arrays vetores
        IF (.NOT. ALLOCATED(state%x)) ALLOCATE(state%x(n))
        IF (.NOT. ALLOCATED(state%m)) ALLOCATE(state%m(n))
        IF (.NOT. ALLOCATED(state%sigma)) ALLOCATE(state%sigma(n))
        IF (.NOT. ALLOCATED(state%epsilo)) ALLOCATE(state%epsilo(n))
        IF (.NOT. ALLOCATED(state%d)) ALLOCATE(state%d(n))
        IF (.NOT. ALLOCATED(state%zeta)) ALLOCATE(state%zeta(4))    ! Fixo em 4 elementos (note que no modelo eh de 0 a 3, aqui de 1 a 4
        IF (.NOT. ALLOCATED(state%am)) ALLOCATE(state%am(7))        ! Fixo em 7 elementos
        IF (.NOT. ALLOCATED(state%bm)) ALLOCATE(state%bm(7))        ! Fixo em 7 elementos
        
        ! Aloca arrays matrizes
        IF (.NOT. ALLOCATED(state%g_ij)) ALLOCATE(state%g_ij(n, n))
        IF (.NOT. ALLOCATED(state%e_ij)) ALLOCATE(state%e_ij(n, n))
        IF (.NOT. ALLOCATED(state%s_ij)) ALLOCATE(state%s_ij(n, n))
        IF (.NOT. ALLOCATED(state%grad_g_ij_rho)) ALLOCATE(state%grad_g_ij_rho(n, n))
    END SUBROUTINE allocate_state_arrays
    
    SUBROUTINE destroy_state(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        !!! Dealoca as arrays vetores
        IF (ALLOCATED(state%x)) DEALLOCATE(state%x)
        IF (ALLOCATED(state%m)) DEALLOCATE(state%m)
        IF (ALLOCATED(state%sigma)) DEALLOCATE(state%sigma)
        IF (ALLOCATED(state%epsilo)) DEALLOCATE(state%epsilo)
        IF (ALLOCATED(state%d)) DEALLOCATE(state%d)
        IF (ALLOCATED(state%zeta)) DEALLOCATE(state%zeta)    
        IF (ALLOCATED(state%am)) DEALLOCATE(state%am)       
        IF (ALLOCATED(state%bm)) DEALLOCATE(state%bm)        
        
        ! Dealoca as arrays matrizes
        IF (ALLOCATED(state%g_ij)) DEALLOCATE(state%g_ij)
        IF (ALLOCATED(state%e_ij)) DEALLOCATE(state%e_ij)
        IF (ALLOCATED(state%s_ij)) DEALLOCATE(state%s_ij)
        IF (ALLOCATED(state%grad_g_ij_rho)) DEALLOCATE(state%grad_g_ij_rho)
        
        ! Zera os escalares do state
        state%T = 0.0d0
        state%P = 0.0d0
        state%rho = 0.0d0
        state%eta = 0.0d0
        state%a_res = 0.0d0
        state%a_hs = 0.0d0
        state%a_hc = 0.0d0
        state%a_disp = 0.0d0
        state%Z = 0.0d0
        state%Z_hs = 0.0d0
        state%Z_hc = 0.0d0
        state%Z_disp = 0.0d0
        state%I_1 = 0.0d0
        state%I_2 = 0.0d0
        state%m2es3 = 0.0d0
        state%m2e2s3 = 0.0d0
        state%C_1 = 0.0d0
        state%C_2 = 0.0d0
        state%prime_I_1 = 0.0d0
        state%prime_I_2 = 0.0d0
        state%mean_m = 0.0d0

    END SUBROUTINE destroy_state
    
    SUBROUTINE init_mixture(state, ncomp, x, m, sigma, epsilo)
        implicit none
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER, INTENT(IN) :: ncomp
        REAL(8), INTENT(IN) :: x(:), m(:), sigma(:), epsilo(:)
        
        state%ncomp = ncomp
        CALL allocate_state_arrays(state)
        
        state%x = x
        state%m = m
        state%sigma = sigma
        state%epsilo = epsilo
    END SUBROUTINE init_mixture
    
    SUBROUTINE update_state(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        !!! Atualiza as variaveis do modelo PC-Saft de acordo com o eta e Temperatura entradas
        CALL calc_diameter_T(state)
        CALL calc_rho(state)
        CALL calc_mean_m(state)
        CALL calc_combining_rules(state)
        CALL calc_zeta(state)
        CALL calc_hard_sphere(state)
        CALL calc_ab_mean_m(state)
        CALL calc_pertubation_integral(state)
        CALL calc_prime_pertubation_integral(state)
        CALL calc_C_12(state)
        CALL calc_abbr_mes(state)
        CALL calc_grad_g_rho(state)
        
        
    END SUBROUTINE update_state
    
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
        
        !!! EQ (A.18) & (A.19)
        aux_1 = (state%mean_m - 1) / state%mean_m
        aux_2 = (state%mean_m -1) * (state%mean_m - 2) / state%mean_m ** 2
        DO i = 1, 7
            state%am(i) = saft_a0(i) + aux_1 * saft_a1(i) + aux_2 * saft_a2(i)
            state%bm(i) = saft_b0(i) + aux_1 * saft_b1(i) + aux_2 * saft_b2(i)
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
    
    SUBROUTINE calc_prime_pertubation_integral(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: j
        
        !!! EQ (A.29) & (A.30)
        state%prime_I_1 = 0.0d0
        state%prime_I_2 = 0.0d0
        DO j = 1, 7
            state%prime_I_1 = state%prime_I_1 + state%am(j) * j * state%eta ** (j - 1)
            state%prime_I_2 = state%prime_I_2 + state%bm(j) * j * state%eta ** (j - 1)
        END DO
    END SUBROUTINE calc_prime_pertubation_integral
    
    SUBROUTINE calc_C_12(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3
        
        !!! EQ (A.11)
        aux_1 = state%mean_m * (8 * state%eta - 2 * state%eta **2) / (1 - state%eta) **4
        aux_2 = (1 - state%mean_m) * (20 * state%eta - 27 * state%eta **2 + 12 * state%eta **3 - 2 * state%eta **4)
        aux_3 = ((1 - state%eta) * (2 - state%eta)) **2
        state%C_1 = 1 / (1 + aux_1 + aux_2 / aux_3)

        !!! EQ (A.31)
        aux_1 = state%mean_m * (-4 * state%eta **2 + 20 * state%eta + 8) / (1 - state%eta) **5 
        aux_2 = (1 - state%mean_m) * (2 * state%eta **3 + 12 * state%eta ** 2 - 48 * state%eta + 40)
        aux_3 = ((1 - state%eta) * (2 - state%eta)) **3
        state%C_2 = - state%C_1 **2 * (aux_1 + aux_2 / aux_3)
    END SUBROUTINE calc_C_12
    
    SUBROUTINE calc_abbr_mes(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        INTEGER :: i, j
        REAL(8) :: x_ij, m_ij, e_ij
         !!! EQ (A.12) & (A.13)
        
        state%m2es3 = 0.0d0
        state%m2e2s3 = 0.0d0
        DO i = 1, SIZE(state%x)
            DO j = 1, SIZE(state%x)
                x_ij = state%x(i) * state%x(j)
                m_ij = state%m(i) * state%m(j)
                e_ij = state%e_ij(i, j) / state%T
                state%m2es3 =   state%m2es3 + (x_ij * m_ij * e_ij * state%s_ij(i, j) ** 3)
                state%m2e2s3 =   state%m2e2s3 + (x_ij * m_ij * e_ij ** 2 * state%s_ij(i, j) ** 3)
            END DO
        END DO
    END SUBROUTINE calc_abbr_mes
    
    SUBROUTINE calc_grad_g_rho(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3, zeta_aux, d_ij
        INTEGER :: i, j
        
        !!! EQ (A.27)
        zeta_aux = 1 - state%zeta(4)
        aux_1 = state%zeta(4) / zeta_aux
        aux_2 = 3*state%zeta(3) / zeta_aux ** 2 + 6 * state%zeta(3) * state%zeta(4) / zeta_aux ** 3
        aux_3 = 4 * state%zeta(3) ** 2 / zeta_aux ** 3 + 6 * state%zeta(3) ** 2 * state%zeta(4) / zeta_aux **4
        
        DO i = 1, SIZE(state%x)
            DO j = 1, SIZE(state%x)
                d_ij = state%d(i) * state%d(j) / (state%d(i) + state%d(j))
                state%grad_g_ij_rho(i, j) = aux_1 + d_ij * aux_2 + d_ij ** 2 * aux_3
            END DO
        END DO
    END SUBROUTINE calc_grad_g_rho
    
    SUBROUTINE calc_rho(state)
        implicit none
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: sum_aux
        INTEGER :: i
        !!! EQ (A.20)
        sum_aux = 0.d0
        DO i = 1, SIZE(state%x)
            sum_aux = sum_aux + state%x(i) * state%m(i) * state%d(i)**3
        end do
        state%rho = (6.0d0 * state%eta / pi) / sum_aux
    END SUBROUTINE calc_rho
    
    SUBROUTINE helmholtz_hard_chain(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3, soma_aux, zeta_aux
        INTEGER :: i
        
        !!! EQ (A.6)
        zeta_aux = 1 - state%zeta(4)
        aux_1 = 3 * state%zeta(2) * state%zeta(3) / zeta_aux
        aux_2 = state%zeta(3) ** 3 / (state%zeta(4) * zeta_aux ** 2)
        aux_3 = (state%zeta(3) ** 3 / state%zeta(4) ** 2 - state%zeta(1)) * LOG(1 - state%zeta(4))
        state%a_hs = (1 / state%zeta(1)) * (aux_1 + aux_2 + aux_3)
        
        soma_aux = 0.0d0
        DO i = 1, SIZE(state%x)
            soma_aux = soma_aux + state%x(i) * (state%m(i) - 1) * LOG(state%g_ij(i, i))
        END DO
        state%a_hc = state%mean_m * state%a_hs - soma_aux
    END SUBROUTINE helmholtz_hard_chain    
    
    SUBROUTINE helmholtz_dispersion(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2
        
        !!! EQ (A.10)
        aux_1 = - 2 * pi * state%rho * state%I_1 * state%m2es3
        aux_2 = - pi * state%rho * state%mean_m * state%C_1 * state%I_2 * state%m2e2s3
        state%a_disp = aux_1 + aux_2
    END SUBROUTINE helmholtz_dispersion
    
    SUBROUTINE helmholtz_residual(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        !!! EQ (A.3)
        CALL helmholtz_hard_chain(state)
        CALL helmholtz_dispersion(state)
        state%a_res = state%a_hc + state%a_disp        
    END SUBROUTINE helmholtz_residual
    
    SUBROUTINE compressibility_hard_chain(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2, aux_3, zeta_aux, sum_aux
        INTEGER :: i
        
        !!! EQ (A.26)
        zeta_aux = 1 - state%zeta(4)
        aux_1 = state%zeta(4) / zeta_aux
        aux_2 = 3 * state%zeta(2) * state%zeta(3) / (state%zeta(1) * zeta_aux ** 2)
        aux_3 = (3 * state%zeta(3) ** 3 - state%zeta(4) * state%zeta(3) ** 3) / (state%zeta(1) * zeta_aux ** 3)
        state%Z_hs = aux_1 + aux_2 + aux_3
        
        !!! EQ (A.25)
        sum_aux = 0.0d0
        DO i = 1, SIZE(state%x)
            sum_aux = sum_aux + state%x(i) * (state%m(i) - 1) * state%g_ij(i, i) ** (-1) * state%grad_g_ij_rho(i, i)
        END DO
        state%Z_hc = state%mean_m * state%Z_hs - sum_aux
    END SUBROUTINE compressibility_hard_chain
    
    SUBROUTINE compressibility_dispersion(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        REAL(8) :: aux_1, aux_2
        
        !!! EQ (A.28)
        aux_1 = - 2 * pi * state%rho * state%prime_I_1 * state%m2es3
        aux_2 = - pi * state%rho * (state%C_1 * state%prime_I_2 + state%C_2 * state%eta * state%I_2) * state%m2e2s3
        state%Z_disp = aux_1 + aux_2
    END SUBROUTINE compressibility_dispersion
    
    SUBROUTINE compressibility(state)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(INOUT) :: state
        !!! EQ (A.24)
        CALL compressibility_hard_chain(state)
        CALL compressibility_dispersion(state)
        state%Z = 1 + state%Z_hc + state%Z_disp
    END SUBROUTINE
    
    
    !!! Funcoes para determinar o eta
    FUNCTION calc_pressao(state) RESULT(P)
        IMPLICIT NONE
        TYPE(PC_saft_state), INTENT(IN) :: state
        REAL(8) :: P
        
        !!! EQ (A.23)
        P = (state%Z * state%rho * state%T * k_boltz) * 1.0d30
    END FUNCTION calc_pressao   
    
    FUNCTION residuo_pressao(eta) RESULT(res)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: eta
        REAL(8) :: res
        TYPE(PC_saft_state) :: state_local
        REAL(8) :: P_calc
        !!! Preenche o estado local
        state_local%eta = eta
        state_local%T = T_global
        !!! inicita
        CALL init_mixture(state_local, ncomp_global, x_global, m_global, sigma_global, epsilo_global)
        CALL update_state(state_local)
        !!! Calcula a compressibilidade
        CALL compressibility(state_local)
        CALL helmholtz_residual(state_local)
        
        P_calc = calc_pressao(state_local)
        res = (P_calc - P_global)**2
        CALL destroy_state(state_local)
    END FUNCTION residuo_pressao
    
    !!! A ROTINA ABAIXO FOI USADA APENAS PARA DETERMINAR O VALOR DE Z E A_RES, SABENDO A RESPOSTA.
     SUBROUTINE teste_Z(eta) 
        IMPLICIT NONE 
        TYPE(PC_saft_state) :: state_local
        REAL(8), INTENT(IN) :: eta
        CALL init_mixture(state_local, ncomp_global, x_global, m_global, sigma_global, epsilo_global)
        state_local%eta = eta
        state_local%T = T_global
        
        !!!
        CALL calc_diameter_T(state_local)
        CALL calc_rho(state_local)
        CALL calc_mean_m(state_local)
        CALL calc_combining_rules(state_local)
        CALL calc_zeta(state_local)
        CALL calc_hard_sphere(state_local)
        CALL calc_ab_mean_m(state_local)
        CALL calc_pertubation_integral(state_local)
        CALL calc_prime_pertubation_integral(state_local)
        CALL calc_C_12(state_local)
        CALL calc_abbr_mes(state_local)
        CALL calc_grad_g_rho(state_local)
        
        !!!
        CALL compressibility(state_local)
        CALL helmholtz_residual(state_local)
        print*, "Rho: ", state_local%rho
        print*, "Z: ", state_local%Z
        print*, "a_res: ", state_local%a_res
        
        CALL destroy_state(state_local)
     END SUBROUTINE teste_Z
     
end module PC_saft_module