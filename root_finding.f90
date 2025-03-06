module root_finding
    implicit none
    interface solve_equation
        module procedure newton_raphson
    end interface
contains

    subroutine swap(a, b)
        implicit none
        real(8), intent(inout) :: a, b
        real(8) :: temp
        temp = a
        a = b
        b = temp
    end subroutine swap
    !!! ABAIXO SAO AS FUNCOES PARA MINIMIZACAO !!!
    function newton_raphson(f, df, x0) result(root)
        implicit none
        real(8) :: x0, tol
        integer :: max_iter
        real(8), external :: f, df
        real(8) :: root, x, fx, dfx
        integer :: iter
        tol = 1e-6
        max_iter = 250
        x = x0
        do iter = 1, max_iter
            fx = f(x)
            dfx = df(x)
            if (ABS(fx) < tol) exit
            x = x - fx / dfx
        end do
    
        root = x
    end function newton_raphson

    FUNCTION simplex(f, x) RESULT(x_root)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x
        REAL(8), EXTERNAL :: f
        REAL(8) :: x_root
        REAL(8) :: alpha, gamma, beta, sigma, delta
        REAL(8) :: x1, x2, xr, xe, xc
        REAL(8) :: f1, f2, fr, fe, fc
        INTEGER :: max_iter, iter
    
        alpha = 1.75d0
        gamma = 1.80d0
        beta = 0.5d0
        sigma = 0.4d0
        max_iter = 250
        
        x1 = x
        x2 = x1 + 0.5d-1

        f1 = f(x1)
        f2 = f(x2)


        DO iter = 1, max_iter
            IF (f2 < f1) THEN
                CALL swap(x1, x2)
                CALL swap(f1, f2)
            END IF
    
            xr = x1 + alpha * (x1 - x2)
            fr = f(xr)
    
            IF (fr < f1) THEN
                xe = x1 + gamma * (xr - x1)
                fe = f(xe)
        
                IF (fe < fr) THEN
                    x2 = xe
                    f2 = fe
                ELSE
                    x2 = xr
                    f2 = fr
                END IF
            ELSE IF (fr < f2) THEN
                x2 = xr
                f2 = fr
            ELSE
                xc = x1 + beta * (x2 - x1)
                fc = f(xc)
        
                IF (fc < f2) THEN
                    x2 = xc
                    f2 = fc
                ELSE
                    x2 = x1 + sigma * (x2 - x1)
                    f2 = f(x2)
                END IF
            END IF
        
            IF (ABS(x2 - x1) < 1.0d-12) EXIT
    
            END DO
        IF (f1 < f2) THEN
            x_root = x1
        ELSE
            x_root = x2
        END IF
    END FUNCTION simplex
end module root_finding
