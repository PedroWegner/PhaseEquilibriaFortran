module root_finding
    implicit none
    interface solve_equation
        module procedure newton_raphson
    end interface
contains
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

end module root_finding
