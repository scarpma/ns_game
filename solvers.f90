module solvers
    use precision
    use bc
    implicit none

    real(sp), parameter :: conviter = 0.005

    interface 
        subroutine set_bnd(b,x)
            use precision    
            integer, intent(in) :: b
            real(sp), intent(inout), dimension(0:,0:) :: x
        end subroutine
    end interface

    contains

    subroutine set_all_bnd(x,u,v,x0,u0,v0,bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), dimension(0:,0:), intent(inout) :: x, u, v, x0, u0, v0
        call bndcnd(0,x)
        call bndcnd(1,u)
        call bndcnd(2,v)
        call bndcnd(0,x0)
        call bndcnd(1,u0)
        call bndcnd(2,v0)
    end subroutine set_all_bnd

    subroutine add_source(x,s)
        real(sp), intent(inout), dimension(0:,0:) :: x
        real(sp), intent(in), dimension(0:,0:) :: s
        integer :: i, j, L, M
        L = size(x,1)-2
        M = size(x,2)-2
        do i=1,L
          do j=1,M
            x(i,j) = x(i,j) + dt*s(i,j)
          end do
        end do
    end subroutine add_source

    subroutine lin_sol_gaussseidel(b,x,x0,a,c,bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: x
        real(sp), intent(in), dimension(0:,0:) :: x0
        real(sp), intent(in) :: a, c
        
        integer :: L, M, k, i, j
   !     real(sp) :: x_av, x_av0
        
        L = size(x,1)-2
        M = size(x,2)-2
        
        do k=1,20
   !         x_av = 0._sp
            do j=1,M
                do i=1,L
                    x(i,j) = ( x0(i,j) + a*(x(i-1,j)+x(i+1,j)+x(i,j-1)+x(i,j+1)) )/c
   !                 x_av = x_av + x(i,j)
                end do
            end do
        call bndcnd(b,x)
   !         if (abs(x_av-x_av0)<conviter) exit
   !         x_av0 = x_av
        end do
        !write(*,'(a1,i1,a1,4(ES15.6),i4)') 'u',b,' ',minval(x),maxval(x),x_av/(L*M),abs(x_av-x_av),k
    end subroutine lin_sol_gaussseidel

    subroutine diffuse(b, x, x0, diff,bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: x
        real(sp), intent(in), dimension(0:,0:) :: x0
        real(sp), intent(in) :: diff

        real(sp) :: a, x_av, x_av0
        integer :: L, M

        L = size(x,1)-2
        M = size(x,2)-2
        a = dt*diff/(h**2._sp)
        call lin_sol_gaussseidel(b,x,x0,a,1.0_sp+4.0_sp*a,bndcnd)  
    end subroutine diffuse
    
    subroutine advect(b, d, d0, u, v, bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: d
        real(sp), intent(in), dimension(0:,0:) :: d0, u, v

        real(sp) :: x, y, s0, t0, s1, t1, dt0, LL, MM
        integer :: i, j, i0, j0, i1, j1, L, M

        L = size(d,1)-2
        M = size(d,2)-2
        i0 = 0; j0 = 0; i1 = 0; j1 = 0
        LL = real(L,sp)
        MM = real(M,sp)
        dt0 = dt/h
        do j=1,M
            do i=1,L
                x = i - dt0*u(i,j)
                y = j - dt0*v(i,j)
                if (x < 0.5_sp)      x = 0.5_sp
                if (x > LL + 0.5_sp) x = LL + 0.5_sp 
                i0 = int(x)
                i1 = i0 + 1
                if (y < 0.5_sp)      y = 0.5_sp
                if (y > MM + 0.5_sp) y = MM + 0.5_sp
                j0 = int(y)
                j1 = j0 + 1
                s1 = x - i0
                s0 = 1._sp - s1
                t1 = y - j0
                t0 = 1._sp - t1
                d(i,j) = s0*(t0*d0(i0,j0)+t1*d0(i0,j1))+s1*(t0*d0(i1,j0)+t1*d0(i1,j1))
            end do
        end do
        call bndcnd(b,d)
    end subroutine advect

    subroutine project(u, v, p, div, bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), intent(inout) :: u(0:,0:), v(0:,0:), p(0:,0:), div(0:,0:)
        integer :: i, j, k, L, M
        real(sp) :: divmax
        L = size(u,1)-2
        M = size(u,2)-2
        do j=1,M
            do i=1,L
                div(i,j) = -0.5_sp*h*(u(i+1,j)-u(i-1,j)+v(i,j+1)-v(i,j-1))
                p(i,j) = 0.0_sp
            end do
        end do
        ! dirichlet conditions on pressure and div(u)
        k = 1
        1020 call bndcnd(2,div); call bndcnd(0,p)
        call lin_sol_gaussseidel(0,p,div,1.0_sp,4.0_sp,bndcnd)
        do j=1,M
            do i=1,L
                u(i,j) = u(i,j) - 0.5*(p(i+1,j)-p(i-1,j))/h
                v(i,j) = v(i,j) - 0.5*(p(i,j+1)-p(i,j-1))/h
            end do
        end do
        ! controllo se la divergenza è nulla
        !do j=1,M
        !    do i=1,L
        !        div(i,j) = -0.5_sp*h*(u(i+1,j)-u(i-1,j)+v(i,j+1)-v(i,j-1))
        !    end do
        !end do
        !divmax = maxval(abs(div(:,:)))
        !if (divmax > 0.01_sp) then
        !    k = k + 1
        !    goto 1020
        !end if
        !print*, k*40
    end subroutine project
    
    subroutine density_step(x, x0, u, v, diff, bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), intent(inout), dimension(0:,0:) :: x, x0, u, v
        real(sp), intent(in) :: diff
        !call add_source(x,x0)
        call diffuse(0,x0,x,diff,bndcnd)
        call advect(0,x,x0,u,v,bndcnd)
    end subroutine density_step

    subroutine vel_step(u,v,u0,v0,p,div,visc,bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), intent(inout), dimension(0:,0:) :: u, v, u0, v0, p, div
        real(sp), intent(in) :: visc
        integer :: L, M
        L = size(u,1)-2
        M = size(u,2)-2
        !!call add_source(u,u0) assuming forces in u0 and v0
        !!call add_source(v,v0)
        !call diffuse(1,u0,u,visc,bndcnd)
        !call diffuse(2,v0,v,visc,bndcnd)
        !!call bnd_cerchio(xc,yc,rc,u0,v0)
        !call project(u0,v0,p,div,bndcnd)
        !!call bnd_cerchio(xc,yc,rc,u0,v0)
        !call advect(1,u,u0,u0,v0,bndcnd)
        !call advect(2,v,v0,u0,v0,bndcnd)
        !!call bnd_cerchio(xc,yc,rc,u,v)
        !call project(u,v,p,div,bndcnd)
        !!call bnd_cerchio(xc,yc,rc,u0,v0)

        !call add_source(u,u0) assuming forces in u0 and v0
        !call add_source(v,v0)
        call advect(1,u0,u,u,v,bndcnd)
        call advect(2,v0,v,u,v,bndcnd)
        !call bnd_cerchio(xc,yc,rc,u,v)
        call diffuse(1,u,u0,visc,bndcnd)
        call diffuse(2,v,v0,visc,bndcnd)
        !call bnd_cerchio(xc,yc,rc,u0,v0)
        call project(u,v,p,div,bndcnd)
        !call bnd_cerchio(xc,yc,rc,u0,v0)
    end subroutine vel_step
    
    function errmax(u,u1)
        real(sp), intent(in), dimension(0:,0:) :: u, u1
        real(sp) :: errmax
        errmax = maxval(abs((u(:,:) - u1(:,:))/u(:,:)))
    end function errmax
    
    subroutine check_uv_maxerr(n,u,v,u1,v1,conv,i,conv_check)
        integer, intent(in) :: n, i
        integer, intent(inout) :: conv_check
        real(sp), intent(in) :: conv
        real(sp), intent(in), dimension(0:,0:) :: u, v
        real(sp), intent(inout), dimension(0:,0:) :: u1, v1
        real(sp) :: errmax_u, errmax_v
        
        if (mod(i,n)==0) then
            errmax_u = errmax(u,u1)
            errmax_v = errmax(v,v1)
            write(*,'(" Iter ", I6, " errmax_u= ", ES7.1, " errmax_v= ", ES7.1)') i,  errmax_u, errmax_v
            if (errmax_u <= conv .and. errmax_v <= conv) then
                write(*,'("Convergenza al ", ES7.1,"% raggiunta. Arresto.")') conv*100
                conv_check = 1
            end if
            u1 = u
            v1 = v
        end if
    end subroutine check_uv_maxerr


end module solvers
