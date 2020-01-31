module solvers
    use precision
    use bc
    use common
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"

    type(c_ptr), public :: planr2c, planc2r
    real(sp), parameter :: conviter = 0.005

    contains

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
        call lin_sol(b,x,x0,a,1.0_sp+4.0_sp*a,bndcnd)  
    end subroutine diffuse
    
    subroutine advect(b, d, d0, u, v, bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: d
        real(sp), intent(in), dimension(0:,0:) :: d0, u, v

        real(sp) :: x, y, LLL, MM
        integer :: i, j, L, M

        L = size(d,1)-2
        M = size(d,2)-2
        LLL = real(L,sp)
        MM = real(M,sp)
        do j=1,M
            do i=1,L
                call particle_tracer(x,y,i,j,u,v)
                if (x < 0.5_sp)      x = 0.5_sp
                if (x > LLL + 0.5_sp) x = LLL + 0.5_sp 
                if (y < 0.5_sp)      y = 0.5_sp
                if (y > MM + 0.5_sp) y = MM + 0.5_sp
                d(i,j) = interpolate(d0,x,y)
            end do
        end do
    end subroutine advect

    subroutine advect_per(b, d, d0, u, v, bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: d
        real(sp), intent(in), dimension(0:,0:) :: d0, u, v

        real(sp) :: x, y, LLL, MM
        integer :: i, j, L, M

        L = size(d,1)
        M = size(d,2)
        !LLL = real(L,sp)
        LLL = real(L-1,sp)
        !MM = real(M,sp)
        MM = real(M-1,sp)
        do j=0,M-1
            do i=0,L-1
                call particle_tracer(x,y,i,j,u,v)
                if (x <= 0._sp)      x = LLL + x!if (x < 0.5_sp)      x = 0.5_sp
                if (x >= LLL)        x = x - LLL!if (x > LLL + 0.5_sp) x = LLL + 0.5_sp 
                if (y <= 0.0_sp)     y = MM + y!if (y < 0.5_sp)      y = 0.5_sp
                if (y >= MM)         y = y - MM!if (y > MM + 0.5_sp) y = MM + 0.5_sp
                d(i,j) = interpolate(d0,x,y)
            end do
        end do
        !call bndcnd(b,d)
    end subroutine advect_per

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
        call lin_sol(0,p,div,1.0_sp,4.0_sp,bndcnd)
        do j=1,M
            do i=1,L
                u(i,j) = u(i,j) - 0.5*(p(i+1,j)-p(i-1,j))/h
                v(i,j) = v(i,j) - 0.5*(p(i,j+1)-p(i,j-1))/h
            end do
        end do
        ! controllo se la divergenza è nulla
        do j=1,M
            do i=1,L
                div(i,j) = -0.5_sp*h*(u(i+1,j)-u(i-1,j)+v(i,j+1)-v(i,j-1))
            end do
        end do
        divmax = maxval(abs(div(:,:)))
        if (divmax > 0.01_sp) then
            k = k + 1
            goto 1020
        end if
        !print*, k*40
    end subroutine project

    subroutine diffuse_and_project(ut,vt,visc,bndcnd)
        procedure(set_bnd) :: bndcnd
        complex(sp), intent(inout), dimension(0:,0:) :: ut, vt
        real(sp), intent(in) :: visc
        
        real(sp) :: k, k1, k2
        complex(sp) :: ww
        integer :: i, j, L, M
        
        L = size(ut,1)
        M = size(ut,2)
        ! DIFFUSE AND PROJECT
        do j=0,M/2
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k==0._sp) cycle
                ww = k1*ut(i,j)+k2*vt(i,j)
                ut(i,j) = (ut(i,j)-ww*k1/k)*exp(-k*visc*dt)! / (1._sp+visc*dt*k)
                vt(i,j) = (vt(i,j)-ww*k2/k)*exp(-k*visc*dt)! / (1._sp+visc*dt*k)
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                ww = k1*ut(i,j)+k2*vt(i,j)
                ut(i,j) = (ut(i,j)-ww*k1/k)*exp(-k*visc*dt)! / (1._sp+visc*dt*k)
                vt(i,j) = (vt(i,j)-ww*k2/k)*exp(-k*visc*dt)! / (1._sp+visc*dt*k)
            end do
        end do

    end subroutine diffuse_and_project

    subroutine density_step(x, x0, u, v, diff, bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), intent(inout), dimension(0:,0:) :: x, x0, u, v
        real(sp), intent(in) :: diff
        !call add_source(x,x0)
        !call diffuse(0,x0,x,diff,bndcnd)
        !call advect(0,x,x0,u,v,bndcnd)
    end subroutine density_step

    subroutine vel_step(u,v,u0,v0,ut,vt,fu,fv,p,div,visc,bndcnd)
        procedure(set_bnd) :: bndcnd
        real(sp), intent(inout), dimension(0:,0:) :: u, v, u0, v0, p, div
        complex(sp), intent(inout), dimension(0:,0:) :: ut, vt, fu, fv
        real(sp), intent(in) :: visc
        integer :: L, M

        L = size(u,1)
        M = size(u,2)

        call advect_per(1,u0,u,u,v,bndcnd)
        call advect_per(2,v0,v,u,v,bndcnd)
        call fftw_execute_dft_r2c(planr2c,u0,ut)
        call fftw_execute_dft_r2c(planr2c,v0,vt)
        call diffuse_and_project(ut,vt,visc,bndcnd)

        call forcing(ut,vt,fu,fv,sigma,TL)

        call fftw_execute_dft_c2r(planc2r,ut,u)
        u = u / (L*M)
        call fftw_execute_dft_c2r(planc2r,vt,v)
        v = v / (L*M)
    end subroutine vel_step
    

end module solvers
