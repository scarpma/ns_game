module solvers
    use precision
    use bc
    use common
    use, intrinsic :: iso_c_binding
    implicit none
    include "fftw3.f03"

    type(c_ptr), public :: planr2c, planc2r

    contains

    subroutine advect_per(d, d0, u, v)
        real(sp), intent(inout), dimension(0:,0:) :: d
        real(sp), intent(in), dimension(0:,0:) :: d0, u, v

        real(sp) :: x, y, LLL, MMM
        integer :: i, j, L, M, r

        r = 0
        L = size(d,1) - 1
        M = size(d,2) - 1
        
        ! B.C.
        d(0:L,M) = d(0:L,0)   !n-s
        d(L,0:M) = d(0,0:M)   !e-w

        LLL = real(L,sp)
        MMM = real(M,sp)
        do j=0,M-1 !il ciclo è fatto escludendo i bordi n-e
            do i=0,L-1
                call particle_tracer(x,y,i,j,u,v)
                !if (x < 0._sp .or. x >= LLL .or. y < 0._sp .or. y >= MMM) then
                !    print*, i,j,x,y
                !    r = 1
                !end if
                if (x < 0._sp)      x = LLL + x
                if (x >= LLL)        x = x - LLL
                if (y < 0.0_sp)     y = MMM + y
                if (y >= MMM)         y = y - MMM
                d(i,j) = interpolate(d0,x,y)
                !if (r==1) then
                !    print*, x,y
                !    r = 0
                !end if
            end do
        end do

        ! B.C.
        d(0:L,M) = d(0:L,0)   !n-s
        d(L,0:M) = d(0,0:M)   !e-w
                
    end subroutine advect_per

    subroutine diffuse_and_project(ut,vt,visc)
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
                ut(i,j) = (ut(i,j)-ww*k1/k) / (1._sp+visc*dt*k)!*exp(-k*visc*dt)
                vt(i,j) = (vt(i,j)-ww*k2/k) / (1._sp+visc*dt*k)!*exp(-k*visc*dt))
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                ww = k1*ut(i,j)+k2*vt(i,j)
                ut(i,j) = (ut(i,j)-ww*k1/k) / (1._sp+visc*dt*k)!*exp(-k*visc*dt))
                vt(i,j) = (vt(i,j)-ww*k2/k) / (1._sp+visc*dt*k)!*exp(-k*visc*dt))
            end do
        end do

    end subroutine diffuse_and_project

    subroutine vel_step(u,v,u0,v0,ut,vt,fu,fv,visc)
        real(sp), intent(inout), dimension(0:,0:) :: u, v, u0, v0
        complex(sp), intent(inout), dimension(0:,0:) :: ut, vt, fu, fv
        real(sp), intent(in) :: visc
        integer :: L, M

        L = size(u,1)
        M = size(u,2)

        call advect_per(u0,u,u,v)
        call advect_per(v0,v,u,v)

        call fftw_execute_dft_r2c(planr2c,u0,ut)
        call fftw_execute_dft_r2c(planr2c,v0,vt)

        call diffuse_and_project(ut,vt,visc)
        call forcing(ut,vt,fu,fv,sigma,TL)

        call fftw_execute_dft_c2r(planc2r,ut,u)
        call fftw_execute_dft_c2r(planc2r,vt,v)
        u = u / (L*M)
        v = v / (L*M)
    end subroutine vel_step

end module solvers
