module common
    use precision
    use bc
    implicit none

    contains

    subroutine lin_sol(b,x,x0,a,c,bndcnd)
        procedure(set_bnd) :: bndcnd
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: x
        real(sp), intent(in), dimension(0:,0:) :: x0
        real(sp), intent(in) :: a, c
        
        integer :: L, M, k, i, j
   !     real(sp) :: x_av, x_av0
        
        L = size(x,1)-2
        M = size(x,2)-2
        
        do k=1,10
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
    end subroutine lin_sol

    subroutine particle_tracer(x,y,i,j,u,v)
        integer, intent(in) :: i, j
        real(sp), intent(in), dimension(0:,0:) :: u, v
        real(sp), intent(inout) :: x, y
        x = i - dt*u(i,j)/h
        y = j - dt*v(i,j)/h
    end subroutine particle_tracer

    function interpolate(d0,x,y)
        real(sp), intent(in), dimension(0:,0:) :: d0
        real(sp), intent(in) :: x, y
        real(sp) :: interpolate, s0, t0, s1, t1
        integer :: i0, j0, i1, j1
        !print*, x,y
        i0 = int(x)
        i1 = i0 + 1
        j0 = int(y)
        j1 = j0 + 1

        s1 = x - i0
        s0 = 1._sp - s1
        t1 = y - j0
        t0 = 1._sp - t1
        !print*, i0,i1,j0,j1

        interpolate = s0*(t0*d0(i0,j0)+t1*d0(i0,j1))+s1*(t0*d0(i1,j0)+t1*d0(i1,j1))
    end function interpolate

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

    subroutine c8_normal(z)
        real(sp) :: u1, u2
        complex(sp), intent(inout) :: z
        call random_number(u1)
        call random_number(u2)
        z%re = sqrt(-2._sp * log(u1)) * cos(2._sp * pi * u2)
        z%im = sqrt(-2._sp * log(u1)) * sin(2._sp * pi * u2)
    end subroutine c8_normal

    subroutine forcing(ut,vt,fu,fv,sigma,TL)
        complex(sp), intent(inout), dimension(0:,0:) :: ut, vt, fu, fv
        real(sp), intent(in) :: sigma, TL

        complex(sp) :: ww
        real(sp) :: k, k1, k2
        integer :: L, M, i, j, seed
        
        L = size(fu,1)
        M = size(fu,2)
        do j=0,M/2
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k==0._sp .or. k > KF2) cycle
                call ou_eulero(TL,sigma,fu(i,j))
                call ou_eulero(TL,sigma,fv(i,j))
                ww = k1*fu(i,j)+k2*fv(i,j)
                ut(i,j) = ut(i,j) + dt*( fu(i,j) - ww*k1/k)
                vt(i,j) = vt(i,j) + dt*( fv(i,j) - ww*k2/k)
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k==0._sp .or. k > KF2) cycle
                ww = k1*ut(i,j)+k2*vt(i,j)
                call ou_eulero(1._sp/TL,sigma,fu(i,j))
                call ou_eulero(1._sp/TL,sigma,fv(i,j))
                ww = k1*fu(i,j)+k2*fv(i,j)
                ut(i,j) = ut(i,j) + dt*( fu(i,j) - ww*k1/k)
                vt(i,j) = vt(i,j) + dt*( fv(i,j) - ww*k2/k)
            end do
        end do
        
    end subroutine forcing
    
    subroutine ou_eulero(theta,sigma,z)
        real(sp), intent(in) :: theta, sigma
        complex(sp), intent(inout) :: z

        complex(sp) :: dw

        call c8_normal(dw)
        dw = dw * sqrt ( dt )
        z = z - dt*theta*z + sigma*dw
    
    end subroutine ou_eulero

    subroutine count_forced_modes(Nf,L,M)
        integer, intent(in) :: L,M
        integer, intent(inout) :: Nf
        
        real(sp) :: k1, k2, k
        integer :: LL, MM, i, j
    
        LL = (L+2)/2
        MM = M + 1
        
        Nf = 0
        do j=0,M/2
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k==0._sp .or. k > KF2) cycle
                Nf = Nf + 1
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k==0._sp .or. k > KF2) cycle
                Nf = Nf + 1
            end do
        end do
    end subroutine count_forced_modes

    subroutine e_spect(ut,vt,ek)
        complex(sp), intent(in), dimension(0:,0:) :: ut, vt
        complex(sp), intent(inout), dimension(0:) :: ek
        
        integer :: i, j, L, M, kk
        real(sp) :: k1, k2, k
               
        L = size(ut,1)
        M = size(ut,2)
        ek = complex(0.0,0.0)
        
        do j=0,M/2
            do i=0,L-1
                kk = nint(sqrt(real(i**2+j**2,sp)))
                ek(kk) = ek(kk) + ut(i,j)*conjg(ut(i,j)) + vt(i,j)*conjg(vt(i,j))
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                kk = nint(sqrt(real(i**2+(M-j)**2,sp)))
                ek(kk) = ek(kk) + ut(i,j)*conjg(ut(i,j)) + vt(i,j)*conjg(vt(i,j))
            end do
        end do
    
        !L = size(ut,1)
        !M = size(ut,2)
        !ek = complex(0.0,0.0)
        !
        !! i=0, j=0
        !kk = nint(sqrt(real(0**2+0**2,sp)))
        !ek(kk) = ek(kk) + ut(0,0)*ut(0,0) + vt(0,0)*vt(0,0)
        !! j=0, i
        !do i=1,L-1
        !    kk = nint(sqrt(real(i**2+0**2,sp)))
        !    ek(kk) = ek(kk) + ut(i,0)*ut(L-i,0) + vt(i,0)*vt(L-i,0)
        !end do

        !do j=1,M/2
        !        kk = nint(sqrt(real(0**2+j**2,sp)))
        !        ek(kk) = ek(kk) + ut(0,j)*ut(0,M-j) + vt(0,j)*vt(0,M-j)
        !    do i=1,L-1
        !        kk = nint(sqrt(real(i**2+j**2,sp)))
        !        ek(kk) = ek(kk) + ut(i,j)*ut(L-i,M-j) + vt(i,j)*vt(L-i,M-j)
        !    end do
        !end do
        !do j=M/2+1,M-1
        !    !i=0, j
        !    kk = nint(sqrt(real(0**2+j**2,sp)))
        !    ek(kk) = ek(kk) + ut(0,j)*ut(0,M-j) + vt(0,j)*vt(0,M-j)
        !    do i=1,L-1
        !        kk = nint(sqrt(real(i**2+(M-j)**2,sp)))
        !        ek(kk) = ek(kk) + ut(i,j)*ut(L-i,M-j) + vt(i,j)*vt(L-i,M-j)
        !    end do
        !end do
    end subroutine e_spect



end module
