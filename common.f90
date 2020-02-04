module common
    use precision
    use bc
    implicit none

    contains

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
        i0 = int(x)
        i1 = i0 + 1
        j0 = int(y)
        j1 = j0 + 1

        s1 = x - i0
        s0 = 1._sp - s1
        t1 = y - j0
        t0 = 1._sp - t1

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
        real(sp) :: k, k1, k2, norm
        integer :: L, M, i, j, seed
        
        L = size(fu,1)
        M = size(fu,2)
        norm = sqrt(real(((L-1)*2)*M,sp))
        do j=0,M/2
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k < KI2 .or. k > KF2) cycle
                call ou_eulero(1._sp/TL,sigma,fu(i,j))
                call ou_eulero(1._sp/TL,sigma,fv(i,j))
                ww = complex(k1,0._sp)*fu(i,j)+complex(k2,0._sp)*fv(i,j)
                !write(*,'(i10,4(ES12.4))') M*j+i,fu(i,j)%re, fv(i,j)%re, fu(i,j)%im, fv(i,j)%im
                !write(234,'(i10,4(ES12.4))') M*j+i,fu(i,j)%re, fv(i,j)%re, fu(i,j)%im, fv(i,j)%im
                ut(i,j) = ut(i,j) + norm*dt*( fu(i,j) - ww*k1/k)
                vt(i,j) = vt(i,j) + norm*dt*( fv(i,j) - ww*k2/k)
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k < KI2 .or. k > KF2) cycle
                call ou_eulero(1._sp/TL,sigma,fu(i,j))
                call ou_eulero(1._sp/TL,sigma,fv(i,j))
                ww = complex(k1,0._sp)*fu(i,j)+complex(k2,0._sp)*fv(i,j)
                !write(*,'(i10,4(ES12.4))') M*j+i,fu(i,j)%re, fv(i,j)%re, fu(i,j)%im, fv(i,j)%im
                !write(234,'(i10,4(ES12.4))') M*j+i,fu(i,j)%re, fv(i,j)%re, fu(i,j)%im, fv(i,j)%im
                ut(i,j) = ut(i,j) + norm*dt*( fu(i,j) - ww*k1/k)
                vt(i,j) = vt(i,j) + norm*dt*( fv(i,j) - ww*k2/k)
            end do
        end do
        
    end subroutine forcing
    
    subroutine ou_eulero(theta,sigma,z)
        real(sp), intent(in) :: theta, sigma
        complex(sp), intent(inout) :: z

        complex(sp) :: dw

        call c8_normal(dw)
        dw = dw * sqrt ( dt )
        z = z - complex(dt,0.0_sp)*complex(theta,0._sp)*z + complex(sigma,0._sp)*dw
        !print*, z
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
                if (k < KI2 .or. k > KF2) cycle
                Nf = Nf + 1
            end do
        end do
        do j=M/2+1,M-1
            do i=0,L-1
                k1 = real(i,sp)*k0
                k2 = real(j-M,sp)*k0
                k = k1**2._sp+k2**2._sp
                if (k < KI2 .or. k > KF2) cycle
                Nf = Nf + 1
            end do
        end do
    end subroutine count_forced_modes

    subroutine create_shells(sh,ut,nn)
        implicit none
        real(sp), intent(inout) :: sh(0:)
        complex(sp), intent(inout) :: ut(0:,0:)
        integer, intent(in) :: nn
        
        integer :: i,j,k,n_modes,check,kkk,kk,L,M
        
        L = size(ut,1)
        M = size(ut,2)
        
        k=1
        kkk=1
        check = 0
        sh(0) = 0.
        n_modes = 0
        do while (1>0)
            do j=0,M/2
                do i=0,L-1
                    kk = nint(sqrt(real(i**2+j**2,sp)))
                    if (kk==k) then
                        n_modes = n_modes + 1
                        if (n_modes >= 20) then
                            check = 1
                        end if
                    end if
                end do
            end do
            do j=M/2+1,M-1
                do i=0,L-1
                    kk = nint(sqrt(real(i**2+(j-M)**2,sp)))
                    if (kk==k) then
                        n_modes = n_modes + 1
                        if (n_modes >= nn) then
                            check = 1
                        end if
                    end if
                end do
            end do
            if (check == 1) then
                sh(kkk) = k
                kkk = kkk + 1
                check = 0
                n_modes = 0
            end if
            k = k + 1
            
            if (k>kmax) exit

        end do

    end subroutine

    function tot_en(u,v)
        real(sp), intent(in), dimension(0:,0:) :: u, v
        real(sp) :: tot_en
        integer :: i, j, L, M
        
        L = size(u,1)
        M = size(u,2)
        
        tot_en = 0._sp
        do j=0,M-1
            do i=0,L-1
                tot_en = tot_en + 0.5_sp*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            end do
        end do
        
        return
        
    end function

    subroutine e_spect(ut,vt,ek)
        complex(sp), intent(in), dimension(0:,0:) :: ut, vt
        complex(sp), intent(inout), dimension(0:) :: ek
        
        real(sp) :: norm
        integer :: i, j, L, M, kk
               
        L = size(ut,1)
        M = size(ut,2)
        norm = real(((L-1)*2)*M,sp)
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
        
        ek = ek / norm
    
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
