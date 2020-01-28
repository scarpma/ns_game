module bc
    use precision
    implicit none

    !interface 
    !    subroutine set_bnd(b,x)
    !        integer, intent(in) :: b
    !        real(sp), intent(inout), dimension(0:,0:) :: x
    !    end subroutine
    !end interface

    contains

    subroutine set_bnd_wt(b,x)
        real(sp), intent(inout), dimension(0:,0:) :: x
        integer, intent(in) :: b
        integer :: L, M
        L = size(x,1) - 2
        M = size(x,2) - 2
        if (b == 0) then ! densità
            x(0,1:M)   = x(1,1:M)
            x(L+1,1:M) = x(L,1:M)
            x(1:L,0)   = x(1:L,1)
            x(1:L,M+1) = x(1:L,M)
            !x(0,N/2-larg:N/2+larg) = 10._sp*dt
            x(0,0)     = (x(0,1)+x(1,0))/2._sp
            x(0,M+1)   = (x(1,M+1)+x(0,M))/2._sp
            x(L+1,0)   = (x(L+1,1)+x(L,0))/2._sp
            x(L+1,M+1) = (x(L,M+1)+x(L+1,M))/2._sp
            
        end if
        if (b == 1) then ! u
            x(0,1:M)   = 0._sp!u_in           ! parete in
            x(L+1,1:M) = 0._sp!u_in           ! parete out
            x(1:L,0)   = 0._sp!x(1:N,1)       ! parete sotto
            x(1:L,M+1) = 0._sp!x(1:N,N)       ! parete sopra
            x(0,0)     = (x(0,1)+x(1,0))/2._sp
            x(0,M+1)   = (x(1,M+1)+x(0,M))/2._sp
            x(L+1,0)   = (x(L+1,1)+x(L,0))/2._sp
            x(L+1,M+1) = (x(L,M+1)+x(L+1,M))/2._sp
        end if
        if (b == 2) then ! v
            x(0,1:M)   = x(1,1:M)       ! parete in
            x(L+1,1:M) = x(L,1:M)       ! parete out
            x(1:L,0)   = -x(1:L,1)      ! parete sotto
            x(1:L,M+1) = -x(1:L,M)      ! parete sopra
            x(0,0)     = (x(0,1)+x(1,0))/2._sp
            x(0,M+1)   = (x(1,M+1)+x(0,M))/2._sp
            x(L+1,0)   = (x(L+1,1)+x(L,0))/2._sp
            x(L+1,M+1) = (x(L,M+1)+x(L+1,M))/2._sp
        end if
    end subroutine set_bnd_wt

    subroutine set_bnd_box(b,x)
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: x

        integer :: i, L, M, j
        integer :: i0, j0, r

        L = size(x,1) - 2
        M = size(x,2) - 2
        
        select case(b)
        case(0)
          ! NEUMANN B.C.
          x(0,1:M) = 4.0_sp/3.0_sp*x(1,1:M)-1.0_sp/3.0_sp*x(2,1:M)
          x(L+1,1:M) = 4.0_sp/3.0_sp*x(L,1:M)-1.0_sp/3.0_sp*x(L-1,1:M)
          x(1:L,0) =  4.0_sp/3.0_sp*x(1:L,1)-1.0_sp/3.0_sp*x(1:L,2)
          x(1:L,M+1) = 4.0_sp/3.0_sp*x(1:L,M)-1.0_sp/3.0_sp*x(1:L,M-1)
          ! fix pressure on one point (does not seem to work)
          !x(L/2,M/2) = 0.0_sp
          !x(:,:) = x(:,:) - x(2,3) !Pressure undefined to a constant. Fix at least 1 value of pressure
        case(1)
          ! DIRICHLET B.C.
          x(0,1:M) = 0._sp!-x(1,1:M)
          x(L+1,1:M) = 0._sp!-x(L,1:M)
          x(1:L,0) = 0._sp!-x(1:L,1)
          x(1:L,M+1) = 1.0_sp!0._sp!-x(1:L,M)
        case(2)
          ! DIRICHLET B.C.
          x(0,1:M) = 0._sp!-x(1,1:M)
          x(L+1,1:M) = 0._sp!-x(L,1:M)
          x(1:L,0) = 0._sp!-x(1:L,1)
          x(1:L,M+1) = 0._sp!-x(1:L,M)
        end select

        x(0,0) = 0.5_sp*(x(1,0) + x(0,1))
        x(0,M+1) = 0.5_sp*(x(1,M+1) + x(0,M+1))
        x(L+1,0) = 0.5_sp*(x(L,0) + x(L+1,1))
        x(L+1,M+1) = 0.5_sp*(x(L,M+1) + x(L+1,M))

    end subroutine set_bnd_box

    subroutine set_bnd_per(b,x)
        integer, intent(in) :: b
        real(sp), intent(inout), dimension(0:,0:) :: x

        integer :: L, M

        L = size(x,1) - 1
        M = size(x,2) - 1
        
        select case(b)
            case (1:2)
                x(0:L,0) = x(0:L,M)   !n-s
                x(0,0:M) = x(L,0:M)     !w-e
            case(0)
                print *, "case 0 not defined"
        end select 
        
    end subroutine set_bnd_per

    !subroutine set_bnd_box(b,x)
    !    real(sp), intent(inout), dimension(0:,0:) :: x
    !    integer, intent(in) :: b
    !    integer :: i, L, M, j
    !    integer :: i0, j0, r
    !    L = size(x,1) - 2
    !    M = size(x,2) - 2
    !    ! CONDIZIONI SU LATI ORIZZONTALI (||u)
    !    if(b==1) then
    !        x(0,1:M) = 0._sp!-x(1,1:M)
    !        x(L+1,1:M) = 0._sp!-x(L,1:M)
    !    else
    !        x(0,1:M) = 0._sp!x(1,1:M)
    !        x(L+1,1:M) = 0._sp!x(L,1:M)
    !    endif

    !    ! CONDIZIONI SU LATI VERTICALI (||v)
    !    if(b==2) then
    !        x(1:L,0) = 0._sp!-x(1:L,1)
    !        x(1:L,M+1) = 0._sp!-x(1:L,M)
    !    else
    !        x(1:L,0) = 0._sp!x(1:L,1)
    !        x(1:L,M+1) = 0._sp!x(1:L,M)
    !    endif

    !    ! SLIP SLIT
    !    if(b==1) then
    !        x(1:L,M+1) = u_in
    !    endif
    !    
    !    !if (b==2) then
    !    !    i0=L/2+1
    !    !    j0=4.0_sp*M/5
    !    !    r=10
    !    !    do j=1,M
    !    !        do i=1,L
    !    !            if((i-i0)**2+(j-j0)**2 <= r**2) then
    !    !                x(i,j)= u_in
    !    !            end if
    !    !        end do
    !    !    end do
    !    !end if

    !    ! vertici
    !    x(0,0) = 0.5*(x(1,0) + x(0,1))
    !    x(0,M+1) = 0.5*(x(1,M+1) + x(0,M+1))
    !    x(L+1,0) = 0.5*(x(L,0) + x(L+1,1))
    !    x(L+1,M+1) = 0.5*(x(L,M+1) + x(L+1,M))
    !end subroutine set_bnd_box

    subroutine init_sources(d,u,v)
        real(sp), dimension(0:,0:), intent(out) :: d, u, v
        integer :: i, j, L, M
        integer :: i0, j0, r

        ! pacchetto circolare
        L = size(d,1)-2
        M = size(d,2)-2
        i0=3*L/5
        j0=M/2
        r=M/10
        !print*, "sources at:", i0, j0, r
        do j=1,M
            do i=1,L
                if((i-i0)**2+(j-j0)**2 <= r**2) then
                    u(i,j) = 1.0_sp
                    v(i,j) = 1.0_sp
                end if
            end do
        end do
    end subroutine init_sources
    
    subroutine bnd_cerchio(i0,j0,r,u,v)
        integer, intent(in) :: i0, j0, r ! centro e raggio cerchio  
        real(sp), intent(inout), dimension(0:,0:) :: u, v
        real(sp) :: xn, a, b, nx, ny, proj
        integer :: i, j, L, M
        L = size(u,1) - 2
        M = size(u,2) - 2
        do j=1,M
            do i=1,L
                if ((i-i0)**2+(j-j0)**2 < r**2) then
                    u(i,j) = 0._sp; v(i,j) = 0._sp
                else if ((i-i0)**2+(j-j0)**2 == r**2) then
                    ! componenti versore normale al cerchio
                    nx = real(i-i0,sp)/real(r,sp) 
                    ny = real(j-j0,sp)/real(r,sp)
                    ! prodotto scalare
                    proj = nx*u(i,j) + ny*v(i,j)
                    u(i,j) = u(i,j) - proj*nx
                    v(i,j) = v(i,j) - proj*ny
                end if
            end do
        end do
    end subroutine bnd_cerchio

end module bc
