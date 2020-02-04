module bc
    use precision
    implicit none

    contains

    subroutine set_all_bnd(u,v,u0,v0)
        real(sp), dimension(0:,0:), intent(inout) :: u, v, u0, v0
        call set_bnd_per(u)
        call set_bnd_per(v)
        call set_bnd_per(u0)
        call set_bnd_per(v0)
    end subroutine set_all_bnd

    subroutine set_bnd_per(x)
        real(sp), intent(inout), dimension(0:,0:) :: x

        integer :: L, M

        L = size(x,1) - 1
        M = size(x,2) - 1
        
        x(0:L,M) = x(0:L,0)   !n-s
        x(L,0:M) = x(0,0:M)   !w-e
        
    end subroutine set_bnd_per

    subroutine init_sources(u,v)
        real(sp), dimension(0:,0:), intent(inout) :: u, v
        integer :: i, j, L, M
        integer :: i0, j0, r

        ! pacchetto circolare
        L = size(u,1)-1
        M = size(v,2)-1
        i0=3*L/5
        j0=M/2
        r=M/10
        do j=0,M
            do i=0,L
                if((i-i0)**2+(j-j0)**2 <= r**2) then
                    u(i,j) = 1.0_sp
                    v(i,j) = 1.0_sp
                end if
            end do
        end do
    end subroutine init_sources
    
end module bc
