! Stable Fluid Solver subroutines
! version: 0.6
! Copyleft 2016 Dipa

module common 
  use precision
  use linear_solver
  use mat_def
  implicit none
  private

  public :: init_sources_2D
  public :: dens_step_2D
  public :: vel_step_2D

  public :: init_sources_3D
  public :: dens_step_3D
  public :: vel_step_3D
  public :: vel_step_3D_cg



  contains

  ! inizializza sorgenti !

  !-----------------------------------------------------------------------------
  subroutine init_sources_2D(d,u,v,L,M)
    real(dp), dimension(:,:), intent(out) :: d, u, v
    integer, intent(in) :: L, M

    integer :: i, j
    integer :: i0, j0, r

    ! pacchetto circolare
    i0=L/2+1
    j0=M/5

    r=5

    do j=2,M+1
      do i=2,L+1
        if((i-i0)**2+(j-j0)**2 <= r**2) then
          d(i,j)= 5.0_dp
          v(i,j)= 5.0_dp
        end if
      end do
    end do

  end subroutine init_sources_2D

  !-----------------------------------------------------------------------------
  subroutine init_sources_3D(d,u,v,w,L,M,N)
    real(dp), dimension(:,:,:), intent(out) :: d, u, v, w
    integer, intent(in) :: L, M, N

    integer :: i, j, k
    integer :: i0, j0, k0, r

    ! pacchetto sferico
    i0=L/2+1
    j0=M/2+1
    k0=N/5

    r=5

    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          if((i-i0)**2+(j-j0)**2+(k-k0)**2 <= r**2) then
            d(i,j,k)= 10.0_dp
            w(i,j,k)= 10.0_dp
          end if
        end do
      end do
    end do

  end subroutine init_sources_3D

  ! termine di sorgente !

  !-----------------------------------------------------------------------------
  subroutine add_source_2D(x,s,L,M,dt)
    real(dp), dimension(:,:), intent(inout) :: x
    real(dp), dimension(:,:), intent(in) :: s
    integer, intent(in) :: L, M
    real(dp), intent(in) :: dt

    integer :: i, j

    do j=1,M+2
      do i=1,L+2
        x(i,j) = x(i,j) + dt*s(i,j)
      end do
    end do
  end subroutine add_source_2D

  !-----------------------------------------------------------------------------
  subroutine add_source_3D(x,s,L,M,N,dt)
    real(dp), dimension(:,:,:), intent(inout) :: x
    real(dp), dimension(:,:,:), intent(in) :: s
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: dt

    integer :: i, j, k

    do k=1,N+2 
      do j=1,M+2
        do i=1,L+2
          x(i,j,k) = x(i,j,k) + dt*s(i,j,k)
        end do
      end do
    end do
  end subroutine add_source_3D

  ! termine di diffusione !

  !-----------------------------------------------------------------------------
  subroutine diffuse_2D(x,x0,diff,b,L,M,h,dt)
    real(dp), dimension(:,:), intent(inout) :: x
    real(dp), dimension(:,:), intent(in) :: x0
    integer, intent(in) :: b, L, M
    real(dp), intent(in) :: diff, h, dt

    integer :: i, j, q
    real(dp) :: a

    a=dt*diff/(h**2)

    do q=1,20
      do j=2,M+1
        do i=2,L+1
          x(i,j) = (x0(i,j) + a*(x(i-1,j)+x(i+1,j) &
                   +x(i,j-1)+x(i,j+1)))/(1+4*a)
        end do
      end do
      call set_bnd_2D(L,M,b,x)
    end do
  end subroutine diffuse_2D

  !-----------------------------------------------------------------------------
  subroutine diffuse_3D(x,x0,diff,b,L,M,N,h,dt)
    real(dp), dimension(:,:,:), intent(inout) :: x
    real(dp), dimension(:,:,:), intent(in) :: x0
    integer, intent(in) :: b, L, M, N
    real(dp), intent(in) :: diff, h, dt

    integer :: i, j, k, q
    real(dp) :: a

    a=dt*diff/(h**2)

    do q=1,50
      do k=2,N+1
        do j=2,M+1
          do i=2,L+1
            x(i,j,k) = (x0(i,j,k) + a*( x(i-1,j,k)+x(i+1,j,k) &
            +x(i,j-1,k)+x(i,j+1,k)+x(i,j,k-1)+x(i,j,k+1) ))/(1+6*a)
          end do
        end do
      end do
      call set_bnd_3D(L,M,N,b,x)
    end do
  end subroutine diffuse_3D

  ! termine di avvezione !

  !-----------------------------------------------------------------------------
  subroutine advect_2D(d,d0,u,v,b,L,M,h,dt)
    real(dp), dimension(:,:), intent(inout) :: d
    real(dp), dimension(:,:), intent(in) :: d0, u, v
    integer, intent(in) :: b, L, M
    real(dp), intent(in) :: h, dt

    integer :: i, j, i0, j0, i1, j1
    real(dp) :: x, y, r0, s0, r1, s1, dt0

    dt0 = dt/h;

    do j=2,M+1
      do i=2,L+1
        x = i-dt0*u(i,j)
        y = j-dt0*v(i,j)

        if (x<1.5) then
          x=1.5;
        endif
        if (x>L+1.5) then
          x=L+1.5
        endif
        i0=int(x)
        i1=i0+1
        if (y<1.5) then
          y=1.5
        endif
        if (y>M+1.5) then
          y=M+1.5
        endif
        j0=int(y)
        j1=j0+1

        r1 = x-i0
        r0 = 1-r1
        s1 = y-j0
        s0 = 1-s1
        d(i,j) = r0*(s0*d0(i0,j0)+s1*d0(i0,j1)) &
                +r1*(s0*d0(i1,j0)+s1*d0(i1,j1))
      end do
    end do
    call set_bnd_2D(L,M,b,d);
  end subroutine advect_2D

  !-----------------------------------------------------------------------------
  subroutine advect_3D(d,d0,u,v,w,b,L,M,N,h,dt)
    real(dp), dimension(:,:,:), intent(inout) :: d
    real(dp), dimension(:,:,:), intent(in) :: d0, u, v, w
    integer, intent(in) :: b, L, M, N
    real(dp), intent(in) :: h, dt

    integer :: i, j, k, i0, j0, k0, i1, j1, k1
    real(dp) :: x, y, z, r0, s0, t0, r1, s1, t1, dt0

    dt0 = dt/h;
    
    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          x = i-dt0*u(i,j,k)
          y = j-dt0*v(i,j,k)
          z = k-dt0*w(i,j,k)

          if (x<1.5) then
            x=1.5;
          endif
          if (x>L+1.5) then
            x=L+1.5
          endif
          i0=int(x)
          i1=i0+1

          if (y<1.5) then
            y=1.5
          endif
          if (y>M+1.5) then
            y=M+1.5
          endif
          j0=int(y)
          j1=j0+1

          if (z<1.5) then
            z=1.5
          endif
          if (z>N+1.5) then
            z=N+1.5
          endif
          k0=int(z)
          k1=k0+1
  
          r1 = x-i0
          r0 = 1-r1
          s1 = y-j0
          s0 = 1-s1
          t1 = z-k0
          t0 = 1-t1

          d(i,j,k) = r0*(s0*(t0*d0(i0,j0,k0)+t1*d0(i0,j0,k1)) &
                        +s1*(t0*d0(i0,j1,k0)+t1*d0(i0,j1,k1))) &
                    +r1*(s0*(t0*d0(i1,j0,k0)+t1*d0(i1,j0,k1)) &
                        +s1*(t0*d0(i1,j1,k0)+t1*d0(i1,j1,k1)))
        end do
      end do
    end do
    call set_bnd_3D(L,M,N,b,d);
  end subroutine advect_3D

  ! pressure Poisson equation !

  !-----------------------------------------------------------------------------
  subroutine project_2D(u,v,p,div,L,M,h)
    real(dp), dimension(:,:), intent(inout) :: u, v, p, div
    integer, intent(in) :: L, M
    real(dp), intent(in) :: h

    integer :: i, j, q

    ! -h*div(u)
    do j=2,M+1
      do i=2,L+1
        div(i,j) = -0.5*h*(u(i+1,j)-u(i-1,j)+v(i,j+1)-v(i,j-1))      
        p(i,j) = 0
      end do
    end do
 
    call set_bnd_2D(L,M,0,div)
    call set_bnd_2D(L,M,0,p)

    ! ppe
    do q=1,20
      do j=2,M+1
        do i=2,L+1
          p(i,j) = (div(i,j)+p(i-1,j)+p(i+1,j) &
                            +p(i,j-1)+p(i,j+1))/4
        end do
      end do
      call set_bnd_2D(L,M,0,p)
    end do
    
    ! u-grad(p)
    do j=2,M+1
      do i=2,L+1
        u(i,j) = u(i,j)-0.5*(p(i+1,j)-p(i-1,j))/h
        v(i,j) = v(i,j)-0.5*(p(i,j+1)-p(i,j-1))/h
      end do
    end do
    
    call set_bnd_2D(L,M,1,u)
    call set_bnd_2D(L,M,2,v)

  end subroutine project_2D

  !-----------------------------------------------------------------------------
  subroutine project_3D(u,v,w,p,div,L,M,N,h)
    real(dp), dimension(:,:,:), intent(inout) :: u, v, w, p, div
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h

    integer :: i, j, k, q

    ! -h*div(u)
    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          div(i,j,k) = -0.5*h*(u(i+1,j,k)-u(i-1,j,k)+v(i,j+1,k) &
                               -v(i,j-1,k)+w(i,j,k+1)-w(i,j,k-1))    
          p(i,j,k) = 0
        end do
      end do
    end do
 
    call set_bnd_3D(L,M,N,0,div)
    call set_bnd_3D(L,M,N,0,p)
 
    ! ppe
    do q=1,200
      do k=2,N+1
        do j=2,M+1
          do i=2,L+1
            p(i,j,k) = (div(i,j,k)+p(i-1,j,k)+p(i+1,j,k)+p(i,j-1,k) &
                               +p(i,j+1,k)+p(i,j,k-1)+p(i,j,k+1))/6
          end do
        end do
        call set_bnd_3D(L,M,N,0,p)
      end do
    end do

    ! u-grad(p)
    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          u(i,j,k) = u(i,j,k)-0.5*(p(i+1,j,k)-p(i-1,j,k))/h
          v(i,j,k) = v(i,j,k)-0.5*(p(i,j+1,k)-p(i,j-1,k))/h
          w(i,j,k) = w(i,j,k)-0.5*(p(i,j,k+1)-p(i,j,k-1))/h
        end do
      end do
    end do

    call set_bnd_3D(L,M,N,1,u)
    call set_bnd_3D(L,M,N,2,v)
    call set_bnd_3D(L,M,N,3,w)

  end subroutine project_3D

  !-----------------------------------------------------------------------------
  subroutine project_3D_cg(u,v,w,A,L,M,N,h)
    real(dp), dimension(:,:,:), intent(inout) :: u, v, w
    type(r_CSR), intent(in) :: A
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h

    real(dp), dimension(L+2,M+2,N+2) :: tmp
    real(dp), dimension(L*M*N) :: b, p
    integer :: i, j, k
    
    ! -V*div(u)
    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          tmp(i,j,k) = -0.5*h*(u(i+1,j,k) &
          -u(i-1,j,k)+v(i,j+1,k)-v(i,j-1,k)+w(i,j,k+1)-w(i,j,k-1))
        end do
      end do
    end do

    do k=1,N
      do j=1,M
        do i=1,L
          b(i+(j-1)*L+(k-1)*L*M) = tmp(i+1,j+1,k+1)
        end do
      end do
    end do

    p=0

    ! ppe
    call conjugate_gradient(A, b, p, real(1.0e-5,8))

    do k=1,N
      do j=1,M
        do i=1,L
          tmp(i+1,j+1,k+1) = p(i+(j-1)*L+(k-1)*L*M)
        end do
      end do
    end do

    call set_bnd_3D(L,M,N,0,tmp)

    ! u-grad(p)    
    do k=2,N+1
      do j=2,M+1
        do i=2,L+1
          u(i,j,k) = u(i,j,k)-0.5*(tmp(i+1,j,k)-tmp(i-1,j,k))/h
          v(i,j,k) = v(i,j,k)-0.5*(tmp(i,j+1,k)-tmp(i,j-1,k))/h
          w(i,j,k) = w(i,j,k)-0.5*(tmp(i,j,k+1)-tmp(i,j,k-1))/h
        end do
      end do
    end do

    call set_bnd_3D(L,M,N,1,u)
    call set_bnd_3D(L,M,N,2,v)
    call set_bnd_3D(L,M,N,3,w)

  end subroutine project_3D_cg

  ! condizioni al contorno !

  !-----------------------------------------------------------------------------
  subroutine set_bnd_2D(L,M,b,x)
    real(dp), dimension(:,:), intent(inout) :: x
    integer, intent(in) :: b, L, M

    integer :: i, j

    ! lati
    do j=2,M+1
      if(b==1) then
        x(1,j) = -x(2,j)
        x(L+2,j) = -x(L+1,j)
      else
        x(1,j) = x(2,j)
        x(L+2,j) = x(L+1,j)
      endif
    end do

    do i=2,L+1
      if(b==2) then
        x(i,1) = -x(i,2)
        x(i,M+2) = -x(i,M+1)
      else
        x(i,1) = x(i,2)
        x(i,M+2) = x(i,M+1)
      endif
    end do

    ! slip slit
    
    if(b==1) then
      do i=2,L+1
        x(i,M+2) = 10.0_dp
      end do
    endif

    ! vertici
    x(1,1) = 0.5*(x(2,1) + x(1,2))
    x(1,M+2) = 0.5*(x(2,M+2) + x(1,M+1))
    x(L+2,1) = 0.5*(x(L+1,1) + x(L+2,2))
    x(L+2,M+2) = 0.5*(x(L+1,M+2) + x(L+2,M+1))

  end subroutine set_bnd_2D

  !-----------------------------------------------------------------------------
  subroutine set_bnd_3D(L,M,N,b,x)
    real(dp), dimension(:,:,:), intent(inout) :: x
    integer, intent(in) :: b, L, M, N

    integer :: i, j, k

    ! facce
    do k=2,N+1  
      do j=2,M+1
        if(b==1) then
          x(1,j,k) = -x(2,j,k)
          x(L+2,j,k) = -x(L+1,j,k)
        else
          x(1,j,k) = x(2,j,k)
          x(L+2,j,k) = x(L+1,j,k)
        endif
      end do
    end do

    do k=2,N+1   
      do i=2,L+1      
        if(b==2) then
          x(i,1,k) = -x(i,2,k)
          x(i,M+2,k) = -x(i,M+1,k)
        else
          x(i,1,k) = x(i,2,k)
          x(i,M+2,k) = x(i,M+1,k)
        endif
      end do
    end do

    do j=2,M+1   
      do i=2,L+1      
        if(b==3) then
          x(i,j,1) = -x(i,j,2)
          x(i,j,N+2) = -x(i,j,N+1)
        else
          x(i,j,1) = x(i,j,2)
          x(i,j,N+2) = x(i,j,N+1)
        endif
      end do
    end do

    ! slip slit
    !
    if(b==1) then
      do j=2,M+1
        do i=2,L+1
          x(i,j,N+2) = 10.0_dp
        end do
      end do
    endif

    ! spigoli
    do i=2,L+1
      x(i,1,1) = 0.5*(x(i,2,1) + x(i,1,2))
      x(i,1,N+2) = 0.5*(x(i,2,N+2) + x(i,1,N+1))
      x(i,M+2,1) = 0.5*(x(i,M+1,1) + x(i,M+2,2))
      x(i,M+2,N+2) = 0.5*(x(i,M+1,N+2) + x(i,M+2,N+1))
    end do

    do j=2,M+1
      x(1,j,1) = 0.5*(x(2,j,1) + x(1,j,2))
      x(1,j,N+2) = 0.5*(x(2,j,N+2) + x(1,j,N+1))
      x(L+2,j,1) = 0.5*(x(L+1,j,1) + x(L+2,j,2))
      x(L+2,j,N+2) = 0.5*(x(L+1,j,N+2) + x(L+2,j,N+1))
    end do

    do k=2,N+1
      x(1,1,k) = 0.5*(x(2,1,k) + x(1,2,k))
      x(1,M+2,k) = 0.5*(x(2,M+2,k) + x(1,M+1,k))
      x(L+2,1,k) = 0.5*(x(L+1,1,k) + x(L+2,2,k))
      x(L+2,M+2,k) = 0.5*(x(L+1,M+2,k) + x(L+2,M+1,k))
    end do

    ! vertici
    x(1,1,1) = (x(2,1,1) + x(1,2,1) + x(1,1,2))/3.0_dp
    x(1,M+2,1) = (x(2,M+2,1) + x(1,M+1,1) + x(1,M+2,2))/3.0_dp
    x(1,1,N+2) = (x(2,1,N+2) + x(1,2,N+2) + x(1,1,N+1))/3.0_dp
    x(1,M+2,N+2) = (x(2,M+2,N+2) + x(1,M+1,N+2) &
                                 + x(1,M+2,N+1))/3.0_dp
    x(L+2,1,1) = (x(L+1,1,1) + x(L+2,2,1) + x(1,1,2))/3.0_dp
    x(L+2,M+2,1) = (x(L+1,M+2,1) + x(L+2,M+1,1) &
                                 + x(L+2,M+2,2))/3.0_dp
    x(L+2,1,N+2) = (x(L+1,1,N+2) + x(L+2,2,N+2) &
                                 + x(L+2,1,N+1))/3.0_dp
    x(L+2,M+2,N+2) = (x(L+1,M+2,N+2) + x(L+2,M+1,N+2) &
                                     + x(L+2,M+2,N+1))/3.0_dp


  end subroutine set_bnd_3D

  ! step densita` !

  !-----------------------------------------------------------------------------
  subroutine dens_step_2D(x,x0,u,v,diff,L,M,h,dt)
    real(dp), dimension(:,:), intent(inout) :: x, x0
    real(dp), dimension(:,:), intent(in) :: u, v
    integer, intent(in) :: L, M
    real(dp), intent(in) :: diff, h, dt

    call add_source_2D(x,x0,L,M,dt)
    ! swap x and x0
    call diffuse_2D(x0,x,diff,0,L,M,h,dt)
    ! swap x and x0
    call advect_2D(x,x0,u,v,0,L,M,h,dt)

  end subroutine dens_step_2D

  !-----------------------------------------------------------------------------
  subroutine dens_step_3D(x,x0,u,v,w,diff,L,M,N,h,dt)
    real(dp), dimension(:,:,:), intent(inout) :: x, x0
    real(dp), dimension(:,:,:), intent(in) :: u, v, w
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: diff, h, dt

    call add_source_3D(x,x0,L,M,N,dt)
    ! swap x and x0
    call diffuse_3D(x0,x,diff,0,L,M,N,h,dt)
    ! swap x and x0
    call advect_3D(x,x0,u,v,w,0,L,M,N,h,dt)

  end subroutine dens_step_3D

  ! step velocita` !

  !-----------------------------------------------------------------------------
  subroutine vel_step_2D(u,v,u0,v0,visc,L,M,h,dt)
    real(dp), dimension(:,:), intent(inout) :: u, v, u0, v0
    integer, intent(in) :: L, M
    real(dp), intent(in) :: visc, h, dt 

    call add_source_2D(u,u0,L,M,dt)
    call add_source_2D(v,v0,L,M,dt)

    ! swap u/u0 and v/v0
    call diffuse_2D(u0,u,visc,1,L,M,h,dt)
    call diffuse_2D(v0,v,visc,2,L,M,h,dt)

    call project_2D(u0,v0,u,v,L,M,h)        ! u e v array temporanei per p e div

    ! swap u/u0 and v/v0
    call advect_2D(u,u0,u0,v0,1,L,M,h,dt)
    call advect_2D(v,v0,u0,v0,2,L,M,h,dt)

    call project_2D(u,v,u0,v0,L,M,h)

  end subroutine vel_step_2D

  !-----------------------------------------------------------------------------
  subroutine vel_step_3D(u,v,w,u0,v0,w0,visc,L,M,N,h,dt)
    real(dp), dimension(:,:,:), intent(inout) :: u, v, w, u0, v0, w0
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: visc, h, dt 

    call add_source_3D(u,u0,L,M,N,dt)
    call add_source_3D(v,v0,L,M,N,dt)
    call add_source_3D(w,w0,L,M,N,dt)

    ! swap u/u0 and v/v0
    call diffuse_3D(u0,u,visc,1,L,M,N,h,dt)
    call diffuse_3D(v0,v,visc,2,L,M,N,h,dt)
    call diffuse_3D(w0,w,visc,3,L,M,N,h,dt)

    call project_3D(u0,v0,w0,u,v,L,M,N,h)  

    ! swap u/u0 and v/v0
    call advect_3D(u,u0,u0,v0,w0,1,L,M,N,h,dt)
    call advect_3D(v,v0,u0,v0,w0,2,L,M,N,h,dt)
    call advect_3D(w,w0,u0,v0,w0,3,L,M,N,h,dt)

    call project_3D(u,v,w,u0,v0,L,M,N,h)

  end subroutine vel_step_3D

  !-----------------------------------------------------------------------------
  subroutine vel_step_3D_cg(u,v,w,u0,v0,w0,visc,L,M,N,h,dt)
    real(dp), dimension(:,:,:), intent(inout) :: u, v, w, u0, v0, w0
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: visc, h, dt

    type(r_CSR) :: A

    call project_3D_initCSR(A,L,M,N)

    call add_source_3D(u,u0,L,M,N,dt)
    call add_source_3D(v,v0,L,M,N,dt)
    call add_source_3D(w,w0,L,M,N,dt)

    ! swap u/u0 and v/v0
    call diffuse_3D(u0,u,visc,1,L,M,N,h,dt)
    call diffuse_3D(v0,v,visc,2,L,M,N,h,dt)
    call diffuse_3D(w0,w,visc,3,L,M,N,h,dt)

    call project_3D_cg(u0,v0,w0,A,L,M,N,h)

    ! swap u/u0 and v/v0
    call advect_3D(u,u0,u0,v0,w0,1,L,M,N,h,dt)
    call advect_3D(v,v0,u0,v0,w0,2,L,M,N,h,dt)
    call advect_3D(w,w0,u0,v0,w0,3,L,M,N,h,dt)

    call project_3D_cg(u,v,w,A,L,M,N,h)

    call destroy(A)

  end subroutine vel_step_3D_cg

end module common
