! Stable Fluid Solver input/output
! version: 0.6
! Copyleft 2016 Dipa

module inout 
  use precision
  implicit none
  private

  public :: out_gnuplot_2D
  public :: out_gnuplot_3D
  public :: out_gnuplot_3D_section

  public :: out_paraview_2D_dens
  public :: out_paraview_2D_uv

  public :: in_paraview_3D_dens
  public :: out_paraview_3D_dens

  public :: in_paraview_3D_uvw
  public :: out_paraview_3D_uvw
  public :: out_paraview_3D_uvw_section

  contains

  ! ------- gnuplot ------- !

  !-----------------------------------------------------------------------------
  subroutine out_gnuplot_2D(d,u,v,L,M,h,path)
    real(dp), dimension(:,:), intent(in) :: d, u, v
    integer, intent(in) :: L, M
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y
    integer :: i, j

    open(23,file=trim(path)) 
    do j=2,M+1
      y=(j-2.0_dp)*h
      do i=2,L+1
        x=(i-2.0_dp)*h
        write(23,*) x, y, d(i,j), u(i,j), v(i,j)   
      end do
    end do
    close(23) 

  end subroutine out_gnuplot_2D

  !-----------------------------------------------------------------------------
  subroutine out_gnuplot_3D(d,u,v,w,L,M,N,h,path)
    real(dp), dimension(:,:,:), intent(in) :: d, u, v, w
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y, z
    integer :: i, j, k

    open(23,file=trim(path)) 
    do k=2,N+1
      z=(k-2.0_dp)*h
      do j=2,M+1
        y=(j-2.0_dp)*h
        do i=2,L+1
          x=(i-2.0_dp)*h
          write(23,*) x, y, z, d(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k) 
        end do
      end do
    end do
    close(23) 

  end subroutine out_gnuplot_3D

  !-----------------------------------------------------------------------------
  subroutine out_gnuplot_3D_section(d,u,v,w,L,M,N,h,path)
    real(dp), dimension(:,:,:), intent(in) :: d, u, v, w
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y, z
    integer :: i, j, k

    i=L/2+1

    open(23,file=trim(path)) 
    do k=2,N+1
      z=(k-2.0_dp)*h
      do j=2,M+1
        y=(j-2.0_dp)*h
        write(23,*) y, z, d(i,j,k), v(i,j,k), w(i,j,k) 
      end do
    end do
    close(23) 

  end subroutine out_gnuplot_3D_section

  ! ------- paraview 2D ------- !

  !-----------------------------------------------------------------------------
  subroutine out_paraview_2D_dens(d,L,M,h,path)
    real(dp), dimension(:,:), intent(in) :: d
    integer, intent(in) :: L, M
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    integer :: i, j
     
    character :: tab 

    tab = char(9) 

    open(122,file=trim(path))
   
    write(122,'(A)') "# vtk DataFile Version 2.0"
    write(122,'(A)') "density"
    write(122,'(A)') "ASCII"
    write(122,'(A)') "DATASET STRUCTURED_POINTS"
    write(122,'(A,3(I5))') "DIMENSIONS",L,M,1
    write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
    write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
    write(122,'(A,I6)') "POINT_DATA ", L*M
    write(122,'(A)') "SCALARS density float 1"
    write(122,'(A)') "LOOKUP_TABLE default" 
    do j = 2, M+1
      do i = 2, L+1
        write(122,'(F15.7)') d(i,j) 
      end do
    end do
    close(122)

  end subroutine out_paraview_2D_dens

  !-----------------------------------------------------------------------------
  subroutine out_paraview_2D_uv(u,v,L,M,h,path)
    real(dp), dimension(:,:), intent(in) :: u, v
    integer, intent(in) :: L, M
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    integer :: i, j
     
    character :: tab 

    tab = char(9) 

    open(122,file=trim(path))
   
    write(122,'(A)') "# vtk DataFile Version 2.0"
    write(122,'(A)') "vw"
    write(122,'(A)') "ASCII"
    write(122,'(A)') "DATASET STRUCTURED_POINTS"
    write(122,'(A,3(I5))') "DIMENSIONS",L,M,1
    write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
    write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
    write(122,'(A,I6)') "POINT_DATA ", L*M
    write(122,'(A)') "VECTORS vectors float"
    do j = 2, M+1
      do i = 2, L+1
        write(122,'(F15.7,A,F15.7,A,F15.7)') u(i,j), tab, v(i,j), tab, 0.0
      end do
    end do
    close(122)

  end subroutine out_paraview_2D_uv

  ! ------- paraview 3D dens ------- !

  !-----------------------------------------------------------------------------
  subroutine in_paraview_3D_dens(d,L,M,N,path)
    real(dp), dimension(:,:,:), intent(out) :: d
    integer, intent(in) :: L, M, N
    character(*), intent(in) :: path

    integer :: i, j, k
    character(64) :: text
    integer:: error

    error=0

    open(101,file=trim(path),action='read',iostat=error)
      if(error==0) then

        do i=1, 10
          read(101, *) text
        end do

        do k=2,N+1
          do j=2,M+1
            do i=2,L+1
              read(101, *) d(i,j,k)
            end do
          end do
        end do

        else 
          print*, "Errore: impossibile aprire il file"
          stop
        end if
    close(101)

  end subroutine in_paraview_3D_dens

  !-----------------------------------------------------------------------------
  subroutine out_paraview_3D_dens(d,L,M,N,h,path)
    real(dp), dimension(:,:,:), intent(in) :: d
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y, z
    integer :: i, j, k
     
    character :: tab 

    tab = char(9) 

    open(122,file=trim(path))
   
    write(122,'(A)') "# vtk DataFile Version 2.0"
    write(122,'(A)') "density"
    write(122,'(A)') "ASCII"
    write(122,'(A)') "DATASET STRUCTURED_POINTS"
    write(122,'(A,3(I5))') "DIMENSIONS",L,M,N
    write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
    write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
    write(122,'(A,I6)') "POINT_DATA ", L*M*N
    write(122,'(A)') "SCALARS density float 1"
    write(122,'(A)') "LOOKUP_TABLE default" 
    do k = 2, N+1
      do j = 2, M+1
        do i = 2, L+1
          write(122,'(F15.7)') d(i,j,k) 
        end do
      end do
    end do
    close(122)

  end subroutine out_paraview_3D_dens

  ! ------- paraview 3D uvw ------- !

  !-----------------------------------------------------------------------------
  subroutine in_paraview_3D_uvw(u,v,w,L,M,N,path)
    real(dp), dimension(:,:,:), intent(out) :: u, v, w
    integer, intent(in) :: L, M, N
    character(*), intent(in) :: path

    integer :: i, j, k
    character(64) :: text
    integer:: error

    error=0

    open(101,file=trim(path),action='read',iostat=error)
      if(error==0) then

        do i=1, 9
          read(101, *) text
        end do

        do k=2,N+1
          do j=2,M+1
            do i=2,L+1
              read(101, *) u(i,j,k), v(i,j,k), w(i,j,k)
            end do
          end do
        end do

        else 
          print*, "Errore: impossibile aprire il file"
          stop
        end if
    close(101)

  end subroutine in_paraview_3D_uvw

  !-----------------------------------------------------------------------------
  subroutine out_paraview_3D_uvw(u,v,w,L,M,N,h,path)
    real(dp), dimension(:,:,:), intent(in) :: u, v, w
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y, z
    integer :: i, j, k
     
    character :: tab 

    tab = char(9) 

    open(122,file=trim(path))
   
    write(122,'(A)') "# vtk DataFile Version 2.0"
    write(122,'(A)') "uvw"
    write(122,'(A)') "ASCII"
    write(122,'(A)') "DATASET STRUCTURED_POINTS"
    write(122,'(A,3(I5))') "DIMENSIONS",L,M,N
    write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
    write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
    write(122,'(A,I6)') "POINT_DATA ", L*M*N
    write(122,'(A)') "VECTORS vectors float"
    do k = 2, N+1
      do j = 2, M+1
        do i = 2, L+1
          write(122,'(F15.7,A,F15.7,A,F15.7)') u(i,j,k), &
                           tab, v(i,j,k), tab, w(i,j,k)
        end do
      end do
    end do
    close(122)

  end subroutine out_paraview_3D_uvw

  !-----------------------------------------------------------------------------
  subroutine out_paraview_3D_uvw_section(u,v,w,L,M,N,h,path)
    real(dp), dimension(:,:,:), intent(in) :: u, v, w
    integer, intent(in) :: L, M, N
    real(dp), intent(in) :: h
    character(*), intent(in) :: path
   
    real(dp) :: x, y, z
    integer :: i, j, k
     
    character :: tab 

    tab = char(9) 

    open(122,file=trim(path))
   
    write(122,'(A)') "# vtk DataFile Version 2.0"
    write(122,'(A)') "vw"
    write(122,'(A)') "ASCII"
    write(122,'(A)') "DATASET STRUCTURED_POINTS"
    write(122,'(A,3(I5))') "DIMENSIONS",1,M,N
    write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
    write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
    write(122,'(A,I6)') "POINT_DATA ", M*N
    write(122,'(A)') "VECTORS vectors float"
    do k = 2, N+1
      do j = 2, M+1
        i = L/2+1
        write(122,'(F15.7,A,F15.7,A,F15.7)') 0.0, tab, v(i,j,k), tab, w(i,j,k)
      end do
    end do
    close(122)

  end subroutine out_paraview_3D_uvw_section

end module
