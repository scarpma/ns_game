module io_tools
    use precision
    use solvers
    implicit none
    real(sp) :: x_exct(17), v_exct(17)
    real(sp) :: y_exct(17), u_exct(17)
    contains
    
    subroutine get_ic_scalar(p,filename)
        real(sp), intent(inout), dimension(0:,0:) :: p
        character(len=*), intent(in) :: filename
        integer :: i, j, L, ierr
        ierr = 0
        open(288,file=trim(filename),status="old",iostat=ierr)
        if (ierr > 0) then
            print*, "Problema lettura file scalare"
            stop
        end if
        L = size(p,1) - 2
        do i=0,L+1
            read(288,*) p(i,:)
        end do
        close(288)
    end subroutine get_ic_scalar 
    
    subroutine get_ic_vec(u,v,filename)
        real(sp), intent(inout), dimension(0:,0:) :: u, v
        character(len=*), intent(in) :: filename
        integer :: i, j, L, ierr
        ierr = 0
        open(287,file=trim(filename),status="old",iostat=ierr)
        if (ierr > 0) then
            print*, "Problema lettura file vettoriale"
            stop
        end if
        L = size(u,1) - 2 
        do i=0,L+1
            read(287,*) u(i,:)
            read(287,*) v(i,:)
        end do
        close(287)
    end subroutine get_ic_vec
    
    subroutine write_scalar_field(p,filename)
        real(sp), intent(inout), dimension(0:,0:) :: p
        character(len=*), intent(in) :: filename
        integer :: i, j, L
        open(194,file=trim(filename),status="replace")
        L = size(p,1) - 2 
        do i=0,L+1
            !do j=0,N+1
                write(194,*) p(i,:)
            !end do
        end do
        close(194)
    end subroutine write_scalar_field
    
    subroutine write_vec_field(u,v,filename)
        real(sp), intent(inout), dimension(0:,0:) :: u, v
        character(len=*), intent(in) :: filename
        integer :: i, j, L
        open(194,file=trim(filename),status="replace")
        L = size(u,1) - 2 
        do i=0,L+1
            write(194,*) u(i,:)
            write(194,*) v(i,:)
        end do
        close(194)
    end subroutine write_vec_field

    subroutine out_paraview_vtk(p,u,v,filename)
        real(sp), intent(inout), dimension(0:,0:) :: u, v, p
        real(sp), dimension(:), allocatable :: x, y
        character(*) :: filename
        integer :: nx, ny, nz, i,j, NXMAX, NYMAX
        character :: tab 
        
        NXMAX = size(u,1) - 2
        NYMAX = size(u,1) - 2
        tab = char(9) 
        nx = NXMAX
        ny = NYMAX
        nz = 1
        
        ! CREO GRIGLIA STUPIDA
        allocate(x(0:NXMAX+1), y(0:NXMAX+1))
        do i=0,NXMAX+1
            x(i) = real(i,sp); y(i) = real(i,sp)
        end do

        open(122, file=trim(filename), status="replace")
        
        write(122,'(A)') "# vtk DataFile Version 2.0"
        write(122,'(A)') "U V velocities"
        write(122,'(A)') "ASCII"
        write(122,'(A)') "DATASET RECTILINEAR_GRID"
        write(122,'(A,3(I5))') "DIMENSIONS",nx,ny,nz
        write(122,'(A,I3,A)') "X_COORDINATES ",nx," float"
        do i = 0, NXMAX-1
          !write(122,'(F15.7)') Xc(i,1) 
          write(122,'(F15.7)') x(i) 
        end do
        write(122,'(A,I3,A)') "Y_COORDINATES ",ny," float"
        do j = 0, NYMAX-1
          !write(122,'(F15.7)') Yc(1,j) 
          write(122,'(F15.7)') y(j) 
        end do
        write(122,'(A,I3,A)') "Z_COORDINATES ",nz," float"
        write(122,'(F15.7)') 0.0 

        write(122,'(A,I6)') "POINT_DATA",nx*ny*nz
        write(122,'(A)') "VECTORS velocity float"
 
        do j = 0, NXMAX-1
          do i = 0, NYMAX-1
             write(122,'(F15.7,A,F15.7,A,F15.7)') u(i,j), tab, v(i,j), tab, 0.0 
          end do
        end do  
    !    write(122,'(A)') "SCALARS density float"
 
    !    do j = 0, NXMAX-1
    !      do i = 0, NYMAX-1
    !         write(122,'(F15.7)') p(i,j)
    !      end do
    !    end do  
        close(122)
        deallocate(x,y)

    end subroutine out_paraview_vtk

    subroutine out_paraview_2D_uv(u,v,path)
        real(sp), dimension(0:,0:), intent(in) :: u, v
        integer :: L, M
        character(*), intent(in) :: path
        integer :: i, j
        character :: tab 
        tab = char(9) 
        open(122,file=trim(path))
        L = size(u,1)
        M = size(u,2)
        write(122,'(A)') "# vtk DataFile Version 2.0"
        write(122,'(A)') "vw"
        write(122,'(A)') "ASCII"
        write(122,'(A)') "DATASET STRUCTURED_POINTS"
        write(122,'(A,3(I5))') "DIMENSIONS",L,M,1
        write(122,'(A,F15.7,F15.7,F15.7)') "ORIGIN",0.0,0.0,0.0
        write(122,'(A,F15.7,F15.7,F15.7)') "SPACING",h,h,h
        write(122,'(A,I6)') "POINT_DATA ", L*M
        write(122,'(A)') "VECTORS velocity float"
        do j=0,M-1
            do i=0,L-1
                write(122,'(F15.7,A,F15.7,A,F15.7)') u(i,j), tab, v(i,j), tab, 0.0
            end do
        end do
        close(122)
    end subroutine out_paraview_2D_uv
    
    subroutine out_paraview_2D_dens(d,path)
        real(sp), dimension(0:,0:), intent(in) :: d
        integer :: L, M
        character(*), intent(in) :: path
        integer :: i, j
        character :: tab 
        L = size(d,1)
        M = size(d,2)
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
        do j=0,M-1
            do i=0,L-1
                write(122,'(F15.7)') d(i,j) 
            end do
        end do
        close(122)
    
    end subroutine out_paraview_2D_dens

    subroutine print_parameters(L, h, dt, simtime, diff, ReL)
        real(sp), intent(in) :: L, h, dt, simtime, diff, ReL
        write(*,"(' Parametri:')")
        write(*,"(' L=',ES9.3)") L
        write(*,"(' h=',ES9.3)") h
        write(*,"(' dt=',ES9.3)") dt
        write(*,"(' simtime=',ES9.3)") simtime
        write(*,"(' diff=',ES9.3)") diff
        write(*,"(' ReL=',ES9.3)") ReL
        !write(*,"(' xc=',ES9.3)") real(xc)*h
        !write(*,"(' yc=',ES9.3)") real(yc)*h
        !write(*,"(' rc=',ES9.3)") real(rc)*h
    end subroutine print_parameters 
    
    subroutine take_n_snapshots(n,x,u,v,ut,vt,ek,i,j,Niter)
        integer, intent(in) :: i, n, Niter
        integer, intent(inout) :: j
        real(sp), intent(in), dimension(0:,0:) :: x, u, v
        complex(sp), intent(in), dimension(0:,0:) :: ut, vt
        complex(sp), intent(inout), dimension(0:) :: ek
        
        character(64) :: path
        integer :: kkk, kkk_dim
    
        kkk_dim = size(ek,1)
        if (mod(i-1,Niter/n)==0) then
            write(path,'(a,i4.4,a)') "data/vel_out.v",j,".vtk"
            call out_paraview_2D_uv(u,v,path)
            call e_spect(ut,vt,ek)
            write(path,'(a,i4.4,a)') "data/spec_out",j,".dat"
            open(456,file=path,status="replace")
            do kkk = 0,kkk_dim-1
                write(456,"(i0,XX,2(ES10.3,XX))") kkk, real(ek(kkk)), imag(ek(kkk))
            end do
            close(456)
            j = j + 1
        end if
    end subroutine take_n_snapshots

    subroutine init_variables(x0,x,u0,u,u1,v0,v,v1, bndcnd)
        real(sp), intent(inout), dimension(0:,0:) :: x0,x,u0,u,u1,v0,v,v1
        procedure(set_bnd) :: bndcnd
        if (command_argument_count()==5) then
            ! CARICO CONDIZIONI INIZIALI DA ULTIMA SIMULAZIONE
            !call get_ic_scalar(x,"data/dens_out.dat")
            call get_ic_vec(u,v,"data/vel_out.dat")
            x0 = x
            u0 = u
            v0 = v
            u1 = u
            v1 = v
            ! SCRIVO CONDIZIONI INIZIALI USATE
            !call write_scalar_field(x,"data/dens_ic.dat")
            call write_vec_field(u,v,"data/vel_ic.dat")
        else if (command_argument_count()==4) then
            ! INIZIALIZZO VARIABILI
            x = 0._sp
            u = 0._sp
            v = 0._sp
            !call init_sources(x,u,v)
            x0 = x
            u0 = u
            v0 = v
            ! IMPONGO CONDIZIONI AL BORDO SULLA C.I.
            call set_all_bnd(x,u,v,x0,u0,v0,bndcnd)
            u1 = u
            v1 = v
            ! SCRIVO CONDIZIONI INIZIALI USATE
            !call write_scalar_field(x,"data/dens_ic.dat")
            call write_vec_field(u,v,"data/vel_ic.dat")
        end if
    end subroutine init_variables

    subroutine progress(j)
        implicit none
        integer(kind=4)::j,k
        character(len=18)::bar="????% |          |"
        write(unit=bar(1:4),fmt="(i4)") 10*j
        do k=1, j
            bar(7+k:7+k)="*"
        enddo
        ! print the progress bar.
        if (j*10 == 100) then
            write(*, fmt="(a1,a1,a18)") " ", char(13), bar
        else 
            write(*, fmt="(a1,a1,a18)", advance="no") " ", char(13), bar
        end if
        return
    end subroutine progress 

    subroutine exact()
        integer :: i
        !Re=1000, u_in = 10
        x_exct(1)=0.00282901; v_exct(1)=-0.07028612
        x_exct(2)=6.42882E-2; v_exct(2)=2.71771E+0
        x_exct(3)=7.35542E-2; v_exct(3)=2.86215E+0
        x_exct(4)=8.09753E-2; v_exct(4)=3.02728E+0
        x_exct(5)=9.76486E-2; v_exct(5)=3.25421E+0
        x_exct(6)=1.58722E-1; v_exct(6)=3.70750E+0
        x_exct(7)=2.28897E-1; v_exct(7)=3.31348E+0
        x_exct(8)=2.34432E-1; v_exct(8)=3.25139E+0
        x_exct(9)=5.01948E-1; v_exct(9)=1.87997E-1
        x_exct(10)=8.04508E-1; v_exct(10)=-3.33066E+0
        x_exct(11)=8.59780E-1; v_exct(11)=-4.42685E+0
        x_exct(12)=9.07688E-1; v_exct(12)=-5.33693E+0
        x_exct(13)=9.44865E-1; v_exct(13)=-4.07737E+0
        x_exct(14)=9.54199E-1; v_exct(14)=-3.51971E+0
        x_exct(15)=9.61692E-1; v_exct(15)=-2.92069E+0
        x_exct(16)=9.69195E-1; v_exct(16)=-2.25968E+0
        x_exct(17)=1.00098E+0; v_exct(17)=-9.09091E-2

        y_exct(1)=-1.85185E-3; u_exct(1)=0.00000E+0
        y_exct(2)=5.00000E-2; u_exct(2)=-1.84615E+0
        y_exct(3)=5.92593E-2; u_exct(3)=-2.03077E+0
        y_exct(4)=6.66667E-2; u_exct(4)=-2.21538E+0
        y_exct(5)=1.00000E-1; u_exct(5)=-2.98462E+0
        y_exct(6)=1.68519E-1; u_exct(6)=-3.81538E+0
        y_exct(7)=2.77778E-1; u_exct(7)=-2.76923E+0
        y_exct(8)=4.48148E-1; u_exct(8)=-1.07692E+0
        y_exct(9)=4.96296E-1; u_exct(9)=-6.15385E-1
        y_exct(10)=6.09259E-1; u_exct(10)=5.84615E-1
        y_exct(11)=7.31481E-1; u_exct(11)=1.84615E+0
        y_exct(12)=8.50000E-1; u_exct(12)=3.32308E+0
        y_exct(13)=9.50000E-1; u_exct(13)=4.64615E+0
        y_exct(14)=9.57407E-1; u_exct(14)=5.07692E+0
        y_exct(15)=9.64815E-1; u_exct(15)=5.72308E+0
        y_exct(16)=9.74074E-1; u_exct(16)=6.58462E+0
        y_exct(17)=9.96296E-1; u_exct(17)=1.00000E+1
        
        open(194,file="data/exact_v.dat",status="replace")
        do i=1,17
            write(194,*) x_exct(i), v_exct(i)
        end do
        close(194)
        open(195,file="data/exact_u.dat",status="replace")
        do i=1,17
            write(195,*) y_exct(i), u_exct(i)
        end do
        close(195)
    end subroutine exact

    subroutine write_profiles(u,v)
        real(sp), intent(in), dimension(0:,0:) :: u, v
        integer :: i, L, M
        L = size(u, 1)-2
        M = size(u, 2)-2
        open(194,file="data/profile_v.dat",status="replace")
        do i=0,L+1
            write(194,*) h*real(i,sp), v(i,M/2)
        end do
        close(194)
        open(195,file="data/profile_u.dat",status="replace")
        do i=0,M+1
            write(195,*) h*real(i,sp), u(L/2,i)
        end do
        close(195)
    end subroutine write_profiles

end module io_tools
