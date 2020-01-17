program ns_game
    use precision
    use solvers
    use io_tools
    use bc
    implicit none
    
    real(sp), allocatable, dimension(:,:) :: u, v, u0, v0, x, x0, u1, v1
    integer :: L, M, Niter, i, err, ierr, j, argn, k
    real(sp) :: LL, diff, visc, simtime, ReL, inizio, fine, u_max_change, v_max_change
    character(64) :: argv, path
    real(sp), parameter :: conv = 0.03
    
    call cpu_time(inizio)
     
    ! CARICO INPUT DA LINEA DI COMANDO
    if (command_argument_count()<4) then
        print*, "ns_game L M Niter in_velocity [continue_last_sim]"
        stop
    end if
    call get_command_argument(1, argv)
    read(argv, *) L
    call get_command_argument(2, argv)
    read(argv, *) M
    call get_command_argument(3, argv)
    read(argv, *) Niter
    call get_command_argument(4, argv)
    read(argv, *) u_in

    print*, "Stable Fluid Simulator 1.0"
    
    ! CANCELLO FILE OUTPUT PRECEDENTI
    print*, "Cancello file output precedenti"
    call execute_command_line("rm data/*.vtk", wait=.true.)
    
    ! INIZIALIZZO PARAMETRI
    LL = 1.0_sp
    h = LL/real(M,sp)
    dt = 0.0001_sp
    simtime = real(Niter,sp)*dt
    diff = 0.001_sp
    visc = 0.01_sp
    ReL = LL*u_in/visc
    !xc = L/5
    !yc = M/2
    !rc = M/5
    if (rc > xc + 2) then
        print*, "Cerchio troppo vicino alle bc"
        stop
    end if
    call print_parameters(LL, h, dt, simtime, diff, visc, ReL)
 
    ! ALLOCO MEMORIA
    print*, "Alloco memoria e inizializzo variabili"
    allocate(u(0:L+1,0:M+1), v(0:L+1,0:M+1), stat=err)
    allocate(u0(0:L+1,0:M+1), v0(0:L+1,0:M+1), stat=err)
    allocate(u1(0:L+1,0:M+1), v1(0:L+1,0:M+1), stat=err)
    allocate(x(0:L+1,0:M+1), x0(0:L+1,0:M+1), stat=err)
    if (err > 0) then
         print*, "allocation error"
         stop
    end if

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
        call set_all_bnd(x,u,v,x0,u0,v0,set_bnd_box)
        u1 = u
        v1 = v
        ! SCRIVO CONDIZIONI INIZIALI USATE
        !call write_scalar_field(x,"data/dens_ic.dat")
        call write_vec_field(u,v,"data/vel_ic.dat")
    end if
    
    i = 0
    write(path,'(a,i4.4,a)') "data/vel_out.v",i,".vtk"
    !write(argv,'(a,i4.4,a)') "data/dens_out.v",i,".vtk"
    call out_paraview_2D_uv(u,v,path)
    !call out_paraview_2D_dens(x,argv)

    ! INTEGRO
    print*, "Integro "
    j = 0
    do i=1,Niter
        !call get_from_UI(x0,u0,v0)
        !print*, "velstep"
        call vel_step(u,v,u0,v0,visc,set_bnd_box)
        !print*, "densstep"
        !call density_step(x,x0,u,v,diff,set_bnd_box)
        !print*, "write"
        if (mod(i,Niter/30)==0) then
            write(path,'(a,i4.4,a)') "data/vel_out.v",j,".vtk"
            !write(argv,'(a,i4.4,a)') "data/dens_out.v",j,".vtk"
            !print*, "p write"
            call out_paraview_2D_uv(u,v,path)
            !call out_paraview_2D_dens(x,argv)
            !write(path,'(a,i4.4,a)') "data/vel_out.v",j,".dat"
            !write(argv,'(a,i4.4,a)') "data/dens_out.v",j,".dat"
            !call write_scalar_field(x,path)
            !call write_vec_field(u,v,argv)
            j = j + 1
        end if
        if (mod(i,500)==0) then
            u_max_change = maxval(abs((u(:,:) - u1(:,:))/u(:,:)))
            v_max_change = maxval(abs((v(:,:) - v1(:,:))/v(:,:)))
            write(*,'(" Iter ", I6, " u_change= ", ES9.2, " v_change= ", ES9.2)') i, u_max_change, v_max_change
            if (u_max_change <= conv .and. v_max_change <= conv) then
                print*, 'Convergenza al', conv,'% raggiunta. Arresto.'
                exit
            end if
            u1 = u
            v1 = v
        end if
        call progress(10*i/Niter)
    end do
    
    ! SCRIVO FILE OUTPUT
    print*, "Scrivo file output"
    !call write_scalar_field(x,"data/dens_out.dat")
    call write_vec_field(u,v,"data/vel_out.dat")
    
    ! SCRIVO SOLUZIONE ESATTA E PROFILI CALCOLATI
    !call exact()
    call write_profiles(u,v)
    
    ! DEALLOCO
    deallocate(x,x0,u,u0,v,v0,stat=err)
    if (err > 0) then
         print*, "deallocation error"
         stop
    end if
    
    call cpu_time(fine)
    fine = fine - inizio
    write(*,"(' Finito: durata totale ',F7.3,' secondi')") fine

end program ns_game
