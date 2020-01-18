program ns_game
    use precision
    use solvers
    use io_tools
    use bc
    implicit none
    
    real(sp), allocatable, dimension(:,:) :: u, v, u0, v0, x, x0, u1, v1
    integer :: L, M, Niter, i, err, ierr, j, argn, k
    real(sp) :: LL, diff, visc, simtime, ReL, inizio, fine, errmax_u, errmax_v
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
    dt = 0.0005
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

    ! INIZIALIZZO VARIABILI
    call init_variables(x0,x,u0,u,u1,v0,v,v1,set_bnd_box)
    
    ! SCRIVO DATI INIZIALI
    i = 0
    write(path,'(a,i4.4,a)') "data/vel_out.v",i,".vtk"
    !write(argv,'(a,i4.4,a)') "data/dens_out.v",i,".vtk"
    call out_paraview_2D_uv(u,v,path)
    !call out_paraview_2D_dens(x,argv)

    ! INTEGRO
    print*, "Integro "
    j = 0
    do i=1,Niter-1
        !call get_from_UI(x0,u0,v0)
        call vel_step(u,v,u0,v0,visc,set_bnd_box)
        !call density_step(x,x0,u,v,diff,set_bnd_box)
        if (mod(i-1,Niter/30)==0) then
            write(path,'(a,i4.4,a)') "data/vel_out.v",j,".vtk"
            !write(argv,'(a,i4.4,a)') "data/dens_out.v",j,".vtk"
            call out_paraview_2D_uv(u,v,path)
            !call out_paraview_2D_dens(x,argv)
            j = j + 1
        end if

        call progress(10*i/Niter)

        if (mod(i-1,500)==0) then
            errmax_u = errmax(u,u1)
            errmax_v = errmax(v,v1)
            write(*,'(" Iter ", I6, " errmax_u= ", ES7.1, " errmax_v= ", ES7.1)') i,  errmax_u, errmax_v
            if (errmax_u <= conv .and. errmax_v <= conv) then
                write(*,'("Convergenza al ", ES7.1,"% raggiunta. Arresto.")') conv*100
                exit
            end if
            u1 = u
            v1 = v
        end if
    end do
    
    ! SCRIVO FILE OUTPUT
    print*, "Scrivo file output"
    !call write_scalar_field(x,"data/dens_out.dat")
    call write_vec_field(u,v,"data/vel_out.dat")
    
    ! SCRIVO PROFILI PER CONFRONTO CON SOL. ESATTA
    call write_profiles(u,v)
    
    ! DEALLOCO
    deallocate(x,x0,u,u0,v,v0,stat=err)
    if (err > 0) then
         print*, "deallocation error"
         stop
    end if
    
    call cpu_time(fine)
    fine = fine - inizio
    write(*,"(' Finito: durata totale ',F7.1,' secondi')") fine

end program ns_game
