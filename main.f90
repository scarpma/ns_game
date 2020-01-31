program ns_game
    use precision
    use solvers
    use common
    use io_tools
    use bc
    implicit none
    
    real(sp), allocatable, dimension(:,:) :: u, v, u0, v0, x, x0, u1, v1, p, div
    complex(sp), allocatable, dimension(:,:) :: ut, vt, fu, fv
    integer :: L, M, Niter, i, err, j, Nf
    real(sp) :: diff, simtime, ReL, Ref, Tad, inizio, fine, eps, e_in, eta
    character(64) :: argv, path
    

    complex(sp), allocatable :: ek(:)
    integer :: kkk_dim


    call cpu_time(inizio)
    pi = 4.0_sp*atan(1.0_sp)
     
    ! CARICO INPUT DA LINEA DI COMANDO
    if (command_argument_count()<4) then
        print*, "ns_game L M Niter ReL [continue_last_sim]"
        stop
    end if
    call get_command_argument(1, argv)
    read(argv, *) L
    call get_command_argument(2, argv)
    read(argv, *) M
    call get_command_argument(3, argv)
    read(argv, *) Niter
    call get_command_argument(4, argv)
    read(argv, *) ReL

    print*, "Stable Fluid Simulator 1.0"
    
    ! CANCELLO FILE OUTPUT PRECEDENTI
    print*, "Cancello file output precedenti"
    call execute_command_line("rm data/*.vtk", wait=.true.)
    
    ! INIZIALIZZO PARAMETRI
    LL = 1.0_sp
    h = LL/real(M,sp)
    dt = 0.0001! 0.00001_sp ! should be adimensional, so changes with reynolds
    simtime = real(Niter,sp)*dt
    diff = 0.001_sp
    k0 = 2._sp*pi/LL
    sigma = 50.0_sp
    !eps = TL * sigma**2.0_sp
    eps = 10.0_sp
    print*, "eps=",eps
    TL = eps/sigma**2._sp!20*dt
    print*, "TL=",TL
    KF2 = (k0*sqrt(2._sp))**2._sp
    print*, "KF=",sqrt(KF2)
    call count_forced_modes(Nf,L,M)
    print*, "Nf=",Nf
    e_in = 4._sp*Nf*eps
    print*, "e_in=",e_in
    eta = (e_in*ReL**3._sp)**(-1._sp/4._sp)
    print*, "eta=", eta
    print*, "kf/k0=", sqrt(KF2)/k0, "eta*k0=", eta*k0!, "kmax/k0=", real(L+2)
    Ref = eps**(1./3.) * k0**(-4./3.) * ReL
    print*, "Ref=",Ref
    Tad = TL * eps * k0**(2./3.)
    print*, "Tad=", Tad
    

    kkk_dim = nint(sqrt(real(((L+2)/2)**2+((M+2)/2)**2,sp)))
    print*, kkk_dim



    !xc = L/5
    !yc = M/2
    !rc = M/5
    !if (rc > xc + 2) then
    !    print*, "Cerchio troppo vicino alle bc"
    !    stop
    !end if
    call print_parameters(LL, h, dt, simtime, diff, ReL)
 
    ! ALLOCO MEMORIA
    print*, "Alloco memoria e inizializzo variabili"
    allocate(u(0:L+1,0:M+1), v(0:L+1,0:M+1), stat=err)
    allocate(u0(0:L+1,0:M+1), v0(0:L+1,0:M+1), stat=err)
    allocate(u1(0:L+1,0:M+1), v1(0:L+1,0:M+1), stat=err)
    allocate(x(0:L+1,0:M+1), x0(0:L+1,0:M+1), stat=err)
    allocate(p(0:L+1,0:M+1), div(0:L+1,0:M+1), stat=err)
    allocate(ut(0:((L+2)/2),0:M+1), vt(0:((L+2)/2),0:M+1), stat=err)
    allocate(fu(0:((L+2)/2),0:M+1), fv(0:((L+2)/2),0:M+1), stat=err)
    allocate(ek(0:kkk_dim), stat=err)




    if (err > 0) then
         print*, "allocation error"
         stop
    end if

    ! INIZIALIZZO VARIABILI
    ! PREPARO LIBRERIA FFTW3
    planr2c = fftw_plan_dft_r2c_2d(M+2,L+2,v,vt,FFTW_ESTIMATE)
    planc2r = fftw_plan_dft_c2r_2d(M+2,L+2,vt,v,FFTW_ESTIMATE)
    call init_variables(x0,x,u0,u,u1,v0,v,v1,set_bnd_per)
    ut = complex(0._sp,0._sp)
    vt = complex(0._sp,0._sp)
    fu = complex(0._sp,0._sp)
    fv = complex(0._sp,0._sp)
    ek = complex(0._sp,0._sp)
    
    
    ! SCRIVO DATI INIZIALI
    i = 0
    write(path,'(a,i4.4,a)') "data/vel_out.v",i,".vtk"
    !write(argv,'(a,i4.4,a)') "data/dens_out.v",i,".vtk"
    call out_paraview_2D_uv(u,v,path)
    !call out_paraview_2D_dens(x,argv)

   



    ! INTEGRO
    print*, "Integro "
    j = 1
    do i=1,Niter-1
        call vel_step(u,v,u0,v0,ut,vt,fu,fv,p,div,1.0_sp/ReL,set_bnd_per)
        call take_n_snapshots(60,x,u,v,ut,vt,ek,i,j,Niter)
        call progress(10*(i+1)/Niter)
    !    call check_uv_maxerr(500,u,v,u1,v1,conv,i,conv_check)
    !    if (conv_check == 1) exit
    end do
    
    ! SCRIVO FILE OUTPUT
    print*, "Scrivo file output"
    !call write_scalar_field(x,"data/dens_out.dat")
    call write_vec_field(u,v,"data/vel_out.dat")
    
    ! SCRIVO PROFILI PER CONFRONTO CON SOL. ESATTA
    !call write_profiles(u,v)

    ! DEALLOCO
    deallocate(x,x0,u,u0,v,v0,ut,vt,fu,fv,ek,stat=err)
    if (err > 0) then
         print*, "deallocation error"
         stop
    end if
    
    call cpu_time(fine)
    fine = fine - inizio
    write(*,"(' Finito: durata totale ',F7.1,' secondi')") fine

end program ns_game
