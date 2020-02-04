program ns_game
    use precision
    use solvers
    use common
    use io_tools
    use bc
    implicit none
    
    real(sp), allocatable, dimension(:,:) :: u, v, u0, v0
    real(sp), allocatable, dimension(:) :: sh
    complex(sp), allocatable, dimension(:,:) :: ut, vt, fu, fv
    complex(sp), allocatable :: ek(:)
    integer :: L, M, Niter, i, err, j, Nf, kkk_dim, n_modes
    real(sp) :: simtime, ReL, Ref, Tad, inizio, fine, eps, e_in, eta
    character(64) :: argv, path

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

    print*, "Stable Fluid Simulator FFT 1.0"
    
    ! CANCELLO FILE OUTPUT PRECEDENTI
    print*, "Cancello file output precedenti"
    call execute_command_line("rm data/*.vtk", wait=.true.)
    
    ! INIZIALIZZO PARAMETRI
    LL = 1.0_sp
    h = LL/real(M,sp)
    dt = 0.0005
    simtime = real(Niter,sp)*dt
    k0 = 2._sp*pi/LL
    sigma = 50._sp
    print*, "sigma=", sigma
    n_modes = 100
    !eps = 10.0_sp
    !print*, "eps=",eps
    !TL = eps/sigma**2._sp!20*dt
    TL = 100._sp*dt
    print*, "TL=",TL
    eps = TL * sigma**2.0_sp
    print*, "eps=",eps
    KI2 = (k0*sqrt(2._sp))**2._sp
    KF2 = 4._sp*(k0*sqrt(2._sp))**2._sp
    print*, "KF=",sqrt(KF2)
    print*, "KI=",sqrt(KI2)
    kmax = sqrt(real((L/2)**2+(M/2)**2))*k0
    print*, "kmax=",kmax
    call count_forced_modes(Nf,L,M)
    print*, "Nf=",Nf
    e_in = 4._sp*Nf*eps
    print*, "e_in=",e_in
    eta = (e_in*ReL**3._sp)**(-1._sp/4._sp)
    print*, "eta=", eta
    print*, "kf/k0=", sqrt(KF2)/k0, "kmax/keta=", kmax/(2*pi/eta), "kmax/k0=", kmax/k0
    Ref = eps**(1./3.) * k0**(-4./3.) * ReL
    print*, "Ref=",Ref
    Tad = TL * eps * k0**(2./3.)
    print*, "Tad=", Tad
    
    call print_parameters(LL, h, dt, simtime, ReL)

 
    ! ALLOCO MEMORIA
    print*, "Alloco memoria e inizializzo variabili"
    kkk_dim = nint(sqrt(real((L/2)**2+(M/2)**2,sp)))
    print*, kkk_dim
    allocate(u(0:L-1,0:M-1), v(0:L-1,0:M-1), stat=err)
    allocate(u0(0:L-1,0:M-1), v0(0:L-1,0:M-1), stat=err)
    allocate(ut(0:(L/2),0:M-1), vt(0:(L/2),0:M-1), stat=err)
    allocate(fu(0:(L/2),0:M-1), fv(0:(L/2),0:M-1), stat=err)
    allocate(ek(0:kkk_dim), stat=err)
    allocate(sh(0:((L/2+1)*(M/2+1)/n_modes)+100), stat=err)
    if (err > 0) then
         print*, "allocation error"
         stop
    end if

    ! PREPARO LIBRERIA FFTW3
    planr2c = fftw_plan_dft_r2c_2d(M,L,v,vt,FFTW_ESTIMATE)
    planc2r = fftw_plan_dft_c2r_2d(M,L,vt,v,FFTW_ESTIMATE)


    ! INIZIALIZZO VARIABILI
    call random_seed()
    call init_variables(u0,u,v0,v)
    ut = complex(0._sp,0._sp)
    vt = complex(0._sp,0._sp)
    fu = complex(0._sp,0._sp)
    fv = complex(0._sp,0._sp)
    ek = complex(0._sp,0._sp)
    
    
    ! SCRIVO DATI INIZIALI
    i = 0
    open(457,file="en_out.dat",status="replace")
    write(457,"(F10.5,XX,ES10.3,XX)") i*dt, tot_en(u,v)
    write(path,'(a,i4.4,a)') "data/vel_out.v",i,".vtk"
    !write(argv,'(a,i4.4,a)') "data/dens_out.v",i,".vtk"
    call out_paraview_2D_uv(u,v,path)
    !call out_paraview_2D_dens(x,argv)
    
    ! INIZIALIZZO SHELLS PER CALCOLARE LO SPETTRO
    call create_shells(sh,ut, n_modes)
    !print*, sh


    ! INTEGRO
    print*, "Integro "
    j = 1
    do i=1,Niter-1
        call vel_step(u,v,u0,v0,ut,vt,fu,fv,1.0_sp/ReL)
        call take_n_snapshots(60,u,v,ut,vt,ek,i,j,Niter)
        call progress(10*(i+1)/Niter)
    end do
    
    ! SCRIVO FILE OUTPUT
    print*, "Scrivo file output"
    call write_vec_field(u,v,"data/vel_out.dat")
    
    ! DEALLOCO
    deallocate(u,u0,v,v0,ut,vt,fu,fv,ek,stat=err)
    if (err > 0) then
         print*, "deallocation error"
         stop
    end if
    close(457)
    
    call cpu_time(fine)
    fine = fine - inizio
    write(*,"(' Finito: durata totale ',F7.1,' secondi')") fine

end program ns_game
