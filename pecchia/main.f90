! Stable Fluid Solver main
! version: 0.6
! Copyleft 2016 Dipa

program main

 use precision
 use common
 use inout
 implicit none

 integer :: i,j,k,q, error

 integer, parameter :: L=40, M=40, N=40
 real(dp), dimension(L+2,M+2) :: u, v, u_prev, v_prev
 real(dp), dimension(L+2,M+2) :: dens, dens_prev

! real(dp), dimension(L+2,M+2,N+2) :: u, v, w, u_prev, v_prev, w_prev
! real(dp), dimension(L+2,M+2,N+2) :: dens, dens_prev
 real(dp) :: visc, diff, h, dt, tot
 character(32) :: path

 ! init variabili !
 
 u=0
 v=0
 !w=0
 dens=0

! ! init da file ! 
! path="uvw.vtk"
! call in_paraview_3D_uvw(u,v,w,L,M,N,path)
! path="dens.vtk"
! call in_paraview_3D_d(dens,L,M,N,path)

 dt=1.0_dp
 visc = 0.000015_dp
 diff = 0.000015_dp
 h=1.0_dp

 ! simulazione !
 do i=1,500

   u_prev=u
   v_prev=v
   !w_prev=w
!   dens_prev=0


  ! step 2D !
   !call init_sources_2D(dens_prev,u_prev,v_prev,L,M)   
   call vel_step_2D(u,v,u_prev,v_prev,visc,L,M,h,dt)
   !call dens_step_2D(dens,dens_prev,u,v,diff,L,M,h,dt)


!   ! step 3D !
!   if(i==1) then
!     call init_sources_3D(dens_prev,u_prev,v_prev,w_prev,L,M,N)
!   end if
!   call vel_step_3D_cg(u,v,w,u_prev,v_prev,w_prev,visc,L,M,h,dt)
!   call dens_step_3D(dens,dens_prev,u,v,w,diff,L,M,N,h,dt)


!   print*, "writing output"
!   write(path,'(a,i4.4,a)') "data/data.d",i,".vtk"
!   call out_paraview_2D_dens(dens,L,M,h,path)
   write(path,'(a,i4.4,a)') "data/data.v",i,".vtk"
   call out_paraview_2D_uv(u,v,L,M,h,path)
   

!    ! controllo densità !   
!    tot=0
!    do k = 2, N+1
!      do j = 2, M+1
!        do q = 2, L+1
!          tot=tot+dens(q,j,k) 
!        end do
!      end do
!    end do

   print*,"step", i      ! , "tot density", tot
 end do

! call out_paraview_3D_dens(dens,L,M,N,h,"dens_out.vtk")
! call out_paraview_3D_uvw(u,v,w,L,M,N,h,"uvw_out.vtk")

end program main
