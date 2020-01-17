module linear_solver
  use precision
  use mat_def
  implicit none
  private

  public :: conjugate_gradient
  public :: project_3D_initCSR

contains


  subroutine project_3D_initCSR(A_CSR,L,M,N)
    type(r_CSR), intent(out) :: A_CSR
    integer, intent(in) :: L, M, N
    integer :: i, j, k, count, nnz
  
    nnz = L*M*N*7 - (L*M+L*N+M*N)*2
    call create(A_CSR,L*M*N,L*M*N,nnz)
  
    count=1
    A_CSR%rowpnt(1) = 1

    do k = 1, N
      do j = 1, M
        do i = 1, L
          
          if(k/=1) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i + (j-1)*L + (k-2)*L*M
            count=count+1
          end if

          if(j/=1) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i + (j-2)*L + (k-1)*L*M
            count=count+1
          end if

          if(i/=1) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i-1 + (j-1)*L + (k-1)*L*M
            count=count+1
          end if


          A_CSR%nzval(count) = 6.0_dp
          if(i==1 .or. i==L) then
            A_CSR%nzval(count) = A_CSR%nzval(count)-1.0_dp           
          end if
          if(j==1 .or. j==M) then
            A_CSR%nzval(count) = A_CSR%nzval(count)-1.0_dp           
          end if
          if(k==1 .or. k==N) then
            A_CSR%nzval(count) = A_CSR%nzval(count)-1.0_dp           
          end if
          A_CSR%colind(count) = i + (j-1)*L + (k-1)*L*M
          count=count+1



          if(i/=L) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i+1 + (j-1)*L + (k-1)*L*M
            count=count+1
          end if

          if(j/=M) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i + j*L + (k-1)*L*M
            count=count+1
          end if

          if(k/=N) then
            A_CSR%nzval(count) = -1.0_dp
            A_CSR%colind(count) = i + (j-1)*L + k*L*M
          count=count+1
          end if

          A_CSR%rowpnt(i+1+(j-1)*L+(k-1)*L*M) = count

        end do
      end do
    end do

  end subroutine project_3D_initCSR

  
  subroutine conjugate_gradient(A,b,x,acc)
    type(r_CSR) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x 
    real(dp), intent(in) :: acc
    
    integer :: i, Nmax, err
    real(dp), dimension(:), allocatable :: r, d, y 
    real(dp) :: alfa, beta, gam
    
    Nmax=size(b)
    allocate(r(Nmax), stat=err)
    if (err /= 0) stop 'ERRORE'   
    allocate(d(Nmax), stat=err)
    if (err /= 0) stop 'ERRORE'   
    allocate(y(Nmax), stat=err)
    if (err /= 0) stop 'ERRORE'   
    
    x = 0.0_dp
    d = b
    r = b
    
    do i = 1, Nmax
       
      call amux(A,d,y)
       
      alfa = dot_product(r,r)/dot_product(d,y)
       
      gam = dot_product(r,r)
      
      
      if (sqrt(gam) < acc) then
        ! write(*,*) 'iter=',i,'acc=',sqrt(gam)
        exit
      end if
        
      x = x + alfa * d

      r = r - alfa * y
        
      beta = dot_product(r,r)/gam
      d = r + beta * d

    end do   
    
    deallocate(r)
    deallocate(d)
    deallocate(y)

  end subroutine conjugate_gradient
   
 
  subroutine amux(A, x, y)
    type(r_CSR), intent(in) :: A
    real(dp), intent(in) ::  x(:)
    real(dp), intent(out) :: y(:)

    integer :: i, k
    real(dp) :: t
     
    do i = 1, A%nrow 

      t = 0.0_dp
      do k = A%rowpnt(i), A%rowpnt(i+1)-1
        t = t + A%nzval(k) * x(A%colind(k))
      end do
        
      y(i) = t
    end do

  end subroutine amux   

end module linear_solver
