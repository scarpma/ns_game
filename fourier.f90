module fourier
    use precision
    use, intrinsic :: ISO_C_BINDING
    implicit none    
    include 'fftw3.f03'

    contains
    
    subroutine FFTp(b,u,ut)
        real(sp), dimension(:,:), intent(inout) :: u
        complex(sp), dimension(:,:), intent(inout) :: ut
        integer, intent(in) :: b
        type(C_PTR) :: plan
        integer :: L, M
        
        L = size(u,1)
        M = size(u,2)
        
        if (b == 1) then
            plan = fftw_plan_dft_r2c_2d(M,L,u,ut,FFTW_ESTIMATE)
            call fftw_execute_dft_r2c(plan,u,ut)
            call fftw_destroy_plan(plan)
        else
            plan = fftw_plan_dft_c2r_2d(M,L,ut,u,FFTW_ESTIMATE)
            call fftw_execute_dft_c2r(plan,ut,u)
            u = u/(L*M)
            call fftw_destroy_plan(plan)
        end if

    end subroutine FFTp
end module fourier
