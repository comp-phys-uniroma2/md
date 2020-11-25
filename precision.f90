module precision
  implicit none
  private

  integer, parameter, public :: dp = selected_real_kind(15)
  public :: mach

  contains
 
  function mach() result(eps)
    real(dp) :: eps

    eps = 1.0_dp

    do while (1.0_dp+eps /= 1.0_dp)
      eps = eps / 2.0_dp
    enddo
        
  end function mach


      
end module precision

