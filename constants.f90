module constants
  use precision    
  real(dp), parameter :: Pi =  3.14159265358979323844_dp 
  real(dp), parameter :: Kb = 8.6173850d-5    ! eV/K
  real(dp), parameter :: m_H= 1822.8966d0     ! m_e
  real(dp), parameter :: c  = 299.792458d0    ! nm/fs
  real(dp), parameter :: me = 0.510999d6      ! eV
  real(dp), parameter :: M2F= m_H*me/(c**2)   ! eV * (fs/nm)^2
  real(dp), parameter :: M2E= m_H*me          ! eV/AMU
end module constants
