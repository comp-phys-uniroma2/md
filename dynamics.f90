module dynamics
  use constants, only : dp
  use parameters
  implicit none
  private

  public :: verlet
  
  interface  
    subroutine Tforces(x,F,U,virial)
      use constants, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), intent(out) :: U
      real(dp), intent(out) :: virial
    end subroutine Tforces
  end interface

  contains

  subroutine verlet(x,v,x1,v1,U,virial,dt,forces)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(Tforces) :: forces

    real(dp), dimension(:,:), allocatable :: F, F1
    integer :: err

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,F,U,virial) ! ALGORITMO VERLET !

    !print*,'F:',minval(F),maxval(F)

    x1 = x + v * dt + 0.5d0*dt*dt*F/Mass

    call forces(x1,F1,U,virial)

    !print*,'F1:',minval(F1),maxval(F1)

    v1 = v + 0.5d0 * (F+F1)/Mass * dt 

    deallocate(F,F1)  

  end subroutine verlet

end module dynamics
