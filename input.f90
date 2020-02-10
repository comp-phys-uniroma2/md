module input
  use parameters
  implicit none
  private

  public :: read_input

  contains

  subroutine read_input()

    read(*,*) Natoms 
    read(*,*) Lx 
    read(*,*) Ly 
    read(*,*) Lz
    read(*,*) Nx, Ny, Nz
    read(*,*) Rc
    read(*,*) Temp
    read(*,*) Mass 
    read(*,*) eps 
    read(*,*) sigma 
    read(*,*) tinit
    read(*,*) tsim
    read(*,*) dt
    read(*,*) scaling
    read(*,*) dr

  end subroutine



end module input
