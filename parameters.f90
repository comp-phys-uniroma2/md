module parameters
 use constants 
 implicit none

 integer :: Natoms       ! Number of atoms 
 integer :: Nx, Ny, Nz   ! Number of atoms in each direction
 real(dp) :: Lx, Ly, Lz  ! Box sizes  [nm]
 real(dp) :: Vol
 real(dp) :: Rc          ! Cutoff radius [nm]
 real(dp) :: Temp        ! Temperature [K]
 real(dp) :: tinit       ! initialization time [fs]
 real(dp) :: tsim        ! simulation time [fs]
 real(dp) :: dt          ! time step [fs] 
 integer :: Nsteps       ! number of time steps

 real(dp) :: eps         ! LJ Energy [eV]
 real(dp) :: sigma       ! LJ sigma [nm]

 real(dp) :: Mass        ! in AMU (1822.886 me)
 real(dp) :: dr          ! step in sampling g(r) 

 logical :: scaling     ! velocity rescaling
 logical :: print_xyz     ! print xyz files
 integer :: xyz_interval  ! print every so many steps 

 real(dp), dimension(:,:), allocatable :: x  ! x(1:3, j) j indice di partic
 real(dp), dimension(:,:), allocatable :: v 

 contains

 subroutine create_xv()
   integer :: err

   allocate(x(3,Natoms),stat=err)
   allocate(v(3,Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v'
 end subroutine create_xv 

 subroutine destroy_xv()
   if (allocated(x)) deallocate(x)
   if (allocated(v)) deallocate(v)
 end subroutine destroy_xv
  
 subroutine transform_units()
   
     ! transform particle mass from AMU to 
     ! something such that
     ! a = F/M  ([F]=eV/nm [a]=nm/fs^2)
     ! [M2F] = eV * (fs/nm)^2 
     Mass = Mass * M2F

     Vol = Lx*Ly*Lz

 end subroutine transform_units


end module parameters
