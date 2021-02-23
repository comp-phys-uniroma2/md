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
 real(dp) :: v_drift     ! velocity of boundary particles
 
 integer :: Ngram       !Gram-S steps
 real(dp):: D0          !initial lenght tangent vectors

 logical :: scaling     ! velocity rescaling
 logical :: print_xyz
 integer :: xyz_interval 

 real(dp), dimension(:,:), allocatable :: x
 real(dp), dimension(:,:), allocatable :: v
 real(dp), dimension(:,:,:), allocatable :: eta ! Nose Hoover
 real(dp), dimension(:,:,:), allocatable :: dx
 real(dp), dimension(:,:,:), allocatable :: dv
 
 contains

 subroutine create_xv() !agg
   integer :: err
   allocate(x(3,Natoms),stat=err)
   allocate(v(3,Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_xv

 subroutine create_eta() !agg
   integer :: err
   allocate(eta(3,Natoms,5), stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_eta

 subroutine destroy_xv()
   if (allocated(x)) deallocate(x)
   if (allocated(v)) deallocate(v)
 end subroutine destroy_xv

 subroutine destroy_eta()
   if (allocated(eta)) deallocate(eta)
 end subroutine destroy_eta

 subroutine create_xvly() !agg
   integer :: err
   allocate(dx(3,Natoms,6*Natoms),stat=err)
   allocate(dv(3,Natoms,6*Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR xl or vl or etal'
 end subroutine create_xvly

 subroutine destroy_xvly()
   if (allocated(dx)) deallocate(dx)
   if (allocated(dv)) deallocate(dv)
 end subroutine destroy_xvly

 subroutine transform_units()
   
     ! transform particle mass from AMU to 
     ! something such that
     ! a = F/M  ([F]=eV/nm [a]=nm/fs^2)
     ! [M2F] = eV * (fs/nm)^2
   
     Mass = Mass * M2F

     Vol = Lx*Ly*Lz

     

 end subroutine transform_units


end module parameters
