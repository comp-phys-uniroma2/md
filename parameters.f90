module parameters
 use constants 
 implicit none

 integer :: Natoms       ! Number of atoms 
 integer :: Nx, Ny   ! Number of atoms in each direction
 real(dp) :: Lx, Ly  ! Box sizes  [nm]
 real(dp) :: Area
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
 

 real(dp):: Q=5.0_dp     ! Nose-Hoover mass
 
 logical :: scaling=.false.    ! velocity rescaling
 logical :: nose_hoover=.false. ! nose-hoover thermostat
 logical :: print_xyz=.false.  ! if xyz should be printed
 integer :: print_interval=1   ! xyz print interval




 real(dp), dimension(:,:), allocatable :: x
 real(dp), dimension(:,:), allocatable :: v
 real(dp), dimension(:,:,:), allocatable :: eta ! Nose Hoover
 
 contains

 subroutine create_xv() !agg
   integer :: err
   allocate(x(2,Natoms),stat=err)
   allocate(v(2,Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_xv

 subroutine create_eta() !agg
   integer :: err
   allocate(eta(2,Natoms,5), stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_eta

 subroutine destroy_xv()
   if (allocated(x)) deallocate(x)
   if (allocated(v)) deallocate(v)
 end subroutine destroy_xv

 subroutine destroy_eta()
   if (allocated(eta)) deallocate(eta)
 end subroutine destroy_eta

 

 subroutine transform_units()
   
     ! transform particle mass from AMU to 
     ! something such that
     ! a = F/M  ([F]=eV/nm [a]=nm/fs^2)
     ! [M2F] = eV * (fs/nm)^2
   
     Mass = Mass * M2F

     Area = Lx*Ly

     

 end subroutine transform_units


end module parameters
