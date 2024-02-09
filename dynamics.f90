module dynamics
  use constants, only : dp
  use parameters
  use boxes, only : update_boxes
  implicit none
  private

  public :: verlet          
  public :: verlet_nh15  

  interface  
    subroutine Tforces(x, v, F,UU,virial, Sxy, xpt)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(in) :: v
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
      real(dp), intent(out) :: Sxy
      real(dp), intent(out) :: xpt
    end subroutine TForces
  end interface

contains

  ! ---------------------------------------------------------------------------------
  ! BASIC VERLET ALGORITHM 
  subroutine verlet(x,v,x1,v1,U,virial,Sxy,xpt,dt, forces, K)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(out) :: Sxy
    real(dp), intent(out) :: xpt
    real(dp), intent(in) :: dt
    procedure(Tforces) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1

    integer :: err

    allocate(F(2,Natoms), stat=err)
    allocate(F1(2,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,v,F,U,virial, Sxy,xpt)
    x1 = x + v * dt + 0.5_dp*dt*dt*F/Mass

    call forces(x1,v1,F1,U,virial, Sxy,xpt)
    v1 = v + 0.5_dp * (F+F1)/Mass * dt

    K=kinetic(v1)

    deallocate(F,F1)  

  end subroutine verlet

 

  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with 5 levels chain
  subroutine verlet_nh15(x,v,U,virial,Sxy,xpt,dt,forces,K,K0,eta,etaf,xf,vf)
 
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: xf, vf
    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf

    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(out) :: Sxy
    real(dp), intent(out) :: xpt
    real(dp), intent(in) :: K0
    real(dp), intent(inout) :: K
    real(dp), intent(in) :: dt

    procedure(Tforces) :: forces


    ! locals

    real(dp) :: K_step,K_first
    
    real(dp), dimension(:,:,:), allocatable :: etaf_half
    real(dp), dimension(:,:), allocatable:: vf_half
   
    real(dp), dimension(:,:), allocatable :: F, F1
    integer :: err,i
   
    allocate(F(2,Natoms), stat=err)
    allocate(F1(2,Natoms), stat=err) 

    allocate(etaf_half(2,Natoms,5),stat=err)
    allocate(vf_half(2,Natoms),stat=err)

    if (err.ne.0) STOP 'ALLOCATION ERROR'

    call forces(x,v,F,U,virial,Sxy,xpt)

    xf=x + v*dt + ((F/Mass)-eta(:,:,1)*v)*dt*dt*0.5_dp
    vf_half=v + dt*0.5_dp*((F/Mass) - eta(:,:,1)*v)

    
    call forces(xf,vf,F1,U,virial,Sxy, xpt)
    

    K_first=kinetic(v)


    K_step= kinetic(vf_half)
        
    etaf_half(:,:,1)=eta(:,:,1)+dt*0.5_dp*(1/Q)*(K_first-K0)
    
    etaf(:,:,1)=etaf_half(:,:,1)+dt*0.5_dp*(1/Q)*(K_step-K0)-eta(:,:,2)*eta(:,:,1)*dt*(1/Q)
    
    etaf_half(:,:,2)=etaf(:,:,1)+dt*0.5_dp*(1/Q)*(K_first/Natoms-K0/Natoms)
    
    etaf(:,:,2)=etaf_half(:,:,2)+dt*0.5_dp*(1/Q)*(K_step/Natoms-K0/Natoms)-eta(:,:,3)*eta(:,:,2)*dt*(1/Q)
    
    etaf_half(:,:,3)=etaf(:,:,2)+dt*0.5_dp*(1/Q)*(K_first/Natoms-K0/Natoms)  

    etaf(:,:,3)=etaf_half(:,:,3)+dt*0.5_dp*(1/Q)*(K_step/Natoms-K0/Natoms)-eta(:,:,4)*eta(:,:,3)*dt*(1/Q)
    
    etaf_half(:,:,4)=etaf(:,:,3)+dt*0.5_dp*(1/Q)*(K_first/Natoms-K0/Natoms)
    
    etaf(:,:,4)=etaf_half(:,:,4)+dt*0.5_dp*(1/Q)*(K_step/Natoms-K0/Natoms)-eta(:,:,5)*eta(:,:,4)*dt*(1/Q)
    
    etaf_half(:,:,5)=etaf(:,:,4)+dt*0.5_dp*(1/Q)*(K_first/Natoms-K0/Natoms)
    
    etaf(:,:,5)=etaf_half(:,:,5)+dt*0.5_dp*(1/Q)*(K_step/Natoms-K0/Natoms)
    
    vf=(vf_half+dt*0.5_dp*F1/Mass)/(1.0_dp+dt*0.5_dp*etaf(:,:,5))

    K=kinetic(vf)
    
    deallocate(F,F1,etaf_half,Vf_half)

  end subroutine verlet_nh15

  
  ! ---------------------------------------------------------------------

  function kinetic(v) result(Ek)
     real(dp), dimension(:,:), intent(in):: v
     real(dp) :: Ek
     integer :: n
     
     Ek =0.d0
     do n = 1, Natoms
        Ek = Ek + dot_product(v(:,n),v(:,n))   
     end do
     Ek = Ek*Mass*0.5_dp

  end function kinetic

end module dynamics
