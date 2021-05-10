module dynamics
  use constants, only : dp
  use parameters, only : Natoms, Mass, Q
  use boxes, only : update_boxes
  implicit none
  private

  public :: verlet
  public :: verlet2
  public :: verlet_nh1
  public :: verlet_nh2
  public :: verlet_nh1_5
  public :: verlet_nh2_5
  
  public :: sym4

  interface  
    subroutine Tforces1(x,F,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(inout) :: F
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
    end subroutine TForces1

    subroutine Tforces2(x,dx,F,lF,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:,:), intent(in) :: dx
      real(dp), dimension(:,:), intent(inout) :: F
      real(dp), dimension(:,:,:), intent(inout) :: lF
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
    end subroutine Tforces2
  end interface

contains

  ! ---------------------------------------------------------------------------------
  ! BASIC VERLET ALGORITHM 
  subroutine verlet(x,v,x1,v1,U,virial,dt,forces,K)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(Tforces1) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1
    integer :: err

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,F,U,virial)
    x1 = x + v * dt + 0.5_dp*dt*dt*F/Mass

    call forces(x1,F1,U,virial) 
    v1 = v + 0.5_dp * (F+F1)/Mass * dt

    K=kinetic(v1)

    deallocate(F,F1)  

  end subroutine verlet

  ! ---------------------------------------------------------------------------------
  ! BASIC VERLET ALGORITHM + Linearized 
  subroutine verlet2(x,v,x1,v1,U,virial,dt,forces,K,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf  
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(Tforces2) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    integer :: err

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,dx,F,lF,U,virial)
    x1 = x + v * dt + 0.5_dp*dt*dt*F/Mass
    dxf = dx + dv*dt + 0.5_dp*dt*dt*lF/Mass

    call forces(x1,dxf,F1,lF1,U,virial) 
    v1 = v + 0.5_dp * (F+F1)/Mass * dt
    dvf= dv+ 0.5_dp*dt*(lF+lF1)/Mass
    
    K=kinetic(v1)

    deallocate(F,F1, lF, lF1)  

  end subroutine verlet2


  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with 5 levels chain
  subroutine verlet_nh1_5(x,v,U,virial,dt,forces,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
 
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf
    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf

    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf
    
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: K0 
    real(dp), intent(inout) :: K 
    real(dp), intent(in) :: dt
    procedure(Tforces2) :: forces

    ! locals
    real(dp), dimension(:,:), allocatable :: F, F1   
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    integer :: err,i
   
 
    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err) 
    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    call forces(x,dx,F,lF,U,virial)
    
    xf = x + v * dt + 0.5_dp*dt*dt*F/Mass

    dxf=dx+dv*dt+ 0.5_dp*dt*dt*lF/Mass

    call forces(xf,dxf,F1,lF1,U,virial)

    do i=1,6*Natoms
      dvf(:,:,i)=dv(:,:,i)+dt*(lF(:,:,i)+lF1(:,:,i))*0.5_dp/Mass-(eta(:,:,1)/(Q*Natoms))*dv(:,:,i)*dt
    end do

    vf = v+ 0.5_dp*((F+F1)/Mass)*dt-eta(:,:,1)*v*dt/(Q*Natoms)
    
    K=kinetic(v)

    etaf(:,:,1)=eta(:,:,1)+ (K-K0)/Q*dt - eta(:,:,1)*eta(:,:,2)*dt

    etaf(:,:,2)=eta(:,:,2)+1/Q*(etaf(:,:,1)*etaf(:,:,1)*Q-K0/Natoms)*dt-eta(:,:,2)*eta(:,:,3)*dt

    etaf(:,:,3)=eta(:,:,3)+1/Q*(etaf(:,:,2)*etaf(:,:,2)*Q-K0/Natoms)*dt-eta(:,:,3)*eta(:,:,4)*dt

    etaf(:,:,4)=eta(:,:,4)+1/Q*(etaf(:,:,3)*etaf(:,:,3)*Q-K0/Natoms)*dt-eta(:,:,4)*eta(:,:,5)*dt
    
    etaf(:,:,5)=eta(:,:,5)+1/Q*(eta(:,:,4)*eta(:,:,4)-K0/Natoms)*dt

    deallocate(F,F1,lF1,lF)

  end subroutine verlet_nh1_5


  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with 5 levels chain
  subroutine verlet_nh2_5(x,v,U,virial,dt,forces,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf

    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf

    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf
    
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(inout) :: K 
    real(dp), intent(in) :: K0 
    real(dp), intent(in) :: dt
    
    procedure(Tforces2) :: forces
    
    
    real(dp) :: K_step,K_first
    !real(dp) :: Q
    real(dp), dimension(:,:), allocatable :: F, F1
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    real(dp), dimension(:,:,:), allocatable :: etaf_half
    real(dp), dimension(:,:), allocatable:: vf_half
    integer :: i, err
   

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)

    allocate(etaf_half(3,Natoms,5),stat=err)
    
    allocate(vf_half(3,Natoms),stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    call forces(x,dx,F,lF,U,virial)
   
    dxf=dx+ dv*dt + 0.5_dp*lF/Mass
    xf=x + v*dt + (F/Mass-eta(:,:,1)*v)*dt*dt*0.5_dp   
    vf_half=v + dt*0.5_dp*F/Mass - eta(:,:,1)*v*dt
    
    call forces(xf,dxf,F1,lF1,U,virial)
    
    do i = 1, 6*Natoms
      dvf(:,:,i) = dv(:,:,i) + dt*0.5_dp*(lF(:,:,i)+lF1(:,:,i))/Mass &
                 & - eta(:,:,1)*dv(:,:,i)*dt 
    end do

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
    
    deallocate(F,F1,etaf_half,lF,lF1,Vf_half)  

  end subroutine verlet_nh2_5

  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with single level
  subroutine verlet_nh1(x,v,U,virial,dt,forces,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf
    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf
    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(inout) :: K 
    real(dp), intent(in) :: K0 
    real(dp), intent(in) :: dt
    !real(dp) :: Q
    procedure(Tforces2) :: forces
   

    real(dp), dimension(:,:), allocatable :: F, F1   
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    real(dp), dimension(:,:,:), allocatable :: etaf_half
    real(dp), dimension(:,:), allocatable:: vf_half
    integer :: err
   
    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err) 

    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)

    allocate(etaf_half(3,Natoms,5),stat=err)

    
    allocate(vf_half(3,Natoms),stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    !Q=Natoms*dt*kb*temp*100000.0_dp

    call forces(x,dx,F,lF,U,virial)
    
    xf = x + v * dt + 0.5_dp*dt*dt*F/Mass

    dxf=dx+dv*dt
 
    call forces(xf,dxf,F1,lF1,U,virial)
    
    dvf=dv+dt*(lF+lF1)*0.5_dp/Mass

    vf = v+ 0.5_dp*((F+F1)/Mass)*dt-eta(:,:,1)*v*dt
    
    K=kinetic(vf)

    etaf(:,:,1)=eta(:,:,1)+(K -K0)*(1/Q)*dt-eta(:,:,1)*dt

    deallocate(F,F1,lF1,lF)

  end subroutine verlet_nh1

  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with single level
  subroutine verlet_nh2(x,v,U,virial,dt,forces,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf
    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf
    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(inout) :: K 
    real(dp), intent(in) :: K0 
    real(dp), intent(in) :: dt
    real(dp) :: K_step,K_first
    !real(dp) :: Q
    procedure(Tforces2) :: forces


    real(dp), dimension(:,:), allocatable :: F, F1
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    real(dp), dimension(:,:,:), allocatable :: etaf_half
    real(dp), dimension(:,:), allocatable:: vf_half
    integer :: err
    integer :: i
   

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
       

    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)

    allocate(etaf_half(3,Natoms,5),stat=err)

    
    allocate(vf_half(3,Natoms),stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    !Q=Natoms*dt*kb*temp*400.0_dp
    
    call forces(x,dx,F,lF,U,virial)
   
    dxf=dx+dv*dt
    xf=x+v*dt+(F/Mass-eta(:,:,1)*v)*dt*dt*0.5_dp   
    vf_half=v+dt*0.5_dp*(F/Mass)-eta(:,:,1)*v*dt
    
    call forces(xf,dxf,F1,lF1,U,virial)
    
    do i = 1, 6*Natoms
      dvf(:,:,i)=dv(:,:,i)+dt*(lF(:,:,i)+lF1(:,:,i))*0.5_dp/Mass &
            & - eta(:,:,1)*dv(:,:,i)*dt
    end do

    K_first=kinetic(v)

    K_step= kinetic(vf_half)
        
    etaf_half(:,:,1)=eta(:,:,1)+dt*0.5_dp*(1/Q)*(K_first-K0)
    
    etaf(:,:,1)=etaf_half(:,:,1)+dt*0.5_dp*(1/Q)*(K_step-K0)-eta(:,:,2)*eta(:,:,1)*dt*(1/Q)
    
    vf=(vf_half+dt*0.5_dp*F1/Mass)/(1+dt*0.5_dp*etaf(:,:,1))

    K=kinetic(vf)
    
    deallocate(F,F1,etaf_half,lF,lF1,Vf_half)  

  end subroutine verlet_nh2

  ! ---------------------------------------------------------------------------------
  ! VERLET + ??
  subroutine verlet_g(x,v,U,virial,dt,forces,K,K0,a,xf,vf,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf

    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf
    
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(inout) :: K 
    real(dp), intent(in) :: K0 
    real(dp), intent(in) :: dt
    procedure(Tforces2) :: forces
    

    real(dp) :: a, K_step,K_first
    real(dp), dimension(:,:), allocatable :: F, F1
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    real(dp),dimension(Natoms):: Fv, vv
    real(dp), dimension(3,Natoms):: vs
    integer :: i, err
   
    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    
    call forces(x,dx,F,lF,U,virial)

    xf=x+v*dt+(F/Mass)*dt*dt*0.5d0  
   
    dxf=dx+dv*dt+(lF/Mass)*dt*dt*0.5d0
 
    call forces(xf,dxf,F1,lF1,U,virial)

    do i=1,Natoms
      FV(i) = dot_product(0.5d0*(F1(:,i)+F(:,i))/Mass,v(:,i))
      VV(i) = dot_product(vf(:,i),vf(:,i))
    end do

    a=0.0_dp
    write(*,*) sum(FV(:)),sum(VV(:))
   
    vf = v+dt*0.5d0*(F1+F)/Mass-a*v*dt
    dvf= dv+dt*0.5d0*(lF1+lF)/Mass-a*dv*dt
  
    K=kinetic(vf)
        
    deallocate(F,F1,lF,lF1)  

  end subroutine verlet_g

  

  ! ---------------------------------------------------------------------------------
  ! Symplectic integration of 4th order
  subroutine sym4(x,v,x1,v1,U,virial,dt,forces,K)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(TForces1) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1
    integer :: err
    real(dp) :: dt1
    real(dp), parameter :: alpha = 1.0_dp/(2.0_dp - 2.0_dp**(1/3))
    real(dp), parameter :: beta = -2.0_dp**(1/3)*alpha

    real(dp), parameter :: a1=0.205177661542290_dp 
    real(dp), parameter :: a2=0.403021281604210_dp 
    real(dp), parameter :: a3=-0.12092087633891_dp 
    real(dp), parameter :: a4=0.512721933192410_dp 
    real(dp), parameter :: a5=0.0_dp 

    real(dp), parameter :: b1=0.061758858135626_dp
    real(dp), parameter :: b2=0.33897802655364_dp
    real(dp), parameter :: b3=0.61479130717558_dp
    real(dp), parameter :: b4=-0.14054801465937_dp 
    real(dp), parameter :: b5=0.12501982279453_dp 

    allocate(F(3,Natoms), stat=err)
    !allocate(F1(3,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    ! See S.K.Gray et al., J.Chem.Phys. 101, 4062 (1994)
    call forces(x,F,U,virial)
    v1 = v + F/Mass * b1 * dt 
    x1 = x + v1 * a1 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b2 * dt 
    x1 = x1 + v1 * a2 * dt
    
    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b3 * dt 
    x1 = x1 + v1 * a3 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b4 * dt 
    x1 = x1 + v1 * a4 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b5 * dt 
    
    K=kinetic(v1)
    deallocate(F)

  end subroutine sym4
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
