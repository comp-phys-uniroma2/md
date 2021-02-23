module dynamics
  use constants, only : dp
  use parameters
  use boxes, only : update_boxes
  use lyap_n
  implicit none
  private

  public :: verlet
  public :: verlet_nh1
  public :: verlet_nh2
  public :: verlet_nh1_5
  public :: verlet_nh2_5
  
  public :: linearizzazione
  public :: sym4

  interface  
    subroutine Tforces1(x,F,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
    end subroutine TForces1

    subroutine Tforces2(x,dx,F,lF,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:,:), intent(in) :: dx
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), dimension(:,:,:), intent(out) :: lF
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
    real(dp), intent(inout) :: K 
    real(dp), intent(inout) :: K0 
    real(dp), intent(in) :: dt
    real(dp) :: Q
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

      
    Q=Natoms*dt*kb*temp*100000.0_dp

    call forces(x,dx,F,lF,U,virial)
    
    xf = x + v * dt + 0.5_dp*dt*dt*F/Mass

    dxf=dx+dv*dt
    
    call forces(xf,dxf,F1,lF1,U,virial)
    
    dvf=dv+dt*(lF+lF1)*0.5_dp/Mass

    vf = v+ 0.5_dp*((F+F1)/Mass)*dt-eta(:,:,1)*v*dt
    
    K=kinetic(vf)

    etaf(:,:,1)=eta(:,:,1)+(K -K0)*(1/Q)*dt-eta(:,:,1)*eta(:,:,2)*dt

    etaf(:,:,2)=eta(:,:,2)+1/Q*((etaf(:,:,1))*(etaf(:,:,1))*Q-K0/Natoms)*dt-eta(:,:,2)*eta(:,:,3)*dt

    etaf(:,:,3)=eta(:,:,3)+1/Q*((etaf(:,:,2))*(etaf(:,:,2))*Q-K0/Natoms)*dt-eta(:,:,3)*eta(:,:,4)*dt

    etaf(:,:,4)=eta(:,:,4)+1/Q*((etaf(:,:,3))*(etaf(:,:,3))*Q-K0/Natoms)*dt-eta(:,:,4)*eta(:,:,5)*dt
    
    etaf(:,:,5)=eta(:,:,5)+1/Q*((eta(:,:,4))*(eta(:,:,4))-K0/Natoms)*dt

    deallocate(F,F1,lF1,lF)

  end subroutine verlet_nh1_5

  ! ---------------------------------------------------------------------------------
  ! VERLET ALGORITHM + Nose-Hoover with 5 levels chain
  subroutine verlet_nh2_5(x,v,U,virial,dt,forces,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
    integer :: i
    real(dp), dimension(:,:), intent(in) :: x, v 
    real(dp), dimension(:,:), intent(out) :: xf, vf

    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf

    real(dp),dimension(:,:,:), intent(in):: eta
    real(dp),dimension(:,:,:), intent(out):: etaf
    
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(inout) :: K 
    real(dp), intent(inout) :: K0 
    real(dp), intent(in) :: dt
    
    real(dp) :: K_step,K_first
    real(dp) :: Q
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

    Q=Natoms*dt*kb*temp*50000.0_dp

    call forces(x,dx,F,lF,U,virial)
   
    dxf=dx+dv*dt
    
    xf=x+v*dt+(F/Mass-eta(:,:,1)*v)*dt*dt*0.5_dp   
    vf_half=v+dt*0.5_dp*(F/Mass)-eta(:,:,1)*v*dt
    
    call forces(xf,dxf,F1,lF1,U,virial)
    
    dvf=dv+dt*(lF+lF1)*0.5_dp/Mass
    
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
    real(dp), intent(inout) :: K0 
    real(dp), intent(in) :: dt
    real(dp) :: Q
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

    Q=Natoms*dt*kb*temp*100000.0_dp

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
    real(dp), intent(inout) :: K0 
    real(dp), intent(in) :: dt
    real(dp) :: K_step,K_first
    real(dp) :: Q
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

    Q=Natoms*dt*kb*temp*100000.0_dp
    
    call forces(x,dx,F,lF,U,virial)
   
    dxf=dx+dv*dt
    xf=x+v*dt+(F/Mass-eta(:,:,1)*v)*dt*dt*0.5_dp   
    vf_half=v+dt*0.5_dp*(F/Mass)-eta(:,:,1)*v*dt
    
    call forces(xf,dxf,F1,lF1,U,virial)
    
    dvf=dv+dt*(lF+lF1)*0.5_dp/Mass

    K_first=kinetic(v)

    K_step= kinetic(vf_half)
        
    etaf_half(:,:,1)=eta(:,:,1)+dt*0.5_dp*(1/Q)*(K_first-K0)
    
    etaf(:,:,1)=etaf_half(:,:,1)+dt*0.5_dp*(1/Q)*(K_step-K0)-eta(:,:,2)*eta(:,:,1)*dt*(1/Q)
    
    vf=(vf_half+dt*0.5_dp*F1/Mass)/(1+dt*0.5_dp*etaf(:,:,1))

    K=kinetic(vf)
    
    deallocate(F,F1,etaf_half,lF,lF1,Vf_half)  

  end subroutine verlet_nh2

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


  subroutine gram(A,B,C,D,E,F,pq,Apq,Bpq,Cpq,Dpq,Epq,Fpq,ngram,nstep1)
    real(dp) :: no0,no1,no2,no3,no4,no5,no6
    real(dp) :: co1,co2,co3,co4,co5,co6
    real(dp) :: BA , BC , AC, BD, CD, AD, AE, BE, CE, DE, AF, BF, CF, DF, EF
    integer :: q,err,pq
    integer :: ngram,nstep1
    real(dp), dimension(:,:),intent(inout) :: A,B,C,D,E,F
    real(dp), dimension(:,:),allocatable :: Estar,Fstar
    real(dp), dimension(:,:,:),intent(out) :: Apq,Bpq,Cpq,Dpq,Epq,Fpq
    real(dp) :: NA,NB,NC,ND,NE,NF

    allocate(Fstar(6,Natoms), stat=err)
    allocate(Estar(6,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'  
   
    no0=1.d0
  
    do q=1,Natoms
      Apq(:,pq,q)=A(:,q)
      Bpq(:,pq,q)=B(:,q)
      Cpq(:,pq,q)=C(:,q)
      Dpq(:,pq,q)=D(:,q)
      Epq(:,pq,q)=E(:,q)
      Fpq(:,pq,q)=F(:,q)
  
      NB=sqrt(dot_product(B(:,q),B(:,q)))
      NC=sqrt(dot_product(C(:,q),C(:,q)))
      ND=sqrt(dot_product(D(:,q),D(:,q)))
      NE=sqrt(dot_product(E(:,q),E(:,q)))
      NF=sqrt(dot_product(F(:,q),F(:,q)))
  
      no1=sqrt(dot_product(A(:,q),A(:,q)))
      co1=no0/no1
  
      BA=dot_product(B(:,q),A(:,q))
   
      no2=sqrt(dot_product(B(:,q)-(BA)*A(:,q),B(:,q)-(BA)*A(:,q)))
      co2=no0/no2
  
      B(:,q)=B(:,q)-(BA)*A(:,q)/NA
     
      BC=dot_product(C(:,q),B(:,q))
      AC=dot_product(C(:,q),A(:,q))
      
      no3=sqrt(dot_product(C(:,q)-(AC)*A(:,q)-(BC)*B(:,q),C(:,q)-(AC)*A(:,q) -(BC)*B(:,q)))
      co3=no0/no3
      C(:,q)=C(:,q)-(AC)*A(:,q)/NA -(BC)*B(:,q)/NB
  
      BD=dot_product(D(:,q),B(:,q))
      AD=dot_product(D(:,q),A(:,q))
      CD=dot_product(D(:,q),C(:,q))
      
      no4=sqrt(dot_product(D(:,q)-(AD)*A(:,q)-(BD)*B(:,q)-(CD)*C(:,q),D(:,q)-(AD)*A(:,q) -(BD)*B(:,q)-(CD)*C(:,q)))
      co4=no0/no4
      D(:,q)=D(:,q)-(AD)*A(:,q)/NA -(BD)*B(:,q)/NB -(CD)*C(:,q)/NC
  
      BE=dot_product(E(:,q),B(:,q))
      AE=dot_product(E(:,q),A(:,q))
      CE=dot_product(E(:,q),C(:,q))
      DE=dot_product(E(:,q),D(:,q))
  
      Estar(:,q)=E(:,q)-(AE)*A(:,q)-(BE)*B(:,q)-(CE)*C(:,q)-(DE)*D(:,q)
  
      no5=sqrt(dot_product(Estar(:,q),Estar(:,q)))
      co5=no0/no5
      E(:,q)=E(:,q)-(AE)*A(:,q)/NA -(BE)*B(:,q)/NB -(CE)*C(:,q)/NC-(DE)*D(:,q)/ND
  
      BF=dot_product(F(:,q),B(:,q))
      AF=dot_product(F(:,q),A(:,q))
      CF=dot_product(F(:,q),C(:,q))
      DF=dot_product(F(:,q),D(:,q))
      EF=dot_product(F(:,q),E(:,q))
  
      
      Fstar(:,q)=F(:,q)-(AF)*A(:,q)-(BF)*B(:,q)-(CF)*C(:,q)-(DF)*D(:,q)-(EF)*E(:,q)
  
      no6=sqrt(dot_product(Fstar(:,q),Fstar(:,q)))
      co6=no0/no6
      F(:,q)=F(:,q)-(AF)*A(:,q)/NA -(BF)*B(:,q)/NB -(CF)*C(:,q)/NC-(DF)*D(:,q)/ND -(EF)*E(:,q)/NE
      
      A(:,q)=A(:,q)*co1
      B(:,q)=B(:,q)*co2
      C(:,q)=C(:,q)*co3
      D(:,q)=D(:,q)*co4
      E(:,q)=E(:,q)*co5
      F(:,q)=F(:,q)*co6
  
    end do
    deallocate(Fstar,Estar)

  end subroutine gram

  subroutine linearizzazione(x,Xl,Yl,Zl,Xlf,Ylf,Zlf,Vxl,Vyl,Vzl,Vxlf,Vylf,Vzlf,linfor,dt)
   real(dp),dimension(:,:),intent(in):: x
   
   integer :: err,i,j
   real(dp),dimension(:,:),intent(in):: Xl,Yl,Zl
   real(dp),dimension(:,:),intent(out):: Xlf,Ylf,Zlf
   real(dp),dimension(:,:),intent(in):: Vxl,Vyl,Vzl
   real(dp),dimension(:,:),intent(out):: Vxlf, Vylf, Vzlf
   real(dp)::dt
   interface 
     subroutine linfor(x,lF) !Agg
        use precision, only : dp
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(:,:,:), intent(out) :: lF
     end subroutine linfor
   end interface
   
   real(dp), dimension(:,:,:), allocatable :: lF,lF2
   real(dp), dimension(:,:), allocatable :: A,B,C,D,E,F
   allocate(lF(3,3,Natoms), stat=err)
   allocate(lF2(6,6,Natoms), stat=err)

   allocate(A(6,Natoms), stat=err)

   A=0.d0
   
   call linfor(x,lF)

   do j=1,6
      if (j .le. 3) then
        do i=1,6
        if (i .le. 3) then
          lF2(i,j,:)=lF(i,j,:)*dt*0.5_dp/Mass
        else
          lF2(i,j,:)=lF(i,j,:)/Mass
        end if
        end do
      else
         do i=1,6
            if (i.ge.3) then
               lF2(i,j,:)=0.d0
            else
               if(i==j) then
                 lF2(i,j,:)=1.d0/Mass
               else
                 lF2(i,j,:)=0.d0
               end if
            end if
         end do
      end if
   end do

   do i=1,6
     A(1,:)=A(1,:)+ (lF2(1,i,:)*Vxl(i,:))
     A(2,:)=A(2,:)+ (lF2(2,i,:)*Vxl(i,:))
     A(3,:)=A(3,:)+ (lF2(3,i,:)*Vxl(i,:))

     A(4,:)=A(4,:)+ (lF2(4,i,:)*Vxl(i,:))
     A(5,:)=A(5,:)+ (lF2(5,i,:)*Vxl(i,:))
     A(6,:)=A(6,:)+ (lF2(6,i,:)*Vxl(i,:))
   end do

   Vxlf(1,:)=Vxl(1,:)+dt*A(1,:)
   Vxlf(2,:)=Vxl(2,:)+dt*A(2,:)
   Vxlf(3,:)=Vxl(3,:)+dt*A(3,:)

   
   Vxlf(4,:)=Vxl(4,:)+dt*A(4,:)
   Vxlf(5,:)=Vyl(5,:)+dt*A(5,:)
   Vxlf(6,:)=Vzl(6,:)+dt*A(6,:)


   do i=1,6
      A(1,:)=A(1,:)+ (lF2(1,i,:)*Vyl(i,:))
      A(2,:)=A(2,:)+ (lF2(2,i,:)*Vyl(i,:))
      A(3,:)=A(3,:)+ (lF2(3,i,:)*Vyl(i,:))

      A(4,:)=A(4,:)+ (lF2(4,i,:)*Vyl(i,:))
      A(5,:)=A(5,:)+ (lF2(5,i,:)*Vyl(i,:))
      A(6,:)=A(6,:)+ (lF2(6,i,:)*Vyl(i,:))
   end do

   Vylf(1,:)=Vyl(1,:)+dt*A(1,:)
   Vylf(2,:)=Vyl(2,:)+dt*A(2,:)
   Vylf(3,:)=Vyl(3,:)+dt*A(3,:)

   
   Vylf(4,:)=Vyl(4,:)+dt*A(4,:)
   Vylf(5,:)=Vyl(5,:)+dt*A(5,:)
   Vylf(6,:)=Vyl(6,:)+dt*A(6,:)

   do i=1,6
     A(1,:)=A(1,:)+ (lF2(1,i,:)*Vzl(i,:))
     A(2,:)=A(2,:)+ (lF2(2,i,:)*Vzl(i,:))
     A(3,:)=A(3,:)+ (lF2(3,i,:)*Vzl(i,:))

     A(4,:)=A(4,:)+ (lF2(4,i,:)*Vzl(i,:))
     A(5,:)=A(5,:)+ (lF2(5,i,:)*Vzl(i,:))
     A(6,:)=A(6,:)+ (lF2(6,i,:)*Vzl(i,:))
   end do
         

   Vzlf(1,:)=Vzl(1,:)+dt*A(1,:)
   Vzlf(2,:)=Vzl(2,:)+dt*A(2,:)
   Vzlf(3,:)=Vzl(3,:)+dt*A(3,:)

   
   Vzlf(4,:)=Vzl(4,:)+dt*A(4,:)
   Vzlf(5,:)=Vzl(5,:)+dt*A(5,:)
   Vzlf(6,:)=Vzl(6,:)+dt*A(6,:)


   do i=1,6
      A(1,:)=A(1,:)+ (lF2(1,i,:)*Xl(i,:))
      A(2,:)=A(2,:)+ (lF2(2,i,:)*Xl(i,:))
      A(3,:)=A(3,:)+ (lF2(3,i,:)*Xl(i,:))

      A(4,:)=A(4,:)+ (lF2(4,i,:)*Xl(i,:))
      A(5,:)=A(5,:)+ (lF2(5,i,:)*Xl(i,:))
      A(6,:)=A(6,:)+ (lF2(6,i,:)*Xl(i,:))
   end do
         

    Xlf(1,:)=Xl(1,:)+dt*A(1,:)
    Xlf(2,:)=Xl(2,:)+dt*A(2,:)
    Xlf(3,:)=Xl(3,:)+dt*A(3,:)
    Xlf(4,:)=Xl(4,:)+dt*A(4,:)
    Xlf(5,:)=Xl(5,:)+dt*A(5,:)
    Xlf(6,:)=Xl(6,:)+dt*A(6,:)

    do i=1,6
      A(1,:)=A(1,:)+ (lF2(1,i,:)*Yl(i,:))
      A(2,:)=A(2,:)+ (lF2(2,i,:)*Yl(i,:))
      A(3,:)=A(3,:)+ (lF2(3,i,:)*Yl(i,:))

      A(4,:)=A(4,:)+ (lF2(4,i,:)*Yl(i,:))
      A(5,:)=A(5,:)+ (lF2(5,i,:)*Yl(i,:))
      A(6,:)=A(6,:)+ (lF2(6,i,:)*Yl(i,:))
    end do

    Ylf(1,:)=Yl(1,:)+dt*A(1,:)
    Ylf(2,:)=Yl(2,:)+dt*A(2,:)
    Ylf(3,:)=Yl(3,:)+dt*A(3,:)
     
    Ylf(4,:)=Yl(4,:)+dt*A(4,:)
    Ylf(5,:)=Yl(5,:)+dt*A(5,:)
    Ylf(6,:)=Yl(6,:)+dt*A(6,:)

    do i=1,6
       A(1,:)=A(1,:)+ (lF2(1,i,:)*Zl(i,:))
       A(2,:)=A(2,:)+ (lF2(2,i,:)*Zl(i,:))
       A(3,:)=A(3,:)+ (lF2(3,i,:)*Zl(i,:))

       A(4,:)=A(4,:)+ (lF2(4,i,:)*Zl(i,:))
       A(5,:)=A(5,:)+ (lF2(5,i,:)*Zl(i,:))
       A(6,:)=A(6,:)+ (lF2(6,i,:)*Zl(i,:))
    end do

    Zlf(1,:)=Zl(1,:)+dt*A(1,:)
    Zlf(2,:)=Zl(2,:)+dt*A(2,:)
    Zlf(3,:)=Zl(3,:)+dt*A(3,:)
    
    Zlf(4,:)=Zl(4,:)+dt*A(4,:)
    Zlf(5,:)=Zl(5,:)+dt*A(5,:)
    Zlf(6,:)=Zl(6,:)+dt*A(6,:)
    
    deallocate(lF,lF2,A)

  end subroutine linearizzazione


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
