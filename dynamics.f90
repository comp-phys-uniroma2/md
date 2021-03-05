module dynamics
  use constants, only : dp
  use parameters
  use boxes, only : update_boxes
  implicit none
  private

  public :: verlet
  public :: reynolds 
  
  interface  
    subroutine Tforces(x,F,U,virial)
      use constants, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), intent(out) :: U
      real(dp), dimension(:,:), intent(out) :: virial
    end subroutine Tforces
  end interface

  contains

  subroutine verlet(x,v,x1,v1,U,virial,dt,forces)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), dimension(:,:), intent(out) :: virial
    real(dp), intent(in) :: dt
   
    procedure(Tforces) :: forces

    real(dp), dimension(:,:), allocatable :: F, F1 
    
    integer :: err

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
   
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    !print*,'x0max=',maxval(x(1,:)),maxval(x(2,:)),maxval(x(3,:))
    !print*,'x0min=',minval(x(1,:)),minval(x(2,:)),minval(x(3,:))
    !print*,'v0max=',maxval(v(1,:)),maxval(v(2,:)),maxval(v(3,:))
    !print*,'v0min=',minval(v(1,:)),minval(v(2,:)),minval(v(3,:))

    call forces(x,F,U,virial) ! ALGORITMO VERLET !

    !print*,'F:',minval(F),maxval(F)
    !F2 = F
    x1 = x + v * dt + 0.5d0*dt*dt*F/Mass
     
    !call update_boxes(x,x1)

    !print*,'x1max=',maxval(x1(1,:)),maxval(x1(2,:)),maxval(x1(3,:))
    !print*,'x1min=',minval(x1(1,:)),minval(x1(2,:)),minval(x1(3,:))
    call forces(x1,F1,U,virial)

    !print*,'F1:',minval(F1),maxval(F1)

    v1 = v + 0.5d0 * (F+F1)/Mass * dt 
    !print*,'v1max=',maxval(v1(1,:)),maxval(v1(2,:)),maxval(v1(3,:))
    !print*,'v1min=',minval(v1(1,:)),minval(v1(2,:)),minval(v1(3,:))
   
    deallocate(F,F1)
    

  end subroutine verlet



  subroutine reynolds(x,v,virial)
    
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(in) :: virial
    real(dp) :: vel, Pxz
    real(dp) :: Vol, density
    integer :: err, i, maxi, mini
    real(dp) :: const, eta, dv, dz

    real(dp), dimension(:), allocatable :: z, vx

    allocate(z(Natoms))
    allocate(vx(Natoms))
    
    Vol = Lx*Ly*Lz
    density = Mass*Natoms/Vol

    !if (err.ne.0) STOP 'ALLOCATION ERROR'

    ! vx(z) = dvx/dz * z = gam * z
    ! vx_k = gam * z_k
    ! Xi^2 = Sum_k (vx_k - gam z_k)^2 
    ! d/d gam Xi^2 = 0
    ! gam = Sum_k (z_k vx_k) / Sum_k (z_k)
    !  

    vel = 0.0_dp
    do i=1,Natoms 
       vel = vel + sqrt(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
       Pxz = Pxz + Mass*v(1,i)*v(3,i) 
       if (x(3,i) .gt. 0.8_dp .and. x(3,i) .lt. 2.5_dp) then   
          z(i) = x(3,i)
          vx(i) = v(1,i)
       else
          z(i) = 0.0_dp
          vx(i) = 0.00001_dp                 
       endif
    enddo
    Pxz = Pxz + virial(1,3)
    const = Pxz/Vol
    write(*,*) 'Media velocità', vel/Natoms
    !viscosity = -Pxz/(maxval(vx) - minval(vx))/ (z(maxloc(vx))-z(minloc(vx)))
!!$    _________________________________________________________________________
!!$    Because the flow profile is linear, the gradient is constant.
!!$    Linear interpolation in the central part provides the gamma coefficient,
!!$    necessary in order to calculate the coefficient of shear viscosity
!!$    _________________________________________________________________________

    write(*,*) 'gamma is',  (maxval(vx) - minval(vx))/(z(maxloc(vx))-z(minloc(vx)))
    write(*,*) maxval(vx), minval(vx),z(maxloc(vx)), z(minloc(vx))
    write(*,*) 'pressure is', const

!!$    ____________________________________________________________________________
!!$    The coefficient of shear viscosity is calculated as Pxz/gamma.
!!$    Pxz must be calculated from equation 13.3.9 once the equilibrium is reached.
!!$    ____________________________________________________________________________

    dv = maxval(vx) - minval(vx)
    mini = minloc(vx,1)
    maxi = maxloc(vx,1) 
    dz = z(maxi)-z(mini)
    eta =  -const/(dv/dz)
    write(*,*) 'Shear viscosity = ', eta*M2F*eV2J*1.0d-15/1.0d-27*1.0d2,'mP' 
    write(*,*) 'Reynold number is', density*(vel/Natoms)/eta

    deallocate(z,vx)
    end subroutine reynolds

    
end module dynamics
