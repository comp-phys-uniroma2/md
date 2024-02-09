module simulations

  use constants 
  use parameters
  use list
  use boxes
  use forces
  use dynamics
  use clock
  implicit none
  private

  public :: transform_units
  public :: init_positions_fcc
  public :: init_velocities
  public :: init_velocities_couette
  public :: init_seed
  public :: nve_sim
  public :: write_coords
  public :: compute_g
  

  contains

  ! ---------------------------------------------------------
  ! Initialize positions of particles as in fcc lattice
  subroutine init_positions_fcc()
    integer :: i,err
    real(dp) :: rndr, ax, ay, aa(2)
    real(dp) :: bb(2,2)
    integer :: ii,jj
    integer :: l,k

    
    ax = Lx/Nx; ay = Ly/Ny

    bb(:,1) = (/        ax/4.0_dp,        ay/4.0_dp       /)
    bb(:,2) = (/ 3.0_dp*ax/4.0_dp,  3.0_dp*ay/4.0_dp        /)

    i = 1

      do k = 1, Ny
         do l= 1, Nx

            aa(1) = (l-1)*Lx/Nx
            aa(2) = (k-1)*Ly/Ny
            ! Initialize positions in the box
            x(:,i) = aa(:) + bb(:,1)
            call boxind(x(:,i),ii,jj)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj
            call add(boxlists(ii,jj),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,2)
            call boxind(x(:,i),ii,jj)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj
            call add(boxlists(ii,jj),i)
            i = i + 1

         enddo
      enddo

  end subroutine init_positions_fcc

  ! ---------------------------------------------------------
  ! according to the MB distribution
  subroutine init_velocities()

    integer :: i

    do i = 1, Natoms

      !print*,'set velocities'
      ! Initialize random MB velocities.
      v(1,i) = maxwell_boltzmann()
      v(2,i) = maxwell_boltzmann()
    end do

  end subroutine init_velocities

  ! ---------------------------------------------------------
  ! Init velocities for the couette BC
subroutine init_velocities_couette()

    integer :: i
    real(dp) :: upperlayers, lowerlayers

    upperlayers = maxval(x(2,:)) - Ly/Ny !first two layers of the fluid
    lowerlayers = minval(x(2,:)) + Ly/Ny

    do i = 1, Natoms

       if(x(2,i) .ge. upperlayers) then
          v(1,i) = v_drift
          v(2,i) = 0.d0

        elseif(x(2,i) .le. lowerlayers) then
          v(1,i) = 0.d0
          v(2,i) = 0.d0          

       else
          v(1,i) = maxwell_boltzmann()
          v(2,i) = maxwell_boltzmann()
       end if
    end do
    print*,'init velocities couette:'
    print*,'vmax=',maxval(v(1,:)),maxval(v(2,:))
    print*,'vmin=',minval(v(1,:)),minval(v(2,:))
    print*,'done'

  end subroutine init_velocities_couette

    

  
  ! ---------------------------------------------------------
  ! In ogni direzione P = (m/2PikT)^1/2 exp(-m v^2 / 2kT)
  ! P0=(m/2PikT)^1/2
  ! P/P0 = exp[-m v^2/ 2kT]
  ! => v^2 = +/- sqrt( - m/2kT * ln(P/P0) )
  function maxwell_boltzmann() result(vv)
    real(dp) :: vv
    real(dp) :: rndr, b, eps

    ! Define machine precision
    ! This is used to change the interval [0,1) into (0,1]
    eps = mach()

    call random_number(rndr)
    b= cos(2.0_dp*Pi*rndr)
    do
       call random_number(rndr)

       vv = b*sqrt( - kb * Temp/Mass * log(rndr+eps) )
       exit
    enddo

  end function maxwell_boltzmann

  ! ---------------------------------------------------------
  subroutine init_seed(seed_in)
    integer, intent(in), optional :: seed_in

    integer, dimension(:),allocatable :: seed
    integer :: j

    call random_seed(size=j)
    allocate(seed(j))

    if (present(seed_in)) then
      seed=seed_in
    else
      call system_clock(j)
      seed=j
    end if
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_seed


  ! ---------------------------------------------------------
  ! Perform an NVE simulation
  subroutine nve_sim()
    real(dp), dimension(:,:), allocatable :: xf,vf
    real(dp), dimension(:,:,:), allocatable :: etaf

    real(dp), dimension(:),allocatable :: kine
    
    integer :: nstep1, nstep2, err
    integer :: n, iter, pq, i, j
    character(3) :: ind

    real(dp) :: U           ! Potential energy
    real(dp) :: K           ! Kinetic energy
    real(dp) :: K0          ! kin ?
    real(dp) :: virial      ! virial = - Sum_i Sum_j (r_ij * F_ij)
    real(dp) :: P           ! Pressure (from virial)
    real(dp), allocatable :: Sxy1(:)        ! Shear stress warm up phase
    real(dp), allocatable :: Sxy2(:)        ! Shear stress simulation phase
    real(dp) :: lambda      ! scaling parameter
    real(dp) :: Kav,Uav,Pav ! Averaged quantities
    real(dp) :: R(2),R0(2)  ! For diffusivity 
    real(dp), allocatable :: dR(:)  ! For diffusivity 
    real(dp) :: R2, R02, Diff, Rcm(2)
    real(dp) :: viscosity, xpt1, xpt

    real(dp), dimension(:), allocatable :: viscosityHelft

    Pav = 0.0_dp
    Uav = 0.0_dp
    Kav = 0.0_dp
   
    Diff= 0.0_dp 
    Rcm = 0.0_dp

    nstep1=nint(tinit/dt)
    nstep2=nint(tsim/dt)


    allocate(dR(2*Natoms),stat=err)
    allocate(etaf(2,Natoms,5),stat=err)
    allocate(Sxy1(nstep1),stat=err)
    allocate(Sxy2(nstep2),stat=err)
    allocate(viscosityHelft(nstep2), stat=err)
    
    allocate(xf(2,Natoms),stat=err)
    allocate(vf(2,Natoms),stat=err)


    allocate(kine(nstep2),stat=err)

    
    if (err /= 0) STOP 'ALLOCATION ERROR'
    
    call init_lj(Natoms,eps,sigma,Rc)
    
    call set_clock()
    !write(ind,'(i3.3)') n
    !open(101,file='coord'//ind//'.xyz')

    if (print_xyz) then
      open(101,file='data/coords.xyz')
      write(101,'(i0)') Natoms
      write(101,*) 'Frame',0
      call write_xyz(101)
    end if
    open(102,file='data/R2.dat')
    open(103,file='data/kin.dat')
    open(113,file='data/e_tot.dat')
    open(123,file='data/U.dat')
    open(133,file='data/P.dat')
    open(134,file='data/Sxy.dat')
    open(135,file='data/viscosityhelft.dat')
    

    ! init time
    write(*,*) 'Warm up phase:',nstep1,'steps'
    ! Target mean Kinetic energy:
    K0=(Natoms*kb*Temp)
    write(*,*) 'Target T=',Temp,'K=',K0


    do n=1,nstep1

       if (print_xyz .and. mod(n,print_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif

       if (nose_hoover) then
         call verlet_nh15(x,v,U,virial,Sxy1(n),xpt,dt,lj,K,K0,eta,etaf,xf,vf)
       else  
         call verlet(x,v,xf,vf,U,virial,Sxy1(n),xpt,dt,lj,K)
       end if

       call update_boxes(x,xf)


       P = (Natoms*kb*Temp + virial/2.d0)/Area

       if (mod(n,print_interval) == 0) then
         write(*,'(i6,a,i6,3x,4(a3,ES14.6,2x))') n,'/',nstep1,'Ek=',K,'U=',U,'E=',K+U,'P=',P
       end if

       if (scaling) then
          lambda = sqrt((Natoms*kb*Temp)/K)
          v = vf * lambda
       else 
          v = vf
       endif
       x = xf
       eta = etaf

       
    end do

    ! simulation time
    R0 = x(:,Natoms/2)
    R02 = dot_product(R0,R0)
    dR = 0.0_dp   
   
    write(*,*) '**********************************************************************' 
    write(*,*) 'Simulation phase:',nstep2,'steps'
    ! init time

    do n=1, nstep2
     
       if (print_xyz .and. mod(n,print_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif
       
       if (nose_hoover) then
         call verlet_nh15(x,v,U,virial,Sxy2(n),xpt,dt,lj,K,K0,eta,etaf,xf,vf)
       else  
         call verlet(x,v,xf,vf,U,virial,Sxy2(n),xpt,dt,lj,K)
       end if

       

       call update_boxes(x,xf)
 
       P = (Natoms*kb*Temp + virial/2.d0)/Area
      
        if (n==1) then
          xpt1=xpt
        end if
       viscosityHelft(n) = (((xpt - xpt1)/Natoms)**2)/(2*kb*Temp*Area)

       if (mod(n,print_interval) == 0) then
         write(*,'(i6,a,i6,3x,4(a3,ES14.6,2x))') n,'/',nstep2,'Ek=',K,'U=',U,'E=',K+U,'P=',P
       end if 
       !write(103,*) n*dt, K
       !write(113,*) n*dt, U+K
       !write(123,*) n*dt, U
       !write(133,*) n*dt, P
       write(134,*) n*dt, Sxy2(n)
       write(135,*) n*dt, viscosityHelft(n)
       kine(n)=K
       Uav = Uav + U/nstep2





       if (scaling) then 
          lambda = sqrt((Natoms*kb*Temp)/K)
          v = vf * lambda
       else 
          v = vf
       endif

       x = xf     
       eta=etaf
      
      dR = dR + reshape(v ,(/ 2*Natoms /))
      dR = dR + reshape(vf,(/ 2*Natoms /))
      R2 = dot_product(dR,dR)*dt*dt/Natoms/4.0_dp
      write(102,*) n*dt, R2

      Pav = Pav + P/nstep2
      Kav = Kav + K/nstep2
      Uav = Uav + U/nstep2
      R = x(:,Natoms/2) - R0
      Rcm = Rcm + R/nstep2
      R2 = dot_product(R,R)
      Diff = Diff + R2/nstep2
   
    end do


    close(101)
    close(102)
    close(103)
    close(113)
    close(123)
    close(133)
    close(104)
    close(105)
    close(115)
    close(125)
    close(106)
    close(134)
    close(135)

    viscosity = gk_viscosity(Sxy2, dt, nstep2)*Area/(kb*Temp)   ! Green Kubo viscosity

    write(*,*)
    write(*,*) 'K0 =', K0
    write(*,*) '<K> =' , sum(kine)/nstep2
    write(*,*) '|<K> - K0| =', abs(K0-sum(kine)/nstep2)
    write(*,*) 'sqrt(<K^2>-<K>^2) =',sqrt( (sum((kine-(sum(kine)/nstep2))**2))/nstep2)
    write(*,*) '<U> =', Uav
    write(*,*) 'Shear Viscosity Coefficient =', viscosity
    write(*,'(a16,3x)',advance='NO') 'Simulation time:'
    call write_clock()


    deallocate(xf,vf,etaf,viscosityhelft,Sxy1,Sxy2,dR)

    open(102,file='data/av.dat')
    write(102,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav
    write(102,*) 'Rcm=',Rcm
    write(102,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
    write(102,*) 'Shear Viscosity Coefficient=',viscosity, 'Pa *fs^2/nm'
    close(102)

  end subroutine nve_sim


  function kinetic() result(Ek)
     real(dp) :: Ek
     integer :: n
     
     Ek =0.0_dp
     do n = 1, Natoms
        Ek = Ek + dot_product(v(:,n),v(:,n))
     end do
     Ek = Ek*Mass*0.5_dp 

  end function kinetic

  function gk_viscosity(Sxy, dt, nstep) result(viscosity)

    integer, intent (in) :: nstep
    integer i, j
    real(dp), intent (in) :: Sxy(nstep),dt
    real(dp) :: viscosity, sum, sacf(nstep)
    
    viscosity = 0.0_dp

    do j=1,nstep-1
    sum =0.0_dp
      do i=1, nstep-j
        sum = sum + Sxy(i+j)*Sxy(i)
      end do
      sacf(j) = sum / (nstep - j)  !stress autocorrelation function
    end do


    do i = 1, nstep - 1
        viscosity = viscosity + (sacf(i) + sacf(i + 1)) * dt / 2.0
    end do


    viscosity = viscosity/nstep

  end function gk_viscosity

  subroutine write_xyz(id)
    integer :: id

    integer :: n, ii, jj
 
    do n = 1, Natoms      
       !call boxind(x(:,n),ii,jj)
      write(id,'(a,4(f12.6),3(ES20.8))') 'He  ', x(:,n)*10.0_dp, 0.0_dp  !,v(:,n)
    enddo

  end subroutine write_xyz

  subroutine write_coords(id)
    integer :: id

    integer :: n
 
    write(id,'(1X,L2,I7,3E23.15)') .true.,Natoms,Lx/sigma,Ly/sigma
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') x(:,n)/sigma
    enddo
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') v(:,n)*sqrt(Mass/eps)
    enddo
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') v(:,n)/dt*(sigma*Mass/eps)
    enddo

  end subroutine write_coords


  
  subroutine compute_g()
    real(dp) :: rij(2), g(2), r, a, b 
    integer :: ii, jj, ci, cj, u,v
    integer :: m, l
    type(TNode), pointer :: it  
     
    real(dp),dimension(:), allocatable :: gg
    integer :: Nk, basket, err
    
    Nk = aint(Rc/dr)
    allocate(gg(Nk),stat=err)
    if(err.ne.0) STOP 'ALLOCATION ERROR gg'

    gg = 0.0_dp

    do m = 1, Natoms

       ! cerca la scatola ci,cj,ck di i
       call boxind(x(:,m),ci,cj)          
    
         do v=-1,1 
           do u=-1,1
              ii = ci + u
              jj = cj + v
                    
              ! g e' vettore supercella
              call folding(ii,jj,g)

              ! Iterates over atoms in box (ii,jj)
              it => boxlists(ii,jj)%start
         
              do while (associated(it))     
             
                  l = it%val

                  if (l .eq. m) then 
                      it => it%next
                      cycle
                  endif

                  ! segno corretto rij = rj - ri
                  rij(:) = x(:,l)-x(:,m)+g(:)
           
                  r = sqrt(dot_product(rij,rij))
                
                  if (r<Rc) then
                    basket = aint(r/dr) + 1
                    gg(basket) = gg(basket) + 1
                  endif

                  it => it%next
                
              end do

            enddo
          enddo


     end do
      

     open(101,file='data/g.dat')
     do m = 1, Nk
        r = m*dr
        write(101,*) r, gg(m)*Area/(2.d0*Natoms*Natoms*Pi*r*r*dr) 
     enddo
 
     close(101)

   end subroutine compute_g




end module simulations

