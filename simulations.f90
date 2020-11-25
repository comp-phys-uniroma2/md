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
  public :: init_particles 
  public :: nve_sim
  public :: write_coords
  public :: compute_g

  contains

  ! ---------------------------------------------------------
  ! Initialize position and velocities of particles
  ! according to the MB distribution
  subroutine init_particles()
    integer :: i,err  
    real(dp) :: rndr, ax, ay, az, aa(3)
    real(dp) :: bb(3,4)
    integer :: ii,jj,kk
    integer :: l,k,m

    if (4*Nx*Ny*Nz .ne. Natoms) then
       STOP 'Error: Nx*Ny*Nz != Natoms'
    endif

    ax = Lx/Nx; ay = Ly/Ny; az = Lz/Nz 
    bb(:,1) = (/        ax/4.0_dp,        ay/4.0_dp, az/4.0_dp /) 
    bb(:,2) = (/ 3.0_dp*ax/4.0_dp, 3.0_dp*ay/4.0_dp, az/4.0_dp /)    
    bb(:,3) = (/ ax/4.0_dp, 3.0_dp*ay/4.0_dp, 3.0_dp*az/4.0_dp /) 
    bb(:,4) = (/ 3.0_dp*ax/4.0_dp, ay/4.0_dp, 3.0_dp*az/4.0_dp /) 

    i = 1

    do m = 1, Nz
      do k = 1, Ny
         do l= 1, Nx

            aa(1) = (l-1)*Lx/Nx
            aa(2) = (k-1)*Ly/Ny
            aa(3) = (m-1)*Lz/Nz
            ! Initialize positions in the box
            x(:,i) = aa(:) + bb(:,1)      
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,2)      
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,3)      
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,4)      
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

         enddo
      enddo
    enddo

    call init_seed()

    do i = 1, Natoms

      !print*,'set velocities'
      ! Initialize random MB velocities.
      v(1,i) = maxwell_boltzmann()
      v(2,i) = maxwell_boltzmann()
      v(3,i) = maxwell_boltzmann()

    end do  

    !open(101,file='data/positions.dat')
    !call write_coords(101)
    !close(101)  

  end subroutine init_particles

  ! ---------------------------------------------------------
  ! In ogni direzione P = (m/2PikT)^1/2 exp(-m v^2 / 2kT)
  ! P0=(m/2PikT)^1/2 
  ! P/P0 = exp[-m v^2/ 2kT] 
  ! => v^2 = +/- sqrt( - m/2kT * ln(P/P0) ) 
  function maxwell_boltzmann() result(v)
    real(dp) :: v
    real(dp) :: rndr, b, eps

    ! Define machine precision 
    ! This is used to change the interval [0,1) into (0,1]
    eps = mach()  

    call random_number(rndr)
    b= cos(2.0_dp*Pi*rndr)
    do 
       call random_number(rndr)
       
       v = b*sqrt( - 2.0_dp * kb * Temp/Mass * log(rndr+eps) )
       exit
    enddo

  end function maxwell_boltzmann 
 
  ! ---------------------------------------------------------
  subroutine init_seed()
    integer, dimension(:),allocatable :: seed
    integer :: j

    call random_seed(size=j)
    allocate(seed(j))
    call system_clock(j)
    seed=j
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_seed


  ! ---------------------------------------------------------
  ! Perform an NVE simulation
  subroutine nve_sim()
    real(dp), dimension(:,:), allocatable :: x1,v1   
    
    integer :: nstep1,nstep2, err
    integer :: n
    character(3) :: ind
    real(dp) :: U      ! Potential energy
    real(dp) :: K      ! Kinetic energy
    real(dp) :: virial ! virial = - Sum_i Sum_j (r_ij * F_ij)
    real(dp) :: P      ! Pressure (from virial)
    real(dp) :: lambda ! scaling parameter 
    real(dp) :: Kav,Uav,Pav ! Averaged quantities
    real(dp) :: R(3),R0(3)  ! For diffusivity 
    real(dp), allocatable :: dR(:)  ! For diffusivity 
    real(dp) :: R2, R02, Diff, Rcm(3)

    Pav = 0.d0
    Uav = 0.d0
    Kav = 0.d0
    Diff = 0.d0 
    Rcm = 0.d0

    allocate(x1(3,Natoms),stat=err)
    allocate(v1(3,Natoms),stat=err)
    allocate(dR(3*Natoms),stat=err)
    if (err /= 0) STOP 'ALLOCATION ERROR'
    
    nstep1=nint(tinit/dt)
    nstep2=nint(tsim/dt)

    call init_lj(Natoms,eps,sigma,Rc)
    
    write(*,*) 'Target T=',Temp,'K=',(3.d0*Natoms*kb*Temp/2.d0)

    call set_clock()
    
    if (print_xyz) then
      open(101,file='data/coord.xyz')
      write(101,'(i0)') Natoms
      write(101,*) 'Frame',0
      call write_xyz(101)
    end if  
    open(102,file='data/R2.dat')

    ! init time
    do n=1,nstep1
       
       if (print_xyz .and. mod(n,xyz_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n 
          call write_xyz(101)
       endif

       !write(*,*) 'propagate verlet'   
       call verlet(x,v,x1,v1,U,virial,dt,lj)

       K = kinetic()  
       P = (Natoms*kb*Temp + virial/3.d0)/Vol

       if (mod(n,xyz_interval) == 0) then
         write(*,'(i6,a,i6,3x,3(a3,ES14.6,2x))') n,'/',nstep1,'Ek=',K,'U=',U,'P=',P 
       end if

       call update_boxes(x,x1)

       x = x1

       if (scaling) then 
          lambda = sqrt((3.d0*Natoms*kb*Temp/2.d0)/K)
          v = v1 * lambda
       else
          v = v1
       endif      

    end do 
        
    ! simulation time
    R0 = x(:,Natoms/2)
    R02 = dot_product(R0,R0)
    dR = 0.0_dp

    do n=1,nstep2
       
       if (print_xyz .and. mod(n,xyz_interval) == 0) then
         write(101,*) Natoms
         write(101,*) 'Frame',n 
         call write_xyz(101)
       endif

       call verlet(x,v,x1,v1,U,virial,dt,lj)

       K = kinetic()  
       P = (Natoms*kb*Temp + virial/3.d0)/Vol

       if (mod(n,xyz_interval) == 0) then
         write(*,'(i6,a,i6,3x,3(a3,ES14.6,2x))') n,'/',nstep2,'Ek=',K,'U=',U,'P=',P 
       end if

       call update_boxes(x,x1)

       x = x1

       dR = dR + reshape(v ,(/ 3*Natoms /)) 
       dR = dR + reshape(v1,(/ 3*Natoms /)) 
       R2 = dot_product(dR,dR)*dt*dt/Natoms/4.0_dp
       write(102,*) n*dt, R2
       
       if (scaling) then 
          lambda = sqrt((3.d0*Natoms*kb*Temp/2.d0)/K)
          v = v1 * lambda
       else
          v = v1
       endif      
   
       Pav = Pav + P/nstep2
       Kav = Kav + K/nstep2
       Uav = Uav + U/nstep2
       R = x(:,Natoms/2) - R0
       Rcm = Rcm + R/nstep2
       R2 = dot_product(R,R) 
       Diff = Diff + R2/nstep2

    end do

    if (print_xyz) then
      close(101)
    end if
    close(102)

    write(*,*)
    write(*,'(a16,3x)',advance='NO') 'Simulation time:'
    call write_clock()

    write(*,*)
    write(*,*) 'Averaged quantities:'
    write(*,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
    write(*,*) 'Rcm=',Rcm
    write(*,*) 'Diffusivity=',Diff/dt,'nm^2/fs'

    deallocate(x1,v1)  

    open(102,file='data/av.dat')
    write(102,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
    write(102,*) 'Rcm=',Rcm
    write(102,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
    close(102)

  end subroutine nve_sim
 
  function kinetic() result(Ek)
     real(dp) :: Ek
     integer :: n
     
     Ek =0.d0
     do n = 1, Natoms
        Ek = Ek + dot_product(v(:,n),v(:,n))   
     end do
     Ek = Ek*Mass*0.5d0 

  end function kinetic

  subroutine write_xyz(id)
    integer :: id

    integer :: n, ii, jj, kk
 
    do n = 1, Natoms      
       !call boxind(x(:,n),ii,jj,kk)
       write(id,'(a,3(f12.6),3(ES20.8))') 'He ', x(:,n)*10.d0  !,v(:,n)
    end do

  end subroutine write_xyz

  subroutine write_coords(id)
    integer :: id

    integer :: n
 
    write(id,'(1X,L2,I7,3E23.15)') .true.,Natoms,Lx/sigma,Ly/sigma,Lz/sigma
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


  
  subroutine compute_g() ! CALCOLA FUNZIONE DI CORRELAZIONE RADIALE
    real(dp) :: rij(3), g(3), r, a, b 
    integer :: ii, jj, kk, ci, cj, ck, u,v,w
    integer :: m, l
    type(TNode), pointer :: it  
     
    real(dp),dimension(:), allocatable :: gg
    integer :: Nk, basket, err
    
    Nk = aint(Rc/dr)
    allocate(gg(Nk),stat=err)
    if(err.ne.0) STOP 'ALLOCATION ERROR gg'

    gg = 0.d0

    do m = 1, Natoms

       ! cerca la scatola ci,cj,ck di i
       call boxind(x(:,m),ci,cj,ck)          
       
       do w=-1,1
         do v=-1,1 
           do u=-1,1
              ii = ci + u
              jj = cj + v     
              kk = ck + w  
                    
              ! g e' vettore supercella
              call folding(ii,jj,kk,g)

              ! Iterates over atoms in box (ii,jj,kk)
              it => boxlists(ii,jj,kk)%start
         
              do while (associated(it))     
             
                  l = it%val

                  if (l .eq. m) then ! evita il caso l=m (stessa particella)
                      it => it%next
                      cycle ! ricomincia il ciclo do dall'inizio
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
        enddo

     end do
      

     open(101,file='data/g.dat')
     do m = 1, Nk
        r = m*dr
        write(101,*) r, gg(m)*Vol/(4.d0*Natoms*Natoms*Pi*r*r*dr)  
     enddo
 
     close(101)

  end subroutine compute_g


end module simulations

