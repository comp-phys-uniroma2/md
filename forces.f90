module forces 
  use constants
  use list
  use boxes
  implicit none
  private

  public :: lj
  public :: init_lj

  type Tpar
     integer :: Natoms
     real(dp) :: eps
     real(dp) :: sigma
     real(dp) :: Rc
  end type

  type(Tpar) :: par 
  real(dp) :: sg2, ra2, rc2, Fc, Fa, Uc, Ua

  contains

  subroutine init_lj(Natoms,eps,sigma,Rc)
     integer :: Natoms
     real(dp) :: eps, sigma, Rc
     real(dp) :: rm2, rm6, rm12

     par%Natoms = Natoms
     par%eps = eps
     par%sigma = sigma
     par%Rc = Rc

     ! Compute LJ forces and energy at cutoff
     sg2 = par%sigma*par%sigma
     rc2 = par%Rc*par%Rc
     rm2 = sg2/rc2
     rm6 = rm2*rm2*rm2
     rm12 = rm6*rm6
     Fc = 24.d0*rm2*(2.d0*rm12-rm6)
     Uc = 4.d0*(rm12-rm6)

     ! Compute LJ forces and energy at Rmin= sigma/2
     ra2 = sg2/4.d0
     rm2 = sg2/ra2
     rm6 = rm2*rm2*rm2
     rm12 = rm6*rm6
     Fa = 24.d0*rm2*(2.d0*rm12-rm6)
     Ua = 4.d0*(rm12-rm6)

     write(*,*) 'Uc=',Uc*par%eps,'Fc=',Fc*par%eps
     write(*,*) 'Um=',Ua*par%eps,'Fm=',Fa*par%eps
     
  end subroutine init_lj

  subroutine lj(x,F,UU,virial)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), intent(out) :: F 
    real(dp), intent(out) :: UU
    real(dp), intent(out) :: virial

    real(dp) :: rij(3), g(3), r2, rm2, rm6, rm12, tmp
    integer :: ii, jj, kk, ci, cj, ck, u,v,w
    integer :: m, l
    type(TNode), pointer :: it   

    ! Virial should be corrected due to cutoff potential

    UU = 0.d0
    virial = 0.d0
    do m = 1, par%Natoms

       ! cerca la scatola ci,cj,ck di i
       call boxind(x(:,m),ci,cj,ck)          

       F(:,m) = 0.d0

       do w=-1, +1 ! CALCOLA FORZE SU PARTICELLA m A PARTIRE DALLE
                   ! PARTICELLE NEI BOX INTORNO E NELLO STESSO.
         do v=-1, +1
           do u=-1, +1

              ii = ci + u
              jj = cj + v
              kk = ck + w
              ! controlla se la scatola IN QUESTION È una copia periodica
              ! ED IN CASO PRENDE I DATI DA QUELLA ORIGINALE
              ! g e' vettore supercella
              call folding(ii,jj,kk,g)

              ! Iterates over atoms in box (ii,jj,kk)
              it => boxlists(ii,jj,kk)%start
         
              do while (associated(it))     
             
                  l = it%val

                  if (l .eq. m) then 
                      it => it%next
                      cycle
                  endif

                  ! segno corretto rij = rj - ri 
                  rij(:) = x(:,l)+g(:) - x(:,m)

                  r2 = dot_product(rij,rij)
                  
                  if (r2 .le. ra2) then
                     F(:,m) = F(:,m) - (Fa-Fc)*rij
                     UU = UU + (Ua-Uc)
                     virial = virial + dot_product(rij,F(:,m))
                     it => it%next
                     cycle
                  endif
   
                  if (r2 .ge. rc2) then
                  else
                     rm2 = sg2/r2
                     rm6 = rm2*rm2*rm2
                     rm12 = rm6*rm6

                     tmp = 24.d0*rm2*(2.d0*rm12-rm6)     
                     F(:,m) = F(:,m) - (tmp-Fc) * rij 
                     UU = UU + (4.d0*(rm12-rm6) - Uc)  
                     virial = virial + dot_product(rij,F(:,m))
                  endif
                  
                  it => it%next

              end do

            enddo
          enddo
        enddo 
        F(:,m) = F(:,m) * par%eps/sg2
     end do

     UU = UU * 0.5d0 * par%eps
     virial = virial * 0.5d0 * par%eps/sg2
      
     !print*,'F(1)=',F(:,1)
     !print*,F(:,200)

  end subroutine lj

end module forces
