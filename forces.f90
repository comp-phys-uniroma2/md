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
     Fc = 24.0_dp*rm2*(2.0_dp*rm12-rm6)
     Uc = 4.0_dp*(rm12-rm6)

     ! Compute LJ forces and energy at Rmin= sigma/2
     ra2 = sg2/4.0_dp
     rm2 = sg2/ra2
     rm6 = rm2*rm2*rm2
     rm12 = rm6*rm6
     Fa = 24.0_dp*rm2*(2.0_dp*rm12-rm6)
     Ua = 4.0_dp*(rm12-rm6)

     write(*,*) 'Uc=',Uc*par%eps,'Fc=',Fc*par%eps
     write(*,*) 'Um=',Ua*par%eps,'Fm=',Fa*par%eps
     
  end subroutine init_lj

  subroutine lj(x,F,UU,virial)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), intent(out) :: F 
    real(dp), intent(out) :: UU
    real(dp), intent(out) :: virial

    real(dp) :: rij(3), g(3), r2, rm2, rm6, rm12, tmp, Fx, Fy, Fz
    integer :: ii, jj, kk, ci, cj, ck, u,v,w
    integer :: m, l, Natoms
    type(TNode), pointer :: it   

    ! Virial should be corrected due to cutoff potential

    UU = 0.0_dp
    virial = 0.0_dp
    Natoms = par%Natoms
    !$OMP PARALLEL DO DEFAULT(PRIVATE), SHARED(map,boxlists,Fa,Fc,Ua,Uc,ra2,rc2,sg2,x,F) & 
    !$OMP&   REDUCTION( + : UU, virial) 
    do m = 1, Natoms

       ! cerca la scatola ci,cj,ck di m 
       call boxind(x(:,m),ci,cj,ck)          

       Fx = 0.0_dp; Fy = 0.0_dp; Fz = 0.0_dp

       ! CALCOLA FORZE SU PARTICELLA m A PARTIRE DALLE
       ! PARTICELLE NEI BOX INTORNO E NELLO STESSO.

       ! !$OMP PARALLEL DO DEFAULT(PRIVATE), SHARED(map,boxlists,m,Fa,Fc,Ua,Uc,ra2,rc2,sg2,ci,cj,ck,x) & 
       ! !$OMP&   REDUCTION( + : Fx, Fy, Fz, UU, virial) 
       do u = 1, 27

         ii = ci + map(1,u)
         jj = cj + map(2,u)
         kk = ck + map(3,u)

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
                Fx = Fx - (Fa-Fc)*rij(1)
                Fy = Fy - (Fa-Fc)*rij(2)
                Fz = Fz - (Fa-Fc)*rij(3)
                UU = UU + (Ua-Uc)
                virial = virial + rij(1)*Fx+rij(2)*Fy+rij(3)*Fz
                it => it%next
                cycle
             endif
   
             if (r2 .ge. rc2) then
             else
                rm2 = sg2/r2
                rm6 = rm2*rm2*rm2
                rm12 = rm6*rm6

                tmp = 24.0_dp*rm2*(2.0_dp*rm12-rm6)     
                Fx = Fx - (tmp-Fc)*rij(1)
                Fy = Fy - (tmp-Fc)*rij(2)
                Fz = Fz - (tmp-Fc)*rij(3)
                UU = UU + (4.0_dp*(rm12-rm6) - Uc)  
                virial = virial + rij(1)*Fx+rij(2)*Fy+rij(3)*Fz
             endif
             
             it => it%next

          end do

        end do 
        ! !$OMP END PARALLEL DO
        F(1,m) = Fx
        F(2,m) = Fy
        F(3,m) = Fz
     end do
     !$OMP END PARALLEL DO

     F = F * par%eps/sg2
     UU = UU * 0.5_dp * par%eps
     virial = virial * 0.5_dp * par%eps/sg2
      
  end subroutine lj

end module forces
