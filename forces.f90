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

  subroutine lj(x,dx,F,lF,UU,virial)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:,:), intent(in) :: dx
    real(dp), dimension(:,:), intent(out) :: F
    real(dp), dimension(:,:,:), intent(out) :: lF
    
    real(dp), intent(out) :: UU
    real(dp), intent(out) :: virial

    real(dp) :: rij(3), g(3), r2, rm1, rm2, rm6, rm12, tmp
    real(dp) :: drij(3), Fm(3), ltmp
    integer :: ii, jj, kk, ci, cj, ck, u,v,w,h
    integer :: m, l, Natoms
    type(TNode), pointer :: it

    Natoms = par%Natoms
    
    ! Virial should be corrected due to cutoff potential

    UU = 0.0_dp
    virial = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(PRIVATE), SHARED(map,boxlists,Fa,Fc,Ua,Uc,ra2,rc2,sg2,x,F,lF) &
    !$OMP&   REDUCTION( + : UU, virial)
    do m = 1, Natoms

       ! cerca la scatola ci,cj,ck di m
       call boxind(x(:,m),ci,cj,ck)

       Fm = 0.d0
       lF(:,m,:)=0.d0
  
       do u = 1, 27
               
         ii = ci + u
         jj = cj + v
         kk = ck + w
         ! controlla se la scatola una copia periodica
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

             rij(:) = x(:,l)+g(:) - x(:,m)
              
             r2 = dot_product(rij,rij)
               
             ! Questo non si verifica mai
             !if (r2 .le. ra2) then
             !   Fm(:) = Fm(:) - (Fa-Fc)*rij(:)
             !   UU = UU + (Ua-Uc)
             !   virial = virial + dot_product(rij,Fm)
             !   it => it%next
             !   cycle
             !endif

             if (r2 .ge. rc2) then
             else
               rm1 = 1.0_dp/sqrt(r2)
               rm2 = sg2/r2
               rm6 = rm2*rm2*rm2
               rm12 = rm6*rm6

               Fm(:) = Fm(:) - (tmp-Fc) * rij(:)
               UU = UU + (4.0_dp*(rm12-rm6)-Uc)  
               virial = virial + dot_product(rij,Fm(:))

               tmp = 24.0_dp*rm2*(2.0_dp*rm12-rm6)
               ltmp = 24.0_dp*(8.0_dp*rm6-28.0_dp*rm12)*rm2*rm1
               do h=1,6*Natoms
                 drij(:) = dx(:,l,h) - dx(:,m,h)
                 lF(:,m,h)=lF(:,m,h)-tmp*drij(:)-ltmp*rij(:)*dot_product(drij(:),rij(:))/sqrt(r2)
               end do
          
             endif
             
             it => it%next

          end do

        end do           
        F(:,m) = Fm(:)          
     end do
     !$OMP END PARALLEL DO

     lF = lF * par%eps/sg2
     F = F * par%eps/sg2
     UU = UU * 0.5_dp * par%eps
     virial = virial * 0.5_dp * par%eps/sg2
      
  end subroutine lj

end module forces
