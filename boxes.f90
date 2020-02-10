module boxes
  use constants 
  use list
  implicit none
  private

  public :: boxlists
  public :: create_boxes
  public :: destroy_boxes
  public :: update_boxes
  public :: boxind
  public :: boxinfo
  public :: folding

   
  type(Tlist), dimension(:,:,:), allocatable :: boxlists

  ! PRIVATE STUFF:
  integer :: Nx, Ny, Nz
  real(dp) :: lx, ly, lz
  real(dp) :: LLx, LLy, LLz

  contains

  subroutine create_boxes(x,y,z,cutoff)
    real(dp) :: x, y, z, cutoff
    integer :: err

   
    LLx=x
    LLy=y
    LLz=z

    Nx = floor(LLx/cutoff) ! CON IL FLOOR MI ASSICURO CHE 
    Ny = floor(LLy/cutoff) ! CUTOFF SIA MAGGIORE DI lx ...
    Nz = floor(LLz/cutoff)

    lx = LLx/Nx
    ly = LLy/Ny
    lz = LLz/Nz

    allocate(boxlists(Nx,Ny,Nz), stat=err)     
  
    if (err.ne.0) STOP 'ERROR ALLOCATION Boxlist'

  end subroutine create_boxes

  ! ---------------------------------------------------- 
  subroutine destroy_boxes()
     integer i,j,k

     do i=1,Nx
      do j=1,Ny
       do k=1,Nz
        !print*,'box:',i,j,k 
        call delete_list(boxlists(i,j,k))
       enddo
      enddo
     enddo

     !print*,'remove boxlists'
     deallocate(boxlists)     
  end subroutine destroy_boxes 
  ! ---------------------------------------------------- 
  
  subroutine update_boxes(x,x1)
    real(dp), dimension(:,:) :: x,x1

    integer :: i, natoms, ii,jj,kk, ii1,jj1,kk1
    integer :: ii2,jj2,kk2
    type(TNode), pointer :: node
    real(dp) :: g(3)

    natoms = size(x,2)

    do i = 1, natoms
       call boxind(x(:,i),ii,jj,kk)  ! IN CHE BOX SI TROVA X
       call boxind(x1(:,i),ii1,jj1,kk1) ! IN CHE BOX SI TROVA X1

       if (ii.ne.ii1 .or. jj.ne.jj1 .or. kk.ne.kk1) then ! ESEGUO SOLO SE
                                                         ! ESCO AD UNA BOX
          call folding(ii1,jj1,kk1,g)
          !write(*,'(a4,i4,a5,3i3,a4,3i3)'),'move',i,' from',ii,jj,kk,' to ',ii1,jj1,kk1     
          x1(:,i) = x1(:,i) - g(:) ! SOTTRAGGO LUNGHEZZA DI CUI DEVO RIENTRARE
          !----- check ----------------------------------------
          !call boxind(x1(:,i),ii2,jj2,kk2)
          !if (ii2.ne.ii1 .or. jj2.ne.jj1 .or. kk2.ne.kk1) then
          !   print*, ii2,jj2,kk2
          !   print*, x1(:,i)
          !   print*, g     
          !   STOP 'error of box'
          !endif 
          ! ---------------------------------------------------  
          call remove(boxlists(ii,jj,kk),node,i) ! RIMUOVE PARTICELLA i DAL BOX
          if (.not.associated(node)) then
              print*, i,'not found in',ii,jj,kk   
              STOP 'error'
          endif     
          call add(boxlists(ii1,jj1,kk1),node) ! AGGIUNGE PARTICELLE I AL BOX NUOVA
       end if
    end do  

  end subroutine update_boxes


  !//////////////////////////////////////////////////////////////
  subroutine boxinfo()
    
    write(*,*) 'Simulation box properties:'
    write(*,*) 'Lx=', LLx,  'Nx=', Nx
    write(*,*) 'Ly=', LLy,  'Ny=', Ny
    write(*,*) 'Lz=', LLz,  'Nz=', Nz
     
  end subroutine boxinfo  

  !//////////////////////////////////////////////////////////////
  subroutine boxind(r,i,j,k) 
     real(dp) :: r(3)

     integer :: i,j,k

     i = floor(r(1)/lx) + 1
     j = floor(r(2)/ly) + 1
     k = floor(r(3)/lz) + 1

  end subroutine boxind

  !//////////////////////////////////////////////////////////////
  subroutine folding(ii,jj,kk,g) ! SERVE PER LA PERIODICITÀ E
                                 ! E PER ASSEGNARE LE BOXES
     integer :: ii, jj, kk, m
     real(dp) :: g(3)
     !integer :: i,j,k
     
     !i=ii; j=jj; k=kk;

     g=0.d0

     if (ii.lt.1) then
        ! Versione che puo' far saltare piu' celle (raro)     
        m=(ii-Nx)/Nx ! DI QUANTE VOLTE STO FUORI I CONFINI PERIODICI
        g(1)=LLx*m ! LUNGHEZZA DI CUI DEVO RIENTRARE
        ii=Nx+mod(ii,Nx)
        ! Versione semplificata 
        !g(1) = -LLx
        !ii = Nx 
     else if (ii.gt.Nx) then
        ! Versione che puo' far saltare piu' celle (raro)        
        m=ii/Nx
        g(1)=LLx*m
        ii=mod(ii,Nx)
        !g(1) = LLx
        !ii = 1
     endif

     if (jj.lt.1) then
        ! Versione che puo' far saltare piu' celle (raro) 
        m=(jj-Ny)/Ny
        g(2)=LLy*m
        jj=Ny+mod(jj,Ny)
        !g(2) = -LLy  
        !jj = Ny    
     else if (jj.gt.Ny) then
        m=jj/Ny
        g(2)=LLy*m
        jj=mod(jj,Ny)
        !g(2) = LLy
        !jj = 1 
     endif

     if (kk.lt.1) then
        m=(kk-Nz)/Nz
        g(3)=LLz*m
        kk=Nz+mod(kk,Nz)
        !g(3) = -LLz
        !kk = Nz     
     else if (kk.gt.Nz) then
        m=kk/Nz
        g(3)=LLz*m
        kk=mod(kk,Nz)
        !g(3) = LLz
        !kk = 1
     endif 
  
  end subroutine folding
 
   


end module boxes
