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
  public :: map
  public :: init_map

  type(Tlist), dimension(:,:,:), allocatable :: boxlists
  integer :: map(3,27)


  ! PRIVATE STUFF:
  integer :: Nx, Ny, Nz
  real(dp) :: lx, ly, lz
  real(dp) :: LLx, LLy, LLz

  contains

  !               x
  ! |---|---|---|---|---|---|---|---|---|---|
  ! 0 lx          4                         LLx

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

    ! Algebra modulare, funziona bene partendo da 0
    allocate(boxlists(0:Nx-1,0:Ny-1,0:Nz-1), stat=err)

    if (err.ne.0) STOP 'ERROR ALLOCATION Boxlist'

  end subroutine create_boxes

  ! ----------------------------------------------------
  subroutine destroy_boxes()
     integer i,j,k

     do i=0,Nx-1
      do j=0,Ny-1
       do k=0,Nz-1
         call delete_list(boxlists(i,j,k))
       enddo
      enddo
     enddo

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
       call boxind(x(:,i),ii,jj,kk)
       call boxind(x1(:,i),ii1,jj1,kk1)

       if (ii.ne.ii1 .or. jj.ne.jj1 .or. kk.ne.kk1) then ! ESEGUO SOLO SE
                                                         ! ESCO AD UNA BOX

          call folding(ii1,jj1,kk1,g)

          x1(:,i) = x1(:,i) - g(:) ! SOTTRAGGO LUNGHEZZA DI CUI DEVO RIENTRARE

          ! ---------------------------------------------------
          ! RIMUOVE PARTICELLA i DAL BOX
          call remove(boxlists(ii,jj,kk),node,i)
          if (.not.associated(node)) then
              print*, i,'not found in',ii,jj,kk
              STOP 'error'
          endif
          ! AGGIUNGE PARTICELLE I AL BOX NUOVA
          call add(boxlists(ii1,jj1,kk1),node)
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

     i = floor(r(1)/lx)
     j = floor(r(2)/ly)
     k = floor(r(3)/lz)

  end subroutine boxind

  !//////////////////////////////////////////////////////////////
  !
  !   ---||---|---|---|---|---||-o-|---|---|---|---|---|
  !      0                   LLx
  !                             x(t+dt)
  !    x(t+dt)=x(t+dt) + g
  subroutine folding(ii,jj,kk,g)

     integer, intent(inout) :: ii, jj, kk
     real(dp), intent(out) :: g(3)

     integer :: m

     g=0.d0

     ! Versione che puo' foldare di piu' celle
     if (ii.lt.0) then
        m=(ii-Nx+1)/Nx
        g(1)=LLx*m
        ii=mod(ii+16*Nx,Nx)
     else if (ii.gt.Nx-1) then
        m=ii/Nx
        g(1)=LLx*m
        ii=mod(ii,Nx)
     endif

     if (jj.lt.0) then
        m=(jj-Ny+1)/Ny
        g(2)=LLy*m
        jj=mod(jj+16*Ny,Ny)
     else if (jj.gt.Ny-1) then
        m=jj/Ny
        g(2)=LLy*m
        jj=mod(jj,Ny)
     endif

     if (kk.lt.0) then
        m=(kk-Nz+1)/Nz
        g(3)=LLz*m
        kk=mod(kk+16*Nz,Nz)
     else if (kk.gt.Nz-1) then
        m=kk/Nz
        g(3)=LLz*m
        kk=mod(kk,Nz)
     endif

     !if (ii<0 .or. jj<0 .or. kk<0 .or. ii>Nx-1 .or. jj>Ny-1 .or. kk>Nz-1) then
     !   print*, "FOLDING ERROR", ii, jj, kk
     !   stop
     !end if

  end subroutine folding

  subroutine init_map()
    integer :: u, v, w, i
    i = 1
    do w=-1, +1
      do v=-1, +1
        do u=-1, +1
          map(:,i) = (/u,v,w/)
          i = i + 1
        end do
      end do
    end do
  end subroutine init_map

end module boxes
