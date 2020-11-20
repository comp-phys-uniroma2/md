module list
  implicit none
  private
  public :: Tnode
  public :: Tlist
  public :: TArray
  public :: add
  public :: remove
  public :: delete 
  public :: write_list
  public :: delete_list

  type Tnode ! INDICE PARTICELLA E PUNTATORE AL TNODE SUCCESSIVO
     integer :: val
     type(Tnode), pointer :: next => null()
  end type Tnode

  Type Tlist ! LISTA DI PARTICELLE (START AND END POINTER TO PARTICLES)
             ! END POINTER NON NECESSARIO, MA PIÙ EFFICIENTE
     type(Tnode), pointer :: start => null()
     type(Tnode), pointer :: last  => null()
  end type Tlist

  Type TArray ! NON UTILUZZATO (CASO CON ARRAY INVECE CHE POINTER)
     integer, dimension(:), allocatable :: list
  end type TArray

  interface add
     module procedure add_int
     module procedure add_node
     module procedure add_array_int
  end interface add     

  interface remove 
     module procedure remove_int
     module procedure remove_node
     module procedure rm_array_int
  end interface remove     

  contains
  
  subroutine write_list(list)
     type(Tlist), intent(in) :: list
  
     type(Tnode), pointer :: node 
     
     node => list%start

     do while (associated(node))  
        write(*,*) node%val   
        node => node%next
     enddo

  end subroutine write_list
  
  !------------------------------------------------------
  ! List%start -> nodo1%next -> nodo2%next -> nodo3%next => node
  ! List%last  ---- ------------------------------------------^
  subroutine add_node(list,node)
     type(Tlist)  :: list
     type(Tnode), pointer :: node

     ! Associated controlla che start /= null
     if (.not.associated(list%start)) then
        list%start => node 
        list%last => node 
     else
        list%last%next => node
        list%last => node
     endif
 
  end subroutine add_node   

  !------------------------------------------------------
  subroutine add_int(list,val)
     type(Tlist), intent(in)  :: list
     integer :: val

     type(Tnode), pointer :: node
     integer :: err
     
     allocate(node, stat= err)

     if (err.ne.0) STOP 'NODE ALLOCATION ERROR' 

     node%val = val

     call add(list,node)

  end subroutine add_int

  !------------------------------------------------------
  ! RIMUOVE DA UNA LISTA IL NODO DELLA
  ! PARTICELLA CON INDICE val. OVVIAMENTE
  ! IL NODO NON VIENE DEALLOCATO, MA
  ! VIENE AGGANCIATO AD UN'ALTRA LISTA.
  !                    |---------------------------|
  ! List%start -> nodo1%next   nodo2%next   nodo3%next => node
  !                             nodo2%val 
  ! List%last  ---------------------------------------------^
  ! 
  !------------------------------------------------------
  subroutine remove_int(list, vnode, val)

     type(Tlist)  :: list
     type(Tnode), pointer :: vnode ! NODE TROVATO DA CONSEGNARE IN USCITA 
     integer, intent(in) :: val
  
     type(Tnode), pointer :: node 
     type(Tnode), pointer :: pnode ! PREVIOUS NODE
     
     if (associated(list%start)) then
         node => list%start
         pnode => list%start ! PNODE È IL NODO PRECEDENTE
     else
         !print*,'empty list'    
         vnode=>null()    
         return
     endif
  
     ! search the node with val and remove it from list   
     if (node%val == val) then  ! SE LA PARTICELLA DA RIMUOVERE È ALL'INIZIO
                                ! DELLA LISTA IN CUI STO CERCANDO
        list%start => node%next
        if (.not.associated(node%next)) then ! SE IL NODO IN QUESTIONE ERA L'ULTIMO
           list%last => null()               ! ALLORA LA LISTA RIMARRÀ VUOTA
        endif     
        node%next => null()
        vnode => node
     else
        do  
           node => node%next
           if (associated(node)) then 
              if (node%val == val) then ! SE HO TROVATO LA PARTICELLA DA RIMUOVERE
                 pnode%next => node%next
                 if (.not.associated(node%next)) then
                    list%last => pnode
                 endif     
                 node%next => null()
                 vnode => node
                 exit
              endif
           else                        ! NON HO TROVATO LA PARTICELLA
              vnode => null()      
              exit
           endif  
           pnode => node
        enddo
     endif

  end subroutine remove_int

  subroutine remove_int2(list, vnode, val)
    type(Tlist)  :: list
    type(Tnode), pointer :: vnode 
    integer, intent(in) :: val
  
    type(Tnode), pointer :: node 
    type(Tnode), pointer :: pnode
 
    pnode => list%start
    node => list%start
!    if (node%val == val) goto 1002

    do
        node => node%next
        if (associated(node)) then
            if (node%val==val) then
                vnode => node
                pnode%next => node%next
                ! CASO list%start%val == val
                if (node%val==pnode%val) then
                    list%start => node%next
                    pnode => null()
                end if
                ! CASO list%last%val == val
                if (.not.associated(node%next)) list%last => pnode

                node%next => null()
                exit
            end if
        else
            vnode => null()
            exit
        end if
        pnode => node
    end do
  end subroutine remove_int2
  !----------------------------------------------------
  
  subroutine remove_node(list, node)
     type(Tlist) :: list
     type(Tnode), pointer :: node 

     type(Tnode), pointer :: it
     type(Tnode), pointer :: itp

     if(associated(node)) then
         ! to be implemented
     end if 

  end subroutine remove_node

  !----------------------------------------------------
  subroutine delete(node)
     type(Tnode), pointer :: node 

     if (associated(node)) then
        !print*,'remove',node%val    
        deallocate(node)
     end if

  end subroutine delete   

  !----------------------------------------------------
  subroutine search(list,node,val)
     type(Tlist), intent(in)  :: list
     type(Tnode), pointer :: node 
     integer, intent(in)  :: val

     type(Tnode), pointer :: it

     if (associated(list%start)) then
         it => list%start
     else
         node => null()    
         return
     endif

     ! search the node with val    
     if (it%val == val) then
        node => it
     else   
        do
           it => it%next
           if (associated(it)) then
              if (it%val == val) then
                 node => it
                 exit
              endif     
           else
              node => null()
           endif        
        enddo      
     endif

  end subroutine

  !----------------------------------------------------
  subroutine delete_list(list)
     type(Tlist), intent(in)  :: list
     type(Tnode), pointer :: node 
     type(Tnode), pointer :: nnode 

     node => list%start

     do while (associated(node))
        nnode => node%next
        call delete(node)
        node => nnode
     enddo

  end subroutine delete_list

  !//////////////////////////////////////////////////////////
  subroutine rm_array_int(list,val)
    integer, dimension(:) :: list 
    integer, intent(in) :: val
    
    integer :: i, maxsz
    maxsz = size(list) 
   
    i =1
    do while (i.le.maxsz)
       if (list(i) .eq. val) exit
       i = i + 1
    enddo

    if (i.lt.maxsz) then
       list(i:maxsz-1) = list(i+1:maxsz)
       list(maxsz) = 0
    elseif (i.eq.maxsz) then   
       list(maxsz) = 0
    endif

  end subroutine rm_array_int

  !//////////////////////////////////////////////////////////
  subroutine add_array_int(list,val)
    integer, dimension(:) :: list 
    integer, intent(in) :: val

    integer :: i, maxsz
    maxsz = size(list) 
  
    i =1
    do while (i.le.maxsz)
       if (list(i) .eq. val) exit
       i = i + 1
    enddo
    
    if (i.lt.maxsz) then
       list(i+1:maxsz) = list(i:maxsz-1) 
       list(i) = val 
    elseif (i.eq.maxsz) then   
       list(maxsz) = val 
    endif

  end subroutine add_array_int

end module list

