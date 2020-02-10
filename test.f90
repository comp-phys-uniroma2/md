program test
  use list
  use clock
  implicit none

  type(Tlist) :: lista
  type(TNode), pointer :: node 
  type(TArray) :: array

  integer :: i, N, maxsz, val
  

  N = 50 
  maxsz= 4*N

  do i =1, N
     call add(lista,i)
  enddo
    
  node => lista%start
  print*, "node%val",node%val 

  call message_clock('remove from list')
  do i = 1, 1 
    do val = 1, N, 3
      print*, "node%val", node%val, "node%next%val",node%next%val
      call remove(lista,node,val)
      call delete(node) 
      node => lista%start
      print*, node%val 
    end do
    do val =1, N, 3
       call add(lista,val)
    enddo
  enddo
  call write_clock()

  !call write_list(lista)


  allocate(array%list(maxsz))
  array%list = 0
  do i = 1, N
    call add(array%list,i)
  enddo

  call message_clock('remove from array')
  do i = 1, 10000 
    do val = 1, N, 3
       call remove(array%list,val)
    end do
    do val =1, N, 3
       call add(lista,val)
    enddo
  enddo
  call write_clock()

  !call write_list(lista)


end program test
