program testp
  implicit none

  integer, pointer :: pint
  integer, dimension(:), pointer :: parr
  integer, target :: i1, i2
  !integer, dimension(:), allocatable, target :: arr
  integer, dimension(:), pointer :: arr

  allocate(arr(10))

  arr=1

  parr => arr

  print*,parr(4)


  parr(4) = 6

  print*, arr(4)

  

  
  

  
  
  i1 = 4
  i2 = 8

  pint => i1
  
  print*,pint

  pint => i2

  print*,pint 
  
end program testp  
