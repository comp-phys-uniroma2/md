program test3
    implicit none
    integer :: i
    
    type mytype
        integer :: ind
        real(8) :: val
    end type mytype
    
    type tarray
        integer, dimension(:), allocatable :: list
    end type tarray
    
    !type(mytype) :: tt(10)
    type(tarray) :: variable_array(5)
    real(8), pointer :: rp
    type(mytype) :: myvar
    !tt(1)%ind = 5
    !tt(1)%val = 8.d0
    !print*, tt
    
    myvar%ind=1
    myvar%val=2.345
    
    rp => myvar%val
    print*, "rp=", rp

    do i=1,5
        allocate(variable_array(i)%list(2*i))
        variable_array(i)%list = i
    end do
    
    do i=1,5 
        print*, variable_array(i)%list
    end do

end program test3
