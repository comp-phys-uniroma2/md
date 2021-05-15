module vector
  use constants
  use parameters
  implicit none

  public :: grams
  public :: qr

contains



  subroutine grams(q,p,v,u)
    real(dp),dimension(3,Natoms,6*Natoms),intent(inout) :: q, p
    real(dp),dimension(6*Natoms,6*Natoms),intent(out):: V
    real(dp),dimension(6*Natoms,6*Natoms),intent(out):: U

    integer:: i,j
    
    !Combina q-p da 2 array (3, Natoms, 6*Natoms) -> (6Natoms, 6Natoms)
    V=newshape2(q,p)
    
    U(:,1) = V(:,1)/sqrt(dot_product(V(:,1),V(:,1)))

    do i = 2, 6*Natoms
      U(:,i) = V(:,i)
      do j = 1,i-1
        U(:,i) = U(:,i) - U(:,j)*dot_product( U(:,j),U(:,i) )/dot_product(U(:,j),U(:,j))
      end do
      U(:,i) = U(:,i)/sqrt(dot_product(U(:,i),U(:,i)))
    end do

    ! Orthonormalization check 
    !do i=1,6*Natoms
    !   do j=1,6*Natoms
    !      if (j.ne.i) then
    !         if (abs(dot_product(u(:,i),u(:,j))).ge.0.0001) then
    !            write(*,*) 'ortornormal error',i,j
    !            stop
    !         end if
    !      end if
    !   end do
    !end do
   
    !riporta u agli array q p 
    call oldshape(U,q,p)

  end subroutine grams

   
  subroutine qr(q,p,lyap_ex)
    real(dp),dimension(3,Natoms,6*Natoms),intent(in) :: q,p
    real(dp),dimension(6*Natoms),intent(inout) :: lyap_ex
    
    integer:: i,j,info,lwork,n
    real(dp),dimension(6*Natoms,6*Natoms)::AB
    real(dp),allocatable:: tau(:),work(:)

    n=6*Natoms
    lwork=n
    allocate(tau(n),work(lwork))

    AB=newshape2(q,p)
   
    call dgeqr2p(n,n,AB,n,tau,work,info)

    do i=1,6*Natoms
      if (AB(i,i)> 0.0_dp) then
        lyap_ex(i)=lyap_ex(i)+log(AB(i,i))/tsim
      end if    
    end do
        
  end subroutine qr

  

  function Projection(a,b) result(c)
    real(dp),dimension(6*Natoms),intent(in) :: a,b
      real(dp),dimension(6*Natoms) :: c
      
      c=dot_product(a,b)*a/dot_product(a,a)
      write(*,*)'oooooooooooooooooooooooo'
      write(*,*) dot_product(a,a),dot_product(a,b)
      write(*,*)'oooooooooooooooooooooooo'

  end function Projection



  

  function newshape(a,b) result(d)
    integer:: i,k,j
    real(dp),dimension(3,Natoms,6*Natoms),intent(in) :: a,b
    real(dp),dimension(6,Natoms,6*Natoms) :: c
    real(dp),dimension(6*Natoms,6*Natoms) :: d
     
    do i=1,6*Natoms
       c(1:3,:,i)=a(1:3,:,i)
       c(4:6,:,i)=b(1:3,:,i)
    end do

    do i=1,6*Natoms
      d(:,i)=reshape(c(:,:,i),(/6*Natoms/))
    end do

  end function newshape

   
  function newshape2(a,b) result(d)
    real(dp),dimension(3,Natoms,6*Natoms),intent(in) :: a,b
    real(dp),dimension(6*Natoms,6*Natoms) :: d
     
    integer:: i,k,j
    real(dp),dimension(3*Natoms,6*Natoms) :: e,f
        
    do i=1,6*Natoms
      e(:,i)=reshape(a(:,:,i),(/3*Natoms/))
      f(:,i)=reshape(b(:,:,i),(/3*Natoms/))
    end do  

    do i=1,6*Natoms
      d(1:3*Natoms,i)=e(1:3*Natoms,i)
      d(3*Natoms+1:6*Natoms,i)=f(1:3*Natoms,i)
    end do

   end function newshape2


   subroutine oldshape(a,b,c)
     real(dp),dimension(6*Natoms,6*Natoms),intent(in)::a
     real(dp),dimension(3,Natoms,6*Natoms),intent(out)::b,c

     real(dp),dimension(3*Natoms,6*Natoms):: e, f
     integer:: i,k

     do i=1,6*Natoms
        e(1:3*Natoms,i) = a(1:3*Natoms,i) 
        f(1:3*Natoms,i) = a(3*Natoms+1:6*Natoms,i) 
        b(:,:,i)=reshape(e(:,i),(/3,Natoms/))
        c(:,:,i)=reshape(f(:,i),(/3,Natoms/))
     end do

   end subroutine oldshape
     

 end module vector

 

   


        



     
  
