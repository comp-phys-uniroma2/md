module lyap_n
  use constants
  use parameters
 
  
  implicit none


  public:: lyap_numbers


contains

  subroutine lyap_numbers(a,c,n,m)
    integer, intent(in) :: n, m 
    real(dp),dimension(6*Natoms,6*Natoms,n),intent(in):: a
    real(dp),dimension(6*Natoms),intent(out):: c
  
    integer:: i,j,pq,h,l
    integer:: err
    real(dp),dimension(:,:,:),allocatable :: G
    real(dp),dimension(:,:),allocatable :: Vol,Vol1
    real(dp),dimension(:),allocatable:: lyapsommati
    
    allocate(G(6*Natoms,6*Natoms,n),stat=err)
    allocate(Vol(6*Natoms,n),stat=err)
    allocate(Vol1(6*Natoms,n),stat=err)
    allocate(lyapsommati(6*Natoms),stat=err)
     if (err /= 0) STOP 'ALLOCATION ERROR'
    
     G=0.d0
     Vol=0.d0
     Vol1=0.d0
     lyapsommati=0.d0
     c=0.d0

     ! Norma vettore
     open(107,file='l_1.dat')
     ! Generalizzazioni a 5 e 10 vettori
     open(108,file='l_5.dat')
     open(109,file='l_10.dat')
     

     do pq=1,n
       write(*,*) pq,n
        
       do i=1,6*Natoms
         if(i==1) then
           Vol(i,pq)=sqrt(dot_product(a(:,i,pq),a(:,i,pq)))
           ! write(*,*) vol(i,pq)
         else
           Vol(i,pq)=sqrt(GramMatrix(a(:,1:i,pq),i,pq))
         end if

         if(vol(i,pq) .eq. 0) then
            write(*,*) "zero volume"
            write(*,*) Vol(i-1,pq)
            write(*,*) Vol(i,pq)
            do h=1,i
               write(*,*)'..............................................',i
               write(*,*) a(:,h,pq)
            end do
      
              do h=1,i
               write(*,*)'..............................................',i
               do l=1,i
                  if(l.eq.h) then
                  write(*,*)'--------norm-------'
                  write(*,*) dot_product(a(:,h,pq),a(:,l,pq))
                  write(*,*)'--------norm-------'
               else
                  write(*,*) dot_product(a(:,h,pq),a(:,l,pq))
               end if
      
            end do
      
            
            end do
            stop
         end if

      end do

   end do

   do i=1,6*Natoms
      Vol1(i,:)=log(Vol(i,:))
   end do


   do i=1,range(vol1(1,:))
      write(107,*) (Vol1(1,i))/(10*range(Vol1(1,:)))
      write(108,*) (Vol1(5,i)-Vol1(4,i))/(10*range(Vol1(1,:)))
      write(109,*) (Vol1(10,i)-Vol1(9,i))/(10*range(Vol1(1,:)))
   end do
   close(107)
   close(109)
   close(108)
   

   do i=1,6*Natoms
     lyapsommati(i)=sum(Vol1(i,:))/(10*range(Vol1(i,:)))
   end do

   do i=1,6*Natoms
     if (i==1) then
       c(i)=lyapsommati(i) 
     else
       c(i)= lyapsommati(i)-lyapsommati(i-1)
     end if   
   end do

 end subroutine lyap_numbers




 FUNCTION GRAMmatrix(Vectors,n,pq) result(f)

   real(dp),dimension(6*natoms,n) :: Vectors
   real(dp),dimension(n,n) :: G_matrix
   real(dp)::f
   integer,intent(in) :: n
   integer::p,q,pq
   !!$    do p=1,n
   !!$       write(*,*)'...............................',n
   !!$       do q=1,6*Natoms
   !!$          write(*,*) vectors(q,p)
   !!$       end do
   !!$       write(*,*)'...............................'
   !!$    end do


    do p=1,n
       do q=1,n
          G_matrix(p,q)=dot_product(Vectors(:,p),Vectors(:,q))
       end do
    end do

    !!$    do p=1,n
    !!$       write(*,*) G_matrix(:,p)
    !!$    end do

    f= finddet(G_matrix,n)

    !!$    write(*,*)'shape', shape(G_matrix),'det',f
    !!$    write(*,*) 'vol',f,n
    !!$    if(pq==) then
    !!$       stop
    !!$    end if

    
    
  end FUNCTION GRAMmatrix

  
  
  FUNCTION FindDet1(matrix, n) result(finddet)
    IMPLICIT NONE
    REAL(dp), DIMENSION(n,n) :: matrix
    real(dp):: finddet
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1

    !Convert to upper triangular form
    !!$    do i=1,n
    !!$       write(*,*) matrix(:,i),n
    !!$    end do
    !!$         if (n==5)then
    !!$       stop
    !!$    end if
    !!$    write(*,*)'---------------'

    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
     END DO

     
   
END FUNCTION FindDet1

REAL FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL(dp), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
  END FUNCTION FindDet




end module lyap_n

