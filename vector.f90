module vector
  use constants
  use parameters
  implicit none

  public:: grams


  public:: oldshape
 


contains



  subroutine grams(a,b,v,u)
    real(dp),dimension(3,Natoms,6*Natoms),intent(inout) :: a,b
    real(dp),dimension(6*Natoms,6*Natoms),intent(out):: v
    real(dp),dimension(6*Natoms,6*Natoms),intent(out):: u

    integer:: i,j,k,p
    
    real(dp),dimension(6*Natoms) :: proj,Norm
    real(dp) :: D0
    proj=0.d0
    D0=1.0d0
   
    

    
    !salvare v per ogni pq(lyap)
    v=newshape2(a,b)
    

    U(:,1) = V(:,1)/sqrt(dot_product (V(:,1),V(:,1)))
    
  

    do i = 2,6*Natoms
      U(:,i) = V(:,i)
      do j = 1,i-1
        U(:,i) = U(:,i) - U(:,j)*dot_product( U(:,j),U(:,i) )/dot_product(U(:,j),U(:,j))
     end do

     
      U(:,i) = U(:,i)/sqrt(dot_product ( U(:,i),U(:,i)))
   end do

   

!!$        do i=1,6*Natoms
!!$             do j=1,6*Natoms
!!$                if(j.ne.i) then
!!$                   if (dot_product(v(:,i),v(:,j)).ge.0.5) then
!!$                      write(*,*) 'ortornormal errorv',i,j
!!$                      write(*,*) dot_product(v(:,i),v(:,j))
!!$                      write(*,*) sqrt(dot_product(v(:,j),v(:,j)))
!!$                      write(*,*)sqrt( dot_product(v(:,i),v(:,i)))
!!$                      do p=1,6*Natoms
!!$                         write(*,*)v(p,i),v(p,j),p
!!$                         if (v(p,i) .ne.0 .and. v(p,j).ne.0) then
!!$                         write(*,*) v(p,i)/v(p,j)
!!$                      end if
!!$
!!$                      
!!$
!!$                      
!!$                      end do
!!$
!!$                      
!!$                      stop
!!$                   end if
!!$
!!$                end if
!!$
!!$             end do
!!$
!!$          end do
  
          do i=1,6*Natoms
             do j=1,6*Natoms
                if(j.ne.i) then
                   if (dot_product(u(:,i),u(:,j)).ge.0.0001) then
                      !write(*,*) 'ortornormal error',i,j
                      !write(*,*) dot_product(u(:,i),u(:,j))
                      !write(*,*) sqrt(dot_product(v(:,j),v(:,j)))
                      !write(*,*)sqrt( dot_product(v(:,i),v(:,i)))
                      !write(*,*) sqrt(dot_product(u(:,j),u(:,j)))
                      !write(*,*)sqrt( dot_product(u(:,i),u(:,i)))
                      do p=1,6*Natoms
                        ! write(*,*)u(p,i),u(p,j),p
                      end do

                      
                      stop
                   end if

                else if(j.eq.i) then
                   if (nint(dot_product(u(:,i),u(:,j))).ne.1) then
                      !write(*,*) 'norm err'
                      !write(*,*) dot_product(u(:,i),u(:,j))
                      do p=1,6*Natoms
                         !write(*,*)u(p,i),p
                      end do
                      

                   end if

                end if

                
                   
                   

             end do

          end do
   
    
  
    !ortogonalizzazione
!!$    do i=1,6*Natoms
!!$       if (i==1) then
!!$          u(:,i)=v(:,i)
!!$         
!!$
!!$       else
!!$         do j=1,i
!!$          proj=proj+Projection(u(:,i-1),v(:,i))
!!$          
!!$          u(:,i)=v(:,i)-proj
!!$          
!!$
!!$       end if
!!$
!!$    end do
!!$
!!$    !normalizzazione
!!$    do i=1,6*Natoms
!!$       Norm(i)=sqrt(dot_product(u(:,i),u(:,i)))
!!$       if (i.le.6)then
!!$          write(*,*) Norm(i)
!!$       end if
!!$
!!$       
!!$
!!$       
!!$    end do
!!$    do i=1,6*Natoms
!!$       u(:,i)=u(:,i)*D0/(Norm(i))
!!$    end do

    
    

    !ricombinare u in formato 3-D
    call oldshape(u,a,b)
   

  
         

  end subroutine grams


  

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


       do k=1,6
          if (k.le.3) then
             
             c(k,:,i)=a(k,:,i)
             
          else
             
             c(k,:,i)=b(k-3,:,i)



          end if
           
        end do

     end do

      
     do i=1,6*Natoms
        d(:,i)=reshape(c(:,:,i),(/6*Natoms/))

     end do

!!$     do i=1,6
!!$        do k=1,6*Natoms
!!$           write(*,*) d(k,i)
!!$           end do
!!$        
!!$            write(*,*) '--------------------------'
!!$     end do
!!$stop
   end function newshape

   
  function newshape2(a,b) result(d)
    integer:: i,k,j
    real(dp),dimension(3,Natoms,6*Natoms),intent(in) :: a,b
    real(dp),dimension(3*Natoms,6*Natoms) :: e,f
    real(dp),dimension(6*Natoms,6*Natoms) :: d
     
        
     do i=1,6*Natoms
        e(:,i)=reshape(a(:,:,i),(/3*Natoms/))
        f(:,i)=reshape(b(:,:,i),(/3*Natoms/))

     end do  


     

        
    do i=1,6*Natoms


       do k=1,6*Natoms
          if (k.le.3*natoms) then
             
             d(k,i)=e(k,i)
             
          else
             
             d(k,i)=f(k-3*natoms,i)



          end if
           
        end do

     end do

!!$     do i=1,3*Natoms
!!$        write(*,*)  e(i,6*Natoms),d(i,6*Natoms)
!!$     end do
!!$
!!$     write(*,*)'-----'
!!$
!!$     do i=1,3*Natoms
!!$        write(*,*)  f(i,6*Natoms),d(3*Natoms+i,6*Natoms)
!!$     end do
!!$stop
!!$     do i=1,3
!!$        write(*,*) a(i,1:3,3*natoms+2)
!!$     end do
!!$   do i=1,3
!!$        write(*,*) b(i,1:3,3*natoms+2)
!!$     end do
!!$
!!$    do i=1,6*natoms
!!$        write(*,*) d(i,6*natoms)
!!$   end do
!!$     stop
        



!!$     do i=1,6
!!$        do k=1,6*Natoms
!!$           write(*,*) d(k,i)
!!$           end do
!!$        
!!$            write(*,*) '--------------------------'
!!$     end do
!!$stop
   end function newshape2

   subroutine oldshape(a,b,c)
     real(dp),dimension(3,Natoms,6*Natoms),intent(out)::b,c
     real(dp),dimension(6,Natoms,6*Natoms)::d
     real(dp),dimension(6*Natoms,6*Natoms),intent(in)::a
     integer:: i,k

     do i=1,6*Natoms
        d(:,:,i)=reshape(a(:,i),(/6,Natoms/))
     end do


     do i=1,6*Natoms


       do k=1,6
          if (k.le.3) then
             
             b(k,:,i)=d(k,:,i)
             
          else
             
             c(k-3,:,i)=d(k,:,i)



          end if
           
        end do

     end do
 

    

          

   end subroutine oldshape
   
        

     

 end module vector

 

   


        



     
  
