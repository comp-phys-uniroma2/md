module input
  use parameters
  implicit none
  private

  public :: parse_input
 
  type(TLyapunovEnum), parameter, public :: LyAlgo = TLyapunovEnum()

  contains
  
  subroutine parse_input()

    character(20) :: token, val
    character(256) :: line
    integer :: idx1, idx2, err

    do  
      read(*,fmt='(A)',iostat=err) line
      if (err/=0) then
         exit
      endif   
      idx1 = index(trim(line),"=")
      idx2 = index(trim(line),"!")
      if (idx1==0 .or. idx2==1) cycle

      token = line(1:idx1-1)
      if (idx2>0) then
        val = line(idx1+1:idx2-1)
      else
        val = line(idx1+1:idx1+20)
      end if 
      write(*,'(A10,3x,A)') trim(adjustl(token)), trim(adjustl(val))


      if (trim(val)=="") then
        write(*,*) trim(token)    
        stop "error empty val"
      end if   

      call readvalue(token, val)
    end do   
 
  end subroutine parse_input


  subroutine readvalue(token, val)
    character(20) :: token, val

    if (trim(token)=="Natoms") then
       read(val,*) Natoms
    else if (trim(token)=="Lx") then
       read(val,*) Lx
    else if (trim(token)=="Ly") then
       read(val,*) Ly
    else if (trim(token)=="Lz") then
       read(val,*) Lz
    else if (trim(token)=="Nx") then
       read(val,*) Nx
    else if (trim(token)=="Ny") then
       read(val,*) Ny
    else if (trim(token)=="Nz") then
       read(val,*) Nz
    else if (trim(token)=="Rc") then
       read(val,*) Rc
    else if (trim(token)=="Rc") then
       read(val,*) Rc
   else if (trim(token)=="Temp") then
       read(val,*) Temp
    else if (trim(token)=="Mass") then
       read(val,*) Mass
    else if (trim(token)=="Eps") then
       read(val,*) eps
    else if (trim(token)=="Sigma") then
       read(val,*) Sigma
    else if (trim(token)=="Tinit") then
       read(val,*) Tinit
    else if (trim(token)=="Tsim") then
       read(val,*) Tsim
    else if (trim(token)=="dt") then
       read(val,*) dt
    else if (trim(token)=="Scaling") then
       read(val,*) Scaling
    else if (trim(token)=="Nose_Hoover") then
       read(val,*) nose_hoover
    else if (trim(token)=="Print_xyz") then
       read(val,*) print_xyz
    else if (trim(token)=="Print_interval") then
       read(val,*) print_interval 
    else if (trim(token)=="NH_Mass") then
       read(val,*) Q
    else if (trim(token)=="dr") then
       read(val,*) dr
    else if (trim(token)=="Lyapunov") then
       read(val,*) do_lyapunov
    else if (trim(token)=="LyapunovAlgo") then
       if (trim(val) == "volumes") then 
         algorithm = LyAlgo%volumes
       else if (trim(val) == "lengths") then    
         algorithm = LyAlgo%lengths
       else if (trim(val) == "QR") then   
         algorithm = LyAlgo%QR
       end if
    else if (trim(token)=="Tortho") then
       read(val,*) tgram
    else if (trim(token)=="D0") then
       read(val,*) D0
    else if (trim(token)=="Dtemp") then
       read(val,*) Dtemp
    else if (trim(token)=="Tempf") then
       read(val,*) Tempf
    else if (trim(token)=="Shnn") then
       read(val,*) shnn
    else if (trim(token)=="Tsh") then
       read(val,*) tsh
    else if (trim(token)=="Fext") then
       read(val,*) Fe
    else 
       write(*,*) "Token '",trim(token),"' not recognized"  
    end if

  end subroutine readvalue


end module input

