program md
  use parameters
  use input
  use boxes
  use simulations
  implicit none

  integer :: i 


  write(*,*) 'MD Simulator 0.1'
  write(*,*) '-------------------------------'
  call read_input()
 
  call transform_units()

  call create_xv()

  write(*,*) 'Set up simulation box'
  call create_boxes(Lx,Ly,Lz,Rc)

  call boxinfo()

  write(*,*) 'Set up particles'
  call init_particles()


  write(*,*) 'Starting MD run'
  call nve_sim() 

  write(*,*) '-------------------------------'
  write(*,*) 'Compute g(r)'
  call compute_g()


  write(*,*) 'delete boxes'
  call destroy_boxes()

  write(*,*) 'delete vectors'
  call destroy_xv()
  

end program md

