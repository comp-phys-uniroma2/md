program md
  use parameters
  use input
  use boxes
  use simulations
  implicit none

  integer :: i 


  write(*,*) 'MD Simulator 0.1 + Lyapunov spectrum'
  write(*,*) '------------------------------------'
  call parse_input()
  write(*,*) '------------------------------------'
  
  write(*,*) 'transform units'
  call transform_units()

  write(*,*) 'create xv, eta, xvly'
  call create_xv()
  call create_eta()
  call create_xvly()

  write(*,*) 'Set up simulation box'
  call create_boxes(Lx,Ly,Lz,Rc)
  call init_map()

  call boxinfo()

  write(*,*) 'Set up particles'
  call init_seed(111111111)
  call init_positions_fcc()
  call init_velocities()
  !call init_particles()


  write(*,*) '------------------------------------'
  write(*,*) 'Starting MD run'
  call nve_sim() 


  write(*,*) '------------------------------------'
  write(*,*) 'Compute g(r)'
  call compute_g()


  write(*,*) 'delete boxes'
  call destroy_boxes()

  write(*,*) 'delete vectors'
  call destroy_xv()
  call destroy_eta()

  call destroy_xvly()
  

end program md

