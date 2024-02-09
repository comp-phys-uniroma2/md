program md
  use parameters
  use input
  use boxes
  use simulations
  implicit none

  integer :: i 


  write(*,*) 'MD Simulator 2D'
  write(*,*) '------------------------------------'
  call parse_input()
  write(*,*) '------------------------------------'
  write(*,*) 'transform units'
  call transform_units()

  write(*,*) 'create xv, eta'
  call create_xv()
  call create_eta()

  write(*,*) 'Set up simulation box'
  call create_boxes(Lx,Ly,Rc)
  call init_map()

  call boxinfo()

  write(*,*) 'Set up particles'
  call init_seed(111111111)
  call init_positions_fcc()
  if (v_drift .ne. 0.0_dp) then
    call init_velocities_couette()
  else
    call init_velocities()
  endif

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

  

end program md

