PROGRAM verlet_algorithm
  use verlet_module
  IMPLICIT none

  character(len=50) :: file_name
  character(len=2), dimension(:), allocatable :: atom_list
  integer(4) :: n_at, n_traj, ns, i, j, k
  real(8) :: t, dt, kinetic
  real(8), dimension(:,:,:), allocatable :: coord, grad, vel, ma
  real(8), dimension(:,:), allocatable :: m
  real(8), dimension(:), allocatable :: mass


  ! Initial variables
  file_name="hcn.xyz"
  n_traj=2
  ns=10000
  dt=10 !a.u
  n_at = 3
  t=0


  allocate (atom_list (n_at), coord(n_at,3,n_traj), mass(n_at), grad(n_at,3,n_traj), vel(n_at,3,n_traj), m(n_at,3))

  ! Read initial input file
  call read_input(file_name,n_at,atom_list,coord(:,:,1),mass)

  ! Convert units
  vel = 0
  grad = 0
  mass = mass * 1836.152
  coord = coord * 1.88971616463207


  do i=2,n_traj
     coord(:,:,i)= coord(:,:,1)
  end do

  do i=1,n_at
     m(i,:) = mass(i)
  end do


  !PRIVATE:every open mpi thread will have it's own private copy of variables and others don't have acces to it
  !SHARED: every open mpi thread will have access to all these variables

  ! For the first step
  !$OMP PARALLEL PRIVATE(i) SHARED(atom_list,coord,vel,grad,n_at,ns,mass)
  !$OMP DO
  do i=1,n_traj
     call grad_first(atom_list,coord(:,:,i),grad(:,:,i), n_at, i)
  end do
  !$OMP END DO 
  !$OMP END PARALLEL

  do j=1,ns
     !$OMP PARALLEL PRIVATE(i) SHARED(atom_list,coord,vel,grad,n_at,ns,mass)
     !$OMP DO
     do i=1,n_traj
        call verlet_step(atom_list,coord(:,:,i),vel(:,:,i),grad(:,:,i), n_at, dt, m, i)

     end do

     !$OMP END DO 
     !$OMP END PARALLEL
     t=t+dt

     kinetic=0
     do k=1,n_at
        kinetic = 0.5*sum(vel(k,:,i)**2)*mass(k)
     end do
     write(6,*)t,kinetic

  end do


  STOP
END PROGRAM verlet_algorithm
  
