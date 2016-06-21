MODULE verlet_module
  IMPLICIT none

CONTAINS
  subroutine read_input(file_name, n_at, atom_list, coord, mass)
    character(len=50), intent(in) :: file_name
    character(len=2), dimension(:), intent(inout) :: atom_list
    integer(4), intent(inout) :: n_at
    integer(4) :: i,j
    real(8), dimension(:,:),  intent(inout) :: coord
    real(8), dimension(:), intent(inout) :: mass

    open(20,file=file_name)
    read(20,*)
    
    do i=1,n_at
       read(20,*) atom_list(i), (coord(i,j),j=1,3), mass(i)
    end do

    close(20)

  end subroutine read_input


  subroutine write_geom(traj, n_at, atom_list, coord, kinetic)
    character(len=50) :: file_name
    character(len=2), dimension(:), intent(in) :: atom_list
    integer(4), intent(in) :: n_at, traj
    integer(4) :: k, m, out_unit
    real(8), intent(in) :: kinetic
    real(8), dimension(:,:), intent(in) :: coord

    out_unit = traj+10
    
    write(file_name,'(A4,I2.2,A4)') 'geom',traj,'.xyz'
    open(out_unit,file=file_name)
    
    write(out_unit,'(I3)') n_at
    write(out_unit,*) "COMMENT LINE", kinetic

   
    do k=1,n_at
       write(out_unit,*) atom_list(k), (coord(k,m)*0.529177249,m=1,3)
    end do
    
    close(out_unit)

  end subroutine write_geom

  subroutine grad_first(atom_list, coord, grad, n_at, traj)
    character(len=2), dimension(:), intent(in) :: atom_list
    character(len=50) :: run
    integer(4), intent(in) :: n_at, traj
    real(8) :: kinetic
    real(8), dimension(:,:), intent(inout) :: grad, coord

    kinetic = 0
    !Write the geometry of the trajectory nj
    call write_geom(traj,n_at,atom_list,coord,kinetic)

    ! Call the system to calculate V(n+1)
    write(run,'(A12,I2.2)') 'bash run.sh ',traj
    call system(run)

    ! Read the gradient V(n+1)
    call read_grad(traj,n_at,grad)

  end subroutine grad_first

  subroutine verlet_step(atom_list,coord,vel,grad,n_at,dt,mass,traj)
    character(len=2), dimension(:), intent(in) :: atom_list
    character(len=50) :: run
    integer(4), intent(in) :: n_at, traj
    integer(4) :: j
    real(8) :: t, kinetic
    real(8), intent(in) :: dt
    real(8), dimension(:,:), intent(in) :: mass
    real(8), dimension(:,:), intent(inout) :: coord, vel, grad


    ! Calculate V(n+0.5) e R(n+1)
    vel = vel + 0.5 * dt * (-grad/mass)
    coord = coord + vel*dt

    !Write the geometry of the trajectory nj
    call write_geom(traj,n_at,atom_list,coord,kinetic)

    ! Call the system to calculate V(n+1)
    write(run,'(A12,I2.2)') 'bash run.sh ',traj
    call system(run)

    ! Read the gradient V(n+1)
    call read_grad(traj,n_at,grad)

    ! Calculate V(n+1)
    vel = vel + 0.5 * dt * (-grad/mass)

   end subroutine verlet_step

   subroutine read_grad(traj,n_at,grad)
    character(len=50) :: file_name
    integer(4), intent(in) :: n_at, traj
    integer(4) :: k, m, out_unit
    real(8), dimension(:,:), intent(inout) :: grad

    out_unit = traj+10
    write(file_name,'(A8,I2.2,A7)') 'orca_job',traj,'.engrad'
    open(out_unit,file=file_name)

    ! Skip lines that don't containt relevant information
    do k=1,11
       read(out_unit,*)
    end do

    ! Read gradient
    do k=1,n_at
       do m=1,3
          read(out_unit,*) grad(k,m)
       end do
    end do
    
  end subroutine read_grad

END MODULE verlet_module
