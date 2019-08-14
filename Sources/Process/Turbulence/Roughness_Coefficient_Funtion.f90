!==============================================================================!
  subroutine Roughness_Coefficient_Funtion(flow) 
!------------------------------------------------------------------------------!
!                                                                              !
!   This subroutine reads values of roughness length in file and interpolate it!
!   on rough walls in the domain.                                              !
!                                                                              !
!------------------------------------------------------------------------------!
  use Grid_Mod
  use Field_Mod
  use Rans_Mod
  use Comm_Mod                       ! parallel stuff
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
!  type(Grid_Type)  :: grid
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: t
  integer             :: k, s, c, c1, c2, n_points, nearest_cell 
  real                :: new_distance, old_distance
  real                :: Distance
  character(len=80)   :: z_o_map_name
  character(len=80)   :: store_name
  real,allocatable    :: x_coord(:), y_coord(:), z_coord(:)
  real,allocatable    :: z_o_map(:), c_o_map(:) 
  integer,allocatable :: id_map(:)
  logical             :: there
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t 

  ! Set the name for coordinate file
  call Name_File(0, z_o_map_name, ".z_o")

  do c = 1, grid % n_cells
    wall_cells(c)   = -1.0
    ground_cells(c) = -1.0
  end do

  !------------------!
  !   Read 1r file   !
  !------------------!
  inquire(file=z_o_map_name, exist=there)
  if(.not. there) then
    if(this_proc < 2) then
      print *, '#=============================================================='
      print *, '# No name.z_o file. Return                                     '
      print *, '#=============================================================='
    end if
    
    return
  end if

  open(9, file=z_o_map_name)

  ! Write the number of searching intervals
  read(9,*) 
  read(9,*) n_points
  allocate(x_coord(n_points)); x_coord = 0.0
  allocate(y_coord(n_points)); y_coord = 0.0
  allocate(z_coord(n_points)); z_coord = 0.0
  allocate(z_o_map(n_points)); z_o_map = 0.0
  allocate(c_o_map(n_points)); c_o_map = 0.0
  allocate(id_map(n_points));  id_map  = 0

  ! Read the z_o map file
  ! Zone Sarajevo mesoscacle sa infrastrukturom:
  ! 1 - Rijeke
  ! 2 - Visoke zgrade 40m + 
  ! 3 - Zgrade 15 - 40m 
  ! 4 - Samostojeće kuće do 15m 
  ! 5 - Livade i nisko rastinje
  ! 6 - Šume 
  ! 7 - Ceste

  do k = 1, n_points
    read(9,*) id_map(k), x_coord(k), y_coord(k), z_coord(k), z_o_map(k)
    x_coord(k) = x_coord(k) * 0.001  
    y_coord(k) = y_coord(k) * 0.001  
    z_coord(k) = z_coord(k) * 0.001  
    z_o_map(k) = z_o_map(k) * 0.001
    if(id_map(k) == 1) then
      c_o_map(k) = 0.000000001
    else if(id_map(k) == 2) then
      c_o_map(k) = 0.0001
    else if(id_map(k) == 3) then
      c_o_map(k) = 0.0001
    else if(id_map(k) == 4) then
      c_o_map(k) = 0.01
    else if(id_map(k) == 5) then
      c_o_map(k) = 0.000000001
    else if(id_map(k) == 6) then
      c_o_map(k) = 0.000000001
    else if(id_map(k) == 7) then
      c_o_map(k) = 0.01
    end if
  end do
  close(9)

  nearest_cell = 0
  old_distance = HUGE
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
          Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        wall_cells(c1) = 1.0
        if(grid % bnd_cond % color(c2) == 1) then  ! This is "Ground" boundary condition
          ground_cells(c1) = 1.0
          old_distance = HUGE
          do k = 1, n_points
            new_distance = Distance(grid % xc(c2), grid % yc(c2), 0.0, &
                                          x_coord(k), y_coord(k), 0.0)
            if(new_distance <= old_distance) then
              nearest_cell =  k
              old_distance = new_distance
            end if
          end do

          z_o_f(c1) = z_o_map(nearest_cell)
          c_o_f(c1) = c_o_map(nearest_cell)
          id_zone(c1) = id_map(nearest_cell)
        else 
          z_o_f(c1) = z_o
          c_o_f(c1) = 0.0
        end if  
      end if  
    end if  
  end do

  deallocate(x_coord)
  deallocate(y_coord)
  deallocate(z_coord)
  deallocate(z_o_map)
  deallocate(c_o_map)
  deallocate(id_map)

  if(this_proc < 2)  write(6, *) '# Finished with Roughness_Coefficient_Funtion.f90'

  end subroutine
