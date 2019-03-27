!==============================================================================!
  subroutine Roughness_Coefficient_Funtion(grid) 
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
  type(Grid_Type)  :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer             :: k, s, c1, c2, n_points, nearest_cell
  real                :: new_distance, old_distance
  real                :: Distance
  character(len=80)   :: z_o_map_name
  character(len=80)   :: store_name
  real,allocatable    :: x_coord(:), y_coord(:), z_coord(:), z_o_map(:)
  logical             :: there
!==============================================================================!

  ! Set the name for coordinate file
  call Name_File(0, z_o_map_name, ".z_o")

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
  read(9,*) n_points
  allocate(x_coord(n_points)); x_coord = 0.0
  allocate(y_coord(n_points)); y_coord = 0.0
  allocate(z_coord(n_points)); z_coord = 0.0
  allocate(z_o_map(n_points)); z_o_map = 0.0

  ! Read the z_o map file
  do k = 1, n_points
    read(9,*) x_coord(k), y_coord(k), z_o_map(k)
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
      end if  
    end if  
  end do

  deallocate(x_coord)
  deallocate(y_coord)
  deallocate(z_coord)
  deallocate(z_o_map)

  if(this_proc < 2)  write(6, *) '# Finished with Roughness_Coefficient_Funtion.f90'

  end subroutine
