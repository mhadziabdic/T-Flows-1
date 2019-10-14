!==============================================================================!
  subroutine Rivers_Interpolation(flow) 
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
  type(Var_Type),  pointer :: u, v, w, t
  integer             :: k, s, c1, c2, n_points, nearest_cell 
  real                :: new_distance, old_distance
  real                :: Distance
  character(len=80)   :: rivers_map_name
  character(len=80)   :: store_name
  real,allocatable    :: x_coord(:), y_coord(:)
  real,allocatable    :: angle(:) 
  logical             :: there
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u 
  v    => flow % v 
  w    => flow % w 
  t    => flow % t 

  ! Set the name for coordinate file
  call Name_File(0, rivers_map_name, ".riv")

  !------------------!
  !   Read 1r file   !
  !------------------!
  inquire(file=rivers_map_name, exist=there)
  if(.not. there) then
    if(this_proc < 2) then
      print *, '#=============================================================='
      print *, '# No name.riv file needed for river interpolation. Return      '
      print *, '#=============================================================='
    end if
    
    return
  end if

  open(9, file=rivers_map_name)

  ! Write the number of searching intervals
  read(9,*) 
  read(9,*) n_points
  allocate(x_coord(n_points)); x_coord = 0.0
  allocate(y_coord(n_points)); y_coord = 0.0
  allocate(angle(n_points))  ;   angle = 0.0

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

    read(9,*) x_coord(k), y_coord(k), angle(k) 

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
        if(id_zone(c1) == 1) then
          do k = 1, n_points
            new_distance = Distance(grid % xc(c2), grid % yc(c2), 0.0, &
                                          x_coord(k), y_coord(k), 0.0)
            if(new_distance <= old_distance) then
              nearest_cell =  k
              old_distance = new_distance
            end if
          end do

          river_angle(c1) = angle(nearest_cell)
          U % n(c2) = 2.0*cos(river_angle(c1))
          V % n(c2) = 2.0*sin(river_angle(c1))

        end if  
      end if  
    end if  
  end do

  deallocate(x_coord)
  deallocate(y_coord)
  deallocate(angle)

  if(this_proc < 2)  write(6, *) '# Finished with Rivers_Interpolation.f90'

  end subroutine
