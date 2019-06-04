!==============================================================================!
  subroutine Initialize_ref_and_initial_temperature(flow, backup, time_step)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod, only: Field_Type, heat_transfer, density
  use Les_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
  use Bulk_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  integer                  :: i, c, c1, c2, m, s, nks, nvs

  integer                  :: n_points, k, n_wall, time_step
  real, allocatable        :: prof(:,:), x(:), y(:), z(:), dist(:)
  real                     :: t_wall
  logical                  :: backup, present

!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  t_wall = 0.0
  n_wall = 0
  
  if(.not. backup) then
    do c = 1, grid % n_cells

      if(heat_transfer) then
        t % n(c)  = 5.0 + 4.0 * grid % wall_dist(c)
        t % o(c)  = t % n(c)
        t % oo(c) = t % n(c)
      end if
    end do ! through cells
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. wallfl ) then
        t_wall = t_wall + t % n(c2)
        n_wall = n_wall + 1 
      end if
    end if
  end do
  
  call Comm_Mod_Global_Sum_Real(t_wall)
  call Comm_Mod_Global_Sum_Int(n_wall)
 
  call Comm_Mod_Wait

  t_wall = t_wall/n_wall

  if(this_proc < 2) write(*,*) 'Twall = ', t_wall
 
  do c = 1, grid % n_cells
    if(heat_transfer) then
      flow % t_ref_f(c) = max(t_wall, 5.0 + 4.0 * grid % wall_dist(c))
    end if
  end do ! through cells

  end subroutine
