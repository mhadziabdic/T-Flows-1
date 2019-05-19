!==============================================================================!
  subroutine Initialize_ref_and_initial_temperature(flow, backup)
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

  integer                  :: n_points, k
  real, allocatable        :: prof(:,:), x(:), y(:), z(:), dist(:)
  logical                  :: backup, present

!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  
  if (this_proc < 2) print *, '# Grid material: ', grid % material % name

  if(.not. backup) then
    do c = 1, grid % n_cells

      if(heat_transfer) then
        t % n(c)  = 5.0 + 4.0 * grid % wall_dist(c)
!        t % n(c)  = 5.0 + 4.0 * grid % wall_dist(c)
!        t % n(c)  = 21.9 + 21.5 * grid % wall_dist(c)
        t % o(c)  = t % n(c)
        t % oo(c) = t % n(c)
      end if
    end do ! through cells
  end if

  do c = 1, grid % n_cells
    if(heat_transfer) then
      if(grid % wall_dist(c) < 1.0) then
        flow % t_ref_f(c) = 5.0 + 4.0 * 1.0 
!        flow % t_ref_f(c) = 21.9 + 21.5 * 1.0 
      else
        flow % t_ref_f(c) = 5.0 + 4.0 * grid % wall_dist(c) 
!        flow % t_ref_f(c) = 21.9 + 21.5 * grid % wall_dist(c) 
      end if        
    end if
  end do ! through cells

  end subroutine
