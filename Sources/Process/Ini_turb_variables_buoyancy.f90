!==============================================================================!
  subroutine Ini_turb_variables_buoyancy(flow, backup, time_step)
!------------------------------------------------------------------------------!
!   Initialize dependent variables for buoyancy flows.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod, only: Field_Type, heat_transfer, density, viscosity
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
  real                     :: vis_t_max
  logical                  :: backup, present

!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t
  u    => flow % u
  v    => flow % v
  w    => flow % w

  vis_t_max = 0.0
  
  do c = 1, grid % n_cells
    if(grid % wall_dist(c) < 1.0) then
      vis_t_max = max(vis_t_max,vis_t(c)/viscosity)
    end if
  end do ! through cells

  call Comm_Mod_Global_Max_Real(vis_t_max)

  if(time_step == 100000920) then
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) > 1.0) then
        kin % n(c) = 0.000001
        kin % o(c) = kin % n(c) 
        kin % oo(c)= kin % n(c) 
        eps % n(c) = 0.00001
        eps % o(c) = eps % n(c) 
        eps % oo(c)= eps % n(c) 
        zeta % n(c) = 0.001 
        zeta % o(c) = zeta % n(c) 
        zeta % oo(c)= zeta % n(c) 
        f22 % n(c) = 0.01
        f22 % o(c) = f22 % n(c) 
        f22 % oo(c)= f22 % n(c) 
        u % n(c) = 0.0
        u % o(c) = u % n(c) 
        u % oo(c)= u % n(c) 
        v % n(c) = 0.0
        v % o(c) = v % n(c) 
        v % oo(c)= v % n(c) 
        w % n(c) = 0.0
        w % o(c) = w % n(c) 
        w % oo(c)= w % n(c) 
        t % n(c)  = 5.0 + 4.0 * grid % wall_dist(c)
        t % o(c)  = t % n(c)
        t % oo(c) = t % n(c)
      end if
    end do
  end if
  if(time_step < 1) then
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) < 0.25) then
        kin % n(c) = 0.003
        kin % o(c) = kin % n(c) 
        kin % oo(c)= kin % n(c) 
        eps % n(c) = 0.0001
        eps % o(c) = eps % n(c) 
        eps % oo(c)= eps % n(c) 
        zeta % n(c) = 0.1 
        zeta % o(c) = zeta % n(c) 
        zeta % oo(c)= zeta % n(c) 
        f22 % n(c) = 0.01
        f22 % o(c) = f22 % n(c) 
        f22 % oo(c)= f22 % n(c) 
      else
        kin % n(c) = 0.0000001
        kin % o(c) = kin % n(c) 
        kin % oo(c)= kin % n(c) 
        eps % n(c) = 0.001
        eps % o(c) = eps % n(c) 
        eps % oo(c)= eps % n(c) 
        zeta % n(c) = 0.0001 
        zeta % o(c) = zeta % n(c) 
        zeta % oo(c)= zeta % n(c) 
        f22 % n(c) = 0.01
        f22 % o(c) = f22 % n(c) 
        f22 % oo(c)= f22 % n(c) 
      end if 
    end do 
  end if

  if(this_proc < 2) write(*,*) 'VISt max below z = 1 is ', vis_t_max

  end subroutine
