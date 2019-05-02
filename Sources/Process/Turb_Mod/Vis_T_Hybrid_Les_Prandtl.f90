!==============================================================================!
  subroutine Turb_Mod_Vis_T_Hybrid_Les_Prandtl(turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.
!   In case that, in parallel executions, the subdomain does not have 
!   any nearwall cells, the near(c) is zero.
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig
!   of a simulation.
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod,  only: Comm_Mod_Exchange_Real
  use Field_Mod, only: Field_Type, density, viscosity, buoyancy,  &
                       grav_x, grav_y, grav_z
  !use Les_Mod
  use Grid_Mod,  only: Grid_Type
  use Grad_Mod,  only: Grad_Mod_Array
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w, t
  integer                   :: c, s, c1, c2
  real                      :: fd    ! damping function
  real                      :: hwn   ! grid step size in wall-normal direction 
  real                      :: hmax
  real                      :: cw, u_ff
  real                      :: dw
  real                      :: lf_wm
  real                      :: kappa
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  t    => flow % t 

  !---------------!                                                             
  !   Constants   !                                                             
  !---------------!                                                             
  cw    = 0.15 ! emperical constant [Shur 2008]                             
  kappa = 0.41

  ! Calculate model's eddy viscosity
  do c = 1, grid % n_cells

    hmax = turb % h_w(c)
    hwn  = turb % h_max(c)
    dw   = grid % wall_dist(c)

    ! Wall-modeled LES length scale
    lf_wm = min(max(cw*dw,cw*hmax,hwn),hmax)

    ! If(nearest_wall_cell(c) .ne. 0) is needed for parallel
    ! version since the subdomains which do not "touch" wall
    ! has nearest_wall_cell(c) = 0. 
    if(nearest_wall_cell(c) .ne. 0) then
      u_ff = sqrt( viscosity  &
                  * sqrt(  u % n(nearest_wall_cell(c)) ** 2   &
                         + v % n(nearest_wall_cell(c)) ** 2   &
                         + w % n(nearest_wall_cell(c)) ** 2)  &
                 / (grid % wall_dist(nearest_wall_cell(c))+TINY) )
      y_plus(c) = grid % wall_dist(c) * u_ff / viscosity

      ! Piomelli damping function
      fd = 1.0 - exp(-(y_plus(c)/25.0)**3)
    else
      fd = 1.0
    end if
    vis_t(c) = min((c_smag*lf_wm)**2, (kappa*dw)**2) * flow % shear(c) * fd

  end do

  !-----------------! 
  !   Wall Region   ! 
  !-----------------+---------------------------!
  !  The procedure below calculates the vis..   !
  !  .. at the wall.                            !
  !----------------.----------------------------! 
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      vis_wall(c1) = viscosity                       &
                   +        grid % fw(s) *vis_t(c1)  &
                   + (1.0 - grid % fw(s))*vis_t(c2)
    end if    ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, vis_t)
  call Comm_Mod_Exchange_Real(grid, vis_wall)

  end subroutine
