!==============================================================================!
  subroutine Source_Kin_K_Eps(flow, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grad_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Roughness_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tot2, u_nor, u_nor2, u_tan
  real                       :: kin_vis  ! [m^2/s]
  real                       :: ebf, p_kin_int, p_kin_wf
  real                       :: nx, ny, nz, qx, qy, qz, g_buoy_wall
  real                       :: ut_log_law, vt_log_law, wt_log_law
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  a    => sol % a
  b    => sol % b % val

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  ! Production source:
  do c = 1, grid % n_cells
    p_kin(c) = vis_t(c) * shear(c)**2
    b(c)     = b(c) + p_kin(c) * grid % vol(c)

    a % val(a % dia(c)) = a % val(a % dia(c)) + &
         density * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)

    if (buoyancy) then
      buoy_beta(c) = 1.0
      g_buoy(c) = -buoy_beta(c) * (grav_x * ut % n(c) +  &
                                   grav_y * vt % n(c) +  &
                                   grav_z * wt % n(c)) * density
!      g_buoy(c) = max(g_buoy(c),0.0)
      b(c) = b(c) + g_buoy(c) * grid % vol(c)
    end if
  end do

  kin_vis = viscosity / density

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tot2 = u % n(c1) **2  &
               + v % n(c1) **2  &
               + w % n(c1) **2
        u_nor = ( u % n(c1) * grid % sx(s)     &
                + v % n(c1) * grid % sy(s)     &
                + w % n(c1) * grid % sz(s) )   &
                     / sqrt(  grid % sx(s)**2  &
                            + grid % sy(s)**2  &
                            + grid % sz(s)**2 )
        u_nor2 = u_nor**2

        if(u_tot2  >  u_nor2) then
          u_tan = sqrt(u_tot2 - u_nor2)
        else
          u_tan = TINY
        end if

        if(rough_walls) then
          z_o = Roughness_Coefficient(grid, z_o_f(c1), c1)      
          z_o = max(grid % wall_dist(c1)/(e_log*y_plus(c1)),z_o)
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = u_tau(c1) * (grid % wall_dist(c1) + z_o) &
                     / kin_vis

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                       / log(((grid % wall_dist(c1)+z_o) / z_o))

          p_kin(c1) = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))   &
                    / (kappa*(grid % wall_dist(c1) + z_o))
          b(c1)     = b(c1) + (p_kin(c1) - vis_t(c1) * shear(c1)**2)   &
                    * grid % vol(c1)
        else
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan   &
                       / log(e_log * max(y_plus(c1), 1.05))

          ebf = 0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1))

          p_kin_wf  = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                    / (grid % wall_dist(c1) * kappa)

          p_kin_int = vis_t(c1) * shear(c1)**2

          p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                    * exp(-1.0 / ebf)

          b(c1) = b(c1) + (p_kin(c1) - p_kin_int) * grid % vol(c1)
        end if  ! rough_walls
        if(buoyancy) then
          nx = grid % sx(s) / grid % s(s)
          ny = grid % sy(s) / grid % s(s)
          nz = grid % sz(s) / grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ut_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nx
          vt_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * ny
          wt_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nz

          ut % n(c1) = ut %n(c1)  * exp(-1.0 * EBF) &
                     + ut_log_law * exp(-1.0 / EBF)
          vt % n(c1) = vt %n(c1)  * exp(-1.0 * EBF) &
                     + vt_log_law * exp(-1.0 / EBF)
          wt % n(c1) = wt %n(c1)  * exp(-1.0 * EBF) &
                     + wt_log_law * exp(-1.0 / EBF)
    
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) &
          t % q(c2) = abs(con_wall(c1)*(t % n(c1) &
                      - t % n(c2))/grid % wall_dist(c1))
          g_buoy_wall = abs(grav_z)*sqrt(t % q(c2)*  &
                        0.35*sqrt(abs(t2 % n(c1) * kin % n(c1)))) 

          ! Clean up b(c) from old values of g_buoy         
          b(c1)      = b(c1) - g_buoy(c1) * grid % vol(c1)

          g_buoy(c1) = g_buoy(c1) * exp(-1.0 * EBF) &
                     + g_buoy_wall * exp(-1.0 / EBF)

          ! Add new values of g_buoy based on wall function approach          
          b(c1)      = b(c1) + g_buoy(c1) * grid % vol(c1)
        end if    
      end if    ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
    end if      ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, kin % n)

  end subroutine
