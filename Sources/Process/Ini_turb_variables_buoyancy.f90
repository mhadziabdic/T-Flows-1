!==============================================================================!
  subroutine Ini_turb_variables_buoyancy(flow, backup, time_step)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
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

  integer                  :: n_points, k, n_wall, time_step, n_ref, n_vis
  real, allocatable        :: prof(:,:), x(:), y(:), z(:), dist(:)
  real                     :: vis_t_max, t_wall, t_ref, vis_t_sum
  logical                  :: backup, present

!==============================================================================!

  t_wall = 0.0
  n_wall = 0

  t_ref = 0.0
  vis_t_sum = 0.0
  n_ref = 0
  n_vis = 0


  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t
  u    => flow % u
  v    => flow % v
  w    => flow % w

  vis_t_max = 0.0
  
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % wall_dist(c) > 0.1.and.grid % wall_dist(c) < 0.25) then
      t_ref = t_ref + t % n(c)
      n_ref  = n_ref + 1     
    end if
  end do ! through cells

  call Comm_Mod_Global_Max_Real(vis_t_max)


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
  call Comm_Mod_Global_Sum_Real(t_ref)
  call Comm_Mod_Global_Sum_Int(n_ref)
 
  call Comm_Mod_Wait

  t_wall = t_wall/n_wall
  t_ref  = t_ref /n_ref 

  z_inv = (t_ref-5.0)/4.0

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % wall_dist(c) < z_inv) then
      vis_t_sum = vis_t_sum + vis_t(c)/viscosity
      n_vis     = n_vis + 1
    end if
  end do ! through cells

  call Comm_Mod_Global_Sum_Real(vis_t_sum)
  call Comm_Mod_Global_Sum_Int(n_vis)

  vis_t_sum  = vis_t_sum /n_vis 

  if(this_proc < 2) write(*,*) 'Twall = ', t_wall
  if(this_proc < 2) write(*,*) 'Tref  = ', (t_ref+t_wall)*0.5
  if(this_proc < 2) write(*,*) 'z_inv = ', z_inv
  if(this_proc < 2) write(*,*) 'vis_t_sum = ', vis_t_sum
 
  do c = 1, grid % n_cells
    if(heat_transfer) then
!      flow % t_ref_f(c) = max((t_ref+t_wall)*0.5,5.0 + 4.0 * grid % wall_dist(c))
      flow % t_ref_f(c) = 5.0 + 4.0 * grid % wall_dist(c)
!      flow % t_ref_f(c) = (t_wall, 5.0 + 4.0 * grid % wall_dist(c))
!      if(grid % wall_dist(c) > 0.32) then
!        w % n(c) = 0.0
!        w % o(c) = w % n(c) 
!        w % oo(c)= w % n(c) 
!      end if
    end if
  end do ! through cells

  if(time_step == 300) then
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) < 0.25) then
        kin % n(c) = 0.003
        kin % o(c) = kin % n(c) 
        kin % oo(c)= kin % n(c) 
        eps % n(c) = 0.0001
        eps % o(c) = eps % n(c) 
        eps % oo(c)= eps % n(c) 
        if(turbulence_model .eq. K_EPS_ZETA_F) then 
          zeta % n(c) = 0.1 
          zeta % o(c) = zeta % n(c) 
          zeta % oo(c)= zeta % n(c) 
          f22 % n(c) = 0.1
          f22 % o(c) = f22 % n(c) 
          f22 % oo(c)= f22 % n(c) 
        end if
      else
        kin % n(c) = 0.00001
        kin % o(c) = kin % n(c) 
        kin % oo(c)= kin % n(c) 
        eps % n(c) = 0.00001
        eps % o(c) = eps % n(c) 
        eps % oo(c)= eps % n(c) 
        if(turbulence_model .eq. K_EPS_ZETA_F) then 
          zeta % n(c) = 0.01 
          zeta % o(c) = zeta % n(c) 
          zeta % oo(c)= zeta % n(c) 
          f22 % n(c) = 0.01
          f22 % o(c) = f22 % n(c) 
          f22 % oo(c)= f22 % n(c) 
        end if
      end if
!       vis_t(c)   = viscosity
!        u % n(c) = 0.0
!        u % o(c) = u % n(c) 
!        u % oo(c)= u % n(c) 
!        v % n(c) = 0.0
!        v % o(c) = v % n(c) 
!        v % oo(c)= v % n(c) 
!        w % n(c) = 0.0
!        w % o(c) = w % n(c) 
!        w % oo(c)= w % n(c) 
!      if(grid % wall_dist(c) > z_inv) then
!        t % n(c)  = 5.0 + 4.0 * grid % wall_dist(c)
!        t % o(c)  = t % n(c)
!        t % oo(c) = t % n(c)
!      end if
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
        if(turbulence_model .eq. K_EPS_ZETA_F) then 
          zeta % n(c) = 0.001 
          zeta % o(c) = zeta % n(c) 
          zeta % oo(c)= zeta % n(c) 
        end if
        f22 % n(c) = 0.01
        f22 % o(c) = f22 % n(c) 
        f22 % oo(c)= f22 % n(c) 
      end if 
    end do 
  end if

!  if(this_proc < 2) write(*,*) 'VISt max below z = 1 is ', vis_t_max

  end subroutine
