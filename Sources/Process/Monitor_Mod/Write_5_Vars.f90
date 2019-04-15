!==============================================================================!
  subroutine Monitor_Mod_Write_5_Vars(n, time, flow)
!------------------------------------------------------------------------------!
!   This is to set up monitoring points.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  integer          :: n
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, p, t
  type(Grid_Type), pointer :: grid
  integer                  :: m
  real                     :: vel_mag
  real                     :: time
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  do m = 1, monitor % n_points
    if(monitor % cell(m) > 0) then
      vel_mag = sqrt(( flow % u % n(monitor % cell(m)) )**2 &
                   + ( flow % v % n(monitor % cell(m)) )**2 &
                   + ( flow % w % n(monitor % cell(m)) )**2 )
      write(10+m,'(i9,9e16.6)')  n, time*1000.0/60.0,                  &
                                     flow % u % n(monitor % cell(m)),  &
                                     flow % v % n(monitor % cell(m)),  &
                                     flow % w % n(monitor % cell(m)),  &
                                     vel_mag,                          &
                                     flow % p % n(monitor % cell(m)),  &
                          180.0*ATAN(flow % v % n(monitor % cell(m)) / & 
                              flow % u % n(monitor % cell(m)))/3.14159,&
                                     flow % t % n(monitor % cell(m))   
                                     

    end if
  end do

  end subroutine
