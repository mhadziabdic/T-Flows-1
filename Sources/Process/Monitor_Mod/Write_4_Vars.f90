!==============================================================================!
  subroutine Monitor_Mod_Write_4_Vars(n, flow)
!------------------------------------------------------------------------------!
!   This is to set up monitoring points.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  integer          :: n
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, p
  integer                  :: m
!==============================================================================!

  do m = 1, monitor % n_points
    if(monitor % cell(m) > 0) then
      write(10+m,'(i9,5e16.6)')  n,  flow % u % n(monitor % cell(m)),  &
                                     flow % v % n(monitor % cell(m)),  &
                                     flow % w % n(monitor % cell(m)),  &
                                     flow % p % n(monitor % cell(m)),  &
                          180.0*ATAN(flow % v % n(monitor % cell(m)) / &
                              flow % u % n(monitor % cell(m)))/3.14159

    end if
  end do

  end subroutine
