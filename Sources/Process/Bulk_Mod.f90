!==============================================================================!
  module Bulk_Mod
!------------------------------------------------------------------------------!
!   Mass fluxes, bulk velocities and pressure drops (for each material)        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Bulk type   !
  !---------------!
  type Bulk_Type

    real :: area_x  ! [m^2]
    real :: area_y  ! [m^2]
    real :: area_z  ! [m^2]

    real :: mass_in   ! [kg/s]
    real :: mass_out  ! [kg/s]

    ! Fluxes in x, y and z direction
    real :: flux_x  ! [kg/s]
    real :: flux_y  ! [kg/s]
    real :: flux_z  ! [kg/s]

    ! Old fluxes in x, y and z direction
    real :: flux_x_o  ! [kg/s]
    real :: flux_y_o  ! [kg/s]
    real :: flux_z_o  ! [kg/s]

    real :: p_drop_x  ! [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_y  ! [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_z  ! [N/m^3] = [kg/m^2/s^2]

    real :: u  ! [m/s]
    real :: v  ! [m/s]
    real :: w  ! [m/s]

    ! Monitoring plane coordinates
    real :: xp
    real :: yp
    real :: zp

  end type

  contains

  include 'Bulk_Mod/Adjust_P_Drops.f90'
  include 'Bulk_Mod/Monitoring_Planes_Areas.f90'
  include 'Bulk_Mod/Print_Areas.f90'

  end module
