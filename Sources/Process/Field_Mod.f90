!==============================================================================!
  module Field_Mod
!------------------------------------------------------------------------------!
!   Module for basic flow field plus temperature.                              !
!   It is a bit of a mumbo-jumbo at this moment, it will furhter have to       !
!   differentiate into numerical and physica parts.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Grid_Mod, only: Grid_Type
  use Bulk_Mod, only: Bulk_Type
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    ! Velocity components
    type(Var_Type) :: u
    type(Var_Type) :: v
    type(Var_Type) :: w

    ! Pressure and pressure correction
    type(Var_Type) :: p
    type(Var_Type) :: pp

    ! Temperature
    type(Var_Type) :: t

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars
    type(Var_Type), allocatable :: scalar(:)

    ! Mass fluxes throught cell faces
    real, allocatable :: flux(:)

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk

    ! Reference temperature
    real :: t_ref

  end type

  ! Variables determining if we are dealing with heat transfer and buoyancy
  logical :: heat_transfer
  logical :: buoyancy

  ! Heat flux to the domain (important for periodic case with heat transfer)
  real :: heat_flux, heated_area, heat

  ! Physical properties
  real :: viscosity, density, conductivity, diffusivity, capacity

  ! Angular velocity 
  real :: omega_x, omega_y, omega_z, omega

  ! Gravity
  real :: grav_x, grav_y, grav_z

  ! Thermal expansion coefficient 
  real :: beta_tec

  contains

  include 'Field_Mod/Allocate.f90'

  end module
