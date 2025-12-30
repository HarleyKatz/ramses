! secondary_ionization_module.f90
module secondary_ionization_module
  use amr_parameters, only: dp
  use safe_math, only: safe_exp
  implicit none

  private  ! everything is private by default
  public :: secondary_ionization_fracs

contains 

FUNCTION secondary_ionization_fracs(E_0,x_e,type) result(frac)
  !
  !See Appendix A of https://iopscience.iop.org/article/10.1086/341255/pdf
  !This is based on ricotti 2002
  !
  !Fraction of secondary ionizations per primary electron
  !type = 0 : heat
  !type = 1 : HI ionization
  !type = 2 : HeI ionization
  !"""
  implicit none
  real(dp), intent(in) :: E_0, x_e
  integer, intent(in) :: type
  real(dp) :: frac
  real(dp) :: loc_xe

  loc_xe = MAX(MIN(x_e,1.d0),0.d0)

  frac = 0.d0

  select case (type)

  case (0)  ! Heat
    if (E_0 .ge. 11.d0) then 
       if (x_e.le.1.d-4) then 
          frac = 0.15d0
       else
          frac = 3.9811d0 * ((11.d0 / E_0)**0.7d0) * (loc_xe**0.4d0) * (1.d0 - (loc_xe**0.34d0))**2.d0 + (1.d0 - (1.d0 - (loc_xe**0.2663d0))**1.3163d0)
       end if
    else
       frac = 1.d0
    end if

  case (1)  ! HI
    if (E_0 .ge. 28.d0) then 
       frac = -0.6941d0 * ((28.d0 / E_0)**0.4d0) * (loc_xe**0.2d0) * (1.d0 - (loc_xe**0.38d0))**2.d0 + 0.3908d0 * (1.d0 - (loc_xe**0.4092d0))**1.7592d0
    end if

  case (2)  ! HeI
    if (E_0 .ge. 28.d0) then 
        frac = -0.0984d0 * ((28.d0 / E_0)**0.4d0) * (loc_xe**0.2d0) * (1.d0 - (loc_xe**0.38d0))**2.d0 + 0.0554d0 * (1.d0 - (loc_xe**0.4614d0))**1.6660d0
    end if
  end select

END FUNCTION secondary_ionization_fracs

end module secondary_ionization_module