! cosmic_ray_ionization.f90
module cosmic_ray_ionization_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: cosmic_ray_ionization_rates, cosmic_ray_ionization_rates_induced_UV, cosmic_ray_ionization_rates_induced_UV_heat
  public :: initialize_cr_rates, secondary_cr_rates
  public :: cosmic_ray_ionization_rates_induced_UV_Ca_plus, cosmic_ray_ionization_rates_induced_UV_heat_Ca_plus

  real(dp) :: cosmic_ray_ionization_rates(27,27)
  real(dp) :: cosmic_ray_ionization_rates_induced_UV(27)
  real(dp) :: cosmic_ray_ionization_rates_induced_UV_heat(27)

  ! Need to deal with Ca+ separately (it's the only one that acts like this)
  real(dp) :: cosmic_ray_ionization_rates_induced_UV_Ca_plus
  real(dp) :: cosmic_ray_ionization_rates_induced_UV_heat_Ca_plus

  ! Data taken from UMIST where possible https://umistdatabase.uk/database
  ! Otherwise calculated similar to CHIMES
  ! Harley computed the heating rates from the Gredel 89 spectrum

  ! Hydrogen
  real(dp), parameter :: cosmic_ray_ionization_rates_hydrogen(27)  = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Helium
  real(dp), parameter :: cosmic_ray_ionization_rates_helium(27)    = (/ 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Carbon
  real(dp), parameter :: cosmic_ray_ionization_rates_carbon(27)    = (/ 3.83, 1.6638695, 0.8312262, 0.4541473, 0.06937006, 0.027755102, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Nitrogen
  real(dp), parameter :: cosmic_ray_ionization_rates_nitrogen(27)  = (/ 4.52, 1.8237852, 0.8574613, 0.5265334, 0.30504152, 0.04926644, 0.020386748, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Oxygen
  real(dp), parameter :: cosmic_ray_ionization_rates_oxygen(27)    = (/ 5.637, 1.9201249, 0.98636097, 0.5261845, 0.36398333, 0.21907325, 0.03679156, 0.015607069, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Neon
  real(dp), parameter :: cosmic_ray_ionization_rates_neon(27)      = (/ 4.8251953, 2.2843935, 1.2726027, 0.6982404, 0.43035796, 0.2581805, 0.204098, 0.12864962, 0.022742474, 0.009985316, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Sodium
  real(dp), parameter :: cosmic_ray_ionization_rates_sodium(27)    = (/ 11.72499, 2.26537, 1.31612, 0.81958, 0.49044, 0.31553, 0.19556, 0.16117, 0.10316, 0.02781, 0.00825, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  ! Magnesium
  real(dp), parameter :: cosmic_ray_ionization_rates_magnesium(27) = (/ 9.016716, 5.037889, 1.345602, 0.86565775, 0.5749281, 0.36411676, 0.24167651, 0.15331018, 0.13045996, 0.084565096, 0.015437003, 0.006928171, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Aluminum
  real(dp), parameter :: cosmic_ray_ionization_rates_aluminum(27)  = (/ 6.48529, 4.2986, 2.83717, 0.9014, 0.61631, 0.42696, 0.28137, 0.19101, 0.12355, 0.1077, 0.07056, 0.01954, 0.0059, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  ! Silicon
  real(dp), parameter :: cosmic_ray_ionization_rates_silicon(27)   = (/ 6.37329, 2.4664342, 2.5934594, 1.8422631, 0.64957255, 0.46263242, 0.3302073, 0.22406998, 0.15485357, 0.10161172, 0.09062556, 0.059769776, 0.011156686, 0.0050879163, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Sulfur
  real(dp), parameter :: cosmic_ray_ionization_rates_sulfur(27)    = (/ 7.5164065, 2.8824131, 1.5510099, 0.8586138, 1.2691398, 0.9713803, 0.38636714, 0.28946567, 0.21488377, 0.15199806, 0.107721165, 0.07223445, 0.066572525, 0.044464245, 0.008436725, 0.003892387, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Argon
  real(dp), parameter :: cosmic_ray_ionization_rates_argon(27)     = (/ 10.0, 3.39221, 1.98146, 1.13255, 0.72275, 0.44717, 0.76151, 0.60413, 0.25708, 0.19859, 0.15121, 0.10993, 0.07926, 0.05397, 0.05097, 0.03438, 0.00989, 0.00307, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  ! Chlorine
  real(dp), parameter :: cosmic_ray_ionization_rates_chlorine(27)  = (/ 6.52, 3.363, 1.70566, 1.0128, 0.59972, 0.96532, 0.75473, 0.31175, 0.23753, 0.17885, 0.1284, 0.09186, 0.06212, 0.058, 0.03893, 0.01114, 0.00345, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  ! Calcium
  real(dp), parameter :: cosmic_ray_ionization_rates_calcium(27)   = (/ 12.434942, 6.8917594, 2.1202607, 1.406548, 0.9606202, 0.6236024, 0.42673957, 0.2766994, 0.5097157, 0.41350347, 0.18359101, 0.14470413, 0.112184905, 0.08312979, 0.06079397, 0.04186227, 0.0402239, 0.02735196, 0.005303178, 0.0024862888, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

  ! Iron
  real(dp), parameter :: cosmic_ray_ionization_rates_iron(27)      = (/ 3.192021, 5.842174, 4.724112, 2.5858068, 1.7483323, 1.2082134, 0.85902447, 0.62585855, 0.46490052, 0.3626024, 0.28070357, 0.20538796, 0.15056798, 0.10395802, 0.21419027, 0.18064684, 0.086154655, 0.07006637, 0.056015324, 0.042974997, 0.032204125, 0.02267751, 0.022580078, 0.015563778, 0.0030807566, 0.0014658332, 0.0 /)

CONTAINS

SUBROUTINE initialize_cr_rates()
  implicit none
  integer :: i,j

  ! Initialize
  do i = 1, 27
     do j = 1, 27
        cosmic_ray_ionization_rates(i,j) = 0.d0
        cosmic_ray_ionization_rates_induced_UV_heat(i) = 0.d0
     end do
  end do

  ! Now populate
  do i = 1, 27
     cosmic_ray_ionization_rates(1,i)  = cosmic_ray_ionization_rates_hydrogen(i)
     cosmic_ray_ionization_rates(2,i)  = cosmic_ray_ionization_rates_helium(i)
     cosmic_ray_ionization_rates(6,i)  = cosmic_ray_ionization_rates_carbon(i)
     cosmic_ray_ionization_rates(7,i)  = cosmic_ray_ionization_rates_nitrogen(i)
     cosmic_ray_ionization_rates(8,i)  = cosmic_ray_ionization_rates_oxygen(i)
     cosmic_ray_ionization_rates(10,i) = cosmic_ray_ionization_rates_neon(i)
     cosmic_ray_ionization_rates(11,i) = cosmic_ray_ionization_rates_sodium(i)
     cosmic_ray_ionization_rates(12,i) = cosmic_ray_ionization_rates_magnesium(i)
     cosmic_ray_ionization_rates(13,i) = cosmic_ray_ionization_rates_aluminum(i)
     cosmic_ray_ionization_rates(14,i) = cosmic_ray_ionization_rates_silicon(i)
     cosmic_ray_ionization_rates(16,i) = cosmic_ray_ionization_rates_sulfur(i)
     cosmic_ray_ionization_rates(17,i) = cosmic_ray_ionization_rates_chlorine(i)
     cosmic_ray_ionization_rates(18,i) = cosmic_ray_ionization_rates_argon(i)
     cosmic_ray_ionization_rates(20,i) = cosmic_ray_ionization_rates_calcium(i)
     cosmic_ray_ionization_rates(26,i) = cosmic_ray_ionization_rates_iron(i)
  end do

  ! Cosmic ray data for ionization from induced UV emission
  ! Data from: https://home.strw.leidenuniv.nl/~ewine/photo/cosmic_ray_rates.html
  ! Note that this assumes 10-16 s-1 H2-1 for the primary H2 cosmic ray dissociation rate
  ! Initialize to 0
  do i = 1, 27
     cosmic_ray_ionization_rates_induced_UV(i) = 0.d0
  end do

  cosmic_ray_ionization_rates_induced_UV(1)  = 0.00d+00 !4.08d-16 ! Hydrogen
  cosmic_ray_ionization_rates_induced_UV(2)  = 0.00d+00 ! Helium  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV(6)  = 2.60d-14 ! Carbon
  cosmic_ray_ionization_rates_induced_UV(7)  = 7.34d-17 ! Nitrogen
  cosmic_ray_ionization_rates_induced_UV(8)  = 2.70d-16 ! Oxygen
  cosmic_ray_ionization_rates_induced_UV(10) = 0.00d+00 ! Neon  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV(11) = 1.29d-15 ! Sodium
  cosmic_ray_ionization_rates_induced_UV(12) = 1.12d-14 ! Magnesium
  cosmic_ray_ionization_rates_induced_UV(13) = 2.50d-13 ! Aluminum
  cosmic_ray_ionization_rates_induced_UV(14) = 4.16d-13 ! Silicon
  cosmic_ray_ionization_rates_induced_UV(16) = 7.91d-14 ! Sulfur
  cosmic_ray_ionization_rates_induced_UV(17) = 4.73e-15 ! Chlorine
  cosmic_ray_ionization_rates_induced_UV(18) = 0.00d+00 ! Argon
  cosmic_ray_ionization_rates_induced_UV(20) = 2.71d-14 ! Calcium
  cosmic_ray_ionization_rates_induced_UV(26) = 4.81d-14 ! Iron

  cosmic_ray_ionization_rates_induced_UV_heat(1)  = 0.00d0 ! Hydrogen
  cosmic_ray_ionization_rates_induced_UV_heat(2)  = 0.00d0 ! Helium  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV_heat(6)  = 0.8697 !0.96d0 ! Carbon
  cosmic_ray_ionization_rates_induced_UV_heat(7)  = 0.2851 !1.17d0 ! Nitrogen
  cosmic_ray_ionization_rates_induced_UV_heat(8)  = 0.4766 !1.72d0 ! Oxygen
  cosmic_ray_ionization_rates_induced_UV_heat(10) = 0.00d0 ! Neon  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV_heat(11) = 4.3457 !0.00d0 ! Sodium  
  cosmic_ray_ionization_rates_induced_UV_heat(12) = 1.9733 !2.55d0 ! Magnesium
  cosmic_ray_ionization_rates_induced_UV_heat(11) = 3.4987 !0.00d0 ! Aluminum 
  cosmic_ray_ionization_rates_induced_UV_heat(14) = 2.1836 !2.08d0 ! Silicon
  cosmic_ray_ionization_rates_induced_UV_heat(16) = 1.1771 !1.32d0 ! Sulfur
  cosmic_ray_ionization_rates_induced_UV_heat(17) = 0.00d0 ! Chlorine (i don't have a cross section)
  cosmic_ray_ionization_rates_induced_UV_heat(18) = 0.00d0 ! Argon 
  cosmic_ray_ionization_rates_induced_UV_heat(20) = 3.3717 !0.00d0 ! Calcium 
  cosmic_ray_ionization_rates_induced_UV_heat(26) = 2.1590 !2.32d0 ! Iron

  cosmic_ray_ionization_rates_induced_UV_Ca_plus = 1.54d-16
  cosmic_ray_ionization_rates_induced_UV_heat_Ca_plus = 0.7971

END SUBROUTINE initialize_cr_rates

FUNCTION secondary_cr_rates(xe) result(phi_s)
    ! Secondary CR ionization rate
    implicit none
    real(dp), intent(in) :: xe
    real(dp) :: phi_s

    phi_s = (1.d0 - (xe / 1.2d0)) * (0.670d0 / (1.d0 + (xe / 0.05d0)));
END FUNCTION secondary_cr_rates

end module cosmic_ray_ionization_module