! metal_yields_module.f90
module metal_yields_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: initialize_SN_yields, get_portinari_ejecta_mass, get_portinari_stellar_lifetime
  public :: getStarAgeMyr, get_pop3_lifetime_myr, get_popIII_sn_energy, get_popIII_ejecta
  public :: get_SNIa_ejecta, get_stellar_lifetime_low_mass, get_nomoto_ejecta

  integer, parameter::portinati_N_metal = 5
  integer, parameter::portinati_N_mass = 14

  real(dp):: portinati_metal(portinati_N_metal)
  real(dp):: portinati_mass(portinati_N_mass)
  real(dp):: portinati_lifetimes(portinati_N_metal,portinati_N_mass)
  ! Yields are loaded as H, He, C, N, O, Ne, Mg, Si, S, Fe, Ca
  ! The last one is out of order...
  real(dp):: portinati_yields(portinati_N_metal,portinati_N_mass,11)

  ! Pop III yields
  real(dp), dimension(1:7) :: table_mass_POP3_SNII = (/ 13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0 /)
  real(dp), dimension(1:4) :: table_mass_POP3_HN   = (/ 20.0, 25.0, 30.0, 40.0 /)
  real(dp), dimension(1:6) :: table_mass_POP3_PISN = (/ 140.0, 150.0, 170.0, 200.0, 270.0, 300.0 /)

  real(dp),dimension(1:10,1:7) :: PopIII_SNII
  real(dp),dimension(1:10,1:4) :: PopIII_HN
  real(dp),dimension(1:10,1:6) :: PopIII_PISN

  ! All Nomoto yields
  real(dp), dimension(1:7) :: table_mass_NOMOTO_SNII = (/ 13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0 /)
  real(dp), dimension(1:4) :: table_mass_NOMOTO_HN   = (/ 20.0, 25.0, 30.0, 40.0 /)
  real(dp), dimension(1:4) :: table_metal_NOMOTO     = (/ 0.0, 0.001, 0.004, 0.02 /)
  real(dp), dimension(1:4) :: table_mass_NOMOTO_PISN = (/ 150.0, 170.0, 200.0, 270.0 /)
  real(dp), dimension(1:27,1:7,1:4) :: NOMOTO_SN_YIELDS
  real(dp), dimension(1:27,1:4,1:4) :: NOMOTO_HN_YIELDS
  real(dp), dimension(1:27,1:4) :: NOMOTO_PISN_YIELDS

CONTAINS

SUBROUTINE initialize_SN_yields()
  use amr_commons, only: myid
  implicit none

  integer:: i, j, k, unit_num, ios

  if(myid.eq.1) write(*,*) 'Initializing Metal Yields'

  !!!
  !!! SNII Yields from Portinari 1998
  !!!

  ! Load Portinari yields file
  open(newunit=unit_num, file='./data/yields/Portinari_yields.tsv', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open Portinari yields file'
      return
  end if

  do i=1,portinati_N_metal ! loop over metallicities

     ! Read in the metallicity
     read(unit_num, *, iostat=ios) portinati_metal(i) 

     do j=1,portinati_N_mass ! Loop over stellar masses

        ! Read in mass, stellar lifetime and 
        read(unit_num, *, iostat=ios) portinati_mass(j), portinati_lifetimes(i,j), portinati_yields(i,j,:)

     end do ! end loop over stellar masses

  end do ! end loop over metallicities

  !!!
  !!! Pop III Yields from Nomoto 2013
  !!!

  ! Load Pop III SNII Yields
  open(newunit=unit_num, file='./data/yields/Nomoto2013_PopIII_SNII.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open Pop III SNII yields file'
      return
  end if

  do i=1,10 ! loop over metals
     read(unit_num, *, iostat=ios) PopIII_SNII(i,:)
  end do ! end loop over metals

  ! Load Pop III HN Yields
  open(newunit=unit_num, file='./data/yields/Nomoto2013_PopIII_HN.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open Pop III SN yields file'
      return
  end if

  do i=1,10 ! loop over metals
     read(unit_num, *, iostat=ios) PopIII_HN(i,:)
  end do ! end loop over metals

  ! Load Pop III PISN Yields
  open(newunit=unit_num, file='./data/yields/Nomoto2013_PopIII_PISN.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open Pop III PISN yields file'
      return
  end if

  do i=1,10 ! loop over metals
     read(unit_num, *, iostat=ios) PopIII_PISN(i,:)
  end do ! end loop over metals


  ! Load Nomoto SN yields
  open(newunit=unit_num, file='./data/yields/nomoto_SN_yields.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open nomoto_SN_yields yields file'
      return
  end if

  do i=1,27 ! loop over metals
     do j=1,7 ! loop over masses
        read(unit_num, *, iostat=ios) NOMOTO_SN_YIELDS(i,j,:)
     end do ! end loop over masses
  end do ! end loop over metals

  ! Load Nomoto HN yields
  open(newunit=unit_num, file='./data/yields/nomoto_HN_yields.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open nomoto_HN_yields yields file'
      return
  end if

  do i=1,27 ! loop over metals
     do j=1,4 ! loop over masses
        read(unit_num, *, iostat=ios) NOMOTO_HN_YIELDS(i,j,:)
     end do ! end loop over masses
  end do ! end loop over metals

  ! Load Nomoto PISN yields
  open(newunit=unit_num, file='./data/yields/nomoto_PISN_yields.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open nomoto_PISN_yields yields file'
      return
  end if

  do i=1,27 ! loop over metals
      read(unit_num, *, iostat=ios) NOMOTO_PISN_YIELDS(i,:)
  end do ! end loop over metals


END SUBROUTINE initialize_SN_yields

FUNCTION get_nomoto_ejecta(mass, metallicity, element_idx, is_HN) result(fq)
  implicit none
  real(dp), intent(in) :: mass, metallicity
  integer, intent(in) :: element_idx
  logical, intent(in) :: is_HN
  real(dp):: xq, yq, frac_high, frac_low
  real(dp) :: fq
  integer :: i, j
  real(dp) :: x1, x2, y1, y2
  real(dp) :: f11, f21, f12, f22
  real(dp) :: dx, dy
  real(dp) :: rescale_factor

  ! Initialize
  fq = 0.d0

  ! No ejecta for direct collapse
  if (mass.gt.40.d0.and.mass.lt.150.d0) then
     return
  end if

  ! No ejecta for direct collapse
  if (mass.gt.270.d0) then
     return
  end if

  ! Start with the PISN case
  if (mass.ge.150.d0 .and. mass.le.270.d0) then

       ! Interpolate
       do i = 1,3
          if (mass.ge.table_mass_NOMOTO_PISN(i) .and. mass.le.table_mass_NOMOTO_PISN(i+1)) then
             j = i
          end if
       end do

       frac_high = (mass - table_mass_NOMOTO_PISN(j)) / (table_mass_NOMOTO_PISN(j+1) - table_mass_NOMOTO_PISN(j))
       frac_low = 1.d0 - frac_high

       fq = (frac_low * NOMOTO_PISN_YIELDS(element_idx,j)) + (frac_high * NOMOTO_PISN_YIELDS(element_idx,j+1))

     return
  end if 

  ! Concinue to lower mass SN
  rescale_factor = 1.0 ! Rescale factor for metallicity and mass

  ! Next is the HN case
  if (is_HN) then
     ! Enforce bounds
     xq = MIN(MAX(metallicity,table_metal_NOMOTO(1)),table_metal_NOMOTO(4))
     yq = MIN(MAX(mass,table_mass_NOMOTO_HN(1)),table_mass_NOMOTO_HN(4))

     ! No need to downweight yields
     ! metallicity lower limit is 0 and we only launch HN in a fixed window of mass

     ! Find i such that x(i) <= xq <= x(i+1)
     do i = 1, 3
        if (xq.ge.table_metal_NOMOTO(i) .and. xq.le.table_metal_NOMOTO(i+1)) exit
     end do

     ! Find j such that y(j) <= yq <= y(j+1)
     do j = 1, 3
        if (yq.ge.table_mass_NOMOTO_HN(j) .and. yq.le.table_mass_NOMOTO_HN(j+1)) exit
     end do

     ! Extract corner values
     x1 = table_metal_NOMOTO(i);   x2 = table_metal_NOMOTO(i+1)
     y1 = table_mass_NOMOTO_HN(j);   y2 = table_mass_NOMOTO_HN(j+1)

     f11 = NOMOTO_HN_YIELDS(element_idx, j, i)
     f21 = NOMOTO_HN_YIELDS(element_idx, j, i+1)
     f12 = NOMOTO_HN_YIELDS(element_idx, j+1, i)
     f22 = NOMOTO_HN_YIELDS(element_idx, j+1, i+1)

  ! Otherwise, normal SN
  else
     ! Enforce bounds
     xq = MIN(MAX(metallicity,table_metal_NOMOTO(1)),table_metal_NOMOTO(4))
     yq = MIN(MAX(mass,table_mass_NOMOTO_SNII(1)),table_mass_NOMOTO_SNII(7))

     ! Downweight the yields if we are below lower limits --> do not extrapolate in the case where we are above!
     ! No need to downweight metallicity --> lower limit is 0
     if (mass.lt.table_mass_NOMOTO_SNII(1)) then
       rescale_factor = rescale_factor * (mass / table_mass_NOMOTO_SNII(1))
     endif

     ! Find i such that x(i) <= xq <= x(i+1)
     do i = 1, 3
        if (xq.ge.table_metal_NOMOTO(i) .and. xq.le.table_metal_NOMOTO(i+1)) exit
     end do

     ! Find j such that y(j) <= yq <= y(j+1)
     do j = 1, 6
        if (yq.ge.table_mass_NOMOTO_SNII(j) .and. yq.le.table_mass_NOMOTO_SNII(j+1)) exit
     end do

     ! Extract corner values
     x1 = table_metal_NOMOTO(i);   x2 = table_metal_NOMOTO(i+1)
     y1 = table_mass_NOMOTO_SNII(j);   y2 = table_mass_NOMOTO_SNII(j+1)

     f11 = NOMOTO_SN_YIELDS(element_idx, j, i)
     f21 = NOMOTO_SN_YIELDS(element_idx, j, i+1)
     f12 = NOMOTO_SN_YIELDS(element_idx, j+1, i)
     f22 = NOMOTO_SN_YIELDS(element_idx, j+1, i+1)

  end if

  dx = x2 - x1
  dy = y2 - y1

  ! Bilinear interpolation
  fq = (1.0 / (dx * dy)) * ( &
        f11 * (x2 - xq) * (y2 - yq) + &
        f21 * (xq - x1) * (y2 - yq) + &
        f12 * (x2 - xq) * (yq - y1) + &
        f22 * (xq - x1) * (yq - y1)   &
       ) * rescale_factor

END FUNCTION get_nomoto_ejecta

FUNCTION get_popIII_ejecta(mass, species, is_HN) result(fq)
  implicit none
  real(dp), intent(in) :: mass
  integer, intent(in) :: species
  logical, intent(in) :: is_HN
  real(dp) :: fq
  integer :: i, j, element_idx
  real(dp) :: frac_high, frac_low, loc_mass

  ! Initialize
  fq = 0.d0

  ! No ejecta for direct collapse
  if (mass.gt.40.d0.and.mass.lt.140.d0) then
     return
  end if

  ! No ejecta for direct collapse
  if (mass.gt.300.d0) then
     return
  end if

  select case (species)
    case (1) ! Hydrogen
      element_idx = 1
    case (2) ! Helium
      element_idx = 2
    case (6) ! Carbon
      element_idx = 3
    case (7) ! Nitrogen
      element_idx = 4
    case (8) ! Oxygen
      element_idx = 5
    case (10) ! Neon
      element_idx = 6
    case (12) ! Magnesium
      element_idx = 7
    case (14) ! Silicon
      element_idx = 8
    case (16) ! Sulfur
      element_idx = 9
    case (26) ! Iron
      element_idx = 10
  end select

    ! PISN
    if (mass.ge.140.d0) then

       ! Interpolate
       do i = 1,5
          if (mass.ge.table_mass_POP3_PISN(i).and.mass.lt.table_mass_POP3_PISN(i+1)) then
             j = i
          end if
       end do

       frac_high = (mass - table_mass_POP3_PISN(j)) / (table_mass_POP3_PISN(j+1) - table_mass_POP3_PISN(j))
       frac_low = 1.d0 - frac_high

       fq = (frac_low * PopIII_PISN(element_idx,j)) + (frac_high * PopIII_PISN(element_idx,j+1))

       return
    end if

    ! SNII or HN
    if (is_HN.and.mass.ge.20.d0) then

       ! Interpolate
       do i = 1,3
          if (mass.ge.table_mass_POP3_HN(i).and.mass.lt.table_mass_POP3_HN(i+1)) then
             j = i
          end if
       end do

       frac_high = (mass - table_mass_POP3_HN(j)) / (table_mass_POP3_HN(j+1) - table_mass_POP3_HN(j))
       frac_low = 1.d0 - frac_high

       fq = (frac_low * PopIII_HN(element_idx,j)) + (frac_high * PopIII_HN(element_idx,j+1))
    else
       loc_mass = MAX(mass,table_mass_POP3_SNII(1))

       ! Interpolate
       j = -999
       do i = 1,6
          if (loc_mass.ge.table_mass_POP3_SNII(i).and.loc_mass.lt.table_mass_POP3_SNII(i+1)) then
             j = i
          end if
       end do

       frac_high = (loc_mass - table_mass_POP3_SNII(j)) / (table_mass_POP3_SNII(j+1) - table_mass_POP3_SNII(j))
       frac_low = 1.d0 - frac_high

       fq = (frac_low * PopIII_SNII(element_idx,j)) + (frac_high * PopIII_SNII(element_idx,j+1))

       ! Scale down if below mass limit
       fq = fq * (mass / loc_mass)
    end if

END FUNCTION get_popIII_ejecta

FUNCTION get_SNIa_ejecta(species) result(ejecta)
  ! Ejecta from Seitenzahl 2013 - N100 model
  implicit none
  integer, intent(in) :: species
  real(dp) :: ejecta 

  ejecta = 0.d0

  select case (species)
    case (1) ! Hydrogen
      ejecta = 0.d0
    case (2) ! Helium
      ejecta = 0.d0
    case (6) ! Carbon
      ejecta = 3.04d-03
    case (7) ! Nitrogen
      ejecta = 3.21d-06
    case (8) ! Oxygen
      ejecta = 1.01d-01
    case (10) ! Neon
      ejecta = 3.57d-03
    case (11) ! Sodium
      ejecta = 3.74d-05
    case (12) ! Magnesium
      ejecta = 1.54d-02
    case (13) ! Aluminum
      ejecta = 6.74d-04
    case (14) ! Silicon
      ejecta = 2.87d-01
    case (16) ! Sulfur
      ejecta = 1.15d-01
    case (17) ! Chlorine
      ejecta = 2.13d-04
    case (18) ! Argon
      ejecta = 1.96d-02
    case (20) ! Calcium
      ejecta = 1.48d-02
    case (26) ! Iron
      ejecta = 7.40d-01
  end select

END FUNCTION get_SNIa_ejecta

FUNCTION get_portinari_ejecta_mass(metallicity, mass, species) result(fq)
  implicit none
  real(dp), intent(in) :: metallicity, mass
  integer, intent(in) :: species
  real(dp):: xq, yq
  real(dp) :: fq
  integer :: i, j, element_idx
  real(dp) :: x1, x2, y1, y2
  real(dp) :: f11, f21, f12, f22
  real(dp) :: dx, dy
  real(dp) :: rescale_factor

  rescale_factor = 1.0 ! Rescale factor for metallicity and mass

  ! Enforce bounds
  xq = MIN(MAX(metallicity,portinati_metal(1)),portinati_metal(portinati_N_metal))
  yq = MIN(MAX(mass,portinati_mass(1)),portinati_mass(portinati_N_mass))

  ! Downweight the yields if we are below lower limits --> do not extrapolate in the case where we are above!
  if (metallicity.lt.portinati_metal(1)) then
     rescale_factor = rescale_factor * (metallicity / portinati_metal(1))
  endif

  if (mass.lt.portinati_mass(1)) then
     rescale_factor = rescale_factor * (mass / portinati_mass(1))
  endif

  ! Find i such that x(i) <= xq <= x(i+1)
  do i = 1, portinati_N_metal-1
    if (xq.ge.portinati_metal(i) .and. xq.le.portinati_metal(i+1)) exit
  end do

  ! Find j such that y(j) <= yq <= y(j+1)
  do j = 1, portinati_N_mass-1
    if (yq.ge.portinati_mass(j) .and. yq.le.portinati_mass(j+1)) exit
  end do

  ! Rescale factors from Appendix A3 of https://articles.adsabs.harvard.edu/pdf/2009MNRAS.399..574W
  select case (species)
    case (1) ! Hydrogen
      element_idx = 1
    case (2) ! Helium
      element_idx = 2
    case (6) ! Carbon
      element_idx = 3
      rescale_factor = rescale_factor * 0.5d0
    case (7) ! Nitrogen
      element_idx = 4
    case (8) ! Oxygen
      element_idx = 5
    case (10) ! Neon
      element_idx = 6
    case (12) ! Magnesium
      element_idx = 7
      rescale_factor = rescale_factor * 2.d0
    case (14) ! Silicon
      element_idx = 8
    case (16) ! Sulfur
      element_idx = 9
    case (20) ! Calcium
      element_idx = 11
    case (26) ! Iron
      element_idx = 10
      rescale_factor = rescale_factor * 0.5d0
    case default
      write(*,*) "Element not available in portinari yields", species
      stop
  end select

  ! Extract corner values
  x1 = portinati_metal(i);   x2 = portinati_metal(i+1)
  y1 = portinati_mass(j);   y2 = portinati_mass(j+1)

  f11 = portinati_yields(i, j, element_idx)
  f21 = portinati_yields(i+1, j, element_idx)
  f12 = portinati_yields(i, j+1, element_idx)
  f22 = portinati_yields(i+1, j+1, element_idx)

  dx = x2 - x1
  dy = y2 - y1

  ! Bilinear interpolation
  fq = (1.0 / (dx * dy)) * ( &
        f11 * (x2 - xq) * (y2 - yq) + &
        f21 * (xq - x1) * (y2 - yq) + &
        f12 * (x2 - xq) * (yq - y1) + &
        f22 * (xq - x1) * (yq - y1)   &
       ) * rescale_factor
END FUNCTION get_portinari_ejecta_mass

FUNCTION get_portinari_stellar_lifetime(metallicity, mass) result(fq)
  implicit none
  real(dp), intent(in) :: metallicity, mass
  real(dp):: xq, yq
  real(dp) :: fq
  integer :: i, j
  real(dp) :: x1, x2, y1, y2
  real(dp) :: f11, f21, f12, f22
  real(dp) :: dx, dy

  ! Enforce bounds
  xq = MIN(MAX(metallicity,portinati_metal(1)),portinati_metal(portinati_N_metal))
  yq = MIN(MAX(mass,portinati_mass(1)),portinati_mass(portinati_N_mass))

  ! Find i such that x(i) <= xq <= x(i+1)
  do i = 1, portinati_N_metal-1
    if (xq.ge.portinati_metal(i) .and. xq.le.portinati_metal(i+1)) exit
  end do

  ! Find j such that y(j) <= yq <= y(j+1)
  do j = 1, portinati_N_mass-1
    if (yq.ge.portinati_mass(j) .and. yq.le.portinati_mass(j+1)) exit
  end do

  ! Extract corner values
  x1 = portinati_metal(i);   x2 = portinati_metal(i+1)
  y1 = portinati_mass(j);   y2 = portinati_mass(j+1)

  f11 = portinati_lifetimes(i, j)
  f21 = portinati_lifetimes(i+1, j)
  f12 = portinati_lifetimes(i, j+1)
  f22 = portinati_lifetimes(i+1, j+1)

  dx = x2 - x1
  dy = y2 - y1

  ! Bilinear interpolation
  fq = (1.0 / (dx * dy)) * ( &
        f11 * (x2 - xq) * (y2 - yq) + &
        f21 * (xq - x1) * (y2 - yq) + &
        f12 * (x2 - xq) * (yq - y1) + &
        f22 * (xq - x1) * (yq - y1)   &
       ) * ((yq/mass)**(-1.4d0)) ! Last factor used to extrapolate beyond 9 solar mass limit --> empirically found by harley to put 8 Msol star at 40 Myr lifetime
END FUNCTION get_portinari_stellar_lifetime

FUNCTION get_stellar_lifetime_low_mass(metallicity, mass) result(fq)
   ! Simple function to get the main-sequence lifetime of low-mass stars
   ! We simply scale the lifetime of an 8 Msol star by M^-2.5
   ! --> TODO(code): eventually replace this with something more accurate
   implicit none
   real(dp), intent(in) :: metallicity, mass
   real(dp) :: fq
   real(dp) :: main_sequence_lifetime_8Msun

   ! Main sequence lifetime of an 8 Msol star of the same metallicity
   main_sequence_lifetime_8Msun = get_portinari_stellar_lifetime(metallicity, 8.d0)

   ! Use homology to scale the lifetime
   fq = main_sequence_lifetime_8Msun * ((8.d0 / mass)**(2.5d0))

END FUNCTION

FUNCTION getStarAgeMyr(birth_time) result(age_star)
   use amr_commons, only: dp, use_proper_time
   implicit none
   real(dp), intent(in):: birth_time
   real(dp):: age_star
   real(dp):: age_star_pt

   if (use_proper_time) then
      call getAgeGyr(birth_time, age_star)
   else
      call getProperTime(birth_time, age_star_pt)
      call getAgeGyr(age_star_pt, age_star)
   end if

   ! Convert from Gyr to Myr
   age_star = age_star * 1.d3

END FUNCTION getStarAgeMyr
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
FUNCTION get_pop3_lifetime_myr(mass) result(age_myr)
   use amr_commons, only: dp
   use constants, only: M_sun
   implicit none
   real(dp), intent(in)::mass
   real(dp)::age_myr
   real(dp):: age_gyr, logm
   real(dp), dimension(0:3)::a_fit
   integer::ndeg = 3, i
   real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v, scale_msun

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_msun = scale_l**3*scale_d/M_sun

   a_fit = (/0.7595398e+00, -3.7303953e+00, 1.4031973e+00, -1.7896967e-01/)

   age_gyr = 0d0
   logm = log10(mass*scale_msun)
   do i = 0, ndeg
      age_gyr = age_gyr + a_fit(i)*logm**dble(i)
   end do
   age_gyr = 10d0**age_gyr
   age_myr = age_gyr * 1.d3

END FUNCTION get_pop3_lifetime_myr
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
FUNCTION get_PISN_energy(mass) result(sn_e)
  ! PISN energy interpolating tables from https://arxiv.org/pdf/1511.03040
  ! Here we use the Case A data
  use amr_commons, only: dp
  real(dp), intent(in):: mass
  real(dp):: sn_e
  real(dp), dimension(1:9):: PISN_e = (/ 11.30d0, 14.05d0, 16.94d0, 18.50d0, 33.46d0, 42.58d0, 53.88d0, 56.08d0, 81.91d0 /)
  real(dp), dimension(1:9):: PISN_m = (/ 145.d0, 150.d0, 155.d0, 160.d0, 180.d0, 200.d0, 220.d0, 240.d0, 260.d0 /) 
  real(dp):: mass_loc, frac_low, frac_high
  integer:: i,j

  ! Enforce bounds
  mass_loc = MAX(MIN(mass,PISN_m(9)),PISN_m(1))

  ! Interpolate
  do i=1,8
     if (i.ge.PISN_m(i).and.i.lt.PISN_m(i+1)) then
        j = i
     end if
  end do

  frac_high = (mass_loc - PISN_m(j)) / (PISN_m(j+1) -  PISN_m(j))
  frac_low = 1.d0 - frac_high

  sn_e = (frac_low * PISN_e(j)) + (frac_high * PISN_e(j+1))

  ! Interpolate above and below bounds
  sn_e = sn_e * (mass / mass_loc)

END FUNCTION get_PISN_energy
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
FUNCTION get_popIII_sn_energy(mass,particle_id) result(sn_e)
  use amr_commons, only: dp
  implicit none
  real(dp), intent(in):: mass
  integer, intent(in):: particle_id
  real(dp):: sn_e

  if (mass.ge.140.d0.and.mass.le.300.d0) then 
     ! PISN
     sn_e = get_PISN_energy(mass) * 1.d51
  else if (mass.ge.10.d0.and.mass.le.40.d0) then
     ! Normal SN
     if (mass.lt.20.d0) then
        sn_e = 1.d51
     ! Randomly sample normal SN and HN
     else
        ! 50% of these should be HN --> just use the ID mod 2
        if (MOD(particle_id,2).eq.0) then 
           if (mass.lt.25.d0) then
              sn_e = 1.d52
           else
              sn_e = 1.3333d0 * mass - 23.3333d0
              sn_e = sn_e * 1.d51
           end if
        else
           sn_e = 1.d51
        end if

     end if
  else
     ! This is the direct collapse case
     sn_e = 0.d0
  end if

END FUNCTION get_popIII_sn_energy


end module metal_yields_module