! cross_sections_module.f90
module cross_sections_module
  use amr_parameters, only: dp
  implicit none

  private
  public::initialize_cross_sections, getCrosssection_rtz

  ! These are the nl-resolved cross sections
  ! They are needed for higher enery photons (e.g. x-rays)
  real(dp),dimension(1:27,1:27,1:7,1:7)::verner_cross_sections_nl
  ! These are the total cross sections
  real(dp),dimension(1:27,1:27,1:9)::verner_cross_sections_tot

CONTAINS

SUBROUTINE initialize_cross_sections()
  ! Photoionization cross sections from
  ! https://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
  implicit none

  integer :: a, b, c, d
  real(dp) :: x1, x2, x3, x4, x5, x6, x7, x8, x9, l_state
  integer :: ios, nl_state
  integer :: i1, i2

  ! Initialize cross sections to large negative numbers
  verner_cross_sections_tot = -999.d0

  ! Load the total cross sections
  open(unit=10, file='./data/verner_cross_sections/phfit_tot.dat', status='old', action='read')
  do
     read(10, *, iostat=ios) i1, i2, x1, x2, x3, x4, x5, x6, x7, x8, x9
     if (ios /= 0) exit
     
     if (i1.lt.27) then 
        verner_cross_sections_tot(i1,i1-i2+1,1) = x1 ! E_th
        verner_cross_sections_tot(i1,i1-i2+1,2) = x2 ! E_max
        verner_cross_sections_tot(i1,i1-i2+1,3) = x3 ! E_0
        verner_cross_sections_tot(i1,i1-i2+1,4) = x4 ! sigma_0
        verner_cross_sections_tot(i1,i1-i2+1,5) = x5 ! y_a
        verner_cross_sections_tot(i1,i1-i2+1,6) = x6 ! P
        verner_cross_sections_tot(i1,i1-i2+1,7) = x7 ! y_w
        verner_cross_sections_tot(i1,i1-i2+1,8) = x8 ! y_0
        verner_cross_sections_tot(i1,i1-i2+1,9) = x9 ! y_1
     end if
  end do

  ! Initialize cross sections to large negative numbers
  verner_cross_sections_nl = -999.d0

  ! Load the n,l - resolved cross sections
  open(unit=10, file='./data/verner_cross_sections/phfit.dat', status='old', action='read')

  do
     read(10, *, iostat=ios) a, b, c, d, x1, x2, x3, x4, x5, x6
     if (ios /= 0) exit

     ! Select only up to iron
     if (a.lt.27) then
         if (c.eq.1 .and. d.eq.0) then      !1s
            nl_state = 1
         else if (c.eq.2 .and. d.eq.0) then !2s
            nl_state = 2
         else if (c.eq.2 .and. d.eq.1) then !2p
            nl_state = 3
         else if (c.eq.3 .and. d.eq.0) then !3s
            nl_state = 4
         else if (c.eq.3 .and. d.eq.1) then !3p
            nl_state = 5
         else if (c.eq.3 .and. d.eq.2) then !3d
            nl_state = 6
         else if (c.eq.4 .and. d.eq.0) then !4s
            nl_state = 7
         endif
         verner_cross_sections_nl(a,a-b+1,nl_state,1) = x1 ! E_th
         verner_cross_sections_nl(a,a-b+1,nl_state,2) = x2 ! E_0
         verner_cross_sections_nl(a,a-b+1,nl_state,3) = x3 ! sigma_0
         verner_cross_sections_nl(a,a-b+1,nl_state,4) = x4 ! y_a
         verner_cross_sections_nl(a,a-b+1,nl_state,5) = x5 ! P
         verner_cross_sections_nl(a,a-b+1,nl_state,6) = x6 ! y_w
         verner_cross_sections_nl(a,a-b+1,nl_state,7) = REAL(d,dp) ! l
     end if
  end do

  close(10)

END SUBROUTINE initialize_cross_sections

FUNCTION F_of_y(y, yw, ya, P, Q) result(Fy)
    implicit none
    real(dp), intent(in):: y, yw, ya, P, Q
    real(dp)::Fy

    Fy = (((y - 1d0)**2d0) + (yw**2.)) * (y**(-Q)) * ((1d0 + SQRT(y/ya))**(-P))
END FUNCTION F_of_y

FUNCTION verner_cs_nl(E, sig_0, E_0, yw, ya, P, l) result(cs_nl)
    implicit none
    real(dp), intent(in):: E, sig_0, E_0, yw, ya, P, l
    real(dp)::cs_nl
    real(dp):: y, Q

    y = E / E_0
    Q = 5.5d0 + l - (0.5d0 * P)
    cs_nl =  sig_0 * F_of_y(y,yw,ya,P,Q)

END FUNCTION verner_cs_nl

FUNCTION getCrosssection_rtz(lambda, element, ion) result(cross_sec)
   use constants, only: eV2erg, c_cgs, hplanck
   implicit none
   real(KIND=8), intent(in)::lambda
   integer, intent(in)::element, ion
   real(KIND=8):: cross_sec
   real(KIND=8) :: x, y, F, E
   integer :: i
   real(KIND=8) :: E_0, sig_0, ya, P, yw, l 

   ! Initialize cross section to 0
   cross_sec = 0.d0

   ! Convert lambda into ev
   E = hplanck * c_cgs/(lambda*1.d-8) / eV2erg         ! photon energy in ev

   ! Deal with molecular hydrogen separately
   if (element.eq.1 .and. ion.eq.3) then
      cross_sec = 0.
      if (E .gt. 11.20 .and. E .lt. 13.59) cross_sec = 2.47d-18
      if (E .ge. 13.59 .and. E .le. 15.21) cross_sec = 0.0d0
      if (E .gt. 15.21 .and. E .le. 15.45) cross_sec = 0.09d-18
      if (E .gt. 15.70 .and. E .le. 15.95) cross_sec = 1.15d-18
      if (E .gt. 15.95 .and. E .le. 16.20) cross_sec = 3.00d-18
      if (E .gt. 16.20 .and. E .le. 16.40) cross_sec = 5.00d-18
      if (E .gt. 16.40 .and. E .le. 16.65) cross_sec = 6.75d-18
      if (E .gt. 16.65 .and. E .le. 16.85) cross_sec = 8.00d-18
      if (E .gt. 16.85 .and. E .le. 17.00) cross_sec = 9.00d-18
      if (E .gt. 17.00 .and. E .le. 17.20) cross_sec = 9.50d-18
      if (E .gt. 17.20 .and. E .le. 17.65) cross_sec = 9.80d-18
      if (E .gt. 17.65 .and. E .le. 18.10) cross_sec = 10.10d-18
      if (E .gt. 18.10) cross_sec = (10.10d-18) * ((18.10d0/E)**3.0d0)
      return
   end if

   x = (E / verner_cross_sections_tot(element,ion,3)) - verner_cross_sections_tot(element,ion,8)
   y = sqrt( (x*x) + (verner_cross_sections_tot(element,ion,9)**2.d0) )

   F = (x - 1.d0)**2
   F = F + (verner_cross_sections_tot(element,ion,7)**2.d0)
   F = F * (y**(0.5d0*verner_cross_sections_tot(element,ion,6) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/verner_cross_sections_tot(element,ion,5)))**(-1.d0*verner_cross_sections_tot(element,ion,6)))

   cross_sec = verner_cross_sections_tot(element,ion,4) * F * 1.d-18
   if (E.lt.verner_cross_sections_tot(element,ion,1)) then
      cross_sec = 0.d0
   endif

   ! Above E_max, use the n,l resolved cross sections
   if (E.gt.verner_cross_sections_tot(element,ion,2)) then
      cross_sec = 0.d0
      do i = 1,7
         if (verner_cross_sections_nl(element,ion,i,7) .ge. 0.d0) then 
            E_0   = verner_cross_sections_nl(element,ion,i,2)
            sig_0 = verner_cross_sections_nl(element,ion,i,3)
            ya    = verner_cross_sections_nl(element,ion,i,4)
            P     = verner_cross_sections_nl(element,ion,i,5)
            yw    = verner_cross_sections_nl(element,ion,i,6)
            l     = verner_cross_sections_nl(element,ion,i,7)
            cross_sec = cross_sec + verner_cs_nl(E, sig_0, E_0, yw, ya, P, l)
         end if
      end do
      cross_sec = cross_sec * 1.d-18 ! Convert to cm^-2
   endif

END FUNCTION getCrosssection_rtz

end module cross_sections_module