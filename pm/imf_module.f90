MODULE imf_module
   ! Module for dealing with IMF calculations
   use amr_parameters, only: dp
   implicit none

   public :: sample_IMF_pop3, sample_IMF_pop2, get_SNIa_prob
CONTAINS

FUNCTION sample_IMF_pop3(seed) result(mass)
   !! Samples pop III IMF
   use pm_parameters, only: p3_mchar
   use pm_commons, ONLY: localseed
   use random
   implicit none
   integer, intent(in):: seed
   real(dp):: mass
   real(dp):: total_mass
   real(dp):: m_min, m_max, d_log10m
   real(dp):: total_probability, local_probability
   real(dp)::RandNum
   integer:: n_samples = 1000
   integer:: i

   ! Min and Max of the IMF
   m_min = 1.d0
   m_max = 1000.d0

   ! Get the total_probability
   d_log10m = (LOG10(m_max) - LOG10(m_min)) / n_samples

   ! Get the total_probability
   total_probability = 0.d0
   do i=1,n_samples
      mass = 10.d0 ** (LOG10(m_min) + d_log10m * real(i-1,kind=dp))
      total_probability = total_probability + (mass**-2.3d0) * exp(-(p3_mchar/mass)**1.6d0)
   end do

   ! Loop again to the mass
   local_probability = 0.0
   mass = 0.0
   RandNum = rand_from_seed(seed)
   i = 1
   do while (local_probability.lt.RandNum)
      mass = 10.d0 ** (LOG10(m_min) + d_log10m * real(i-1,kind=dp))
      local_probability = local_probability + ((mass**-2.3d0) * exp(-(p3_mchar/mass)**1.6d0))/total_probability
      i = i + 1
   end do

END FUNCTION sample_IMF_pop3

FUNCTION sample_IMF_pop2(m0,m1,m2,a1,a2,group_mass,lp_mass,seed) result(my_mass)
   ! Samples the IMF mass
   ! Note that the group mass is the mass below which we group all stars
   use random
   use pm_commons, only:localseed
   implicit none
   real(dp), intent(in)::m0,m1,m2,a1,a2,lp_mass,group_mass
   integer, intent(in)::seed
   real(dp)::my_mass
   real(dp)::smooth_factor,total_mass,mean_mass_above_break,N_draw_single
   real(dp)::frac_above_break,frac_below_break,total_Nstars_above_break
   integer::i
   real(dp):: loc_rand,loc_rand2

   if (m1 .gt. group_mass) then
      write(*,*) "M1 cannot be > group mass"
      call clean_stop()
      return
   end if

   ! Get the power-law transition between the two regimes
   smooth_factor = (m1**a1) / (m1**a2)

   ! Calculate the total mass in the IMF
   total_mass = (1.d0/(a1 + 2.d0)) * ((m1**(a1 + 2.d0)) - (m0**(a1 + 2.d0)))
   total_mass = total_mass + ( (smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (m1**(a2 + 2.d0))) )

   ! Calculate the mass fraction above and below the break
   frac_above_break = ( (smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (group_mass**(a2 + 2.d0))) )
   total_Nstars_above_break = ( (smooth_factor/(a2 + 1.d0)) * ((m2**(a2 + 1.d0)) - (group_mass**(a2 + 1.d0))) )

   frac_above_break = frac_above_break/total_mass
   frac_below_break = 1.d0 - frac_above_break

   ! Calculate the mean mass above the break
   mean_mass_above_break = (frac_above_break * total_mass) / total_Nstars_above_break
   ! Calculate the ratio of the lp mass to the mean_mass above break

   ! frac_below_break represents the fraction of the mass in star particles that 
   ! is grouped into the low mass regime. However since we put a fixed mass at 
   ! lp_mass we form many at once. so the break fraction  
   N_draw_single = (lp_mass * frac_above_break) / (frac_below_break * mean_mass_above_break)

   loc_rand = rand_from_seed(seed)
   if (loc_rand <= 1.d0/(N_draw_single+1.d0)) then
      my_mass = -999.d0 ! Set to a large negative value --> needed later to differentiate
   else
      loc_rand2 = rand_from_seed(-seed)
      my_mass = ((loc_rand2 * total_Nstars_above_break / (smooth_factor/(a2 + 1.d0))) + (group_mass**(a2 + 1.d0)))**(1./(a2 + 1.d0))
   endif

END FUNCTION sample_IMF_pop2

FUNCTION get_SNIa_prob(m0,m1,m2,a1,a2,group_mass,star_age,star_msl,seed) result(is_SNIa)
   ! Following https://iopscience.iop.org/article/10.3847/1538-4357/aa8b6e/pdf
   ! where they get N_SNIa/M = 1.3e-3 with a delay time distribution of
   ! t^-1.07
   use amr_parameters, only: dp, h0
   use constants
   implicit none
   real(dp), intent(in)::m0,m1,m2,a1,a2,group_mass
   real(dp), intent(in)::star_age,star_msl
   integer, intent(in)::seed
   logical::is_SNIa
   real(dp)::smooth_factor,total_mass,m_low
   real(dp)::loc_rand,loc_rand2,prob,total_snia_prob
   real(dp)::t0,tnow,tH
   real(dp)::un_normalized_prob,normalized_prob,current_prob

   is_SNIa = .false.

   ! If we don't track low mass stars, then no SNIa
   if (group_mass.gt.8.d0) then
      return
   end if

   ! Get the power-law transition between the two regimes
   smooth_factor = (m1**a1) / (m1**a2)   

   ! Calculate the total mass in the IMF
   total_mass = (1.d0/(a1 + 2.d0)) * ((m1**(a1 + 2.d0)) - (m0**(a1 + 2.d0)))
   total_mass = total_mass + ( (smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (m1**(a2 + 2.d0))) )

   ! Calculate the total number of stars within the SNIa mass range which we assume is 3 Msol - 8 Msol
   m_low = MAX(group_mass,3.d0) ! However, we don't let degenerate star particles explode
   total_SNIa_prob = (smooth_factor/(a2 + 1.d0)) * ((8.d0**(a2 + 1.d0)) - (m_low**(a2 + 1.d0)))

   ! Number of possible SNIa progenitors per solar mass of stars
   prob = (1.d0 / total_mass) * total_SNIa_prob

   ! Check if this specific star will go SNIa
   loc_rand = rand_from_seed(-2 * seed) ! Modify the seed so it's different from before

   ! These stars will go SNIa --> now just need to sample the delay time distribution
   if (loc_rand.le.(1.3d-3/prob)) then
      t0 = star_msl   ! Main sequence lifetime of the star Myr
      tnow = star_age ! Age of the star in Myr
      tH = 1.0 / (h0 / (Mpc2cm / 1e5)) ! Hubble time in Myr
      tH = tH / Myr2sec

      ! Draw a random number --> always the same for all CPUs at all times
      ! Offset so not all random numbers are the same
      loc_rand2 = rand_from_seed(seed + 2345432)

      ! Get the probability that at a given time the s
      un_normalized_prob = (1.d0/(-0.07d0)) * ((tH**(-0.07d0)) - (t0**(-0.07d0)))
      current_prob = (1.d0/(-0.07d0)) * ((tnow**(-0.07d0)) - (t0**(-0.07d0)))
      normalized_prob = current_prob / un_normalized_prob

      if (loc_rand2.lt.normalized_prob) then
         is_SNIa = .true.
      end if

   end if 

END FUNCTION get_SNIa_prob

FUNCTION rand_from_seed(seed) result(r)
    use iso_fortran_env, only: int64, real64
    use pm_parameters, only: uniform_rand_seed
    implicit none
    integer, intent(in) :: seed
    real(real64) :: r
    integer(int64) :: state
    integer(int64) :: word
    integer :: rot
    integer :: result
    
    ! PCG constants
    integer(int64), parameter :: PCG_MULTIPLIER = 6364136223846793005_int64
    integer(int64), parameter :: PCG_INCREMENT = 1442695040888963407_int64
    
    ! Initialize state from seed
    state = int(seed+uniform_rand_seed, int64)
    
    ! Advance the LCG state
    state = state * PCG_MULTIPLIER + PCG_INCREMENT
    
    ! Apply the PCG output function (XSH-RR variant)
    ! XSH: XOR the high and low bits of state, then shift
    word = ieor(ishft(state, -18), state)
    word = ishft(word, -27)
    
    ! RR: Rotate right by the top 5 bits of the original state
    rot = int(ishft(state, -59))
    
    ! Perform the rotation (emulate rotate-right)
    result = int(ior(ishft(int(word), -rot), ishft(int(word), 32-rot)))
    
    ! Convert to [0,1) treating as unsigned 32-bit
    if (result >= 0) then
        r = real(result, kind=real64) / 4294967296.0_real64
    else
        r = (real(result, kind=real64) + 4294967296.0_real64) / 4294967296.0_real64
    endif
END FUNCTION rand_from_seed

END MODULE imf_module