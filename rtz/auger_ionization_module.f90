! auger_ionization_module.f90
module auger_ionization_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: load_auger_yields, get_auger_yields

  real(dp):: auger_yield(26,26,10,300)
  real(dp):: auger_energies(300)

CONTAINS

SUBROUTINE load_auger_yields()
  ! Loads the auger yields and energies from files into the arrays
  use amr_commons, only: myid
  implicit none

  integer:: i, j, k, unit_num, ios
  character(len=50) :: filename

  if(myid.eq.1) write(*,*) 'Initializing auger yields'

  ! Load charge exchange file
  open(newunit=unit_num, file='./data/auger_data/energies.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open auger energies file'
      return
  end if

  ! Reading data from the file into the array
  do i = 1, 300
     read(unit_num, *, iostat=ios) auger_energies(i)
     if (ios /= 0) exit
  end do
  close(unit_num)

  ! Inialize the auger yields to 0
  auger_yield = 0.d0
  do i = 3, 26 ! Loop over elements
     do j = 1, 26 ! Loop over ionization states
        if (j .le. i) then 
           ! Try to open the file
           write(filename, '("./data/auger_data/",I0,"_",I0,"_yields.dat")') i, j
           open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=ios)

           if (ios /= 0) then
              write(*,*) 'Error: Could not open auger file',filename
              call clean_stop
           end if

           ! Read data into the array
           do k = 1, 300
              read(unit_num, *, iostat=ios) auger_yield(i,j,:,k)
           end do

           close(unit_num)
        end if
     end do
  end do

END SUBROUTINE load_auger_yields

SUBROUTINE get_auger_yields(iElement,n_tot_ions,group_egy,electron_probabilities)
  ! For a given energy, computes the probability of ejection of N electrons
  ! Requires photon energy in eV
  implicit none

  integer,intent(in):: iElement,n_tot_ions
  real(dp),dimension(1:NGROUPS),intent(in):: group_egy
  real(dp),dimension(1:27,1:10,1:NGROUPS),intent(inout):: electron_probabilities

  real(dp):: photon_energy, log_photon_energy, u_frac, cum_prob
  integer:: i, ig, energy_idx, idx_ignore

  ! Reset the electron probabilities to default to 1 electron
  electron_probabilities = 0.d0

  do ig = 1,NGROUPS
     ! Set the photon energy
     photon_energy = group_egy(ig)

     do i = 1,27
        electron_probabilities(i,1,ig) = 1.d0
     end do

     ! Make sure that photon energy is within the correct bounds
     if (photon_energy.lt.auger_energies(1) .or. photon_energy.gt.auger_energies(300)) then
        cycle
     end if

     ! Interpolate the photon energy in log space
     log_photon_energy = LOG10(photon_energy)

     ! Get the lower energy bin for the auger yields
     energy_idx = 1 + FLOOR((log_photon_energy - LOG10(auger_energies(1))) / (LOG10(auger_energies(2))-LOG10(auger_energies(1))))

     ! Interpolate the probabilities
     u_frac = (log_photon_energy - LOG10(auger_energies(energy_idx))) / (LOG10(auger_energies(energy_idx+1)) - LOG10(auger_energies(energy_idx)))

     ! Now loop over ionization states and fill out the table
     do i = 1,n_tot_ions
        electron_probabilities(i,:,ig) = auger_yield(iElement,i,:,energy_idx+1) * u_frac
        electron_probabilities(i,:,ig) = electron_probabilities(i,:,ig) + (auger_yield(iElement,i,:,energy_idx) * (1.d0 - u_frac))
     end do

     ! In the case where we don't follow every ion we need to update these tables accordingly or else the yields will be wrong
     ! So we just move the probability into the highest index we follow and set the rest to 0
     do i = 1,n_tot_ions
        idx_ignore = n_tot_ions - i
        if (idx_ignore.gt.0 .and. idx_ignore.lt.10) then 
           electron_probabilities(i,idx_ignore,ig) = SUM(electron_probabilities(i,idx_ignore:,ig))
           electron_probabilities(i,idx_ignore+1:,ig) = 0.d0
        end if
     end do
  end do

END SUBROUTINE get_auger_yields

end module auger_ionization_module