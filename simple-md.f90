! NOTE: We use atomic units throughout the program
program simple_md
   implicit none

   integer, parameter :: N_PARTICLES = 2
   integer, parameter :: N_STEP = 100
   integer, parameter :: N_PRINT = 1
   real(8), parameter :: DT = 20.0D0
   real(8), parameter :: AMU = 1823.0D0
   real(8), dimension(N_PARTICLES) :: mass
   real(8), dimension(N_PARTICLES) :: x, y
   real(8), dimension(N_PARTICLES) :: vx, vy
   real(8), dimension(N_PARTICLES) :: fx, fy

   character(len=20) :: potential = 'harmonic'
   ! Bond legth of N2 molecule
   ! https://cccbdb.nist.gov/exp2x.asp?casno=7727379
   real(8) :: r_eq = 1.098D0
   ! harmonic force constant in a.u. equivalent to 2359 cm^-1 vibration
   real(8) :: force_constant = 1.47497D0
   ! TODO:
   real(8) :: dissociation_energy = 1.0

   real(8) :: potential_energy, kinetic_energy
   real(8) :: time 
   real(8) :: calculate_kinetic_energy
   integer :: istep
   ! File units for output
   integer :: uenergy, ucoords
   integer :: i
   
   ! Defining inital conditions
   ! ================================
   ! Let's start with two particles bound by harmonic potential
   ! Specifically, let's simulate N2 molecule.
   x = 0.0D0
   y = (/ 0.0D0, 1.1D0 /)
   mass(1) = 14.0067 * AMU
   mass(2) = 14.0067 * AMU

   fx = 0.0D0
   fy = 0.0D0
   vx = 0.0D0
   vy = 0.0D0

   potential_energy = 0.0D0

   ! ================================
   open (newunit=uenergy, file='energy.dat', action='write')
   open (newunit=ucoords, file='coords.xyz', action='write')

   ! MD LOOP, using Velocity Verlet integrator
   do istep = 1, N_STEP
      
      ! Implicit iteration here
      ! We can do simple matrix math with arrays as long as they have the same size
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      x = x + vx * DT
      y = y + vy * DT

      if (potential == 'harmonic') then
         call harmonic_forces(force_constant, r_eq, x, y, fx, fy, potential_energy)
      else if (potential == 'morse') then
         call morse_forces(force_constant, dissociation_energy, r_eq, x, y, fx, fy, potential_energy)
      else
         print *, 'Invalid potential'
         stop 1
      end if
      
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      if (modulo(istep, N_PRINT) == 0) then
         kinetic_energy = calculate_kinetic_energy(N_PARTICLES, vx, vy, mass)
         time = DT * (istep - 1)
         write (uenergy, *) time, potential_energy, potential_energy + kinetic_energy
         call write_coords(N_PARTICLES, x, y, ucoords)
      end if

   end do

   close(uenergy)
   close(ucoords)
end program

real(8) function calculate_kinetic_energy(n_particles, vx, vy, mass) result(e_k)
   implicit none
   integer, intent(in) :: n_particles
   real(8), dimension(n_particles), intent(in) :: vx, vy, mass
   integer :: i
   e_k = 0
   do i = 1, N_PARTICLES
      e_k = e_k + 0.5D0 * mass(i) * (vx(i)**2 + vy(i)**2)
   end do
end function

subroutine write_coords(natoms, x, y, channel)
   implicit none
   integer, intent(in) :: natoms
   real(8), dimension(natoms), intent(in) :: x, y
   integer, intent(in) :: channel
   integer :: i

   write (channel, *) natoms
   write (channel, *) 'comment line'
   do i = 1, natoms
      write (channel, '(A,3F10.5)') "N", x(i), y(i), 0.0D0
   end do
end subroutine write_coords

subroutine morse_forces(k, de, r_eq, x, y, fx, fy, potential_energy)
   implicit none
   ! Parameters of the Morse potential
   real(8), intent(in) :: k, de, r_eq
   real(8), dimension(2), intent(in) :: x, y
   real(8), dimension(2), intent(out) :: fx, fy
   real(8), intent(out) :: potential_energy
   write (*, *) 'ERROR: Morse potential not implemented'
   stop 1
end subroutine morse_forces

subroutine harmonic_forces(k, r_eq, x, y, fx, fy, potential_energy)
   implicit none
   real(8), intent(in) :: k, r_eq
   real(8), dimension(2), intent(in) :: x, y
   real(8), dimension(2), intent(out) :: fx, fy
   real(8), intent(out) :: potential_energy
   real(8) :: r, delta_x, delta_y

   delta_x = x(1) - x(2)
   delta_y = y(1) - y(2)
   r = delta_x**2 + delta_y**2
   r = dsqrt(r)

   fx(1) = - k * delta_x * (r - r_eq) / r
   fx(2) = - fx(1)
   fy(1) = - k * delta_y * (r - r_eq) / r
   fy(2) = - fy(1)

   potential_energy = 0.5D0 * k * (r - r_eq)**2
end subroutine harmonic_forces

function get_distance(x, y) result(r)
   implicit none
   real(8), dimension(2), intent(in) :: x, y
   real(8) :: r
   real(8) :: dx, dy

   dx = x(1) - x(2)
   dy = y(1) - y(2)
   r = dx**2 + dy**2
   r = dsqrt(r)
end function get_distance
