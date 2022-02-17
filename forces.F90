
module forces_mod
   implicit none
   character(len=20) :: potential = 'harmonic'
   ! Bond legth of N2 molecule
   ! https://cccbdb.nist.gov/exp2x.asp?casno=7727379
   real(8) :: r_eq = 1.098D0
   ! harmonic force constant in a.u. equivalent to 2359 cm^-1 vibration
   real(8) :: force_constant = 1.47497D0
contains

   subroutine calculate_forces(x, y, fx, fy, potential_energy)
      real(8), dimension(:), intent(in) :: x, y
      real(8), dimension(:), intent(out) :: fx, fy
      real(8), intent(out) :: potential_energy
      real(8) :: r

      if (potential == 'harmonic') then
         call harmonic_forces(x, y, fx, fy, potential_energy)
      else
         call lennard_jones_forces(x, y, fx, fy, potential_energy)
      end if
   end subroutine calculate_forces

   subroutine lennard_jones_forces(x, y, fx, fy, potential_energy)
      real(8), dimension(:), intent(in) :: x, y
      real(8), dimension(:), intent(out) :: fx, fy
      real(8), intent(out) :: potential_energy
      write (*, *) 'ERROR: LJ potential not implemented'
      stop 1
   end subroutine lennard_jones_forces

   subroutine harmonic_forces(x, y, fx, fy, potential_energy)
      real(8), dimension(:), intent(in) :: x, y
      real(8), dimension(:), intent(out) :: fx, fy
      real(8), intent(out) :: potential_energy
      real(8) :: r, delta_x, delta_y

      delta_x = x(1) - x(2)
      delta_y = y(1) - y(2)
      r = delta_x**2 + delta_y**2
      r = dsqrt(r)

      fx(1) = - force_constant * delta_x * (r - r_eq) / r
      fx(2) = - fx(1)
      fy(1) = - force_constant * delta_y * (r - r_eq) / r
      fy(2) = - fy(2)

      potential_energy = 0.5D0 * force_constant * (r - r_eq)**2
   end subroutine harmonic_forces

   function get_distance(x, y, particle1, particle2) result(r)
      real(8), dimension(:), intent(in) :: x, y
      integer, intent(in) :: particle1, particle2
      real(8) :: r
      real(8) :: dx, dy

      dx = x(particle1) - x(particle2)
      dy = y(particle1) - y(particle2)
      r = dx**2 + dy**2
      r = dsqrt(r)
   end function get_distance

end module forces_mod
