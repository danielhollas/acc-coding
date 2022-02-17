program simple_md
   use forces_mod
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

   real(8) :: potential_energy, kinetic_energy
   real(8) :: time 
   real(8) :: calculate_kinetic_energy
   integer :: istep
   integer :: i, funit
   
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
   open(newunit=funit, file='energy.dat', action='write', access='sequential')

   ! MD LOOP, using Velocity Verlet integrator
   do istep = 1, N_STEP
      
      ! Implicit iteration here, FORTRAN is smart
      ! We can do math on arrays as long as they have the same size
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      x = x + vx * DT
      y = y + vy * DT

      call calculate_forces(x, y, fx, fy, potential_energy)
      
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      if (modulo(istep, N_PRINT) == 0) then
         kinetic_energy = calculate_kinetic_energy(vx, vy, mass, N_PARTICLES)
         time = DT * (istep - 1)
         write (funit, *) time, potential_energy, potential_energy + kinetic_energy
      end if

   end do

   close(funit)
end program

real(8) function calculate_kinetic_energy(vx, vy, mass, n_particles) result(e_k)
   real(8), dimension(n_particles), intent(in) :: vx, vy, mass
   integer, intent(in) :: n_particles
   e_k = 0
   do i=1, N_PARTICLES
      e_k = e_k + 0.5D0 * mass(i) * (vx(i)**2 + vy(i)**2)
   end do
end function

subroutine print_trajectory(x, y, istep)
   ! Homework
end subroutine print_trajectory
