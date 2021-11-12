
program simple_md
   use forces_mod
   implicit none

   integer, parameter :: N_PARTICLES = 2
   integer, parameter :: N_STEP = 1000
   integer, parameter :: N_PRINT = 2
   real(8), parameter :: DT = 20.0D0
   ! Add a logical variable

   real(8), parameter :: AMU = 1823.0D0
   real(8), dimension(N_PARTICLES) :: mass = 1.008D0 * AMU 
   real(8), dimension(N_PARTICLES) :: x, y
   real(8), dimension(N_PARTICLES) :: vx, vy
   real(8), dimension(N_PARTICLES) :: fx, fy

   real(8) :: potential_energy, kinetic_energy
   integer :: istep
   integer :: i
   
   ! Defining inital conditions
   ! ================================
   ! Let's start with two particles bound by harmonic potential
   x = 0.0D0
   y = (/ 0.0D0, 0.9D0 /)

   fx = 0.0D0
   fy = 0.0D0
   vx = 0.0D0
   vy = 0.0D0

   potential_energy = 0.0D0

   ! ================================

   ! MD LOOP, using Velocity Verlet integrator
   do istep = 1, N_STEP
      
      ! Implicit iteration here, FORTRAN is smart :-)
      ! We can do math on arrays as long as they have the same dimension
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      x = x + vx * DT
      y = y + vy * DT

      call calculate_forces(x, y, fx, fy, potential_energy)
      
      vx = vx + fx / mass * DT * 0.5D0
      vy = vy + fy / mass * DT * 0.5D0

      kinetic_energy = 0.0D0
      do i=1, N_PARTICLES
         kinetic_energy = kinetic_energy + 0.5D0 * mass(i) * (vx(i)**2 + vy(i)**2)
      end do

      !write (*, *) istep, potential_energy, potential_energy + kinetic_energy
      write (*, *) istep, get_distance(x, y, 1, 2)

   end do

end program
