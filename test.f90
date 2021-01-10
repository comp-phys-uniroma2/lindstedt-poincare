program test
  use precision
  use constants
  use initialization
  use functions
  use solvers
  
  implicit none

  real(dp), dimension(2) :: u, u0
  real(dp) :: t, dt
  integer :: i, Nstep, funit, iargc
  character(50) :: arg, solname

  call init(solname, Nstep, dt, t, u0, u)
  
  open (newunit=funit, file="sol.dat")
  write(funit,*) "time  ", "x(t)  ", "v(t)  " 

  select case(trim(solname))
     case("rk4")
        do i = 1, Nstep
           call rk4(harmonic, t, dt, u0, u)
           u0 = u
           t = t + dt
           write(funit,*) t, u(1), u(2)
        end do

     case("ab4")
        call ab4(harmonic, Nstep, t, dt, u0, u, funit)

     case default
        stop
     end select
     
  close(funit)
  
end program test
