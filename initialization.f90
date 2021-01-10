module initialization
  use precision
  implicit none
  
contains
  subroutine init(solname, Nstep, dt, t, u0, u)
  use precision

  character(50), intent(out) :: solname
  integer, intent(out) :: Nstep
  real(dp), intent(out) :: dt, t
  real(dp), intent(out) :: u0(:), u(:)
  
  integer :: iargc
  character(50) :: arg

  
  if (iargc() < 5) then
     stop "solver, Nstep, dt, x0, v0"
  end if

  call getarg(1,arg)
  read(arg,*) solname

  call getarg(2,arg)
  read(arg,*) Nstep
    
  call getarg(3,arg)
  read(arg,*) dt

  call getarg(4,arg)
  read(arg,*) u0(1)

  call getarg(5,arg)
  read(arg,*) u0(2)
  
  t = 0.0_dp

end subroutine init
end module initialization
