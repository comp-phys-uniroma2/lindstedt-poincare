module functions
  use precision
  use interpolators, only : poly, poly1
  implicit none
  private

  public :: func
  public :: harmonic
  public :: duffing
  public :: linear_duffing
  public :: sys1
  public :: sys2
  public :: sol0
  public :: sol1

  real(dp), public :: eps
  real(dp), public :: qq 
  real(dp), public :: w0
  real(dp), public :: dw

  interface 
    function func(t,u) result(up)   
      use precision    
      real(dp), intent(in) :: t    
      real(dp), intent(in) :: u(:)    
      real(dp), allocatable :: up(:)    
    end function
  end interface
  

  contains

  function harmonic(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))

    up(1) = u(2)/w0 
    up(2) = -u(1)/w0 

  end function harmonic

      
  function duffing(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))

    up(1) = u(2)/w0 
    up(2) = (-1.0_dp - eps*u(1)**2) * u(1)/w0 

  end function duffing

  
  function linear_duffing(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)/w0
    up(2) = (-1.0_dp - 3.0_dp*eps*qq*qq) * u(1)/w0

  end function linear_duffing

  !  dy1/dt = A/w0 y1 + r(x0(t))/w0   
  ! r(t)/w0 = f(x0)/w0 - d/dt x0 
  function sys1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
   
    real(dp), allocatable :: u0(:) 
    allocate(up(size(u)))
    allocate(u0(size(u)))
    ! polynomial interpolation of u0 
    u0 = poly(t)
    ! Set qq also for linear_duffing
    qq = u0(1)
    !print*,'t=',t,'r=',duffing(t,u0) - poly1(t)
    up(:) = linear_duffing(t,u) + duffing(t, u0) - poly1(t) 

  end function sys1

  !  d/dt y2 = A/w0 y2 - 1/w0 d/dt x0
  !
  !  In this way    y = y1 + dw * y2 is the solution to
  !  dy/dt = A/w0 y + r/w0 - dw/w0 d/dt x0  
  function sys2(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    real(dp), allocatable :: u0(:) 
    allocate(up(size(u)))
    allocate(u0(size(u)))
    ! interpolation of u0 
    u0 = poly(t)
    ! Set qq also for linear_duffing
    qq = u0(1)

    up = linear_duffing(t,u) - poly1(t)/w0
    !up = linear_duffing(t,u) - duffing(t,u0)/(w0*w0)
    ! Add residual and dw term with interpolation of derivative 
    !up(:) = up(:) + duffing(t, u0) - (1.0_dp + dw/w0)*poly1(t)
    ! Alternative derivation:
    !up(:) = up(:) +  (1.0_dp - dw/w0) * duffing(t, u0) - poly1(t)
    
  end function sys2

  
  function sol0(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    allocate(u0(2))
    
    u0(1) =  cos(t) 
    u0(2) = -sin(t)

  end function sol0 

  function sol1(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    allocate(u0(2))
    
    u0(1) = -sin(t) 
    u0(2) = -cos(t)

  end function sol1 

end module functions
