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
    up(2) = (-1.d0+eps*u(1)**2)/w0*u(1) 

  end function duffing

  
  function linear_duffing(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)/w0
    up(2) = (-1.0_dp + 3.0_dp*eps*qq*qq) * u(1)/w0

  end function linear_duffing


  function sys1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
   
    real(dp), allocatable :: u0(:) 
    allocate(up(size(u)))
    allocate(u0(size(u)))
    ! cubic interpolation of u0 
    u0 = poly(t)
    ! Set qq also for linear_duffing
    qq = u0(1)

    up = linear_duffing(t,u)
    
    ! Add residual
    up(:) = up(:) + duffing(t, u) - poly(t) 

  end function sys1

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

    up = linear_duffing(t,u)
    ! Add residual and dw term with interpolation of derivative 
    up(:) = up(:) + duffing(t, u0) - (1.0_dp + dw/w0)*poly1(t)

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
