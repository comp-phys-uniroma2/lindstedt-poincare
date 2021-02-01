module functions
  use precision
  use interpolators, only : poly, poly1
  implicit none
  private

  public :: func
  public :: duffing
  public :: linear_duffing
  public :: variant
  public :: linear_variant
  public :: sys0
  public :: sys1
  public :: sys2
  public :: sol0
  public :: sol1

  real(dp), public :: eps
  real(dp), public :: qq, u1, u2 
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
    up(2) = (-1.0_dp - 3.0_dp*eps*u1*u1) * u(1)/w0

  end function linear_duffing

  ! dx/dt = y - y^2 - x * (x^2 - y^2 + 2/3 y^3 + c )
  ! dy/dt = x + (y - y^2)*(x^2 - y^2 + 2/3 y^3 + c )
  ! -1/12+c=0
  function variant(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    real(dp) :: p, qp  
    allocate(up(size(u)))
    
    p = u(1)*u(1) - u(2)*u(2) + 2.0_dp/3.0_dp*u(2)**3 + 1.0_dp/12.0_dp !0.07_dp
    qp = u(2) - u(2)*u(2)

    up(1) = (qp + u(1)*p)/w0
    up(2) = (u(1) + qp * p)/w0

  end function variant

  ! p(x,y) = (x^2 - y^2 + 2/3 y^3 + c )
  ! dp/dx = 2x;  dp/dy = -2y + 2y^2
  !
  ! dfx/dx = - p(x,y) + 2x^2;   dfx/dy = 1 - 2y + 2*x*(y - y^2)
  ! dfy/dx = 1 + 2*x*(y - y^2); dfy/dy = (1 - 2y)*p - 2*(y - y^2)^2
  function linear_variant(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
   
    real(dp) :: p, qp  
    allocate(up(size(u)))
   
    p = u1*u1 - u2*u2 + 2.0_dp/3.0_dp*u2**3 + 1.0_dp/12.0_dp !0.07_dp
    qp = u2 - u2*u2

    up(1) = (2.0_dp*u1*u1 - p) * u(1)/w0 + (1.0_dp - 2.0_dp*u2 + 2*u1*qp) * u(2)/w0
    up(1) = (1.0_dp + 2.0_dp*u1*qp) * u(1)/w0 + ((1.0_dp - 2.0_dp*u2)*p - 2.0_dp*qp*qp) * u(2)/w0

    !up(1) = (-3.0_dp*u1**2+(1.0_dp-0.66666666_dp*u2)*u2**2-0.07_dp)*u(1)/w0 &
    !        + (1.0_dp+2.0_dp*(u1-u1*u2-1.0_dp)*u2)*u(2)/w0
    !up(2) = (1.0_dp+2.0_dp*(u2-u2**2)*u1)*u(1)/w0 &
    !        + (u1**2+(-3.33333_dp*u2**3+6.6666666_dp*u2**2 &
    !        -3.0_dp*u2-2.0_dp*u1**2-0.14_dp)*u2+0.07_dp)*u(2)/w0
    
  end function linear_variant

  ! dM/dt = A/w0 M
  function sys0(t,u) result(up)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    real(dp), allocatable :: u0(:)
    allocate(up(size(u)))
    allocate(u0(size(u)))

    u0 = poly(t)
    u1 = u0(1) 
    u2 = u0(2)

    up(:) = linear_variant(t,u)

  end function


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
    u1 = u0(1)
    u2 = u0(2)
    
    up(:) = linear_variant(t,u) + variant(t,u0) - poly1(t) 
    !up(:) = linear_duffing(t,u) + duffing(t,u0) - poly1(t) 

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
    u1 = u0(1)
    u2 = u0(2)

    up = linear_variant(t,u) + poly1(t)/w0
    !up = linear_duffing(t,u) + poly1(t)/w0
    
  end function sys2

  
  function sol0(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    allocate(u0(2))
    
    u0(1) = 0.5_dp*cos(t)
    u0(2) = 1.0_dp+0.5_dp*sin(t)

    !u0(1) =  cos(t) 
    !u0(2) = -sin(t)

  end function sol0 

  function sol1(t) result(u0)
    real(dp), intent(in) :: t    
    real(dp), allocatable :: u0(:)
    
    allocate(u0(2))

    u0(1) = -0.5_dp*sin(t)
    u0(2) =  0.5_dp*cos(t)
    
    !u0(1) = -sin(t) 
    !u0(2) = -cos(t)

  end function sol1 

end module functions
