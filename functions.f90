module functions
  use precision
  implicit none
  private

  public :: harmonic
  public :: func

  interface 
     function func(t,u) result(up)   
       use precision    
       real(dp), intent(in) :: t    
       real(dp), intent(in) :: u(:)    
       real(dp), allocatable :: up(:)    
     end function func
  end interface

contains

  function harmonic(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = -u(1)
    
  end function harmonic
end module functions


