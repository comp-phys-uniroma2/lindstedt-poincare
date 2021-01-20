module interpolators
  use precision
  implicit none
  private
  public :: set_points
  public :: clean_points 
  public :: poly
  public :: poly1 

  real(dp), allocatable :: c(:,:)

  contains

  ! Performs cubic interpolation
  !        t
  ! |---|--x|---|
  ! t0  t1  t2  t3
  ! 
  subroutine set_points(tt, yy)
    real(dp), intent(in) :: tt(:)
    real(dp), intent(in) :: yy(:,:)
    
    real(dp), allocatable :: A(:,:)
    integer :: ii, jj, nrhs, nn, err
    integer, allocatable :: ipv(:)
    
    if (size(tt) /= size(yy,2)) then
       stop "ERROR: set_points size mismatch"
    end if

    nrhs = size(yy,1)
    nn = size(tt)
    
    if (allocated(c)) then
      deallocate(c)
    end if   
    allocate(c(nn,nrhs))
    allocate(A(nn,nn))   !matrice 5x5
    allocate(ipv(nn))
   
    ! Setup linear system:  A c = y
    ! c1 * t1^0 + c2 * t1^1 + c3 * t1^2 + ... = y1
    ! c1 * t2^0 + c2 * t2^1 + c3 * t2^2 + ... = y2
    ! c1 * t3^0 + c2 * t3^1 + c3 * t3^2 + ... = y3
    ! ...

    do ii = 1, nn
      A(:,ii) = 1.0_dp
      do jj = 2, ii
        A(:,ii) = A(:,ii)*tt(:)  !tt(:)**(ii-1)
      end do  
    end do
        
    c = transpose(yy)

    ! Solve A c = y using LAPACK
    call dgesv(nn,nrhs,A,nn,ipv,c,nn,err)

    ! err=5
    ! dalla documentazione lapack (dove err è chiamato INFO): "if INFO = i, U(i,i) is exactly zero.
    ! The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.
    if (err /= 0) then
       print*, err
       stop "ERROR in dgesv"
    end if

    deallocate(A)
    deallocate(ipv)
  end subroutine set_points


  subroutine clean_points()
    if (allocated(c)) then
      deallocate(c)    
    end if
  end subroutine clean_points


  ! Perform actual polynomial interpolation
  ! f = c1 + c2*t + c3*t*t + c4*t*t*t + ...
  function poly(t) result(f)
    real(dp), intent(in) :: t
    real(dp), allocatable :: f(:)

    integer :: ii
    
    allocate(f(size(c,2)))
    f(:) = 0.0_dp
    
    do ii = size(c,1), 2, -1
      f(:) = (f(:) + c(ii,:))*t
    end do  
      
    f(:) = f(:) + c(1,:)

  end function poly
  
  ! Perform polynomial interpolation of derivative
  ! f = c2 + 2*c3*t + 3*c4*t*t + ...
  function poly1(t) result(f)
    real(dp), intent(in) :: t
    real(dp), allocatable :: f(:)

    integer :: ii
    
    allocate(f(size(c,2)))
    f(:) = 0.0_dp
    
    do ii = size(c,1), 3, -1
      f(:) = (f(:) + (ii-1)*c(ii,:))*t
    end do 

    f(:) = f(:) + c(2,:) 
       
  end function poly1



end module interpolators
