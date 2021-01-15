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
    integer :: ii, nrhs, nn, err
    integer, allocatable :: ipv(:)
    
    if (size(tt) /= size(yy,2)) then
       stop "ERROR: set_points size mismatch"
    end if

    nrhs = size(yy,1)
    nn = size(tt)
    

    allocate(c(nn,nrhs)) !matrice 5x2
    allocate(A(nn,nn))   !matrice 5x5
    allocate(ipv(nn))
    
    do ii = 1, nn
       A(:,ii) = tt(:)**(ii-1)
    end do
        
    c = transpose(yy)

    ! visualizzando i valori di c a schermo, le cose sembrano funzionare fin qui
    !print*, c
    
    ! Solve A c = y using LAPACK
    call dgesv(nn,nrhs,A,nn,ipv,c,nn,err)
    
    ! controllando c adesso si vede come alcuni valori crescano fino ad esplodere
    ! dopo tre iterazioni del ciclo do in lp.f90 
    !print*, c

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
    
    do ii = size(c,1), 1, -1
      f(:) = (f(:) + c(ii,:))*t
    end do  
       
  end function poly
  
  ! Perform polynomial interpolation of derivative
  ! f = c2 + 2*c3*t + 3*c4*t*t + ...
  function poly1(t) result(f)
    real(dp), intent(in) :: t
    real(dp), allocatable :: f(:)

    integer :: ii
    
    allocate(f(size(c,2)))
    f(:) = 0.0_dp
    
    do ii = size(c,1), 2, -1
      f(:) = (f(:) + (ii-1)*c(ii,:))*t
    end do  
       
  end function poly1



end module interpolators
