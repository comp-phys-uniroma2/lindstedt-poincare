program lp
  use precision
  use functions
  use interpolators
  use solvers
  implicit none
  real(dp), parameter :: Pi = 4.0_dp*atan(1.0_dp)
  real(dp), allocatable :: x0(:,:), y(:,:), yy2(:,:), tt(:)
  real(dp), dimension(2) :: y10, y20, y30, y1, y2, y3
  real(dp), allocatable :: M(:,:,:)
  real(dp) :: t, dt, error
  integer :: N, ii, iter1, iter2, max_iter1, max_iter2
  real(dp) :: dw_old, tol_dw, max_error
  logical :: is_on_orbit
  
  !
  !  0                                          2pi
  !  |---|---|---|---|---|---|---|---|---|---|---| 
  !  0   1   2   3                               N
  !
  !  dt = 2*pi/N
  N = 1000
  dt = 2.0_dp * Pi / N
  w0 = 1.0_dp
  eps = 0.1_dp
  tol_dw = 1e-4
  max_iter1 = 5
  max_iter2 = 1
  is_on_orbit = .true.

  ! creiamo un array con punti ridondanti per tenere
  ! facilmente conto della periodicita' 
  allocate(tt(-4:N+4))
  allocate(x0(2,-4:N+4))
  allocate(yy2(2,0:N))
  allocate(y(2,0:N))
  allocate(M(2,2,N))

  ! Assumiamo x0 sia nota e periodica. sol0 in functions.f90
  do ii = -4, N+4
     t = ii*dt
     tt(ii) = t
     x0(:,ii) = sol0(t)
  end do

  ! Primo guess di y(0).
  y00 = 0.001_dp  
  M = 0.0_dp

  ! LOOP ESTERNO
  ! x[i] = x[i-1] + y
  ! w = w0 + dw
  do iter1 = 1, max_iter1
     max_error = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
        !y10 = sol1(t); y20=sol1(t+dt*0.5_dp); y30= sol1(t+dt) 
        !y1 = poly1(t); y2 = poly1(t+dt*0.5_dp); y3 = poly1(t+dt) 
        !print*, abs(y1(1)-y10(1)), abs(y2(1)-y20(1)), abs(y3(1)-y30(1))
        
        y1 = duffing(t,x0(:,ii)) - poly1(t)
        error = max(y1(1), y1(2))
        !if (mod(ii,100)==0) print*,'r=',y1,'err=',sqrt(dot_product(y1, y1))
        if (mod(ii,100)==0) print*,'r=',y1,'err=',error 
        if (error>max_error) then
           max_error = error
        end if
        call clean_points()
     end do
     write(*,*) 'iter1:',iter1, 'w=',w0, 'error=',error

     ! --------------------------------------------------------
     ! Risolvere  w0 dY/dt = A Y    Y(0) = I
     !
     ! A = A(x0) = linear_system @ x0(t)
     !
     ! Nel caso del duffin non serve risolvere per Y(t)
     ! perche' in ogni caso y(0) = 0
     if (.not.is_on_orbit) then 
       !end do
       y10 = 1.0_dp
       y10 = 0.0_dp
       !do ii = 0, N
          M(:,1,ii)=y10
          !t = ii*dt
          call dopri54(sys1, t, dt, y10, y1, error)
          !call rk4(linear_duffin, t, dt, u0, u)

       M(:,1,N)=y10

       y10 = 0.0_dp
       y10 = 1.0_dp
       !do ii = 0, N
          !t = ii*dt
          M(:,2,ii)=y10
          call dopri54(sys1, t, dt, y10, y1, error)
          !call rk4(linear_duffin, t, dt, u0, u)

       M(:,2,N)=y10
     end if
     ! --------------------------------------------------------
     ! Risolvere:  w0 dy2/dt = A y2 - d/dt x0   
     !                 y2(0) = y20 = 0
     !
     ! Soluzione y = y1 + dw * y2  risolve: 
     !           w0 dy/dt = A y + r - dw d/dt x0   
     y20 = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        !if (mod(ii,100)==0) print*, '      t',t,' y2',y20
        yy2(:,ii) = y20(:)
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3)) 
        call dopri54(sys2, t, dt, y20, y2, error)
        !call rk4(sys2, t, dt, y20, y2)
        y20 = y2
        call clean_points()
     end do
     !print*, '      t',N*dt,' y2',y2
     yy2(:,N) = y20(:)
 
     do iter2 = 1, max_iter2
       ! --------------------------------------------------------
       ! Risolvere:  w0 dy1/dt = A y1 + r(x0(t))    
       !                 y1(0) = y10 
       !
       ! r(t) = f(x0) - w0 * d/dt x0 
       ! 
       y10 = y00
       do ii = 0, N-1
          t = ii*dt
          !if (mod(ii,100)==0) print*, '      t',t,' y1',y10
          y(:,ii) = y10(:)
          call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
          call dopri54(sys1, t, dt, y10, y1, error)
          !call rk4(sys1, t, dt, y10, y1)
          y10 = y1
          call clean_points()
       end do
       y(:,N) = y10(:)
       !print*, '      t',N*dt,' y1',y1
 
       ! f1 is last y1 (@ 2pi) 
       ! f2 is last y2 (@ 2pi) 
       if (is_on_orbit) then
          dw = - dot_product(y2,y1)/dot_product(y2,y2)
       else
          ! LINEAR SOLVER
       end if
     
       print*, '      iter2',iter2,'dw=',dw
 
     end do       
     
     ! copia soluzione y(:,ii) su x0(:)
     do ii = 0, N-1
        y(:,ii) = matmul(M(:,:,ii),y00) + y(:,ii) + dw*yy2(:,ii)
        !if (mod(ii,10)==0) print*,'      t',ii*dt,' y',y(:,ii)
        x0(:,ii) = x0(:,ii) + y(:,ii)
     end do

     ! estensioni periodiche
     x0(:,-1) = x0(:,-1) + y(:,N-1) 
     x0(:,-2) = x0(:,-2) + y(:,N-2)
     x0(:,-3) = x0(:,-3) + y(:,N-3)
     x0(:,-4) = x0(:,-4) + y(:,N-4)

     x0(:,N+1) = x0(:,N+1) + y(:,1)     
     x0(:,N+2) = x0(:,N+2) + y(:,2)
     x0(:,N+3) = x0(:,N+3) + y(:,3)
     x0(:,N+4) = x0(:,N+4) + y(:,4)
     
     w0 = w0 + dw
   end do


end program lp
