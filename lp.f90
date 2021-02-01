program lp
  use precision
  use functions
  use interpolators
  use solvers
  implicit none
  real(dp), parameter :: Pi = 4.0_dp*atan(1.0_dp)
  integer, parameter :: neqs = 2
  real(dp), allocatable :: x0(:,:), y(:,:), z2(:,:), tt(:)
  real(dp), dimension(neqs) :: M1, M2, M10, M20, y10, y20, y30, y1, y2, y3, y00, fx0, u0, u
  real(dp), dimension(neqs+1) :: B
  real(dp), allocatable :: M(:,:,:)
  real(dp), dimension(neqs+1, neqs+1) :: A
  real(dp) :: t, dt, error
  integer :: N, ii, iter1, iter2, max_iter1, max_iter2, info, funit
  integer, dimension(neqs+1) :: ipv
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
  w0 = 1.00_dp
  eps = 0.1_dp
  tol_dw = 1e-4
  max_iter1 = 6 
  max_iter2 = 1 
  is_on_orbit = .true. !.false.

  open(newunit=funit, file="solution.dat")

  ! creiamo un array con punti ridondanti per tenere
  ! facilmente conto della periodicita' 
  allocate(tt(-4:N+4))
  allocate(x0(neqs,-4:N+4))
  allocate(z2(neqs,0:N))
  allocate(y(neqs,0:N))
  allocate(M(neqs,neqs,0:N))
  M = 0.0_dp

  ! Assumiamo x0 sia nota e periodica. sol0 in functions.f90
  ! t'=w0 t
  do ii = -4, N+4
     t = ii*dt
     tt(ii) = t
     x0(:,ii) = sol0(t)
  end do

  ! Primo guess di y(0).
  y00 = 0.0_dp  
  
  ! -----------------------------------------
  ! PLOT ORBITA 
  ! ----------------------------------------- 
  u0 = sol0(0.0_dp)
  write(funit,*) u0
  do ii = 1, 9*N-5 
     t = ii*dt
     call dopri54(variant,t,dt,u0,u,error)
     !call dopri54(duffing,t,dt,u0,u,error)
     u0 = u
     write(funit,*) u0
  end do
  do ii = -4, N+4
     t = ii*dt
     tt(ii) = t
     x0(:,ii) = u0
     call dopri54(variant,t,dt,u0,u,error)
     u0 = u 
  end do


  close(funit)
  ! ----------------------------------------
  
  
  ! LOOP ESTERNO
  do iter1 = 1, max_iter1
     max_error = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
        
        y1 = variant(t,x0(:,ii)) - poly1(t)
        !y1 = duffing(t,x0(:,ii)) - poly1(t)
        error = max(abs(y1(1)), abs(y1(2)))
        if (error>max_error) then
           max_error = error
        end if
        call clean_points()
     end do
     write(*,*) 'iter1:',iter1, 'w=',w0, 'error=',error
 
     !if (iter1==2) stop
     ! --------------------------------------------------------
     ! Risolvere  w0 dM/dt = A M    M(0) = I
     !
     ! A = A(x0) = linear_system @ x0(t)
     !
     if (.not.is_on_orbit) then 
        M10(1) = 1.0_dp; M20(1) = 0.0_dp
        M10(2) = 0.0_dp; M20(2) = 1.0_dp
        do ii = 0, N-1
           t = ii * dt
           M(:,1,ii) = M10(:)
           M(:,2,ii) = M20(:)
           call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
           call dopri54(sys0, t, dt, M10, M1, error)
           call dopri54(sys0, t, dt, M20, M2, error)
           M10 = M1
           M20 = M2
           call clean_points()
        end do
        M(:,1,N) = M10(:)
        M(:,2,N) = M20(:) 
     end if

     ! --------------------------------------------------------
     ! Risolvere:  w0 dy2/dt = A y2 - d/dt x0   
     !                 y2(0) = y20 = 0
     !
     ! Soluzione y = y1 + dw * y2  risolve: 
     !           w0 dy/dt = A y + r - dw d/dt x0
     !   
     y20 = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        !if (mod(ii,100)==0) print*, '      t',t,' y2',y20
        z2(:,ii) = y20(:)
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3)) 
        call dopri54(sys2, t, dt, y20, y2, error)
        y20 = y2
        call clean_points()
     end do
     !print*, '      t',N*dt,' y2',y2
     z2(:,N) = y20(:)

     ! LOOP Interno per trovare y(0) e dw
     do iter2 = 1, max_iter2

       ! --------------------------------------------------------
       ! Risolvere:  w0 dy1/dt = A y1 + r(x0(t))    
       !                 y1(0) = y00 
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
          y10 = y1
          call clean_points()
       end do
       y(:,N) = y10(:)
       !print*, '      t',N*dt,' y1',y1

       if (is_on_orbit) then
          dw = dot_product(y2,y1)/dot_product(y2,y2)
       else
          fx0 = variant(0.0_dp, x0(:,0))
          A(1:neqs, 1:neqs) = M(:,:,0)-M(:,:,N)
          A(1:neqs,3) = y2
          A(3,1:neqs) = fx0
          A(3,3) = 0.0_dp
          B(1:neqs)   = y1
          B(3)   = 0.0_dp
          
          call dgesv(3,1,A,3,ipv,B,3,info)
          if (info /= 0) then
             print*, info
             stop "Error in dgesv"
          end if

          y00 = B(1:neqs)
          dw = B(3)
          print*,'iter2=',iter2
          print*,'y(0):',y00
          print*,'dw:',dw
       end if
        
     end do
     

     ! copia soluzione y(:,ii) su x0(:)
     do ii = 0, N-1
        y(:,ii) = y(:,ii) - dw*z2(:,ii)
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
