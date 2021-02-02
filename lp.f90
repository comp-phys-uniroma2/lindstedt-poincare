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
  integer :: N, ii, iter1, max_iter1, info, funit
  integer, dimension(neqs+1) :: ipv
  real(dp) :: dw_old, tol_dw, max_error
  logical :: is_on_orbit
  character(2) :: citer 
  character(10) :: arg

  if (command_argument_count() < 1) then
    print*,'lp niter1'    
    stop
  end if

  call get_command_argument(1, arg)
  read(arg,*) max_iter1

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
  !max_iter1 = 30 
  is_on_orbit = .false.


  ! creiamo un array con punti ridondanti per tenere
  ! facilmente conto della periodicita' 
  allocate(tt(-4:N+4))
  allocate(x0(neqs,-4:N+4))
  allocate(z2(neqs,0:N))
  allocate(y(neqs,0:N))
  allocate(M(neqs,neqs,0:N))
  M = 0.0_dp

  ! Assumiamo x0 sia nota e periodica. sol0 in functions.f90
  do ii = -4, N+4
     t = ii*dt
     tt(ii) = t
     x0(:,ii) = sol0(t)
  end do
  open(newunit=funit, file="solution0.dat")
  do ii = -4, N+4
     write(funit,*) x0(:,ii) 
  end do
  close(funit)

  ! Primo guess di y(0).
  y00 = 0.0_dp  
  
  ! -----------------------------------------
  ! PLOT ORBITA 
  ! ----------------------------------------- 
  ! u0 = x0(:,ii) 
  !do ii = 1, 10*N 
  !   t = ii*dt
  !   call dopri54(variant,t,dt,u0,u,error)
  !   u0 = u
  !end do
  ! ----------------------------------------
  
  ! ALGORITMO DI Lindsted-Poincare
  ! risolve: 
  !          w0 dy/dt = A y + r - dw d/dt x0
  !              y(0) = y0
  !
  ! soluzione generale: y(t) = Y(t)y0 + f1(t) - dw f2(t)
  ! 
  ! poniamo: y = y1 - y2
  !
  ! risolviamo separatamente:
  ! 0) sistema omogeneo
  !    w0 dY/dt = A Y    Y(0) = I 
  ! 1)
  !    w0 dz2/dt = A z2 + d/dt x0   
  !        z2(0) = 0
  ! 2)  
  !    w0 dy1/dt = A y1 + r(x0(t))    
  !        y1(0) = 0
  !
  ! y1(2*pi) = f1(2*pi)
  ! y2(2*pi) = dw f2(2*pi)
  !
  ! y(2*pi) = y0 = M*y0 + f1 - dw f2
  !
  ! => Condizione 1: (I-M) y0 + f2 dw = f1 
  !    Condizione 2: f(x(0)) y0       = 0
  !
  ! M = Y(2pi); f2 = z2(2pi); f1 = y1(2pi)

  do iter1 = 1, max_iter1
     max_error = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
        
        !y1 = variant(t,x0(:,ii)) - poly1(t)
        !y1 = duffing(t,x0(:,ii)) - poly1(t)
        !error = max(abs(y1(1)), abs(y1(2)))
        y1 = poly(t)
        error = abs(y1(1)**2 - y1(2)**2 + 2.0_dp/3.0_dp*y1(2)**3 + cc)
        if (error>max_error) then
           max_error = error
        end if
        call clean_points()
     end do
     write(*,*) 'iter1:',iter1, 'w=',w0, 'error=',error
 
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
     ! Risolvere:  w0 dz2/dt = A z2 + d/dt x0   
     !                 z2(0) = 0
     !
     ! Soluzione y2 = dw * z2 
     !   
     y20 = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        z2(:,ii) = y20(:)
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3)) 
        call dopri54(sys2, t, dt, y20, y2, error)
        y20 = y2
        call clean_points()
     end do
     z2(:,N) = y20(:)

     ! --------------------------------------------------------
     ! Risolvere:  w0 dy1/dt = A y1 + r(x0(t))    
     !                 y1(0) = 0 
     !
     ! r(t) = f(x0) - w0 * d/dt x0 
     !
     y10 = 0.0_dp
     do ii = 0, N-1
        t = ii*dt
        y(:,ii) = y10(:)
        call set_points(tt(ii-2:ii+3), x0(:,ii-2:ii+3))
        call dopri54(sys1, t, dt, y10, y1, error)
        y10 = y1
        call clean_points()
     end do
     y(:,N) = y10(:)


     if (is_on_orbit) then
        dw = dot_product(y2,y1)/dot_product(y2,y2)
     else
        ! Sistema lineare per trovare y(0) e dw
        fx0 = variant(0.0_dp, x0(:,0))
        A(1:neqs, 1:neqs) = M(:,:,0) - M(:,:,N)
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
        !print*,'      y(0):',y00
     end if
        
     ! copia soluzione y(:,ii) su x0(:)
     write(citer,'(I2.2)') iter1
     open(newunit=funit, file="solution"//citer//".dat")
     do ii = 0, N
        y(:,ii) = matmul(M(:,:,ii),y00) + y(:,ii) - dw*z2(:,ii)
        x0(:,ii) = x0(:,ii) + y(:,ii)
        write(funit,*) x0(:,ii)
     end do
     close(funit)

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
