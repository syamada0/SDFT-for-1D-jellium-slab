!
! 1D Stochastic Density Functional Theory
!
! ### Convention ###
!
! # of grid points = 2*Nz+1
! grid points iz = -Nz:Nz
! grid spacing dz = width_z/(2*Nz)
!
! 2nd order differential operator d^2/dz^2 = (2nd order difference operator) / (dz)^2
!
! jellium slab (periodic in the x-y plane)
! # of electrons = (thickness of slab) * (area in the x-y plane) * (jellium density)
!
! The jellium surfaces must be located at the middle point between grid points.
! Example: 
!   Nz = 1000, width_z = 1000 --> dz = 0.5
!   thickness_slab must be [integer + dz] such as "80.5" or "100.5".
!
program main
  use solver
  implicit none
  integer :: Nz,N_so,N_ch,N_iter,i_mode
  real(8) :: width_z,thickness_slab,density_jellium,E_min,E_max,alpha,V0,eps_residual

  namelist/input/width_z &          ! system size
                ,thickness_slab &   ! thickness of the jellium slab
                ,density_jellium &  ! jellium density
                ,E_min,E_max &      ! eigenenergies interval of the Hamiltonian
                ,Nz &               ! # of the grid points
                ,N_so &             ! # of the stochastic orbitals
                ,N_ch &             ! length of the Chebyshev expansion
                ,i_mode &           ! initial potential mode (=0 --> model potential, =1 --> read file)
                ,N_iter &           ! # of SCF iteration
                ,alpha &            ! potential mixing parameter
                ,V0 &               ! parameter for the initial model potential
                ,eps_residual       ! convergence threshold

  integer :: i,iz,iter
  real(8) :: dz,e_fermi,delta_E,E_ave,r,residual
  real(8)   ,allocatable :: V(:),V_ext(:),rho(:),rho_jellium(:) &
                           ,V_H(:),V_xc(:),rand(:,:) 
  complex(8),allocatable :: xi(:,:),chi(:,:)

! input
  i_mode = 0
  eps_residual = 0.0001
  read(10,input)
  write(*,input)

! preparation
  dz = width_z/dble(2*Nz) ! iz=-Nz:Nz --> [# of iz] = 2*Nz+1, width_z <--> 2*Nz

  E_ave = 0.5d0 * ( E_max + E_min )
  delta_E = 0.5d0 * ( E_max - E_min )

  allocate(V(-Nz:Nz),V_ext(-Nz:Nz),V_H(-Nz:Nz),V_xc(-Nz:Nz) &
          ,rho(-Nz:Nz),rho_jellium(-Nz:Nz) &
          ,xi(-Nz:Nz,N_so),chi(-Nz:Nz,N_so))

  do iz=-Nz,Nz
    r = abs( dble(iz) * dz )
    if(r < thickness_slab/2d0) then
      rho_jellium(iz) = - density_jellium ! jellium background charge
    else
      rho_jellium(iz) = 0.d0
    end if
  end do

! stochastic orbitals (SO)
  allocate(rand(-Nz:Nz,N_so))
  call random_number(rand)
  do i=1,N_so
    do iz=-Nz,Nz
      chi(iz,i) = exp( 2.d0 * pi * zI * rand(iz,i) ) / sqrt(dz)
    end do
  end do
  deallocate(rand)

! test  
  call test(N_ch)
  call test2(Nz,N_so,dz,chi)

! initial potential
  if(i_mode==0) then
    do iz=-Nz,Nz
      r = abs(dble(iz) * dz)
      r = r - thickness_slab/2d0
      V(iz) = jellium_surface_potential(r,density_jellium,V0)
    end do
    write(*,*) "set initial V"
  else
    rewind(555)
    read(555,*) iz
    if(iz==Nz) then
      read(555,*) V(-Nz:Nz)
      write(*,*) "read V"
    end if
  end if

  write(*,*) "start Ground State Calculation"

  do iter=1,N_iter

!   potential (V) --> projected SO (xi) & Fermi energy (e_fermi)
    call projection_occupied(Nz,dz,N_so,N_ch,V,delta_E,E_ave,E_max,E_min &
        ,chi,density_jellium*thickness_slab,xi,e_fermi)

!   xi --> rho
    call electron_density_stochastic(Nz,N_so,xi,rho)

!   rho + rho_jellium --> V_H
    call Hartree_potential(Nz,dz,rho+rho_jellium,V_H)

!   rho --> V_xc
    call exchange_correlation_potential(Nz,rho,V_xc)

!   convergence check
    residual = 0d0
    rewind(333)
    do iz=-Nz,Nz
      residual = residual + abs( V(iz) - V_H(iz) - V_xc(iz) )
      r = dble(iz) * dz
      write(333,'(6f15.10)') r,rho(iz),rho_jellium(iz),V_H(iz),V_xc(iz),V(iz)
    enddo
    residual = residual * dz / width_z
    write(*,'(a,f25.17)') "residual(V)= ",residual
    write(*,*) "integral(rho+rho_jellium)= ", sum(rho(:)+rho_jellium(:)) * dz
    if(residual < eps_residual) then
      write(*,*) "GSC converged",iter,residual
      exit
    end if

!   potential mixing
    V = (1.d0 - alpha)* V + alpha* ( V_H + V_xc )
    write(*,*) "GS",iter

  end do

  rewind(555)
  write(555,*) Nz
  write(555,*) V(-Nz:Nz)

  write(*,*) "end"
  deallocate(V,V_ext,V_H,V_xc,rho,rho_jellium,xi,chi)

! test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains

  subroutine test2(Nz,N_so,dz,chi) ! for representability of the stochastic orbitals
    implicit none
    integer   ,intent(in) :: Nz,N_so
    real(8)   ,intent(in) :: dz
    complex(8),intent(in) :: chi(-Nz:Nz,N_so)
    !
    integer :: i,iz
    real(8) :: rz
    complex(8) :: wrk
    complex(8),allocatable :: f(:),g(:)

    allocate(f(-Nz:Nz),g(-Nz:Nz))

    do iz=-Nz,Nz
      rz = dble(iz) * dz
      f(iz) = dcmplx( func(rz) )
    end do

    g = 0d0
    do i=1,N_so
      wrk = dot_product(chi(-Nz:Nz,i),f(-Nz:Nz)) * dz
      g = g + chi(:,i) * wrk
    end do
    g = g / dble(N_so)

    do iz=-Nz,Nz
      rz = dble(iz) * dz
      write(222,*) rz,dble(f(iz)),dble(g(iz))
    end do

    deallocate(f,g)
    return
  end subroutine test2

  function func(x)
    implicit none
    real(8) :: x,func
    func = sin(x/30d0) + cos(x/700d0)
  end function func

  subroutine test(N_ch) ! for Chebyshev expansion
    implicit none
    integer,intent(in) :: N_ch
    !
    integer :: Nz,i
    real(8) :: dz,rz
    real(8),allocatable :: T_ch(:,:),C(:),f(:)

    Nz = 1000
    dz = 1d0/dble(Nz)

    allocate(T_ch(-Nz:Nz,0:N_ch-1),C(0:N_ch-1),f(-Nz:Nz))
    call test_chebyshev_polynomial(Nz,dz,N_ch,T_ch)
    call chebyshev_coefficient(N_ch,0d0,1d0,0d0,C)

    f = 0d0
    do i=0,N_ch-1
      f = f + C(i) * T_ch(:,i) ! Chebyshev expansion
    end do

    do i=-Nz,Nz
      rz = dble(i) * dz
      write(111,*) rz,theta_sqrt(rz,0d0,1d0,0d0),f(i)
    end do

    deallocate(T_ch,C,f)
    return
  end subroutine test

  subroutine test_chebyshev_polynomial(Nz,dz,N_ch,T_ch)
    implicit none
    integer,intent(in)  :: Nz,N_ch
    real(8),intent(in)  :: dz
    real(8),intent(out) :: T_ch(-Nz:Nz,0:N_ch-1)
    !
    integer :: n,iz
    real(8) :: rz
    complex(8),allocatable :: x0(:),x1(:),x_wrk(:)
    !
    allocate(x0(-Nz:Nz),x1(-Nz:Nz),x_wrk(-Nz:Nz))
    x0 = 1
    do iz=-Nz,Nz
      rz = dble(iz) * dz
      x1(iz) = rz
    end do
    T_ch(:,0) = x0
    T_ch(:,1) = x1
    do n=2,N_ch-1
      do iz=-Nz,Nz
        rz = dble(iz) * dz
        x_wrk(iz) = 2d0 * rz * x1(iz) - x0(iz)
      end do
      T_ch(:,n) = x_wrk
      x0 = x1
      x1 = x_wrk
    end do
    deallocate(x1,x0,x_wrk)
    return
  end subroutine test_chebyshev_polynomial
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program main

