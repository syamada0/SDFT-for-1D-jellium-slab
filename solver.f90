module solver
  implicit none
  real(8)   ,parameter :: pi =  3.1415926535897932d0
  complex(8),parameter :: zI = (0.d0,1.d0)
  real(8)   ,parameter :: lap(0:4) = (/ (-205.d0/72.d0) , (8.d0/5.d0), &
                          (-1.d0/5.d0) , (8.d0/315.d0) , (-1.d0/560.d0) /)

contains

  subroutine projection_occupied(Nz,dz,N_so,N_ch,V &
            ,delta_E,E_ave,E_max,E_min,chi,integral0,xi,e_fermi)
    implicit none
    integer   ,intent(in)  :: Nz,N_so,N_ch
    real(8)   ,intent(in)  :: dz,V(-Nz:Nz),delta_E,E_ave,E_max,E_min,integral0
    complex(8),intent(in)  :: chi(-Nz:Nz,N_so)
    complex(8),intent(out) :: xi(-Nz:Nz,N_so)
    real(8)   ,intent(out) :: e_fermi
    !
    integer :: i,iter
    real(8) :: e0,el,eu,diff,integral
    real(8),allocatable :: C(:)
    complex(8),allocatable :: eta(:,:,:)
    !
    integer,parameter :: N_iter = 1000
    real(8),parameter :: eps = 1.d-8

    allocate(C(0:N_ch-1),eta(-Nz:Nz,0:N_ch-1,N_so))

!   Chebyshev polynomial of the Hamiltonian
    do i=1,N_so
      call chebyshev_polynomial(Nz,dz,N_ch,V,delta_E,E_ave,chi(:,i),eta(:,:,i))
    end do

!   bisection method
    eu = E_max
    el = E_min
    do iter=1,N_iter
      e_fermi = (eu + el)*0.5d0
      call density_integral(e_fermi,Nz,dz,N_so,N_ch,delta_E,E_ave,eta,C,xi,integral)
      diff = integral - integral0
      if(diff > eps) then
        eu = e_fermi
      else if(diff < -eps) then
        el = e_fermi
      else
        write(*,*) "Fermi energy converged"
        exit
      end if
      if(iter==N_iter) write(*,*) "Fermi energy : not converged",diff
    end do
    write(*,'(a,2f15.10)') "Fermi energy =",e_fermi

    deallocate(C,eta)
    return
  end subroutine projection_occupied

  subroutine chebyshev_polynomial(Nz,dz,N_ch,V,delta_E,E_ave,chi,eta)
    implicit none
    integer   ,intent(in)  :: Nz,N_ch
    real(8)   ,intent(in)  :: dz,V(-Nz,Nz),delta_E,E_ave
    complex(8),intent(in)  :: chi(-Nz:Nz)
    complex(8),intent(out) :: eta(-Nz:Nz,0:N_ch-1)
    !
    integer :: n
    complex(8),allocatable :: x0(:),x1(:),x_wrk(:)
    !
    allocate(x0(-Nz:Nz),x1(-Nz:Nz),x_wrk(-Nz:Nz))
    x0 = chi
    call hpsi(Nz,dz,V,delta_E,E_ave,x0,x1) ! x1 = H * x0
    eta(:,0) = x0
    eta(:,1) = x1
    do n=2,N_ch-1
      call hpsi(Nz,dz,V,delta_E,E_ave,x1,x_wrk) ! x_wrk = H * x1
      x_wrk = 2d0 * x_wrk - x0
      eta(:,n) = x_wrk
      x0 = x1
      x1 = x_wrk
    end do
    deallocate(x1,x0,x_wrk)
    return
  end subroutine chebyshev_polynomial

  subroutine density_integral(mu,Nz,dz,N_so,N_ch,delta_E,E_ave,eta,C,xi,integral)
    implicit none
    integer    ,intent(in)  :: Nz,N_so,N_ch
    real(8)    ,intent(in)  :: mu,dz,delta_E,E_ave
    complex(8) ,intent(in)  :: eta(-Nz:Nz,0:N_ch-1,N_so)
    real(8)                 :: C(0:N_ch-1)
    complex(8) ,intent(out) :: xi(-Nz:Nz,N_so)
    real(8)    ,intent(out) :: integral
    !
    integer :: i,n
    real(8) :: wrk
    !
    call chebyshev_coefficient(N_ch,mu,delta_E,E_ave,C)
    wrk = 0d0
    do i=1,N_so
      xi(:,i) = 0d0
      do n=0,N_ch-1
        xi(:,i) = xi(:,i) + C(n) * eta(:,n,i)
      end do
      wrk = wrk + DOT_PRODUCT(xi(:,i),xi(:,i)) * dz
    end do
    integral = wrk / dble(N_so)
    return
  end subroutine density_integral

! cf. Ref [3]
  subroutine chebyshev_coefficient(N_ch,mu,delta_E,E_ave,C)
    implicit none
    integer   ,intent(in)  :: N_ch
    real(8)   ,intent(in)  :: mu,delta_E,E_ave
    real(8)   ,intent(out) :: C(0:N_ch-1)
    !
    integer :: j,n
    real(8) :: x1,x2,wrk
    !
    do j=0,N_ch-1
      wrk = 0d0
      do n=1,N_ch
        x1 = cos( pi * ( dble(n) - 0.5d0 ) / dble(N_ch) )
        x2 = cos( pi * dble(j) * ( dble(n) - 0.5d0 ) / dble(N_ch) )
        wrk = wrk + theta_sqrt(x1,mu,delta_E,E_ave) * x2
      end do
      if(j==0) then
        C(j) = wrk / dble(N_ch)
      else
        C(j) = wrk * 2d0 / dble(N_ch)
      end if
    end do
    return
  end subroutine chebyshev_coefficient

  function theta_sqrt(x,mu,delta_E,E_ave) ! x = -1 ~ 1
    implicit none
    real(8) :: x,mu,delta_E,E_ave
    real(8) :: theta_sqrt
    !
    real(8) :: theta,wrk
    wrk = mu - x * delta_E - E_ave
    theta = 0d0
    if(wrk > 0d0) then
      theta = wrk / pi ! cf. Ref [4]
    end if
    theta_sqrt = sqrt( theta )
  end function theta_sqrt

  subroutine electron_density_stochastic(Nz,N_so,xi,rho)
    implicit none
    integer   ,intent(in)  :: Nz,N_so
    complex(8),intent(in)  :: xi(-Nz:Nz,N_so)
    real(8)   ,intent(out) :: rho(-Nz:Nz)
    !
    integer :: i
    rho = 0d0
    do i=1,N_so
      rho(:) = rho(:) + abs(xi(:,i))**2 
    end do
    rho = rho/dble(N_so)
    return
  end subroutine electron_density_stochastic

! conjugate gradient method : (d^2/dz^2) V_H(z) = - 4*pi* rho(z)
  subroutine Hartree_potential(Nz,dz,rho,V_H)
    implicit none
    integer,intent(in)    :: Nz
    real(8),intent(in)    :: dz
    real(8),intent(in)    :: rho(-Nz:Nz)
    real(8),intent(inout) :: V_H(-Nz:Nz)
    !
    integer,parameter :: N_iter = 10000000
    real(8),parameter :: eps = 1.d-5
    !
    integer :: iz,iter
    real(8) :: rr,pDp,bb,alpha,beta,err
    real(8),allocatable :: f(:),r(:),p(:),Dp(:)
    allocate(f(-Nz:Nz),r(-Nz:Nz),p(-Nz:Nz),Dp(-Nz:Nz))

    f = V_H ! initial values

    call second_derivative(Nz,dz,f,Dp)
    rr = 0.d0
    bb = 0.d0
    do iz=-Nz,Nz
      r(iz) = (- 4.d0*pi * rho(iz)) - Dp(iz)
      rr = rr + r(iz) **2
      bb = bb + (- 4.d0*pi * rho(iz)) **2
    end do
    err = sqrt(rr/bb)

    p = 0.d0
    beta = 0.d0
    do iter=1,N_iter
      p = r + beta*p
      call second_derivative(Nz,dz,p,Dp)
      pDp = 0.d0
      do iz=-Nz,Nz
        pDp = pDp + p(iz) * Dp(iz)
      end do
      alpha = rr/pDp
      f = f + alpha * p
      r = r - alpha * Dp
      rr = 0.d0
      do iz=-Nz,Nz
        rr = rr + r(iz)**2
      end do
      err = sqrt(rr/bb)
      if(err < eps) exit
      beta = rr/ (alpha*pDp)
      if(iter==N_iter) write(*,*) "Poisson eq : not converged"
    end do
    V_H = f
    deallocate(f,r,p,Dp)
    return
  end subroutine Hartree_potential

! cf. Exc_Cor_PZ() in ARTED/common/Exc_Cor.f90
  subroutine exchange_correlation_potential(Nz,rho,V_xc)
    implicit none
    integer   ,intent(in)  :: Nz
    real(8)   ,intent(in)  :: rho(-Nz:Nz)
    real(8)   ,intent(out) :: V_xc(-Nz:Nz)
    !
    integer :: iz
    real(8) :: trho,e_xc,de_xc_drho
    !
    do iz=-Nz,Nz
      trho = rho(iz)
      call PZxc(trho,e_xc,de_xc_drho)
      V_xc(iz) = e_xc+trho*de_xc_drho
    enddo
    return
  end subroutine exchange_correlation_potential

! copied from ARTED/common/Exc_Cor.f90
  Subroutine PZxc(trho,exc,dexc_drho)
    implicit none
    real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
    real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
    real(8),parameter :: CU=0.002d0,DU=-0.0116d0
    real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
    real(8) :: ttrho,rs,rssq,rsln
    real(8),intent(in) :: trho
    real(8),intent(out) :: exc
    real(8),intent(out) :: dexc_drho


    ttrho=trho+1d-10
    rs=(3d0/(4*Pi*ttrho))**(1d0/3d0)
    exc=-const/rs
    dexc_drho=exc/(3*ttrho)
    if (rs>1d0) then
      rssq=sqrt(rs)
      exc=exc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
      dexc_drho=dexc_drho+gammaU*(0.5d0*beta1U*rssq+beta2U*rs) &
                /(3*ttrho)/(1+beta1U*rssq+beta2U*rs)**2
    else
      rsln=log(rs)
      exc=exc+AU*rsln+BU+CU*rs*rsln+DU*rs
      dexc_drho=dexc_drho-rs/(3*ttrho)*(AU/rs+CU*(rsln+1)+DU)
    endif
    return
  End Subroutine PZxc

  subroutine hpsi(Nz,dz,V,delta_E,E_ave,psi,H_psi) ! H_psi = H*psi
    implicit none
    integer   ,intent(in)  :: Nz
    real(8)   ,intent(in)  :: dz,V(-Nz:Nz),delta_E,E_ave
    complex(8),intent(in)  :: psi(-Nz:Nz)
    complex(8),intent(out) :: H_psi(-Nz:Nz)
    !
    integer :: iz
  ! kinetic energy part
    call second_derivative_C(Nz,dz,psi*(-0.5d0/delta_E),H_psi)
  ! potential energy part
    do iz=-Nz,Nz
      H_psi(iz) = H_psi(iz) + ( V(iz) - E_ave ) * psi(iz) / delta_E
    end do
    return
  end subroutine hpsi

  subroutine second_derivative(Nz,dz,f,D2f)
    implicit none
    integer,intent(in)  :: Nz
    real(8),intent(in)  :: dz
    real(8),intent(in)  :: f(-Nz:Nz)
    real(8),intent(out) :: D2f(-Nz:Nz)
    !
    integer :: iz,sgn
    !
    D2f = 0d0
    do iz=-Nz+4,Nz-4
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz+1)+f(iz-1)) &
                             + lap(2)*(f(iz+2)+f(iz-2)) &
                             + lap(3)*(f(iz+3)+f(iz-3)) &
                             + lap(4)*(f(iz+4)+f(iz-4))
    end do
    do sgn=-1,1,2 ! sgn=+,-
      iz = Nz * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-1) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-2) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)+f(iz+2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-3) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)+f(iz+2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)+f(iz+3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
    end do
    D2f = D2f * (1d0/dz**2)
    return
  end subroutine second_derivative

  subroutine second_derivative_C(Nz,dz,f,D2f)
    implicit none
    integer   ,intent(in)  :: Nz
    real(8)   ,intent(in)  :: dz
    complex(8),intent(in)  :: f(-Nz:Nz)
    complex(8),intent(out) :: D2f(-Nz:Nz)
    !
    integer :: iz,sgn
    !
    D2f = 0d0
    do iz=-Nz+4,Nz-4
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz+1)+f(iz-1)) &
                             + lap(2)*(f(iz+2)+f(iz-2)) &
                             + lap(3)*(f(iz+3)+f(iz-3)) &
                             + lap(4)*(f(iz+4)+f(iz-4))
    end do
    do sgn=-1,1,2 ! sgn=+,-
      iz = Nz * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-1) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-2) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)+f(iz+2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
      iz = (Nz-3) * sgn
      D2f(iz) = lap(0)*f(iz) + lap(1)*(f(iz-1*sgn)+f(iz+1*sgn)) &
                             + lap(2)*(f(iz-2*sgn)+f(iz+2*sgn)) &
                             + lap(3)*(f(iz-3*sgn)+f(iz+3*sgn)) &
                             + lap(4)*(f(iz-4*sgn))
    end do
    D2f = D2f * (1d0/dz**2)
    return
  end subroutine second_derivative_C

  function jellium_surface_potential(r,rho,V0)
    implicit none
    real(8) :: r,rho,V0,jellium_surface_potential
    !
    real(8) :: kf,rs,A,B,zim,V
    !
    kf = ( 3d0 * pi**2 * rho )**(1d0/3d0)
    rs = ( 3d0/( 4d0 * Pi * rho ) )**(1d0/3d0)
    A = 4d0 * v0 / kf - 1d0
    B = V0 / ( 4d0 * V0 / kf - 1d0 )
    zim = - 0.2d0 * rs + 1.25d0
    if(r<zim) then
      V = - V0 / ( A * exp( B * (r-zim) ) + 1d0 )
    else
      V = - ( 1d0 - exp( -kf * (r-zim) ) ) / ( 4d0 * (r-zim) )
    end if
    jellium_surface_potential = V
  end function jellium_surface_potential

end module solver
