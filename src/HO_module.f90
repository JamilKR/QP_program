Module Harmonic_oscillator
  !
  !This module contains Hermite polinomials and generate HO eigenfunctions
  !
  implicit none
  !
  Double precision, parameter:: pi=3.141592653589793238d0
  !
contains
  !
  Function Hermite (N,x_vec)
    !
    !This routine return a matrix with N first Hermite polynomials
    !evaluated in x_vec points
    !
    ! INPUTs:   N ---------> Maximun order of Hermite polynomials
    !           x_vec -----> x points were Her will be evaluate
    !
    ! OUTPUT:   Hermite ---> Matrix with dimension ( dim(x_vec), 0:N )
    !
    implicit none
    !
    integer,intent(in):: N
    double precision,dimension(:),intent(in):: x_vec
    integer:: xdim  !size(x_vec)
    double precision,allocatable:: Hermite(:,:)
    integer i
    !
    xdim = size(x_vec)
    allocate( Hermite(xdim,0:N))
    !
    Hermite = 0.0d0
    !Define first and second H pol.
    Hermite(:,0) = 1.0d0 
    Hermite(:,1) = 2.0d0*x_vec(:)
    !
    do i=1,N-1
       Hermite(:,i+1) = 2.0d0*x_vec(:)*Hermite(:,i) &
            - 2.0d0*dble(i)*Hermite(:,i-1)
    enddo
    !
  end Function Hermite
  !
  !*****************************************************************************
  !
  Function HO_WF(N,x_vec)
    !
    !This routine return the HO Wavefunctions evaluated in x_vec points
    !
    ! INPUTs:   N ---------> Maximun order of Hermite polynomials
    !           x_vec -----> x points were Her will be evaluate
    !
    ! OUTPUT:   Hermite ---> Matrix with dimension ( dim(x_vec), 0:N )
    !
    implicit none
    !
    integer,intent(in):: N
    double precision,dimension(:),intent(in):: x_vec
    integer:: xdim  !size(x_vec)
    double precision,allocatable:: HO_WF(:,:)
    integer i
    !
    xdim = size(x_vec)
    allocate( HO_WF(xdim,0:N))
    !
    HO_WF = Hermite (N,x_vec)
    !
    do i=0,N
       HO_WF(:,i) = ( 1.0d0 / (dsqrt( (2.0d0)**(dble(i)) * dGAMMA(dble(i+1)) ) &
            * pi**0.25d0 )) * dexp(- 0.5d0* (x_vec(:))**2.0d0)  * &
            HO_WF(:,i)       
    enddo
    !
  end Function HO_WF
  !
end Module Harmonic_oscillator

