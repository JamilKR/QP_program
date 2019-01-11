program HO_test
  !
  use quadpack
  !
  use Harmonic_oscillator
  !
  implicit none
  !
  double precision,allocatable:: HO(:,:),x_vec(:)
  integer,parameter:: xdim=40000
  integer,parameter:: N=40
  integer:: i,j
  double precision:: norma
  !
  allocate(HO(xdim,0:N),x_vec(xdim))
  !
  x_vec(1)=-20.0d0
  do i=2,xdim
     X_vec(i)=X_vec(i-1)+0.01d0
  enddo
  !
  HO = HO_WF(N,x_vec)
  !
  do i=1,xdim
     write(*,*) HO(i,:)
  enddo
  !
  call qags ( F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL,
  !
end program HO_test

