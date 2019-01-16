Module QP_tmp_ev
  !
  ! Temporal evolution of a wave packet
  !
  use QP_module
  use Harmonic_oscillator
  !
  implicit none
  !
  Character(len=1):: temp !input: if = 'g' ---> gaussian WF
  Double complex, allocatable::   &
       Psi(:)  !Wave packet
  Double precision:: &
       sigma,  & !Distribution amplitude
       mu,     & !Distribution center
       x0,     & !Lowest x value
       xstep     ! abs( x_vec(1) - x_vec(2))
  integer:: xdim !Total x-points
  Double precision, allocatable:: &
       Prob(:),   & !Prob density
       HO(:,:),   & !HO Wave Funcions
       x_vec(:),  & !X-axis
       Coef(:),   & !WP coefficients
       x_mtx(:,:)   !< psi_i | x | psi_j >
  !
contains
  !
  Subroutine Coefficients
    !
    ! This routine computes WF coefficients using a Gaussian Distribution
    !
    implicit none
    !
    integer:: i,j,xx
    double precision:: norma, aux
    !
    allocate(coef(0:Nval))
    coef=0.0d0
    !
    do i=0,Nval
       do j=0,Ncon
          aux=0.0d0
          do xx=1,xdim
             aux = aux + dexp( -((x_vec(xx)-mu)**2.0d0)/(4.0d0* sigma**2.0d0)) &
                  *HO(xx,j)
          enddo
          coef(i) = coef(i) + aux * H(j,i) * dabs(x_vec(1)-x_vec(2))
       enddo
    enddo
    !
    coef = coef/dsqrt( sigma*dsqrt(2.0d0*pi) )
    !
    !Norma test
    !
    norma=0.0d0
    do i=0,Nval
       norma = norma + coef(i)*coef(i)
    enddo
    write(*,*)
    write(*,*) norma, 'WavePacket Norm'
    write(*,*)
    !
    if (norma .lt. 0.9999d0 .or. norma .gt. 1.0001d0)  then
       !
       write(*,*)
       Write(*,'(A,F7.4,A)') "Wave Packets Error: 0.9999 < Norm < 1.0001 ! Norm=", norma
       Write(*,'(A)') "Try:"
       Write(*,'(A)') " 1) Choose a higher Nval."
       Write(*,'(A)') " 2) Increase X-axis range."
       write(*,*)
       !
       STOP
       !
    endif
    !
  end Subroutine Coefficients
  !  
  !*****************************************************************************
  !
  Subroutine Trans_x_matrix
    !
    !This routine compute x_mtx(i,j)= < psi_i | x | psi_j > with Hamiltonian's
    !eigenstates.
    !
    implicit none
    !
    integer:: n,k,l
    double precision:: aux
    !
    allocate(x_mtx(0:Nval,0:Nval))
    !
    x_mtx=0.0d0
    !
    do l=0,Nval !column fixed
       do k=0,Nval !row fixed
          !
          aux = H(1,k) * H(0,l) 
          !
          do n = 1,Ncon-1 ! Last addend excluded
             !
             aux = aux + H(n,l) * ( &
                  dble(n+1)*H(n+1,k) + dble(n)*H(n-1,k) )
             !
          enddo
          !
          aux = aux + dble(Ncon) * H(Ncon-1,k) * H(Ncon,l)
          !
          x_mtx(k,l) = aux / dsqrt(2.0d0)         
          !
       enddo
    enddo
    !
  end Subroutine Trans_x_matrix
  !  
  !*****************************************************************************
  !
  Subroutine Prob_dens
    !
    !This routine builds the Wave Packet and the Prob Dens Function
    !
    implicit none
    !
  end Subroutine Prob_dens
  !
end Module QP_tmp_ev
