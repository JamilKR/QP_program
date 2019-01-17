Module QP_tmp_ev
  !
  ! Temporal evolution of a wave packet
  !
  use Harmonic_oscillator
  use QP_module
  !
  implicit none
  !
  Character(len=1):: temp !input: if = 'g' ---> gaussian WF
  Character(len=60):: outfile_tmp !Temporal evolution outfile
  Double complex, allocatable::   &
       Psi(:)  !Wave packet
  Double precision:: &
       sigma,  & !Distribution amplitude
       mu        !Distribution center
  Double precision, allocatable:: &
       Prob(:),   & !Prob density
       Coef(:),   & !WP coefficients
       x_mtx(:,:)   !< psi_i | x | psi_j >
  Double precision, allocatable:: t_vec(:) !Temporal axis points
  Double precision:: &
       tmax,   & !Maximun point
       tstep     !abs( t_vec(1) - t_vec(2))
  Integer:: tdim !Total t-points
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
    ! do i=0,Nval
    !    do j=0,Ncon
    !       aux=0.0d0
    !       do xx=1,xdim
    !          aux = aux + dexp( -((x_vec(xx)-mu)**2.0d0)/(4.0d0* sigma**2.0d0)) &
    !               *HO(xx,j)
    !       enddo
    !       coef(i) = coef(i) + aux * H(j,i) * dabs(x_vec(1)-x_vec(2))
    !    enddo
    ! enddo
    !
    do i = 0,Nval
       !
       aux=0.0d0
       !
       do xx = 1,xdim
          !
          aux = aux + dexp( -((x_vec(xx)-mu)**2.0d0)/(4.0d0* sigma**2.0d0)) * &
               QP_mtx(xx,i)
          !
       enddo
       !
       coef(i) = aux
       !
    enddo
    !         
    coef = (coef*xstep)/dsqrt( sigma*dsqrt(2.0d0*pi) )
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
  ! Subroutine Prob_dens
!     !
!     !This routine builds the Wave Packet and the Prob Dens Function
!     !
!     implicit none
!     !
!     integer:: t,k
!     !
!     allocate(Psi(1:xdim),Prob(1:xdim))
!     !
!     open(unit=21,file=trim(outfile),status='new')
!     !
!     do 10 t = 1,tdim
!        !
!        Psi  = 0.0d0
!        Prob = (0.0d0,0.0d0)
!        !
!        do 20 k = 0,Nval
!           !
!           psi(:) = psi(:) + coef(k) * zexp( (0.0d0,-E(k)*t_vec(t)) ) *
       
! 10  enddo
    
!     !
!     close(21)
!     !
!   end Subroutine Prob_dens
  !
end Module QP_tmp_ev
