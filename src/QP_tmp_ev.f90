Module QP_tmp_ev
  !
  ! Temporal evolution of a wave packet
  !
  use QP_module                 
  !
  implicit none
  !
  Double complex, allocatable::   Psi(:)  !Wave packet
  Double precision, allocatable:: Coef(:) !WP coefficients
  Double precision:: &
       Sigma  & !Distribution amplitude
       Mu       !Distribution center
  Double precision, allocatable:: Prob(:) !Prob density
  !
contains
  !
  Subroutine Coefficients
    !
    ! This routine computes WF coefficients using a Gaussian Distribution
    !
    implicit none
    !
    integer:: i
    double precision:: norma
    !
    coef=0.0d0
    
    !
  end Subroutine Coefficients
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
