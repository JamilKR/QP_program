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
  Double Precision:: Emean !< Psi | H | Psi >
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
  Subroutine Prob_dens
    !
    !This routine builds the Wave Packet and the Prob Dens Function
    !
    implicit none
    !
    Integer:: t,k,i
    Double precision:: x_aux
    Double precision,allocatable::auxX(:)
    !
    allocate(Psi(1:xdim),Prob(1:xdim),auxX(1:xdim))
    !
    open(unit=21,file=trim(outfile_tmp)//".dat",status='replace')
    inquire(21)
    !
    open(unit=31,file=trim(outfile_tmp)//"_Xmean.dat",status='replace')
    inquire(31)
    !
    write(21,11) 
    write(21,12) xdim,x0,xstep
    write(21,13) tmax,tstep
    write(21,14)
    !
    write(31,21) 
    !
11  format("# Rows: Psi(x,t=cte); Columns: Psi(x=cte,t)")
12  format("# X-axis: Total points:", I5, " ; x0:",F7.2" ; x-step:",F7.5)
13  format("# T-evol: Initial t0: 0.0; T-max:",F7.2," ; t-step:",F7.5)
14  format("# The Prob density is displaced the value < Psi | H | Psi >")
    !
21  format("# 1 column: < Psi | x | Psi >")
    !
    do 10 t = 1,tdim
       !
       Psi   = (0.0d0,0.0d0)
       Prob  =  0.0d0
       x_aux =  0.0d0
       auxX  =  0.0d0
       !
       do 20 k = 0,Nval
          !
          ! psi(:) = psi(:) + dcmplx(coef(k),0.0d0)   * &
          !      zexp( dcmplx(0.0d0,-E(k)*t_vec(t)) ) * &
          !      dcmplx(QP_mtx(:,k),0.0d0)
          psi(:) = psi(:) + dcmplx ( &
               + coef(k) * dcos(E(k)*t_vec(t)) * QP_mtx(:,k) , &
               - coef(k) * dsin(E(k)*t_vec(t)) * QP_mtx(:,k)   )
          !
20     enddo
       !
       !WavePackets computed
       !Now we must obtain P = |psi|**2
       !
       Prob(:) = abs(psi(:))**2.0d0
       auxX(:) = Prob(:)*x_vec(:) 
       x_aux   = xstep*sum(auxX)
       !
       write(21,*) ( Prob(i) + Emean, i=1,xdim )
       write(31,*) dble(x_aux)
       !
10  enddo
    !
    close(21)
    close(31)
    !
  end Subroutine Prob_dens
  !  
  !*****************************************************************************
  !
  Subroutine E_mean_WP
    !
    ! This routine computes the E mean value: < Psi | H | Psi >
    !
    implicit none
    !
    Integer:: i
    !
    Emean=0.0d0
    !
    do i = 0,Nval
       !
       Emean = Emean + coef(i)*coef(i)*E(i)
       !
    enddo
    !
  end Subroutine E_mean_WP
  !
end Module QP_tmp_ev
