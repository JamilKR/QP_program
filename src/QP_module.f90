module QP_module
  !
  ! Global definitions
  !
  implicit none
  !
  Double precision:: a,b,c,d !Potential parameters
  Integer:: Ncon, Nval !Convergence dimension and studied levels
  Double precision,allocatable:: &
       H(:,:),&    !Hamiltonian
       E(:),&      !Energies
       IPR_vec(:)  !IPR
  Double precision:: Ezero !Zero energy
  Logical:: CONV,IPR 
  !
Contains
  Subroutine Hamiltonian(Nmax)
    !
    ! This routine builds the diagonal and upper elements of
    ! the Hamiltonian matrix
    !
    implicit none
    !
    Integer:: n,Nmax
    !
    H=0.0d0 !Initialized
    !
    do n=0,Nmax !Diagonal
       H(n,n) = 0.25d0 * (3.0d0*a+2.0d0*c+1.0d0) + (0.5d0+3.0d0*a+c)*dble(n) &
            +1.5d0*a*dble(n)*dble(n-1)
    enddo
    !
    do n=0,Nmax-1 !1st Diagonal
       H(n,n+1) = ( (1.5d0*b+d) + 1.5d0*b*dble(n)) * sqrt(dble(n)+1.0d0) &
            /sqrt(2.0d0)
    enddo
    !
    do n=0,Nmax-2 !2nd Diagonal
       H(n,n+2) = sqrt(dble(n)+1.0d0) *  &
            (0.25d0*(-1.0d0+6.0d0*a+2.0d0*c)*sqrt(dble(n)) + &
            a*dble(n)*sqrt(dble(n)+2.0d0) )
    enddo
    !
    do n=0,Nmax-3 !3rd Diagonal
       H(n,n+3) = b*(0.5d0/sqrt(2.0d0))*sqrt( dble(n+1)*dble(n+2)*dble(n+3))
    enddo
    !
    do n=0,Nmax-4
       H(n,n+4) = a*0.25d0*sqrt( dble(n+1)*dble(n+2)*dble(n+3)*dble(n+4) )
    enddo
    !
  end Subroutine Hamiltonian
  !
  Function IPR_fun(V)
    !
    !This function compute the IPR
    !
    Double precision:: V(0:Ncon,0:Ncon)
    Double precision:: IPR_fun(0:Nval)
    Integer:: i,j
    !
    IPR_fun=0.0d0
    !
    do j = 0,Nval
       do i = 0,Ncon
          IPR_fun(j)=IPR_fun(j)+(V(i,j))**4.0d0
       enddo
       IPR_fun(j)=1.0d0/IPR_fun(j)
    end do
    !
  end Function IPR_fun
  
end module QP_module

  