module QP_module
  !
  ! Global definitions
  !
  use f95_lapack, only : la_syevr
  !
  implicit none
  !
  Double precision:: a,b,c,d,tol !Potential parameters
  Integer:: Ncon, Nval !Convergence dimension and studied levels
  Double precision,allocatable:: &
       H(:,:),&     !Hamiltonian
       E(:),&       !Energies
       IPR_vec(:),& !IPR
       x_prod(:)    !X mean value
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
  !*****************************************************************************
  !
  Function IPR_fun(V)
    !
    !This function compute the IPR
    !
    implicit none
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
  !
  !*****************************************************************************
  !
  Subroutine Convergence
    !
    !This routine studies the energies convergence.
    !The program will increasses the truncation dim with step of 10 and look for
    !convergence in the three eigenstates E(Nval), E(Nval-1), E(Nval-2)
    !
    implicit none
    !
    Double precision:: Eold(1:3) !Last eigenvalues
    Integer:: k
    !
    Eold=E(Nval-2:Nval)
    !
    deallocate(H,E)
    !
    Ncon=Ncon+10 !Increasing HS dim
    !
    allocate(H(0:Ncon,0:Ncon),E(0:Ncon)) !New allocation
    !
    call Hamiltonian(Ncon) !A bigger system
    !
    E=0.0d0 !Initialized
    !
    call la_syevr(A=H, W=E, JOBZ='V', UPLO='U')
    !
    Ezero=minval(E)
    !
    E=E-Ezero
    !
    do while ( & !Convergence's criteria
         abs(Eold(1)-E(Nval-2)) .gt. tol .and. &
         abs(Eold(2)-E(Nval-1)) .gt. tol .and. &
         abs(Eold(3)-E(Nval  )) .gt. tol        )
       !
       Eold=E(Nval-2:Nval)
       !
       Ncon=Ncon+10 !Increasing HS dim
       !
       deallocate(H,E)
       !
       allocate(H(0:Ncon,0:Ncon),E(0:Ncon)) !New allocation
       !
       call Hamiltonian(Ncon) !A bigger system
       !
       E=0.0d0 !Initialized
       !
       call la_syevr(A=H, W=E, JOBZ='V', UPLO='U')
       !
       Ezero=minval(E)
       !
       E=E-Ezero
       !
    enddo
    !
  end Subroutine Convergence
  !   
  !*****************************************************************************
  !
  Subroutine mean_x
    !
    !This subroutine compute <k|x|k>
    !
    implicit none
    !
    integer i,n
    !
    x_prod=0.0d0
    !
    do i=0,Nval
       !
       x_prod(i)= H(0,i) *  H(1,i)
       !
       do n=1,Ncon
          !
          x_prod(i)=x_prod(i) +  H(n,i) * &
               ( sqrt(dble(n)+1.0d0) * H(n+1,i) + sqrt(dble(n)) * H(n-1,i) )
          !
       enddo
       !
       x_prod(i)=x_prod(i)/sqrt(2.0d0)
       !
    enddo
    !
  end Subroutine mean_x
  !
end module QP_module

  
