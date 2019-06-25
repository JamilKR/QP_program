program QP_program
  !
  use Harmonic_oscillator
  use QP_module
  use QP_tmp_ev
  use f95_lapack, only : la_syevr
  !
  implicit none
  !
  Character(len=100):: nm_file
  integer i,j
  !
  Namelist/HIL/ Ncon,Nval
  Namelist/POT/ a,b,c,d
  Namelist/OUT/CONV,IPR,tol,temp
  Namelist/XAXIS/x0,xdim,xstep
  Namelist/GAUSS/sigma, mu
  Namelist/TMP/tmax,tstep,outfile_tmp
  !
  read(*,*) nm_file
  !
  open(11,file=trim(nm_file),status='old',action='read')
  !
  read(11,HIL)
  read(11,POT)
  read(11,OUT)
  read(11,XAXIS)
  read(11,GAUSS)
  read(11,TMP)
  !
  close(11)
  !
  !
  allocate(H(0:Ncon,0:Ncon),E(0:Ncon),x_prod(0:Nval))
  !
  call Hamiltonian(Ncon)
  !
  E=0.0d0 !Initialized
  !
  call la_syevr(A=H, W=E, JOBZ='V', UPLO='U')
  !
  Ezero=minval(E)
  !
  E=E-Ezero
  !
  call mean_x(Ncon) !<x> criteria added to convergence subroutine
  !
  if (CONV) then
     !
     call Convergence
     !
  endif  
  !
  write(*,*) 
  write(*,*) 
  write(*,14)
  write(*,15) Ncon,Nval
  write(*,16) a,b,c,d
  write(*,17) Ezero
  write(*,*) 
  write(*,*) 
  !
  if (IPR) then
     !
     allocate(IPR_vec(0:Nval))
     !
     IPR_vec = IPR_fun(H)
     !
     !call mean_x Not necessary, called beafore
     !
     write(*,'(A)') '# k - E_k - IPR(k) - <k|x|k>'
     !
     do i=0,Nval
        write(*,*) i,E(i),IPR_vec(i),x_prod(i)
     enddo
     !
     open(unit=111,file="EF_HO_basis.dat",status='replace')
     !
     write(111,14)  
     write(111,15) Ncon,Nval
     write(111,16) a,b,c,d
     write(111,17) Ezero
     write(111,'(A)') "# Each line = eigenstate in HO basis"
     !
     do i=1,Nval
        write(111,*) (H(j,i),j=1,Ncon)
     enddo
     !
     close(111)
     !
  else
     !
     write(*,'(A)') '# k - E_k'

     !
     do i=0,Nval
        write(*,*) i,E(i)
     enddo
     !
  endif
  !
14 format("# Potential V(x)=ax**4+bx**3+cx**2+dx")
15 format('# Hilbert Space truncated:',I5,';  States considered: g.s. + ',I5)
16 format('# a = ',F10.4,'; b= ',F10.4,'; c= ',F10.4,'; d= ',F10.4)
17 format('# Zero Point Energy:',F10.4 )
  !
  Select case (temp)
     !
  Case ('g') 
     !
     allocate(HO(xdim,0:Ncon))
     allocate(x_vec(xdim))
     !
     x_vec(1)=x0
     !
     do i=2,xdim
        x_vec(i)=x_vec(i-1)+xstep
     enddo
     !        
     HO = HO_WF(Ncon,x_vec) ! Harmonic Oscillator WFs
     !
     call QP_WF_routine !Quartic Potential eigenfunctions > QP_mtx(x,psi_k)
     !
     call Coefficients    ! WP Coefficients in eigenvectors basis
     !
     call Trans_x_matrix  ! x_mtx(i,j)= < psi_i | x | psi_j >
     !
     call E_mean_WP ! Energy mean value computed.
     !
     !Temporal evolution
     !
     tdim = int(tmax/tstep)
     !
     allocate(t_vec(1:tdim))
     !
     t_vec(1)=0.0d0
     !
     do i=2,tdim
        t_vec(i) = t_vec(i-1) + tstep
     enddo
     !
     call Prob_dens
     !
  Case Default
     !
     write(*,*)
     write(*,'(A)') "Temporal evolution not required"
     write(*,*)
     STOP
     !
  end Select
  !
     
  !
end program QP_program
