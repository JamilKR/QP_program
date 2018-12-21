program QP_program
  !
  use QP_module
  use f95_lapack, only : la_syevr
  !
  implicit none
  !
  Character(len=100):: nm_file
  integer i
  !
  Namelist/HIL/ Ncon,Nval
  Namelist/POT/ a,b,c,d
  Namelist/OUT/CONV,IPR,tol
  !
  read(*,*) nm_file
  !
  open(11,file=trim(nm_file),status='old',action='read')
  !
  read(11,HIL)
  read(11,POT)
  read(11,OUT)
  !
  close(11)
  !
  !
  allocate(H(0:Ncon,0:Ncon),E(0:Ncon))
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
14 format("# Potential V(x)=ax**4+bx**3+cx**2+dx")
15 format('# Hilbert Space truncated:',I5,';  States considered: g.s. + ',I5)
16 format('# a = ',F6.2,'; b= ',F6.2,'; c= ',F6.2,'; d= ',F6.2)
17 format('# Zero Point Energy:',F10.4 )
  !
  if (IPR) then
     !
     allocate(IPR_vec(0:Nval),x_prod(0:Nval))
     !
     IPR_vec = IPR_fun(H)
     !
     call mean_x
     !
     write(*,'(A)') '# k - E_k - IPR(k) - <k|x|k>'
     !
     do i=0,Nval
        write(*,*) i,E(i),IPR_vec(i),x_prod(i)
     enddo
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
end program QP_program
