SUBROUTINE fermi
  USE PW  
  USE PSEUDOPOTENTIAL
  USE gvector
  IMPLICIT NONE 
  double precision,  allocatable ::  eigen_value(:)  
  double precision,  allocatable ::  eigen_value_0(:) 
  double precision,  allocatable ::  fik(:)   
  INTEGER          ,  allocatable ::  ind(:)
  double precision,  parameter   ::  au=2.0*13.6058,eps=1.d-07         
  INTEGER          ,  parameter   ::  N_fermi_iter=200
  INTEGER                         ::  ii,k,ik,ic,nfermi,N_states,i
  double precision               ::  d0,d1,N_e,x,f,min_fermi,max_fermi

  N_states=N_orbitals*N_K_points
  allocate(eigen_value(N_states),eigen_value_0(N_states),fik(N_states),ind(N_states))

  do ik=1,N_K_points
    do i=1,N_orbitals
      ii=(ik-1)*N_orbitals+i
      eigen_value(ii)=eigenvalue(i,ik)*au
    ENDDO 
  ENDDO 
  call indexx(N_states,eigen_value,ind)
  N_e=0.0
  do ii=1,N_states
    eigen_value_0(ii)=eigen_value(ind(ii))
  ENDDO 

  do ii=1,N_states
    d0=N_electrons-N_e
    k=(ind(ii)-1)/N_orbitals+1
    N_e=N_e+W_K_Point(k)*N_sym
    d1=N_electrons-N_e
    if (abs(d1).gt.abs(d0)) then
      nfermi=ii-1
      if (N_electrons/2.-int(N_electrons/2.) < 0.1) then
        E_fermi=(eigen_value_0(nfermi)+eigen_value_0(nfermi+1))/2.0
      else
        E_fermi=eigen_value_0(nfermi)
      endif
      goto 2000
    endif
  ENDDO 
 2000   continue
  write(16,*) 'Fermi level:', E_fermi
  ic=0
  max_fermi =  40000.0    
  min_fermi = -40000.0
 1000 continue

  do ii=1,N_states
    x=(eigen_value(ii)-E_fermi)/kt
    if (x.le.-30.0) then
      fik(ii)=2.0
    else if (x.lt.30.0) then
      fik(ii)=2.0/(1.0+exp(x))
    else
      fik(ii)=0.0
    endif
  ENDDO 
  f = 0.0d0
  do  ii=1,N_states
    k=(ii-1)/N_orbitals+1
    f=f+fik(ii)*W_K_Point(k)*N_sym     
  ENDDO 
  if (abs(f-N_electrons).ge.eps) then
    ic=ic+1
    if (ic.gt.n_fermi_iter) then
       write(16,*) 'Fermi energy not found'
    else
      if (f-N_electrons .gt. 0.0 .and. E_fermi .lt. max_fermi) max_fermi = E_fermi
      if (f-N_electrons .lt. 0.0 .and. E_fermi .gt. min_fermi) min_fermi = E_fermi
      E_fermi = (max_fermi+min_fermi)/2.0
      goto 1000
    endif
  endif

  do ii=1,N_states
    ik=(ii-1)/N_orbitals+1
    i=ii-(ik-1)*N_orbitals
    occupation(i,ik)=fik(ii)
  ENDDO 
  write(16,*) 'Fermi level:', E_fermi

end
