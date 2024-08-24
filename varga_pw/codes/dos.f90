SUBROUTINE calculate_dos
  USE gvector
  USE PW
  USE PSEUDOPOTENTIAL
  IMPLICIT NONE  
  double precision           :: min_ev,max_ev,e_int,dos,energy
  INTEGER                     :: i,k,j
  INTEGER ,parameter          :: N_e_points=1000
  double precision,parameter :: width=0.15d0
  double precision,external  :: broadening_function


! energy window for DOS
  min_ev=eigenvalue(1,1)*2*rydberg
  max_ev=eigenvalue(1,1)*2*rydberg
  do k=1,N_k_points      
    do  i=1,N_orbitals
      if(eigenvalue(i,k).gt.max_ev) max_ev=eigenvalue(i,k)*2*rydberg
      if(eigenvalue(i,k).lt.min_ev) min_ev=eigenvalue(i,k)*2*rydberg
    ENDDO 
  ENDDO 
  e_int = max_ev -  min_ev
! construct dos
  do i=1,N_e_points  
    dos=0.0
    energy=min_ev+((i*1.5)/N_e_points-0.1)*e_int
    do k=1,N_k_points
      do j=1,N_orbitals
        dos=dos+broadening_function(energy,eigenvalue(j,k)*2*rydberg,width)*w_k_point(k)
      ENDDO 
    ENDDO 
    write(20,*)energy-e_fermi,dos    
  ENDDO 

end
      
double precision function broadening_function(energy,eig,width)
  IMPLICIT NONE 
  double precision            :: energy,eig,width
  double precision, parameter :: Pi = 3.141592654 
      
! gaussian
!        broadening_function=2.d0/Sqrt(Pi)/width*exp(-((energy-eig)/width)**2.0)   
! lorentz
        broadening_function=2.d0/Pi*width/((energy-eig)**2.0+width**2.0)
        
end


