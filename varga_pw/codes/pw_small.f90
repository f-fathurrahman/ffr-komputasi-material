MODULE PW_SMALL
  IMPLICIT NONE 
! plane wave cutoff
  double precision                :: E_cut_small
! cutoff for potential and wf
  double precision                :: G_cut_pot,G_cut_wf
! Number of G points for the potential and wf
  INTEGER                          :: N_G_pot
  INTEGER ,allocatable             :: N_G_wf(:)
  INTEGER                          :: N_G_wf_max
! table of g_vectors for the wave function
  INTEGER ,allocatable             :: wf_tab(:,:),g_map(:,:)

CONTAINS  

SUBROUTINE pw_small_init
USE GVECTOR
  IMPLICIT NONE 
  real*8                           :: t,pi
  INTEGER                           :: k,i,ii,ma,jj
  write(6,*)'?'  
  pi=4.d0*atan(1.d0)
  t=(2.0*pi/Lattice_constant)**2
  g_cut_wf=e_cut_small/t
  g_cut_pot=4.0d0*g_cut_wf

! determine the number of potential components
  N_G_pot=0
  do i=1,N_G_vector_max
    if(G_vector_length(i).lt.g_cut_pot) then
      n_g_pot=n_g_pot+1
    endif    
  ENDDO 

! tabulate the pw basis for each k point
  allocate(wf_tab(N_g_pot,N_k_points),N_G_wf(N_k_points))
  ma=0
  wf_tab=0
  do k=1,N_K_points
    ii=0
    jj=0
    do i=1,n_g_vector(k)
      if(gplusk(i,k).lt.g_cut_wf) then
        ii=ii+1
        wf_tab(ii,k)=i
      endif
    ENDDO 
    n_g_wf(k)=ii
    if(ii.gt.ma) ma=ii
  ENDDO 
! create a map to find the matrix elements of the basis functions
  N_G_wf_max=ma
  allocate(g_map(ma,ma))

end SUBROUTINE pw_small_init



SUBROUTINE mapping(k)
USE GVECTOR
! generate a map of  G_i-G_j onto G      
  IMPLICIT NONE 
  INTEGER            :: i1,i2,i3,j1,j2,j3,k1,k2,k3,ii,jj,k,ig,jg,i,j
  INTEGER            :: table(-N_L(1):N_L(1),-N_L(2):N_L(2),-N_L(3):N_L(3))
  INTEGER ,external  :: iflip
!  N_G_Pot+1 --> (0,0)
   do i1=-N_L(1)+1,N_L(1)-1
     do i2=-N_L(2)+1,N_L(2)-1
       do i3=-N_L(3)+1,N_L(3)-1
         table(i1,i2,i3)=N_G_pot+1
       ENDDO 
     ENDDO 
   ENDDO 
   do i=1,N_G_pot
     i1=iflip(G_vector(1,i),N_L(1))
     i2=iflip(G_vector(2,i),N_L(2))
     i3=iflip(G_vector(3,i),N_L(3))
     table(i1,i2,i3)=i
   ENDDO 
   do i=1,N_G_wf(k)
     ig=G_index(wf_tab(i,k),k)
     i1=iflip(G_vector(1,ig),N_L(1))
     i2=iflip(G_vector(2,ig),N_L(2))
     i3=iflip(G_vector(3,ig),N_L(3))
     do j=1,i
       jg=G_index(wf_tab(j,k),k)
       j1=iflip(G_vector(1,jg),N_L(1))
       j2=iflip(G_vector(2,jg),N_L(2))
       j3=iflip(G_vector(3,jg),N_L(3))
       G_map(i,j)=table(i1-j1,i2-j2,i3-j3)
     ENDDO 
   ENDDO 

end SUBROUTINE mapping


END MODULE PW_SMALL



