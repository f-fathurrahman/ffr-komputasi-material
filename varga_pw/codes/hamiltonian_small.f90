SUBROUTINE hamiltonian_small(potential)
!
! solution on a small PW basis by direct diagonalization of the PW Hamiltonian
!  
  use fft_data
  USE GVECTOR
  USE PW_SMALL
  USE PW
  USE PSEUDOPOTENTIAL
  IMPLICIT NONE 

  INTEGER                       :: jj, jg, ijg,i, ik, is,ia,ig,ii,jn,lm
  complex*16                   :: w,ww
  complex*16                   :: c_fft(N_L(1)+fftinc1,N_L(2),N_L(3))
  double precision             :: potential(N_L(1)+fftinc1,N_L(2),N_L(3))
  complex*16,allocatable       :: hm(:,:),eval(:),evec(:,:),vg(:),om(:,:)

  allocate (vg(N_G_pot+1))

!  transform the potential to G-space        
   c_fft=(0.0,0.0)
   potential(N_L(1)+1,:,:)=0.0d0
   c_fft(:,:,:) = cmplx( potential(:,:,:), aimag(c_fft(:,:,:)) )
   call fft(c_fft,N_L(1),N_L(2),N_L(3),.false.)

   do ig=1,n_G_pot
     vg(ig)=c_fft(G_vector(1,ig),G_vector(2,ig),G_vector(3,ig))
   ENDDO 
   vg(N_G_pot+1)=(0.0,0.0)

!  loop over K points
do  ik=1,n_k_points

  allocate(hm(N_G_wf(ik),N_G_wf(ik)),om(N_G_wf(ik),N_G_wf(ik)), &
&         evec(N_G_wf(ik),N_G_wf(ik)),eval(N_G_wf(ik)))

!  setup the Hamiltonian
        call mapping(ik)
        hm=(0.,0.)
        do   ig=1,N_G_wf(ik)
          ii=wf_tab(ig,ik)
          hm(ig,ig)=0.5*gplusk(ii,ik)*tpiba2
          do is=1,n_species
            do ia=1,n_atom(is)
              w=conjg(eigr(ii,ia,is,ik))
              do lm=1,l_max(is)**2+1-2*l_loc(is)
                 ww = pkg_a(lm,ii,is,ik)*wnl(is,lm)*w
                 do jg=1,ig
                    jj=wf_tab(jg,ik)
                    hm(jg,ig)=hm(jg,ig)+ww*pkg_a(lm,jj,is,ik)*eigr(jj,ia,is,ik)
                 ENDDO 
              ENDDO 
            ENDDO 
          ENDDO 
       ENDDO 

   om=(0.d0,0.d0)

  do ig=1,n_G_wf(ik)
    om(ig,ig)=(1.d0,0.d0)
    do jg=1,ig
      ijg=g_map(ig,jg)
      hm(jg,ig) = hm(jg,ig)+ conjg(vg(ijg))
      hm(ig,jg)=conjg(hm(jg,ig))
    ENDDO 
  ENDDO 


  if(N_orbitals.gt.N_G_wf(ik)) then
    write(6,*)'increase initial energy cut off '
    write(6,*)N_orbitals,N_G_wf(ik)
    stop
  endif

  write(500,*)N_G_wf(ik)
  do ig=1,N_G_wf(ik)
    do jg=1,N_G_wf(ik)
      write(500,*)hm(ig,jg)
    ENDDO 
  ENDDO 
   
  call diag_lapack_complex_generalized(hm,om,N_G_wf(ik),eval,evec)

  do i=1,N_G_wf(ik)
    write(501,*)i,eval(i)
  ENDDO 
  

!  do i=1,N_orbitals
!    read(101,*)(evec(ig,i),ig=1,n_g_wf(ik))
!  ENDDO 
  

  do i=1,N_orbitals
    eigenvalue(i,ik)=eval(i)
    wave_function_c(:,i,ik)=(0.0,0.0)        
    do ig=1,n_g_wf(ik)
      wave_function_c(wf_tab(ig,ik),i,ik)=evec(ig,i)
    ENDDO 
    write(149,*) i,ik,real(eigenvalue(i,ik))
  ENDDO 

  deallocate(hm,om,evec,eval)
ENDDO 

deallocate(vg) 

      end





