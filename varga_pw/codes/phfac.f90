SUBROUTINE phase_factor()
  USE Gvector
  USE PSEUDOPOTENTIAL
  USE PW
  IMPLICIT NONE 

  INTEGER :: i, i1, i2, i3, is, ia, j, k, ik, ig, igp, iai,ii1,ii2,ii3
  REAL(8) :: ss, phi1, phi2, phi3, xmat(3,3), xmati(3,3)
  REAL(8) :: taup(3)
  INTEGER, EXTERNAL :: iflip

  xmat = Lattice_vector
  CALL inv_r(xmat,3,xmati)
       
  DO is=1,n_species

    DO ia=1,n_atom(is)
      
      DO i=1,3
        ss =      xmati(i,1)*atpos(1,ia,is)
        ss = ss + xmati(i,2)*atpos(2,ia,is)
        ss = ss + xmati(i,3)*atpos(3,ia,is)
        taup(i) = ss
      ENDDO
      
      phi1 = -2.d0*pi*taup(1)
      phi2 = -2.d0*pi*taup(2)
      phi3 = -2.d0*pi*taup(3)
      
      ! G=0
      ei1(0,ia,is) = (1.d0, 0.d0)
      ei2(0,ia,is) = (1.d0, 0.d0)
      ei3(0,ia,is) = (1.d0, 0.d0)

      DO i = 1,(N_L(1)+1)/2
        ei1(i,ia,is)  = cmplx( cos(i*phi1), sin(i*phi1) )
        ei1(-i,ia,is) = conjg( ei1(i,ia,is) )
      ENDDO
      
      DO j = 1,(N_L(2)+1)/2
        ei2(j,ia,is)  = cmplx( cos(j*phi2), sin(j*phi2) )
        ei2(-j,ia,is) = conjg( ei2(j,ia,is) )
      ENDDO 
          
      DO k = 1,(N_L(3)+1)/2
        ei3(k,ia,is)  = cmplx( cos(k*phi3), sin(k*phi3) )
        ei3(-k,ia,is) = conjg( ei3(k,ia,is) )
      ENDDO

    ENDDO ! loop over all atoms of species

  ENDDO  ! loop over species species 
      
  DO ik=1,n_k_points
    DO  is=1,n_species
      DO  ia=1,n_atom(is)
        DO  ig=1,n_g_vector(ik)
          igp = G_index(ig,ik)
          ii1 = iflip(G_vector(1,igp),N_L(1))
          ii2 = iflip(G_vector(2,igp),N_L(2))
          ii3 = iflip(G_vector(3,igp),N_L(3))
          eigr(ig,ia,is,ik) = ei1(ii1,ia,is)*ei2(ii2,ia,is)*ei3(ii3,ia,is)
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
      
END SUBROUTINE
