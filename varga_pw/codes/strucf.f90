SUBROUTINE structure_factor()
  USE GVECTOR
  USE PSEUDOPOTENTIAL
  IMPLICIT NONE 

  COMPLEX(8) :: stf
  INTEGER :: ia, is, ig,ii1,ii2,ii3
  INTEGER, EXTERNAL :: iflip
  
  DO is = 1,n_species
    DO ig = 1,n_G_vector_max
      ii1 = iflip(G_vector(1,ig),N_L(1))
      ii2 = iflip(G_vector(2,ig),N_L(2))
      ii3 = iflip(G_Vector(3,ig),N_L(3))
      stf = (0.d0, 0.d0)
      DO ia = 1,n_atom(is)
        stf = stf + ei1(ii1,ia,is)*ei2(ii2,ia,is)*ei3(ii3,ia,is)
      ENDDO 
      sfac(is,ig) = stf
    ENDDO 
  ENDDO 
  
END SUBROUTINE
