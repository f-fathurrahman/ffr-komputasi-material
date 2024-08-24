SUBROUTINE form_factor()
  USE Gvector
  USE PSEUDOPOTENTIAL
  USE PW

  IMPLICIT NONE 

  INTEGER :: is, ig, ir,ll
  REAL(8) :: facc, vv, rr, rr2, v, integ, g_old, g_new, arg
  REAL(8) :: vloc(N_PP_max),psi(N_PP_max)

  ! the self-energy E_self of the ionic pseudocharges
  E_self=0.0
  DO is = 1,n_species
    E_self = E_self + charge_pp(is)*charge_pp(is)/beta*n_atom(is)
  ENDDO 
  E_self = E_self/sqrt(2.d0*pi)

  ! form factors of the ionic pseudocharge (rhops)
  DO is = 1,n_species
    facc = 0.25 * tpiba2*beta**2
    DO ig = 2,N_G_vector_max
      rhops(is,ig) = -charge_pp(is)/volume*exp(-facc*g_vector_length(ig))
    ENDDO 
  ENDDO 

  ! the form factors of pseudopotential (vps)
  !  = local part of ps-pot plus field of gaussian pseudo-charge

  DO is=1,n_species
    ll = l_loc(is)
    DO ir=1,N_PP(is)
      vv = -charge_pp(is)/r(ir,is)*derf(r(ir,is)/beta)
      vloc(ir) = vion(ir,is,ll) - vv
    ENDDO 

    ! g = 0 component of pseudo-pot
    DO ir=1,N_PP(is)
      psi(ir) = vloc(ir)*r(ir,is)**3
    ENDDO 
    integ = sum(psi(1:N_PP(is)))*cclog(is)
    vps(is,1) = integ*4.d0*pi/volume


    g_old = -1.d0


    DO ig=2,N_G_vector_max
      g_new = g_vector_length(ig)
      IF( g_new /= g_old) THEN 
        g_old = g_new
        DO ir = 1,N_PP(is)
          arg = r(ir,is)*sqrt(g_new)*tpiba
          psi(ir) = vloc(ir)*sin(arg)/arg*r(ir,is)**3
        ENDDO 
        integ = sum(psi(1:N_PP(is)))*cclog(is)
      ENDIF
      vps(is,ig) = integ*4.d0*pi/volume
    ENDDO 
  ENDDO 

END SUBROUTINE 


SUBROUTINE nl_pp_form_factor()
  
  USE GVECTOR
  USE PSEUDOPOTENTIAL
  USE PW
  
  IMPLICIT NONE 
  INTEGER :: is, ir, ll, i_l, m_l, ik, ig, igp, i, j
  REAL(8) :: fpibo, rr2, v, rr, dv, integ, t, arg, fac,ag
  REAL(8) :: psi(N_PP_max,N_species),psi1(N_PP_max)
  REAL(8) :: cosx,cosy,cosz,gg(3)
  INTEGER :: i_m,j_m,mlmax
  INTEGER, EXTERNAL :: iflip
  INTEGER :: ii1,ii2,ii3
      
  DO is = 1,n_species
    
    m_l = 0
    
    DO i_l = 1,l_max(is)
      
      IF( i_l /= l_loc(is) ) THEN 
        
        DO ir=1,N_PP(is)
          rr=r(ir,is)
          rr2=rr*rr
          dv=vion(ir,is,i_l)-vion(ir,is,l_loc(is))
          psi(ir,is)=dv*p(ir,is,i_l)*rr2
          psi1(ir)=dv*p(ir,is,i_l)**2*rr 
        ENDDO 
        
        integ = sum(psi1(1:N_PP(is)))*cclog(is)
        
        DO j = 1,2*i_l-1
          wnl(is,m_l+j) = (2.0*i_l-1.0)*4.d0*pi/volume/integ
        ENDDO 
        
        DO ik=1,N_k_points
        DO ig=1,n_g_vector(ik)

          ii1 = iflip(G_vector(1,G_index(ig,ik)),N_L(1))
          ii2 = iflip(G_vector(2,G_index(ig,ik)),N_L(2))
          ii3 = iflip(G_vector(3,G_index(ig,ik)),N_L(3))
          
          gg(:) = ii1*R_lattice_vector(:,1) + &
                  ii2*R_lattice_vector(:,2) + &
                  ii3*R_lattice_vector(:,3) + k_point(:,ik)
    
          t = sqrt(gplusk(ig,ik))
          
          arg = t*tpiba

          IF( i_l == 1 ) THEN 
            IF( t < 1.d-4 ) THEN 
              DO ir = 1,N_PP(is)
                psi1(ir) = psi(ir,is)
              ENDDO 
            ELSE 
              DO ir = 1,N_PP(is)
                fac = arg*r(ir,is)
                fac = sin(fac)/fac
                psi1(ir) = fac*psi(ir,is)
              ENDDO 
            ENDIF
            integ = sum(psi1(1:N_PP(is)))*cclog(is)
            pkg_a(m_l+1,ig,is,ik) = integ
          ENDIF

          IF( i_l == 2 ) THEN 
            IF( t < 1.d-4 ) THEN    
              DO j=1,3
                pkg_a(m_l+j,ig,is,ik) = 0.d0
              ENDDO 
            ELSE 
              DO ir = 1,N_PP(is)
                fac = arg*r(ir,is)
                fac = (sin(fac)/fac-cos(fac))/fac
                psi1(ir) = fac*psi(ir,is)
              ENDDO 
              integ = sum(psi1(1:N_PP(is)))*cclog(is)
              pkg_a(m_l+1,ig,is,ik) = integ*gg(1)/t
              pkg_a(m_l+2,ig,is,ik) = integ*gg(2)/t
              pkg_a(m_l+3,ig,is,ik) = integ*gg(3)/t
            ENDIF 
          ENDIF 

          IF( i_l == 3 ) THEN 
            IF( t < 1.d-4 ) THEN 
              DO j = 1,5
                pkg_a(ig,is,ik,m_l+j) = 0.d0
              ENDDO 
            ELSE 
              DO ir=1,N_PP(is)
                ag = arg*r(ir,is)
                fac = (3.0/(ag**2.0)-1.0)*sin(ag)-3.0/ag*cos(ag)
                fac = fac/ag
                psi1(ir) = fac*psi(ir,is)
              ENDDO 
              integ = sum(psi1(1:N_PP(is)))*cclog(is)
              cosx = gg(1)/t
              cosy = gg(2)/t
              cosz = gg(3)/t
              pkg_a(m_l+1,ig,is,ik)=integ*0.5*(3.0*cosz**2.0-1.0)
              pkg_a(m_l+2,ig,is,ik)=integ*sqrt(3.d0)*cosz*cosx
              pkg_a(m_l+3,ig,is,ik)=integ*sqrt(3.d0)*cosz*cosy
              pkg_a(m_l+4,ig,is,ik)=integ*sqrt(3.d0)*cosx*cosy
              pkg_a(m_l+5,ig,is,ik)=integ*sqrt(3.d0)/2.0*(cosx**2.0-cosy**2.0)                    
            ENDIF 
          ENDIF 
        ENDDO 
        ENDDO 

        m_l = m_l + 2*i_l - 1
      
      ENDIF ! if l is not l_loc
    
    ENDDO ! loop over l

  ENDDO ! loop over species
 
END SUBROUTINE 




SUBROUTINE nl_pp()

  USE GVECTOR
  USE PSEUDOPOTENTIAL
  USE PW

  IMPLICIT NONE 

  INTEGER :: k, is, ia, mlmax, iik, ik, m, i, ii, ig, igp
  INTEGER :: i_lm, lm_end
  REAL(8) :: e_nl
  COMPLEX(8) :: t_eigr_pkg(N_G_K_vector_max,(L_PP_max+1)**2)
  COMPLEX(8) :: sf0, tt1

  INTEGER :: i_c, iani, iai
  COMPLEX(8) :: ct,sf

  E_non_local = 0.d0

  DO is = 1,n_species

    DO ia = 1,n_atom(is)
            
      DO iik = 1,N_k_points
        
        ik = iik
        lm_end = l_max(is)**2 + 1 - 2*l_loc(is)
        
        DO i_lm = 1,lm_end
          
          DO ig = 1,n_g_vector(ik)
            t_eigr_pkg(ig,i_lm) = conjg(eigr(ig,ia,is,ik))*pkg_a(i_lm,ig,is,ik)
          ENDDO 
          
          DO i = 1,n_orbitals
            ct = (0.d0,0.d0)
            DO ig = 1,n_g_vector(ik)
              ct = ct + wave_function_c(ig,i,ik)*t_eigr_pkg(ig,i_lm)
            ENDDO 
            fnl(i,ik,is,ia,i_lm) = ct                   
          ENDDO                  
          
          e_nl = 0.0
          
          DO i=1,n_orbitals
            sf = fnl(i,ik,is,ia,i_lm)
            e_nl = e_nl + occupation(i,ik) * ( dble(sf)**2 + aimag(sf)**2 )
          ENDDO 
          E_non_local = E_non_local + N_sym*w_k_point(ik)*e_nl * wnl(is,i_lm)
        ENDDO  
      
      ENDDO ! end loop over kpts  
    
    ENDDO ! loop over Na(is)
  
  ENDDO ! loop over species  

END SUBROUTINE 

