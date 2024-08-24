      SUBROUTINE calculate_potential(rho)
!     input density output potential
      use fft_data
      USE GVECTOR
      USE PSEUDOPOTENTIAL
      USE PW
      IMPLICIT NONE 


      INTEGER      :: k, is, ig, i, ind,ik, ir1,ir2,ir3
      real*8      :: con,V_XC,wkk


      complex*16  ::  rho_e,rhog,vcg,rho_es,rhogs
      complex*16  ::  E_el,E_ps,vp_loc,ps_charge,E_H,E_G,rps
      complex*8   ::  v_rho_es(N_G_vector_max), vpoten(N_G_vector_max)
      complex*16  ::  c_fft(N_L(1)+fftinc1,N_L(2),N_L(3))
      real*8      ::  rho(N_L(1)+fftinc1,N_L(2),N_L(3))
      real*8      ::  vxc_pot(N_L(1)*N_L(2)*N_L(3))
      real*8      ::  rv(N_L(1)*N_L(2)*N_L(3))

!     FFT of the electron density           
      c_fft=(0.,0.)
      c_fft=rho
      call fft(c_fft,N_L(1),N_L(2),N_L(3),.false.)   

      E_el=(0.0,0.0)
      E_H=(0.0,0.0)
      E_G=(0.0,0.0)
      vp_loc=(0.0,0.0)

      do is=1,n_species
        vp_loc=vp_loc+sfac(is,1)*vps(is,1)
      ENDDO 

!   E_ps is g=0 coulomb energy of the charge density and the
!   local pseudopotential-field of gaussian pseudocharges

      E_ps=vp_loc*conjg(c_fft(1,1,1))
      vpoten(1)=vp_loc



      do  ig=2,N_G_vector_max

!    vp_loc = local part of ps-pot + potential of gaussian pseudo charge, 
!         in fourier space (accurate up to 4*Ecut)
!    ps_charge = pseudo charge of pseudo atom in fourier space

        vp_loc=(0.0,0.0)
        ps_charge=(0.0,0.0)

        do is=1,n_species
          vp_loc=vp_loc+sfac(is,ig)*vps(is,ig)
          ps_charge=ps_charge+sfac(is,ig)*rhops(is,ig)
        ENDDO 
        rho_e=c_fft(G_vector(1,ig),G_vector(2,ig),G_vector(3,ig))
        rhog=rho_e+ps_charge
        rho_es=conjg(rho_e)
        v_rho_es(ig) = rho_es
        rhogs=conjg(rhog)
        rps=conjg(ps_charge)
        con=4.d0*pi/(tpiba2*g_vector_length(ig))
        vcg=con*rhog
        E_el=E_el+0.5*vcg*rhogs
        E_H=E_H+0.5*con*rho_e*rho_es
        E_G=E_G+0.5*con*ps_charge*rps
        E_ps=E_ps+rho_es*vp_loc
        vpoten(ig) = vcg+vp_loc
      ENDDO 
      
      do ig=1,N_G_vector_max
        c_fft(G_vector(1,ig),G_vector(2,ig),G_vector(3,ig)) = vpoten(ig)
      ENDDO 
      call fft(c_fft,N_L(1),N_L(2),N_L(3),.true.)   

      ind = 0
      do ir3=1,N_L(3)
        do ir2=1,N_L(2)
          do ir1=1,N_L(1)
            ind=ind+1
            rv(ind) = rho(ir1,ir2,ir3)
          ENDDO 
        ENDDO 
      ENDDO 
      call lda_xc(rv,vxc_pot,E_exchange,N_L(1)*N_L(2)*N_L(3))

      V_XC = 0.0
      ind = 0
      do ir3=1,N_L(3)
        do ir2=1,N_L(2)
          do ir1=1,N_L(1)
            ind=ind+1
            V_XC = V_XC + rho(ir1,ir2,ir3)*vxc_pot(ind)

            write(66,*)ind,rho(ir1,ir2,ir3),real(c_fft(ir1,ir2,ir3)),vxc_pot(ind)
            rho(ir1,ir2,ir3)=real(c_fft(ir1,ir2,ir3))+vxc_pot(ind)
          ENDDO 
        ENDDO 
      ENDDO 



      V_XC=V_XC*volume/(N_L(1)*N_L(2)*N_L(3))
      E_exchange=E_exchange*volume/(N_L(1)*N_L(2)*N_L(3))
      E_es=real(E_el)*volume-E_self
      E_Hartree=real(E_H)*volume
      E_Gauss=real(E_G)*volume
      E_Pseudo=real(E_ps)*volume
      E_total=E_kinetic+E_es+E_Pseudo+E_non_local+E_exchange
      E_eigen=0.0
      do ik=1,N_k_points
        wkk=W_k_point(ik)*N_sym
        do i=1,n_orbitals
          E_eigen=E_eigen+wkk*occupation(i,ik)*eigenvalue(i,ik)
        ENDDO 
      ENDDO 

          write(16,100) E_total,E_kinetic,E_es,E_Hartree,E_Pseudo,E_non_local,E_exchange,V_XC,E_eigen,E_self,E_Gauss
100    format(//&
     &'                         total energy = ',f12.6,' a.u.'/&
     &'                       kinetic energy = ',f12.6,' a.u.'/&         
     &'                 electrostatic energy = ',f12.6,' a.u.'/&              
     &'                  real hartree energy = ',f12.6,' a.u.'/&               
     &'               pseudopotential energy = ',f12.6,' a.u.'/&
     &'           n-l pseudopotential energy = ',f12.6,' a.u.'/&
     &'          exchange-correlation energy = ',f12.6,' a.u.'/&
     &'exchange-correlation potential energy = ',f12.6,' a.u.'/&
     &'          kohn - sham  orbital energy = ',f12.6,' a.u.'/&
     &'                          self energy = ',f12.6,' a.u.'/&
     &'                      gaussian energy = ',f12.6,' a.u.'//)



      end
