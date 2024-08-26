MODULE m_molecular_dynamics
  implicit none
  integer            :: use_periodic_boundary_conditions,thermostat_type,output_frequency
  integer            :: MD_steps,N_s,N_Particles,g,N_term,N_histogram,potential_type,using_sw_three_body
  real(8)             :: K_E,P_E,total_energy,T_av,T_0,L(3),delta,dt,rho,q,s,xi
  real(8)             :: nu,box,r_cutoff,pressure,virial,lj_epsilon,lj_sigma
  real(8),allocatable :: r(:,:),v(:,:),f(:,:),r0(:,:),rdf(:)
  character*255      :: input_filename,coordinate_input_filename

  ! input units
  integer, parameter :: movie_file=7, coordinate_output_file=8, msd_file=9
  integer, parameter :: rdf_file=10, coordinate_input_file=11, stat_file=12

  real(8), parameter :: pi=3.14159265358979d0, pi2=2.d0*3.14159265358979d0
  character*255,parameter :: msd_filename="msd.dat"
  character*255,parameter :: rdf_filename="rdf.dat"
  character*255,parameter :: stat_filename="stat.dat"
  character*255,parameter :: movie_filename="movie.xyz"
  character*255,parameter :: coordinate_output_filename="md.out"

CONTAINS

!----------------------
SUBROUTINE initialize()
!----------------------
  integer :: i,j,coordinate_input_file
  real(8)  :: lm

  ! Silicon
  lj_epsilon = 2.1678d0 ! eV
  lj_sigma = 2.0951d0   ! Angstroms

  ! Argon
!  lj_epsilon=0.010d0 ! eV
!  lj_sigma=3.405d0   ! Angstroms

  ! Read parameters
  open(1,file=trim(adjustl(input_filename)))
  read(1,*) dt                ! time step
  read(1,*) T_0               ! desired temperature
  read(1,*) thermostat_type   ! 1=velocity scaling, 2=Andersen, 3=Nose-Hoover
  read(1,*) potential_type    ! 1=Lennard-Jones, 2=Stillinger-Weber
  read(1,*) use_periodic_boundary_conditions ! 0=no PBC, 1=use PBC
  read(1,*) md_steps          ! number of MD steps
  read(1,*) output_frequency  ! number of MD steps between xyz movie snapshots
  read(1,*) N_term            ! number of equilibration steps
  read(1,*) N_histogram       ! number of histogram bins for RDF
  read(1,'(a)') coordinate_input_filename  ! input coordinate file

  if(potential_type==1) then
    read(1,*) r_cutoff ! Cutoff radius for the Lennard-Jones potential
  elseif(potential_type==2) then
    read(1,*) using_sw_three_body
  endif

  if(thermostat_type==2) then
    read(1,*) nu ! collision frequency for the Andersen method
  elseif(thermostat_type==3) then
    read(1,*) q  ! effective mass for the Nose-Hoover method
  endif
  close(1)

  ! Read the volume size, number of particles, particle positions, and particle velocities
  open(coordinate_input_file,file=trim(adjustl(coordinate_input_filename)))
  read(coordinate_input_file,*) L(1), L(2), L(3)
  read(coordinate_input_file,*) N_particles
  allocate(r(3,N_particles),v(3,N_particles),F(3,N_particles),r0(3,N_particles))
  do i=1,N_particles
    read(coordinate_input_file,*) (r(j,i), j=1,3), (v(j,i), j=1,3)
    ! Convert to reduced units
    r(:,i)=r(:,i)/lj_sigma
    v(:,i)=v(:,i)/lj_sigma
  enddo
  L(:)=L(:)/lj_sigma
  close(coordinate_input_file)

  S = 0
  xi = 0
  G = 3*N_particles
  N_s = 0
  Lm = min(L(1),L(2),L(3))
  delta = Lm/(2.d0*N_histogram)
  allocate(rdf(0:N_histogram))
  rdf=0.d0

  ! Calculate the kinetic energy
  K_E=0.d0
  do i=1,N_particles
    K_E=K_E+0.5d0*(v(1,i)**2+v(2,i)**2+v(3,i)**2)
  enddo

  ! Calculate the particle density
  rho=dfloat(N_particles)/product(L)
  write(6,*)"Particle density: ",rho

  ! Open files
  open(stat_file,file=trim(adjustl(stat_filename)))
  open(msd_file,file=trim(adjustl(msd_filename)))
  open(rdf_file,file=trim(adjustl(rdf_filename)))
  open(movie_file,file=trim(adjustl(movie_filename)))
END SUBROUTINE initialize

SUBROUTINE nose_hoover_simulation
  integer            :: i,iter
  real(8)             :: xio,xin,v2,H,ri,J_ii,delta_x_n,d1,d2,gm_n
  real(8),allocatable :: v_n(:,:),v_o(:,:),gm(:,:)
  logical test
  real(8),parameter   :: tol=1.d-14

  allocate(v_n(3,N_particles),v_o(3,N_particles),gm(3,N_particles))

  r = r + dt*v + 0.5d0*dt**2*(F-xi*v)
  v = v + 0.5d0*dt*(F-xi*v)
  call periodic_boundary_conditions

  v2=sum(v**2)
  S = S + xi*dt + (v2-G*T_0)*0.5d0*dt**2/Q
  xi = xi + (v2-G*T_0)*0.5d0*dt/Q
  call force
  v_n=v
  v2=sum(v_n**2)
  xin = xi
  test = .FALSE.
  iter = 0
  do while (.not.test.and.iter.lt.50)
    iter = iter + 1
    xio = xin
    delta_x_n = 0
    v_o=v_n
    gm = -0.5d0*dt*(F-xio*v_o) - (v-v_o)
    J_ii = -(xio*0.5d0*dt+1)
    gm_n=(-v2+G*T_0)*0.5d0*dt/Q-(xi-xio)
    d1 =sum(v_o*gm)*dt/Q-J_ii*gm_n
    d2=-dt*0.5d0*dt*v2/Q+J_ii
    delta_x_n=d1/d2
    v_n = v_n + (gm+0.5d0*dt*v_o*delta_x_n)/J_ii
    v2=sum(v_n**2)
    xin = xio + delta_x_n
    test = .TRUE.
    i=0
    do while (i.LE.N_Particles.AND.test)
    i=i+1
      if (i.LE.N_Particles) THEN
        if(ABS((v_n(1,i)-v_o(1,i))/v_n(1,i)).GT.tol) test = .FALSE.
        if(ABS((v_n(2,i)-v_o(2,i))/v_n(2,i)).GT.tol) test = .FALSE.
        if(ABS((v_n(3,i)-v_o(3,i))/v_n(3,i)).GT.tol) test = .FALSE.
       else
          if(ABS((xin-xio)/xin).GT.tol) test = .FALSE.
       endif
     enddo
  enddo
  v=v_n
  xi = xin
  K_E = v2/2
  TOTAL_ENERGY = K_E + P_E + (xi**2*Q)/2.d0 + G*T_0*S
  write(6,'(2ES16.8)')K_E/(3.d0*dfloat(N_Particles))*2.d0,Total_energy
  deallocate(v_n,v_o,gm)
END SUBROUTINE nose_hoover_simulation

SUBROUTINE output
  ! Read the volume size, number of particles, particle positions, and particle velocities
  integer :: i,j
  open(coordinate_output_file,file=trim(adjustl(coordinate_output_filename)))
  write(coordinate_output_file,'(3ES16.8)')L(1),L(2),L(3)
  write(coordinate_output_file,*)N_particles
  do i=1,N_particles
    write(coordinate_output_file,'(6ES16.8)')(r(j,i),j=1,3),(v(j,i),j=1,3)
  enddo
  close(coordinate_output_file)
END SUBROUTINE output

SUBROUTINE cleanup
  ! Deallocate memory and close files
  deallocate(r,v,F,r0,rdf)
  close(stat_file)
  close(msd_file)
  close(rdf_file)
  close(movie_file)
END SUBROUTINE cleanup

SUBROUTINE lennard_jones_potential
  integer :: i,j
  real(8)  :: xi,yi,zi,En,dx,dy,dz,r2,Vir,virij,enij,fr,r2i,r6i,rc2,ecut,sigma2

  sigma2=lj_sigma*lj_sigma
  rc2=r_cutoff*r_cutoff
  ECUT=4.d0*lj_epsilon*((sigma2/rc2)**6-(sigma2/rc2)**3)

  P_E=0.d0
  En=0.d0
  Vir=0.d0
  F=0.d0
  do i=1,N_particles-1
    xi=r(1,i)
    yi=r(2,i)
    zi=r(3,i)
    do j=i+1,N_particles
      dx = xi - r(1,j)
      dy = yi - r(2,j)
      dz = zi - r(3,j)
      dx = dx - L(1)*ANINT(dx/L(1))
      dy = dy - L(2)*ANINT(dy/L(2))
      dz = dz - L(3)*ANINT(dz/L(3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2<=RC2) then
        r2i=1.d0/r2
        r6i=r2i*r2i*r2i
        enij=4.d0*(r6i*r6i-r6i)-ECUT
        virij=48.d0*(r6i*r6i-0.5d0*r6i)
        En=En+enij
        Vir=Vir+virij
        fr=virij/r2
        F(1,i)=F(1,i)+fr*dx
        F(2,i)=F(2,i)+fr*dy
        F(3,i)=F(3,i)+fr*dz
        F(1,j)=F(1,j)-fr*dx
        F(2,j)=F(2,j)-fr*dy
        F(3,j)=F(3,j)-fr*dz
        P_E=P_E+4.d0*lj_epsilon*((sigma2/r2)**6-(sigma2/r2)**3)
      endif
    enddo
  enddo
END SUBROUTINE lennard_jones_potential

SUBROUTINE sw_potential
  integer :: i,j,k
  real(8)  :: dxij, dyij, dzij, r2ij, rij,dxik, dyik, dzik, r2ik, rik, dxjk, dyjk, &
&   dzjk, r2jk, rjk, vcuij, evcuij, vlj, res, vv2, C2,rhoij, rij_a, rik_a, rjk_a,vcuik, rhoik, &
&   vcujk, rhojk, vcui, vcuj, vcuk,hjik, hijk, hikj,hij, hik, hjk,ri, rj, rk, si, sj, sk,ei, ej, ek, &
&   cosjik, cosijk, cosikj,cosjik_3, cosijk_3, cosikj_3,eiNjik_x, eiNjik_y, eiNjik_z,eiNkij_x, &
&   eiNkij_y, eiNkij_z, ejNijk_x, ejNijk_y, ejNijk_z,ejNkji_x, ejNkji_y, ejNkji_z,ekNikj_x, ekNikj_y, ekNikj_z, &
&   ekNjki_x, ekNjki_y, ekNjki_z,hRij_x, hRij_y, hRij_z,hRik_x, hRik_y, hRik_z,hRjk_x, hRjk_y, hRjk_z, &
&   dF3i_x, dF3i_y, dF3i_z,dF3j_x, dF3j_y, dF3j_z,dF3k_x, dF3k_y, dF3k_z,  hLx, hLy, hLz

  real(8),parameter   :: AA_sw=7.0495562770d0, B_sw=0.6022245584d0, LAM_sw=21.0d0, &
& GAM_sw=1.2d0,a_sw=1.8d0, BP_sw=2.408898232d0,a2_sw=3.24d0
  integer                      ::  P_sw=4, Q_sw=0, NP_sw=-4, NP_1_sw=-5, NQ_sw=0,NQ_1_sw=-1

  P_E=0.d0
  hLx=L(1)/2.0; hLy=L(2)/2.0; hLz=L(3)/2.0
  F=0.d0

   do i=1,N_particles-1
     do j=i+1,N_particles
       rij=0.d0
       vcuij=0.d0
       dxij  = (r(1,i)-r(1,j))
       dyij  = (r(2,i)-r(2,j))
       dzij  = (r(3,i)-r(3,j))
! periodic B.C.
       if (dxij>hLx)       dxij=dxij-L(1)
       if (dxij<-hLx)      dxij=dxij+L(1)
       if (dyij>hLy)       dyij=dyij-L(2)
       if (dyij<-hLy)      dyij=dyij+L(2)
       if (dzij>hLz)       dzij=dzij-L(3)
       if (dzij<-hLz)      dzij=dzij+L(3)
       r2ij = dxij*dxij + dyij*dyij + dzij*dzij
       rij = sqrt(r2ij)
       if (r2ij<a2_sw) then
   vcuij = 1.0/(rij - a_sw)
   evcuij = 0.0
   if (vcuij > -30.0) evcuij = exp(vcuij)
   res = 1.0
   if (NQ_sw.ne.0) res = rij**NQ_sw
         vlj = B_sw * rij**NP_sw - res
         vv2=AA_sw*vlj*evcuij
   C2  = BP_sw*rij**NP_1_sw
   if (Q_sw.ne.0.d0) C2 = C2-Q_sw*rij**NQ_1_sw
   C2 =C2*AA_sw*evcuij
   C2 = C2+vv2*vcuij*vcuij
   C2 = C2/rij
   F(1,i) =F(1,i)+ dxij*C2
   F(1,j) =F(1,j)- dxij*C2
   F(2,i) =F(2,i)+ dyij*C2
   F(2,j) =F(2,j)- dyij*C2
   F(3,i) =F(3,i)+ dzij*C2
   F(3,j) =F(3,j)- dzij*C2
   P_E=P_E+vv2
       endif
if(using_sw_three_body==1) then
       rij_a=0.d0; rhoij=0.d0; vcuij=0.d0
       if (rij<a_sw) then
   rij_a = rij - a_sw;
   rhoij = rij*rij_a*rij_a;
   vcuij = GAM_sw/rij_a;
         if (vcuij > -30) then
           vcuij = exp(vcuij)
         else
           vcuij = 0.0
         endif
       endif
       do k=j+1,N_particles
         hijk=0.0d0; hjik=0.0d0; hikj=0.d0
         hij=0.0d0; hik=0.0d0; hjk=0.d0
   ri = 0.d0;  rj = 0.d0; rk = 0.0d0
   ei = 0.d0;  ej = 0.d0; ek = 0.0d0
   cosjik = 0.d0; cosijk = 0.d0; cosikj = 0.0d0
   cosjik_3 = 0.d0;  cosijk_3 = 0.d0; cosikj_3 = 0.0d0

   eiNjik_x = 0.0; eiNjik_y = 0.0; eiNjik_z = 0.0
   eiNkij_x = 0.0; eiNkij_y = 0.0; eiNkij_z = 0.0
   ejNijk_x = 0.0; ejNijk_y = 0.0; ejNijk_z = 0.0
   ejNkji_x = 0.0; ejNkji_y = 0.0; ejNkji_z = 0.0
   ekNikj_x = 0.0; ekNikj_y = 0.0; ekNikj_z = 0.0
   ekNjki_x = 0.0; ekNjki_y = 0.0; ekNjki_z = 0.0

   rik_a =0.d0;  vcuik = 0.d0; rhoik = 0.d0
   rjk_a =0.d0;  vcujk = 0.d0; rhojk = 0.d0

   dxik  = (r(1,i)-r(1,k))
   dyik  = (r(2,i)-r(2,k))
   dzik  = (r(3,i)-r(3,k))
         if (dxij>hLx)       dxij=dxij-L(1)
         if (dxij<-hLx)      dxij=dxij+L(1)
         if (dyij>hLy)       dyij=dyij-L(2)
         if (dyij<-hLy)      dyij=dyij+L(2)
         if (dzij>hLz)       dzij=dzij-L(3)
         if (dzij<-hLz)      dzij=dzij+L(3)
   r2ik = dxik*dxik + dyik*dyik + dzik*dzik
         rik = sqrt(r2ik)
   dxjk  = (r(1,j)-r(1,k))
   dyjk  = (r(2,j)-r(2,k))
   dzjk  = (r(3,j)-r(3,k))
         if (dxjk>hLx)       dxjk=dxjk-L(1)
         if (dxjk<-hLx)      dxjk=dxjk+L(1)
         if (dyjk>hLy)       dyjk=dyjk-L(2)
         if (dyjk<-hLy)      dyjk=dyjk+L(2)
         if (dzjk>hLz)       dzjk=dzjk-L(3)
         if (dzjk<-hLz)      dzjk=dzjk+L(3)
   r2jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
         rjk = sqrt(r2jk)
   if (r2ik<a2_sw) then
     rik_a = rik - a_sw
     rhoik = rik*rik_a*rik_a
     vcuik = GAM_sw/rik_a
           if (vcuik > -30.0) then
             vcuik = exp(vcuik)
           else
             vcuik = 0.0
           endif
   endif
   if (r2jk<a2_sw) then
     rjk_a = rjk - a_sw
     rhojk = rjk*rjk_a*rjk_a
     vcujk = GAM_sw/rjk_a
           if (vcujk > -30.0) then
             vcujk = exp(vcujk)
           else
             vcujk = 0.0
           endif
         end if
   vcui = vcuij*vcuik
   if (vcui.ne.0.d0) then
     ri = rij*rik
     cosjik = (dxij*dxik+dyij*dyik+dzij*dzik)/ri
     cosjik_3 = cosjik + 1.d0/3.d0
     si = LAM_sw*vcui*cosjik_3
     ei = 2*si/ri
           hjik = si*cosjik_3
         endif
   vcuj = vcuij*vcujk
   if (vcuj.ne.0.d0) then
     rj = rij*rjk
     cosijk = -(dxij*dxjk+dyij*dyjk+dzij*dzjk)/rj
     cosijk_3 = cosijk + 1.d0/3.d0
     sj = LAM_sw*vcuj*cosijk_3
     ej = 2*sj/rj
     hijk = sj*cosijk_3
   endif
   vcuk = vcuik*vcujk
   if (vcuk.ne.0.d0) then
     rk = rik*rjk
     cosikj = (dxjk*dxik+dyjk*dyik+dzjk*dzik)/rk
     cosikj_3 = cosikj + 1.d0/3.d0
     sk = LAM_sw*vcuk*cosikj_3
     ek = 2*sk/rk
           hikj = sk*cosikj_3
         endif

   eiNjik_x = ei*(dxik - rik/rij*cosjik*dxij)
   eiNjik_y = ei*(dyik - rik/rij*cosjik*dyij)
   eiNjik_z = ei*(dzik - rik/rij*cosjik*dzij)

   eiNkij_x = ei*(dxij - rij/rik*cosjik*dxik)
   eiNkij_y = ei*(dyij - rij/rik*cosjik*dyik)
   eiNkij_z = ei*(dzij - rij/rik*cosjik*dzik)

   ejNijk_x = ej*(dxjk + rjk/rij*cosijk*dxij)
   ejNijk_y = ej*(dyjk + rjk/rij*cosijk*dyij)
   ejNijk_z = ej*(dzjk + rjk/rij*cosijk*dzij)

   ejNkji_x = ej*(-dxij - rij/rjk*cosijk*dxjk)
   ejNkji_y = ej*(-dyij - rij/rjk*cosijk*dyjk)
   ejNkji_z = ej*(-dzij - rij/rjk*cosijk*dzjk)

   ekNikj_x = ek*(-dxjk + rjk/rik*cosikj*dxik)
   ekNikj_y = ek*(-dyjk + rjk/rik*cosikj*dyik)
   ekNikj_z = ek*(-dzjk + rjk/rik*cosikj*dzik)

   ekNjki_x = ek*(-dxik + rik/rjk*cosikj*dxjk)
   ekNjki_y = ek*(-dyik + rik/rjk*cosikj*dyjk)
   ekNjki_z = ek*(-dzik + rik/rjk*cosikj*dzjk)

   if(rhoij.ne.0.d0) hij = GAM_sw*(hijk+hjik)/rhoij
   if(rhoik.ne.0.d0) hik = GAM_sw*(hjik+hikj)/rhoik
   if(rhojk.ne.0.d0) hjk = GAM_sw*(hikj+hijk)/rhojk

      hRij_x = hij * dxij
      hRij_y = hij * dyij
      hRij_z = hij * dzij
      hRik_x = hik * dxik
      hRik_y = hik * dyik
      hRik_z = hik * dzik
      hRjk_x = hjk * dxjk
      hRjk_y = hjk * dyjk
      hRjk_z = hjk * dzjk

      dF3i_x = hRij_x + hRik_x - eiNjik_x - eiNkij_x + ejNijk_x + ekNikj_x
      dF3i_y = hRij_y + hRik_y - eiNjik_y - eiNkij_y + ejNijk_y + ekNikj_y
      dF3i_z = hRij_z + hRik_z - eiNjik_z - eiNkij_z + ejNijk_z + ekNikj_z

      dF3j_x = -hRij_x + hRjk_x - ejNijk_x - ejNkji_x + eiNjik_x + ekNjki_x
      dF3j_y = -hRij_y + hRjk_y - ejNijk_y - ejNkji_y + eiNjik_y + ekNjki_y
      dF3j_z = -hRij_z + hRjk_z - ejNijk_z - ejNkji_z + eiNjik_z + ekNjki_z

      dF3k_x = -hRjk_x - hRik_x - ekNikj_x - ekNjki_x + eiNkij_x + ejNkji_x
      dF3k_y = -hRjk_y - hRik_y - ekNikj_y - ekNjki_y + eiNkij_y + ejNkji_y
      dF3k_z = -hRjk_z - hRik_z - ekNikj_z - ekNjki_z + eiNkij_z + ejNkji_z

      f(1,i) = f(1,i)+dF3i_x
      f(2,i) = f(2,i)+dF3i_y
      f(3,i) = f(3,i)+dF3i_z
      f(1,j) = f(1,j)+dF3j_x
      f(2,j) = f(2,j)+dF3j_y
      f(3,j) = f(3,j)+dF3j_z
      f(1,k) = f(1,k)+dF3k_x
      f(2,k) = f(2,k)+dF3k_y
      f(3,k) = f(3,k)+dF3k_z

      P_E=P_E+(hijk + hjik + hikj)
      enddo
endif
    enddo
  enddo
END SUBROUTINE sw_potential

SUBROUTINE force
  if(potential_type==1) then
    call lennard_jones_potential
  else
    call sw_potential
  endif
  pressure=2.d0*K_E/3.d0/product(L)-virial
END SUBROUTINE force

SUBROUTINE periodic_boundary_conditions
  integer :: i,j
  if(use_periodic_boundary_conditions==0) return
  do i=1,N_particles
    do j=1,3
      do while(r(j,i)<0.d0)
        r(j,i)=r(j,i)+L(j)
      enddo
      do while(r(j,i)>L(j))
        r(j,i)=r(j,i)-L(j)
      enddo
    enddo
  enddo
END SUBROUTINE periodic_boundary_conditions

FUNCTION random_gaussian(sigma)
  ! Calculate a random number from a Gaussian distribution
  real(8) :: random_gaussian,sigma,x,y
  call random_number(x)
  call random_number(y)
  random_gaussian=dsqrt(sigma)*dsqrt(-2.d0*dlog(x))*dcos(pi2*y)
END FUNCTION random_gaussian

SUBROUTINE andersen_simulation
  integer :: i
  real(8)  :: rn,su

  ! Update the particle positions and velocities using the velocity Verlet method
  r=r+v*dt+0.5d0*dt**2*F
  v=v+0.5d0*dt*F
  call periodic_boundary_conditions
  call force
  v=v+0.5d0*dt*F

  if(abs(nu)>0.d0) then
    ! Andersen thermostat
    do i=1,N_particles
      call random_number(rn)
      if(rn<(nu*dt)) then
        v(1,i)=random_gaussian(T_0)
        v(2,i)=random_gaussian(T_0)
        v(3,i)=random_gaussian(T_0)
      endif
    enddo
  endif
END SUBROUTINE andersen_simulation

SUBROUTINE velocity_scaling
  real(8) :: temperature

  ! Update the particle positions and velocities using the velocity Verlet method
  r=r+v*dt+0.5d0*dt**2*F
  v=v+0.5d0*dt*F
  call periodic_boundary_conditions
  call force
  v=v+0.5d0*dt*F

  ! Velocity scaling thermostat
  temperature=2.d0*K_E/(3.d0*N_particles)
  v=v*sqrt(T_0/temperature)
END SUBROUTINE velocity_scaling

SUBROUTINE stats(k)
  integer :: i,k
  real(8)  :: su,rn

  ! Calculate the kinetic energy
  K_E=0.d0
  do i=1,N_particles
    K_E=K_E+0.5d0*(v(1,i)**2+v(2,i)**2+v(3,i)**2)
  enddo

  ! Calculate the mean square displacement
  if(k==N_term) then
    r0=r
  else if(k>N_term) then
    su=0.d0
    do i=1,N_particles
      su=su+(r(1,i)-r0(1,i))**2+(r(2,i)-r0(2,i))**2+(r(3,i)-r0(3,i))**2
    enddo
    write(msd_file,*)k-N_term,su/N_particles
    call radial_distribution
  endif

  total_energy=K_E+P_E
  write(6,'(i10,5E16.8)')k,K_E/dfloat(N_particles),total_energy/dfloat(N_particles), &
  2.d0*K_E/3.d0/dfloat(N_particles),T_av/dfloat(k),pressure
  write(stat_file,'(i10,5E16.8)')k,K_E/dfloat(N_particles),total_energy/dfloat(N_particles), &
  2.d0*K_E/3.d0/dfloat(N_particles),T_av/dfloat(k),pressure
END SUBROUTINE stats

SUBROUTINE radial_distribution
  integer :: i,j,k
  real(8)  :: del(3),dd

  N_s=N_s+1
  do i=1,N_particles-1
    do j=i+1,N_particles
      del(:)=r(:,i)-r(:,j)
      del(:)=del(:)-L(:)*nint(del(:)/L(:))
      dd=sqrt(sum(del**2))
      k=int(dd/delta)
      if(k<N_histogram) rdf(k)=rdf(k)+2
    enddo
  enddo
END SUBROUTINE radial_distribution

SUBROUTINE output_xyz_coordinates
  integer       :: i
  character*255 :: temp

  write(temp,*)N_particles
  write(movie_file,'(a)')trim(adjustl(temp))
  write(movie_file,'(a)')"Atoms"
  ! Convert from reduced units
  do i=1,N_particles
    write(movie_file,'(a,(3ES16.8))')"Si",r(1,i)*lj_sigma,r(2,i)*lj_sigma,r(3,i)*lj_sigma
  enddo
END SUBROUTINE output_xyz_coordinates

END MODULE m_molecular_dynamics

PROGRAM md
  USE m_molecular_dynamics
  implicit none
  integer          :: i,k,num_args
  real(8)           :: del(3), dd, shell_volume, np

  ! Process command-line arguments
  num_args = iargc()
  if(num_args/=1) then
    write(*,*)
    write(*,*) "Usage: md <input_filename>"; write(6,*)
    stop
  else
    call getarg(1,input_filename)
  endif

  call initialize()

  ! Perform the simulation
  T_av = 0.d0
  call force()
  call output_xyz_coordinates()
  do i=1,md_steps
    if(thermostat_type==1) then
      call velocity_scaling()
    elseif(thermostat_type==2) then
      call andersen_simulation()
    else
      call nose_hoover_simulation()
    endif
    call stats(i)
    T_av = T_av + 2.d0*K_E/(3.d0*N_particles)
    if(mod(i,output_frequency)==0) call output_xyz_coordinates()
  enddo

  ! Calculate the RDF
  do k = 0,N_histogram
    shell_volume = (4.d0/3.d0)*pi*((k+1)**3 - k**3)*delta**3
    np = shell_volume*rho
    rdf(k) = rdf(k)/(N_s*N_particles*Np)
    write(rdf_file,*) (k+0.5d0)*delta, rdf(k)
  enddo

  call output()  ! Output final positions and velocities to a file
  call output_xyz_coordinates()
  call cleanup() ! Deallocate memory
END PROGRAM

