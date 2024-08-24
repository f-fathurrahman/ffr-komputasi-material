MODULE GVECTOR

  IMPLICIT NONE 
  ! Lattice dimension
  INTEGER :: N_L(3)
  
  ! cut off in G_space
  REAL(8) :: E_cut,G_cut
  
  ! max number of G_vectors
  INTEGER :: N_G_vector_max
  
  ! G_vectors (must be multiplied by the B_Lattice vector
  INTEGER, ALLOCATABLE :: G_vector(:,:)
  REAL(8), ALLOCATABLE :: G_vector_length(:)
  
  ! K points
  INTEGER :: N_K_points
  REAL(8),ALLOCATABLE :: K_point(:,:),W_K_point(:)
  
  ! number of G vectors for a given K
  INTEGER, ALLOCATABLE :: N_G_vector(:)
  
  ! max number of G_vectors for a given K
  INTEGER :: N_G_K_vector_max
  
  ! index pointing to the position of the G vector 
  INTEGER, ALLOCATABLE :: G_index(:,:)
  
  ! (G+k)**2
  REAL(8), ALLOCATABLE :: Gplusk(:,:)
  
  ! fft index
  INTEGER, ALLOCATABLE :: fft_ind(:,:),ind3f1(:,:),ind3f2(:,:),ind2f1(:),ind2f2(:)    
  
  ! Lattice  vector
  REAL(8) :: Lattice_vector(3,3)
  REAL(8) :: Lattice_constant
  
  ! Lattice vector length            
  REAL(8) :: Lattice_vector_length(3)

  !
  ! Reciprocal space lattice
  !
  
  ! Lattice  vector
  REAL(8) :: R_Lattice_vector(3,3),Volume
  
  ! Lattice vector length            
  REAL(8) :: R_Lattice_vector_length(3)
  
  ! symmetry
  INTEGER :: N_sym,sym_mat(3,3,48)

  ! Units
  REAL(8) :: tpiba,tpiba2
  
  CONTAINS 


SUBROUTINE initialize_lattice()

  IMPLICIT NONE 
  INTEGER :: i,j
  REAL(8) :: g(3),a(3,3),b(3,3),pi
  
  pi = 4.d0*atan(1.d0)

  OPEN(1,file='lattice.inp')
  READ(1,*) Lattice_constant
  DO i=1,3
    read(1,*)(Lattice_vector(j,i),j=1,3)
  ENDDO 
  CLOSE(1)

  OPEN(1,file='kpoints.inp')
  READ(1,*) N_k_points
  ALLOCATE( K_point(3,N_K_Points), W_K_point(N_K_points) )
  DO i=1,N_k_points
    READ(1,*) ( K_point(j,i),j=1,3 ), W_K_point(i)
  ENDDO 
  CLOSE(1)
  
  !
  ! rescale lattice vectors
  !  Lattice_vector=Lattice_vector*Lattice_constant
  Lattice_vector=Lattice_vector
  
  ! length of lattice vectors
  DO i=1,3
    Lattice_vector_length(i) = SQRT(dot_product(Lattice_vector(:,i),Lattice_vector(:,i)))
  ENDDO 
  
  !
  ! reciprocal lattice
  !
  a = Lattice_vector
  CALL inv_r(a,3,b)
  R_lattice_vector = Transpose(b)*Lattice_constant
  
  ! length of lattice vectors
  DO i=1,3
    R_Lattice_vector_length(i) = SQRT(dot_product(R_Lattice_vector(:,i),R_Lattice_vector(:,i)))
  ENDDO 
  
  g(1) = Lattice_vector(2,2)*Lattice_vector(3,3) - Lattice_vector(3,2)*Lattice_vector(2,3)
  g(2) = Lattice_vector(3,2)*Lattice_vector(1,3) - Lattice_vector(1,2)*Lattice_vector(3,3)
  g(3) = Lattice_vector(1,2)*Lattice_vector(2,3) - Lattice_vector(2,2)*Lattice_vector(1,3)
  
  Volume = sum( Lattice_vector(:,1)*g(:) )

  N_L = 2*int( sqrt(e_cut)*Lattice_vector_length/Pi + 1 )
  
  WRITE(33,*) R_Lattice_vector
  WRITE(33,*) 'N_L',N_L
  G_cut = e_cut*Lattice_constant**2/(2.d0*Pi)**2
  WRITE(33,*) G_cut


  tpiba = 2.0*pi/Lattice_constant
  tpiba2 = tpiba**2


END SUBROUTINE 


SUBROUTINE initialize_symmetry()

  IMPLICIT NONE 
  INTEGER :: i,j,k,i1,i2,i3
  INTEGER :: s(3,3,48)
  
  OPEN(1,file='symmetry.inp')
  READ(1,*) N_Sym
  DO i=1,N_sym
    READ(1,*)
    DO k=1,3
      read(1,*) ( s(j,k,i), j=1,3 )
    ENDDO 
  ENDDO 
  CLOSE(1)
  
  DO i=1,N_sym
    DO i1=1,3  
      DO i2=1,3 
        if( mod(s(i1,i2,i)*N_L(i2),N_L(i1)) /= 0 ) write(6,*) 'wrong symmetry'
        sym_mat(i1,i2,i) = dble(s(i1,i2,i)*N_L(i2))/dble(N_L(i1))
      ENDDO 
    ENDDO 
  ENDDO 
  
  DO i=1,N_k_points
    W_K_point(i) = W_K_point(i)/dble(N_sym)
  ENDDO 

END SUBROUTINE 



SUBROUTINE Generate_G_vectors()
  IMPLICIT NONE 
  INTEGER :: i1,i2,i3,n1,n2,n3,ng,ii,i,j,ig(3),k
  INTEGER, ALLOCATABLE :: igv(:,:),inc(:)
  REAL(8), ALLOCATABLE :: gv(:)
  REAL(8) :: g(3), gg
  LOGICAL, ALLOCATABLE :: ic(:,:)
  LOGICAL :: con
  INTEGER, EXTERNAL :: iflip_inv,iflip
    
  !   two cycles: first count and allocate then store 
  DO ii = 1,2
  
    ! upper half: count the g vectors with g(1).ge.0
    ng = 1

    DO i3 = 0,N_L(3)/2
      
      n2 = -N_L(2)/2
      
      IF( i3 == 0) n2 = 0
      
      DO i2 = n2,N_L(2)/2
        
        n1 = -N_L(1)/2
        IF( i3==0 .AND. i2==0 ) n1 = 1
        
        DO i1=n1,N_L(1)/2
          
          g(:) = i1*R_Lattice_vector(:,1) + i2*R_Lattice_vector(:,2) + i3*R_Lattice_vector(:,3)
          
          gg = sum(g(:)**2)
          
          IF( gg <= 4.d0*G_cut ) THEN
            ng = ng + 1  
            IF( ii==2 ) THEN 
              gv(ng) = gg
              igv(1,ng) = i1
              igv(2,ng) = i2
              igv(3,ng) = i3
            ENDIF
          ENDIF
        ENDDO 
      ENDDO 
    ENDDO 
    ! upperhalf +lowerhalf+ g=0    
    N_G_vector_max = 2*ng-1
    ! now allocate arrays to store this
    IF( ii==1 ) THEN 
      ALLOCATE( gv(ng), igv(3,ng), inc(ng) )
      ! G=0 component
      gv(1) = 0.d0
      igv(:,1) = 0
    ENDIF
  ENDDO 
  !
  ! reorder with increasing magnitude
  !
  CALL ordering(ng, gv, inc)
  ! call indexx(ng,gv,inc)
  ALLOCATE( G_vector(3,N_G_vector_max), G_vector_length(N_G_vector_max) )
  ! G=0
  G_vector(:,1)=1
  G_vector_length(1) = 0.d0
  ii = 2
  DO i = 2,ng
    DO k = 1,3
      G_vector(k,ii+1) = iflip_inv(-igv(k,inc(i)),N_L(k))
      G_vector(k,ii) = iflip_inv(igv(k,inc(i)),N_L(k))
    ENDDO 
    G_vector_length(ii) = gv(inc(i))
    G_vector_length(ii+1) = gv(inc(i))
    ii = ii + 2
  ENDDO 
!
! G vectors for a given K point
!
  ALLOCATE( ic(N_L(1),N_L(2)) )
  ALLOCATE( N_G_vector(N_K_points),G_index(N_G_vector_max,N_K_points) )
  ALLOCATE( GplusK(N_G_vector_max,N_K_points),fft_ind(N_G_vector_max,N_K_points) )
  ALLOCATE( ind3f1(N_L(1),N_K_points),ind3f2(N_L(1),N_K_points) )
  ALLOCATE( ind2f1(N_K_points),ind2f2(N_K_points) )

  N_G_K_vector_max=0
  
  DO k=1,N_K_Points
    ic = .true.
    ii = 0
    DO i = 1,N_G_vector_max
      ig(1) = iflip(G_vector(1,i),N_L(1))
      ig(2) = iflip(G_vector(2,i),N_L(2))
      ig(3) = iflip(G_vector(3,i),N_L(3))
      g(:) = ig(1)*R_Lattice_vector(:,1) + ig(2)*R_Lattice_vector(:,2) + ig(3)*R_Lattice_vector(:,3) + K_point(:,k)
      gg = sum(g(:)**2)
      IF(gg <= G_cut) THEN 
        ii = ii + 1
        G_index(ii,k) = i 
        Gplusk(ii,k) = gg
        fft_ind(ii,k) = G_vector(1,i) + (N_L(1)+1)*((G_vector(2,i)-1) + N_L(2)*(G_vector(3,i)-1))
        ic(G_vector(1,i),G_vector(2,i)) = .FALSE.
      ENDIF
    ENDDO 
    N_G_vector(k) = ii
    IF( N_G_vector(k) > N_G_K_vector_max) N_G_K_vector_max=N_G_vector(k)
    
    !
    ! index limits for the FFT for the wf and potential
    !
    ind2f1(k) = 1
    ind2f2(k) = N_L(1) + 1
    DO i = 1,N_L(1)
      con = .TRUE.
      ind3f1(i,k) = 0
      ind3f2(i,k) = N_L(2) + 1
      DO j = N_L(2)/2,1,-1
        IF( .NOT. ic(i,j) .AND. con) THEN 
          ind3f1(i,k) = j
          con = .FALSE.
        ENDIF
      ENDDO 
      con = .TRUE.
      DO j=N_L(2)/2,N_L(2)
        IF( .NOT. ic(i,j) .AND. con) THEN 
          ind3f2(i,k) = j
          con = .FALSE.
        ENDIF
      ENDDO 
      IF( (ind3f1(i,k) /= 0) .AND. (ind3f2(i,k) /= N_L(2)+1) .AND. ( ind2f2(k) == N_L(1) + 1) ) ind2f1(k)=i+1
      IF( (ind3f1(i,k) == 0) .AND. (ind3f2(i,k) == N_L(2)+1) ) ind2f2(k) = i + 1
    ENDDO 
    
  ENDDO 
  
END SUBROUTINE 


END MODULE GVECTOR
