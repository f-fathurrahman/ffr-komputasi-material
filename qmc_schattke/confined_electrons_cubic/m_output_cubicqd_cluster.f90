!-----------------------------------------------------------------------
module m_output
  use m_highlevel
  use m_midlevel
  use m_orbital,only:ALPHA
  use m_jastrow, only:CJAS
  use m_observables
  implicit none
  integer,public     :: ios
  public  ::  INITOUT,OUTWRITE,OUTLOG
  
  contains


subroutine INITOUT
!  Read variables if not an input file is present 
       open(unit=11,file='INPUT_CUBICQD_CLUSTER.DAT', iostat=ios,status='old')
       if (ios == 0) then 
        read(11,*) STEPMAX
        read(11,*) CJAS
        read(11,*) ALPHA
       else
        write(*,*)'STEPMAX=?'
        read(*,*) STEPMAX
        write(*,*)'STEPMAX = ',STEPMAX
        write (*,*)'CJAS=?'
        read (*,*) CJAS
        write (*,*)'CJAS = ',CJAS
        write (*,*) 'ALPHA=?'
        read (*,*) ALPHA
        write (*,*)'ALPHA = ',ALPHA
       end if
       close (11)
       open(unit=31,file='AVERAGE_ENERGY.DAT',status='unknown', position='append')
!       open(unit=35,file='RAD_DENSITY_EL1.DAT',status='unknown')
!       open(unit=36,file='RAD_DENSITY_EL2.DAT',status='unknown')
!       open(unit=37,file='RAD_DENSITY_EL3.DAT',status='unknown')
      end subroutine INITOUT


subroutine OUTWRITE
!       do jj=1,NE
!        do j=1,NRHO
!         write(34+jj,*) j*DRHO,AVRHORAD(j,jj)*(j*DRHO)**2
!        end do
!       end do
       write(31,'(3f7.3,2f12.6)') CJAS, ALPHA, AVTOTALL, VARTOTALL
end subroutine OUTWRITE


subroutine OUTLOG
       write(*,*)'Maximum step in QMC jump = ',STEPMAX
       write(*,*)'Number of Monte Carlo steps after thermalization= ', &
     &            MCMAX-MCPRE, &
     &           'Ratio of accepted steps = ', &
     &            100._dp*DBLE(MCOUNT)/DBLE(NE*(MCMAX-MCPRE+1)),' %  '
       write(*,*) ' '
       write(*,*)'Total energy averaged over electrons=',AVTOTALL
       write(*,*)'Variance=',VARTOTALL
       write(*,*) ' '
       write(*,*)'Total energy per electron=',AVELEN
       write(*,*)'Variance=',VARELEN
       write(*,*)'Sum=',sum(AVELEN(1:NE))
       write(*,*)'Normalized sum=',sum(AVELEN(1:NE))/dble(NE)
       write(*,*) ' '
       write(*,*)'Block average of total energy per electron=',AVBLOCKEL
       write(*,*)'Variance=',VARBLOCKEL
       write(*,*)'Sum=',sum(AVBLOCKEL(1:NE))
       write(*,*)'Normalized sum=',sum(AVBLOCKEL(1:NE))/dble(NE)
       write(*,*) ' '
       write(*,*)'Kinetic energy per electron=',AVKINEL
       write(*,*)'Variance=',VARKINEL
       write(*,*)'Sum=',sum(AVKINEL(1:NE))
       write(*,*)'Normalized sum=',sum(AVKINEL(1:NE))/dble(NE)
       write(*,*) ' '
       write(*,*)'Velocity form of the kinetic energy per electron=', AVVELEL
       write(*,*)'Variance=',VARVELEL
       write(*,*)'Sum=',sum(AVVELEL(1:NE))
       write(*,*)'Normalized sum=',sum(AVVELEL(1:NE))/dble(NE)
       write(*,*) ' '
       write(*,*)'Potential energy per electron=',AVPOTEL
       write(*,*)'Variance=',VARPOTEL
       write(*,*)'Sum=',sum(AVPOTEL(1:NE))
       write(*,*)'Normalized sum=',sum(AVPOTEL(1:NE))/dble(NE)
       write(*,*) ' '
       write(*,*)'Coulomb energy per electron=',AVINTEL
       write(*,*)'Variance=',VARINTEL
       write(*,*)'Sum=',sum(AVINTEL(1:NE))
       write(*,*)'Normalized sum=',sum(AVINTEL(1:NE))/dble(NE)
!
       close(31)
!       close(35)
!       close(36)
!       close(37)
      end subroutine OUTLOG
!-----------------------------------------------------------------------
      end module
