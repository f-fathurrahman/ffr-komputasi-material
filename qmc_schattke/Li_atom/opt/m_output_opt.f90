C-----------------------------------------------------------------------
      module output
       use highlevel
       use midlevel
       use orbital,only:SLAP
       use jastrow, only:CJAS
       use observables
       implicit none
       public  ::  INITOUT,OUTWRITE,OUTLOG
      contains
C-----------------------------------------------------------------------
      subroutine INITOUT
       open(unit=31,file='AVERAGE_ENERGY.DAT',status='unknown',
     &  position='append')
C       open(unit=35,file='RAD_DENSITY_EL1.DAT',status='unknown')
C       open(unit=36,file='RAD_DENSITY_EL2.DAT',status='unknown')
C       open(unit=37,file='RAD_DENSITY_EL3.DAT',status='unknown')
      end subroutine INITOUT
C-----------------------------------------------------------------------
      subroutine OUTWRITE
       integer       :: j,jj
C       do jj=1,NE
C        do j=1,NRHO
C         write(34+jj,*) j*DRHO,AVRHORAD(j,jj)*(j*DRHO)**2
C        end do
C       end do
       write(31,'(3f7.3,4f12.6)')SLAP(1),SLAP(2),CJAS,AVTOTALL,EGUESS,
     &                           OPTVAR,MINOPT
      end subroutine OUTWRITE
C-----------------------------------------------------------------------
      subroutine OUTLOG
       write(*,*)'STEPMAX = ',STEPMAX
       write(*,*)'run: MCMAX-MCPRE= ',MCMAX-MCPRE,' acc. ratio = ',
     &           100._dp*DBLE(MCOUNT)/DBLE(NE*(MCMAX-MCPRE+1)),' %  '
       write(*,*)'AVELEN=',AVELEN
       write(*,*)'VARELEN=',VARELEN
       write(*,*)'AVTOTALL=',AVTOTALL
       write(*,*)'VARTOTALL=',VARTOTALL
       write(*,*)'AVBLOCKEL=',AVBLOCKEL
       write(*,*)'VARBLOCKEL=',VARBLOCKEL
       write(*,*)'AVKINEL=',AVKINEL
       write(*,*)'VARKINEL=',VARKINEL
       write(*,*)'AVVELEL=',AVVELEL
       write(*,*)'VARVELEL=',VARVELEL
       write(*,*)'AVPOTEL=',AVPOTEL
       write(*,*)'VARPOTEL=',VARPOTEL
       write(*,*)'AVINTEL=',AVINTEL
       write(*,*)'VARINTEL=',VARINTEL
       write(*,*)
       write(*,*)
C       close(35)
C       close(36)
C       close(37)
      end subroutine OUTLOG
C-----------------------------------------------------------------------
      end module output
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
