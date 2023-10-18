module m_output
  use m_highlevel
  use m_midlevel
  use m_orbital_Li_HF, only:SLAP
  use m_jastrow, only:CJAS
  use m_observables
  implicit none
  public  ::  INITOUT,OUTWRITE,OUTLOG
      contains


subroutine INITOUT
  open(unit=31,file='AVERAGE_ENERGY.DAT',status='unknown', position='append')
! open(unit=35,file='RAD_DENSITY_EL1.DAT',status='unknown')
! open(unit=36,file='RAD_DENSITY_EL2.DAT',status='unknown')
! open(unit=37,file='RAD_DENSITY_EL3.DAT',status='unknown')
end subroutine INITOUT


subroutine OUTWRITE
  write(31,'(3f7.3,3f12.6)') SLAP(1), SLAP(2), CJAS, AVTOTALL, VARTOTALL
  ! OPTVAR (defined in vars Limc opt)
end subroutine OUTWRITE



subroutine OUTLOG
  write(*,*) 'STEPMAX = ',STEPMAX
  write(*,*) 'run: MCMAX-MCPRE= ',MCMAX-MCPRE
  write(*,*) 'acc. ratio = ', 100._dp*DBLE(MCOUNT)/DBLE(NE*(MCMAX-MCPRE+1)),' %  '
  write(*,*) 'AVELEN=',AVELEN
  write(*,*) 'VARELEN=',VARELEN
  write(*,*) 'AVTOTALL=',AVTOTALL
  write(*,*) 'VARTOTALL=',VARTOTALL
  write(*,*) 'AVBLOCKEL=',AVBLOCKEL
  write(*,*) 'VARBLOCKEL=',VARBLOCKEL
  write(*,*) 'AVKINEL=',AVKINEL
  write(*,*) 'VARKINEL=',VARKINEL
  write(*,*) 'AVVELEL=',AVVELEL
  write(*,*) 'VARVELEL=',VARVELEL
  write(*,*) 'AVPOTEL=',AVPOTEL
  write(*,*) 'VARPOTEL=',VARPOTEL
  write(*,*) 'AVINTEL=',AVINTEL
  write(*,*) 'VARINTEL=',VARINTEL
  write(*,*)
  write(*,*)
end subroutine OUTLOG

end module

