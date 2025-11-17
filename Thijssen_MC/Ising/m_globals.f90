module m_globals

IMPLICIT NONE
INTEGER :: Size, MCSteps, Magnetisation
integer, PARAMETER :: MaxSize=256
INTEGER :: Spins(0:MaxSize, 0:MaxSize)

REAl(8) :: J, H, Energy, MCWeights(0:4)

end module
