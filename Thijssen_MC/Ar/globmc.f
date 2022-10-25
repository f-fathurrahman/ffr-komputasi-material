      IMPLICIT NONE

      INTEGER PartNum, MaxCorr

      PARAMETER (PartNum = 108, MaxCorr=1000)

      DOUBLE PRECISION Qx(PartNum), Qy(PartNum), Qz(PartNum), 
     .       Temperature, Volume, Beta, Potential, Virial, 
     .       VolSize, Density, XVolSize, MaxRange, CorrStep,
     .       ROuter, RInner, Sep, Displace

      INTEGER CorrArray(MaxCorr), PairList(PartNum, PartNum),
     .        DispInt, InitStep, SimStep, ScaleStep,
     .        NeighNum(PartNum), UpdateStep, OutSteps


      COMMON /MainCommon/ Temperature, Volume, Qx, Qy, Qz, Displace, 
     .       VolSize, Density, XVolSize, MaxRange, CorrStep, Beta, 
     .       ROuter, RInner, Sep, CorrArray, Potential, 
     .       Virial, NeighNum, OutSteps, UpdateStep,
     .       PairList, DispInt, SimStep, InitStep, ScaleStep


