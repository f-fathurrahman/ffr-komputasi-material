complex(8) FUNCTION inner_prod(Arr1, Arr2)
  use globals
  IMPLICIT NONE

  complex(8) :: Arr1(0:GridSize-1,0:GridSize-1,0:GridSize-1), &
                    Arr2(0:GridSize-1,0:GridSize-1,0:GridSize-1)
  Inner_Prod = SUM(CONJG(Arr1)*Arr2)
END FUNCTION
