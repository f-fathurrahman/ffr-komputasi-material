SUBROUTINE ordering(n, order, indexing)
  IMPLICIT NONE 
  INTEGER  :: j, i, n, indexing(n), ir, l, indxt
  real*8  :: order(n), q

      do  j=1,n
        indexing(j)=j
      ENDDO 

      l=n/2+1
      ir=n

10    continue
        if (l.gt.1)then
         l=l-1
         indxt=indexing(l)
         q=order(indxt)
        else
         indxt=indexing(ir)
         q=order(indxt)
         indexing(ir)=indexing(1)
         ir=ir-1
         if (ir.eq.1)then
           indexing(1)=indxt
           goto 30
         endif
        endif
        i=l
        j=l+l
20      if (j.le.ir)then
         if (j.lt.ir)then
           if (order(indexing(j)).lt.order(indexing(j+1)))j=j+1
         endif
         if (q.lt.order(indexing(j)))then
           indexing(i)=indexing(j)
           i=j
           j=j+j
         else
           j=ir+1
         endif
        go to 20
        endif
        indexing(i)=indxt
      go to 10
 30   continue
end SUBROUTINE ordering     
