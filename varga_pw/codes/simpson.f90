      real*8  function simpson_integral(n,f,h) 
       IMPLICIT NONE 
       INTEGER  n, i
       real*8  su, f, h
       dimension f(n)

       su=0.d0
       do i=1,n
         su=su+f(i)
       ENDDO 
       
       simpson_integral=su*h
       end

