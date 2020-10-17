! =====================================================================  
! Objective function value
SUBROUTINE  funobj(p,z) 
USE datos 
!.. Array Arguments ..
INTEGER, INTENT(IN)::p(n)
!..Outputs..
DOUBLE PRECISION, INTENT(OUT)::z
!.. Local Scalars ..
INTEGER i,j 

z = 0.d0 
DO i = 1,n
   DO j = 1,n  
     z = z + A(i,j) * B(p(i),p(j)) 
   END DO
END DO
END SUBROUTINE funobj 
