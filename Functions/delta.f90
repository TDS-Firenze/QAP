! =====================================================================  
! delta function 
DOUBLE PRECISION FUNCTION delta(p,i1,i2) !variacion de la funcion objectivo para (i1,i2)-->(i2,i1) 
USE datos 
!.. Scalar Arguments ..     
INTEGER, INTENT(IN):: i1,i2
!.. Array Arguments ..
INTEGER, INTENT(IN):: p(n)  
!.. Local Scalars ..
INTEGER j

delta = 0.d0
IF (sim == 0)  THEN  
   delta = (A(i1,i1)-A(i2,i2)) * (B(p(i2),p(i2))-B(p(i1),p(i1))) &
         + (A(i1,i2)-A(i2,i1)) * (B(p(i2),p(i1))-B(p(i1),p(i2))) 
   DO j = 1,n   
      IF((j /= i1) .AND. (j /= i2) ) THEN
         delta = delta + (A(i1,j)-A(i2,j)) * (B(p(i2),p(j))-B(p(i1),p(j)))& 
                       + (A(j,i1)-A(j,i2)) * (B(p(j),p(i2))-B(p(j),p(i1))) 
      END IF
   END DO         
ELSE IF (sim == 1) THEN  
  DO j = 1,n   
     IF((j /= i1) .AND. (j /= i2) ) THEN
       delta = delta + (A(j,i1)-A(j,i2)) * (B(p(j),p(i2))-B(p(j),p(i1)))  
     END IF
  END DO
  delta = delta * 2 
END IF
     
END FUNCTION  delta
