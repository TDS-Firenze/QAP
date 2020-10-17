! =====================================================================  
! 2-optimum: Best Improvement  
SUBROUTINE opt2best(p,zp,m,zm) 
USE datos
!.. Scalar Arguments ..   
DOUBLE PRECISION, INTENT(IN):: zp
!.. Array Arguments ..   
INTEGER, INTENT(IN):: p(n)
!.. Local Scalars ..
INTEGER i1, i2 , j1,j2, aux
DOUBLE PRECISION:: temp, dmin
!.. Output ..
INTEGER, INTENT(OUT):: m(n) 
DOUBLE PRECISION, INTENT(OUT)::zm
!.. External functions ..
DOUBLE PRECISION, EXTERNAL:: DELTA

m  = p 
zm = zp 
55 dmin = 0.d0 
DO i1 = 1,n-1
   DO i2 = i1+1, n 
      temp = delta(m,i1,i2)
      IF(temp < dmin)  THEN
        j1   = i1
        j2   = i2
        dmin = temp
      END IF
   END DO
END DO
IF (dmin < 0.d0) THEN 
   aux   = m(j1)
   m(j1) = m(j2)
   m(j2) = aux  
   zm    = zm + dmin 
   GOTO 55
ENDIF
ENDSUBROUTINE opt2best
                                               
! =====================================================================  
! 3-optimum: First Improvement  
SUBROUTINE opt3first(p,zp,m,zm) 
USE datos
!.. Scalar arguments ..
DOUBLE PRECISION, INTENT(IN):: zp
!.. Array arguments ..
INTEGER, INTENT(IN):: p(n)
!.. Local scalars ..
INTEGER i1, i2 ,i3,aux
DOUBLE PRECISION temp
!.. Outputs ..
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm  
!.. External functions..
DOUBLE PRECISION, EXTERNAL:: delta1, delta2

  
m  = p 
zm = zp
  
60 DO i1 = 1,n-2 
    DO i2 = i1+1,n-1
      DO i3 = i2+1,n 
           temp = delta1(m,i1,i2,i3)  
          IF (temp < 0.d0) THEN
            aux   = m(i1)
            m(i1) = m(i2) 
            m(i2) = m(i3)
            m(i3) = aux        
            zm    = zm + temp
            GO TO 60 
          END IF

           temp = delta2(m,i1,i2,i3)
           IF (temp < 0.d0) THEN
             aux   = m(i1)
             m(i1) = m(i3)
             m(i3) = m(i2)
             m(i2) = aux        
             zm    = zm + temp
             GO TO 60
           END IF
      END DO
   END DO
END DO 
END SUBROUTINE  opt3first  

! =====================================================================  
! 3-optimum: Best Improvement 
SUBROUTINE opt3best(p,zp,m,zm)  
USE datos                           
!.. Scalar arguments ..
DOUBLE PRECISION, INTENT(IN):: zp
!.. Array arguments ..
INTEGER, INTENT(IN):: p(n)
!.. Local scalars ..
INTEGER :: i1, i2 , i3, j1, j2, j3, l, aux
DOUBLE PRECISION :: temp, dmin
!.. Outputs ..
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm  
!.. External functions..
DOUBLE PRECISION, EXTERNAL :: delta1, delta2
 
m  = p 
zm = zp 
65 dmin = 0.d0 
DO i1 = 1,n-2
   DO i2 = i1+1,n-1
      DO i3 = i2+1,n
         temp = delta1(m,i1,i2,i3)
         IF (temp < dmin)  THEN
           j1 = i1
           j2 = i2 
           j3 = i3
           l = 1
           dmin = temp
         END IF  
         temp = DELTA2(m,i1,i2,i3)
         IF(temp < dmin)  THEN
           j1 = i1
           j2 = i2 
           j3 = i3
           l = 2
           dmin = temp
         END IF 
      END DO
   END DO
END DO

IF (dmin < 0.d0) THEN
   aux = m(j1)
   IF (l == 1) THEN
       m(j1) = m(j2)
       m(j2) = m(j3)
       m(j3) = aux        
   ELSE IF (l == 2) THEN 
       m(j1) = m(j3)
       m(j3) = m(j2)
       m(j2) = aux
   END IF 
   zm = zm + dmin 
   GOTO 65
END IF

END SUBROUTINE opt3best      
