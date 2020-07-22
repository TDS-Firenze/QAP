
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

! =====================================================================  
! 2-optimum: First Improvement
SUBROUTINE opt2first(p,zp,m,zm) 
USE datos
!.. Scalar Arguments ..   
DOUBLE PRECISION, INTENT(IN):: zp
!.. Array Arguments ..   
INTEGER, INTENT(IN):: p(n)
!.. Local scalars ..
INTEGER i1, i2 ,aux 
DOUBLE PRECISION temp
!.. Outputs ..
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm   
!.. External functions ..
DOUBLE PRECISION, EXTERNAL:: delta
  
m = p 
zm = zp
50 DO i1 = 1, n-1 
     DO i2 = i1+1, n             
       temp = delta(m,i1,i2)
       IF (temp < 0.d0 ) THEN 
         zm    = zm + temp
         aux   = m(i1)
         m(i1) = m(i2)
         m(i2) = aux 
         GOTO 50 
       END IF
     END DO 
   END DO
END SUBROUTINE opt2first  

! =====================================================================  
! 2-optimum: Best Improvement  
SUBROUTINE opt2best(p,zp,m,zm) ! mayor mejora 2-optimo 
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
   DO i2 = i1+1,n 
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
! Delta^1
  
DOUBLE PRECISION FUNCTION delta1(p,i1,i2,i3)  ! for (i1,i2,i3) --> (i2,i3,i1)
USE datos 
!.. Scalar arguments ..
INTEGER, INTENT(IN) :: i1,i2,i3
!.. Array arguments ..
INTEGER, INTENT(IN) :: p(n)
!.. Local scalars ..  
INTEGER :: j 

 
IF (sim==0) THEN ! si el problema no es simetrico 
   DELTA1=A(i1,i1)*(B(p(i2),p(i2))-B(p(i1),p(i1)))+A(i1,i2)*(B(p(i2),p(i3))-B(p(i1),p(i2)))+&
          A(i1,i3)*(B(p(i2),p(i1))-B(p(i1),p(i3)))+A(i2,i1)*(B(p(i3),p(i2))-B(p(i2),p(i1)))+&
          A(i2,i2)*(B(p(i3),p(i3))-B(p(i2),p(i2)))+A(i2,i3)*(B(p(i3),p(i1))-B(p(i2),p(i3)))+ &
          A(i3,i1)*(B(p(i1),p(i2))-B(p(i3),p(i1)))+A(i3,i2)*(B(p(i1),p(i3))-B(p(i3),p(i2)))+ &
          A(i3,i3)*(B(p(i1),p(i1))-B(p(i3),p(i3)))
          
    DO j=1,n   
        IF((j .NE. i1) .AND. (j .NE. i2). AND. (j.NE.i3) ) THEN
          DELTA1=DELTA1+  A(i1,j)*(B(p(i2),p(j))-B(p(i1),p(j)))+A(i2,j)*(B(p(i3),p(j))-B(p(i2),p(j)))+&
                          A(i3,j)*(B(p(i1),p(j))-B(p(i3),p(j)))+A(j,i1)*(B(p(j),p(i2))-B(p(j),p(i1)))+&
                          A(j,i2)*(B(p(j),p(i3))-B(p(j),p(i2)))+A(j,i3)*(B(p(j),p(i1))-B(p(j),p(i3)))
        ENDIF
    ENDDO   
ELSEIF (sim==1) THEN ! si el problema es simetrico 
    DELTA1=A(i1,i2)*(B(p(i2),p(i3))-B(p(i1),p(i2)))+A(i1,i3)*(B(p(i1),p(i2))-B(p(i1),p(i3)))+&
           A(i2,i3)*(B(p(i1),p(i3))-B(p(i2),p(i3)))
    DO j=1,n  
          IF((j .NE. i1) .AND. (j .NE. i2). AND. (j.NE.i3) ) THEN
             DELTA1=DELTA1+  A(j,i1)*(B(p(j),p(i2))-B(p(j),p(i1)))+A(j,i2)*(B(p(j),p(i3))-B(p(j),p(i2)))+&
                             A(j,i3)*(B(p(j),p(i1))-B(p(j),p(i3)))
           ENDIF
    ENDDO  
    DELTA1= DELTA1*2
ENDIF      

ENDFUNCTION  DELTA1 

! =====================================================================  
! Delta^2 
DOUBLE PRECISION FUNCTION delta2(p,i1,i2,i3) ! for (i1,i2,i3) --> (i3,i1,i2) 
USE datos       
!.. Scalar arguments ..
INTEGER, INTENT(IN) :: i1,i2,i3
!.. Array arguments ..
INTEGER, INTENT(IN) :: p(n)
!.. Local scalars ..  
INTEGER :: j  

 
IF (sim == 0) THEN ! si el problema no es simetrico 
   delta2 = A(i1,i1)*(B(p(i3),p(i3))-B(p(i1),p(i1))) + A(i1,i2)*(B(p(i3),p(i1))-B(p(i1),p(i2))) &
          + A(i1,i3)*(B(p(i3),p(i2))-B(p(i1),p(i3))) + A(i2,i1)*(B(p(i1),p(i3))-B(p(i2),p(i1))) &
          + A(i2,i2)*(B(p(i1),p(i1))-B(p(i2),p(i2))) + A(i2,i3)*(B(p(i1),p(i2))-B(p(i2),p(i3))) &
          + A(i3,i1)*(B(p(i2),p(i3))-B(p(i3),p(i1))) + A(i3,i2)*(B(p(i2),p(i1))-B(p(i3),p(i2))) &
          + A(i3,i3)*(B(p(i2),p(i2))-B(p(i3),p(i3)))          
    DO j = 1,n   
        IF((j /= i1) .AND. (j /= i2). AND. (j /= i3)) THEN
          delta2 = delta2 + A(i1,j)*(B(p(i3),p(j))-B(p(i1),p(j))) + A(i2,j)*(B(p(i1),p(j))-B(p(i2),p(j)))&
                          + A(i3,j)*(B(p(i2),p(j))-B(p(i3),p(j))) + A(j,i1)*(B(p(j),p(i3))-B(p(j),p(i1)))&
                          + A(j,i2)*(B(p(j),p(i1))-B(p(j),p(i2))) + A(j,i3)*(B(p(j),p(i2))-B(p(j),p(i3)))
        END IF
    END DO   
ELSE IF (sim == 1) THEN ! if problem is symmetric 
    delta2 = A(i1,i2)*(B(p(i1),p(i3))-B(p(i1),p(i2))) + A(i1,i3)*(B(p(i2),p(i3))-B(p(i1),p(i3)))&
           + A(i2,i3)*(B(p(i1),p(i2))-B(p(i2),p(i3)))
    DO j = 1,n  
          IF((j /= i1) .AND. (j /= i2). AND. (j /= i3)) THEN
             delta2 = delta2+  A(j,i1)*(B(p(j),p(i3))-B(p(j),p(i1))) + A(j,i2)*(B(p(j),p(i1))-B(p(j),p(i2)))&
                            +  A(j,i3)*(B(p(j),p(i2))-B(p(j),p(i3)))
           ENDIF
    ENDDO  
    delta2 = delta2*2
END IF      
END FUNCTION  delta2


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
            GOTO 60 
          END IF
           temp = delta2(m,i1,i2,i3)
           IF (temp < 0.d0) THEN
             aux   = m(i1)
             m(i1) = m(i3)
             m(i3) = m(i2)
             m(i2) = aux        
             zm    = zm + temp
             GOTO 60
           END IF
      END DO
   END DO
END DO 
END SUBROUTINE  opt3first  
! =====================================================================  
! 3-optimum: Best Improvement 
SUBROUTINE opt3best(p,zp,m,zm)  ! mayor mejora 3-optimo
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
         temp=DELTA2(m,i1,i2,i3)
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

IF (dmin <0.d0) THEN
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
