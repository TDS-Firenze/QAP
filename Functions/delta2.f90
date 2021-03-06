! =====================================================================  
! Delta^2 
DOUBLE PRECISION FUNCTION delta2(p,i1,i2,i3) ! for {i1,i2,i3} --> {i3,i1,i2} 
USE datos       
!.. Scalar arguments ..
INTEGER, INTENT(IN) :: i1,i2,i3
!.. Array arguments ..
INTEGER, INTENT(IN) :: p(n)
!.. Local scalars ..  
INTEGER :: j  

 
IF (sim == 0) THEN ! if problem is not symmetric 
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
    delta2 = 2.d0 * delta2
END IF      
END FUNCTION  delta2
