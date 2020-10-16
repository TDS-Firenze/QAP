MODULE datos
  INTEGER  n , sim ! sim = 1 if A and B are symmetric, sim = 0 otherwise
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:)     
  DOUBLE PRECISION, PARAMETER :: t_max = 1.d0 
  DOUBLE PRECISION :: BKS
ENDMODULE datos 

!=====================================================================
 PROGRAM QAP   
 USE datos   
 IMPLICIT NONE                         
!.. Local Scalars ..    
 DOUBLE PRECISION tcpu1, tcpu2!, BKS 
 CHARACTER(LEN=20) nomfile, nomfile2
 INTEGER m,i,j,useless,r,s,t,itOpt  
!.. External Functions ..
 DOUBLE PRECISION, EXTERNAL:: DELTA, DELTA1, DELTA2   
 INTEGER, EXTERNAL:: MODDISTANCIA 
!.. External Subroutines .. 
 EXTERNAL funobj
 EXTERNAL greedy, greedy2, greedy3  
 EXTERNAL opt2first, opt2best, opt3first, opt3best
 EXTERNAL TEMPLE

OPEN (10,FILE='todos.txt')  
OPEN (90,FILE='todos_sln.txt')
OPEN (12,file='solucion.sol')       

CALL CPU_TIME (tcpu1)
DO
  READ(10,'(A)',END=10) nomfile
  OPEN(11, FILE=nomfile)
  READ(11,*) sim, n                       
  WRITE(*,'(A)') nomfile        
!  WRITE(*,9998) 'Dimension = ', n        
  WRITE(12,'(/A,3X,I0)') nomfile, n   
  READ(90,'(A)',END=10) nomfile2   
  OPEN(91, FILE=nomfile2)  
  READ(91,*) useless, BKS        
!  WRITE(*,'(/A,3X,I10)')'Best known solution', INT(BKS)
  ALLOCATE(A(n,n),B(n,n))           
  DO i=1,n
    READ(11,*) A(i,:)
  ENDDO      
  DO i=1,n
    READ(11,*) B(i,:)
  ENDDO             
  CALL resolucion     
  DEALLOCATE(A,B)           
ENDDO  
10 CALL CPU_TIME (tcpu2)
WRITE (*,'(A,F10.3)') 'total time = ', tcpu2-tcpu1
9998 FORMAT(T5,(A),T20,i5)
CONTAINS                            

!=====================================================================
SUBROUTINE resolucion  
USE datos  
!.. Local Scalars .. 
 INTEGER itOpt, it  
 INTEGER  i,j ,ii, tabu_size, n1, n2, q, i1, i2 ,i3 ,i4
 DOUBLE PRECISION mu,k          
 DOUBLE PRECISION z_greedy1, z_greedy2, z_greedy3     
 DOUBLE PRECISION z,ztemp,z_random     
!.. Local Arrays ..              
DOUBLE PRECISION s(2) , q_ant, alpha_ant
INTEGER :: m_ant , S_ant , R_ant
INTEGER :: pop(n,m)
INTEGER x(n), p(n),  xstar(n)  , p_ant(n)         
INTEGER p_greedy1(n),p_greedy2(n),p_greedy3(n), p_random(n)    
DOUBLE PRECISION :: SOLUTIONS(4)
DOUBLE PRECISION  z_localsearch(4), z_ant
!Greedy
Double precision:: TEMPO(6)    

q = 0
CALL greedy3 (q , n1 , n2 , p_greedy3 , z_greedy3)  


itOpt = 1
SOLUTIONS(:) = 1.d30

!-------------!
! Tabu Search !
!-------------!  


DO it = 1 , itopt
  CALL TABU(p, z)        
  IF (z < SOLUTIONS(1)) SOLUTIONS(1) = z
ENDDO
 SOLUTIONS(1) =          (SOLUTIONS(1) - BKS) / BKS * 100.d0   


!-----!
! ACO ! 
!-----!    

DO it = 1 , itopt                                 
  CALL ACO(p , z_ant)   
  IF (z_ant < SOLUTIONS(2)) SOLUTIONS(2) = z_ant
END DO  
SOLUTIONS(2) = (SOLUTIONS(2) - BKS)*100d0/BKS                   


!------!
! GVNS ! 
!------!

DO it = 1 , itopt        
   CALL GVNSfirst(p_greedy3, z_greedy3, x, z)   
   IF (z < SOLUTIONS(3)) SOLUTIONS(3) = z
END DO
SOLUTIONS(3) = (SOLUTIONS(3) - BKS) / BKS * 100.d0   

DO it = 1 , itopt        
   CALL GVNSbest(p_greedy3, z_greedy3, x, z)   
   IF (z < SOLUTIONS(4)) SOLUTIONS(4) = z
END DO
   SOLUTIONS(4) =          (SOLUTIONS(4) - BKS) / BKS * 100.d0   



WRITE(*,1235) nomfile, SOLUTIONS



949  FORMAT(t4,A,1x,1(' & ',f6.2,1x),'\\')         
998  FORMAT(t4,44('-'))  
999  FORMAT(t6, 'Algorithm',t34,'Solution')
9999 FORMAT(T6,(A),T30,F13.2)
1235  FORMAT(A,1x,4('&'1x,f5.2,1x),'\\')
1234 FORMAT('\textt{'A'}',1x'&'1x,A,1x'&'1x,4(f4.2,1x'&'1x),'\\')

END SUBROUTINE resolucion    
END PROGRAM QAP 



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

! =====================================================================  
! Greedy1
  
 SUBROUTINE greedy1(q,n1,n2,p_best,z_best)        
 USE datos       
!.. Scalar Arguments ..           
 INTEGER, INTENT(IN):: q,n1,n2
!.. Local Scalars .. 
 INTEGER i,j,k,r,s
 DOUBLE PRECISION v0,v1 , z
!.. Local Arrays ..                  
 INTEGER f(n),l(n), aux(n),vec(n) , p(n)
 DOUBLE PRECISION aa(n,n)
!.. Outputs .. 
 INTEGER, INTENT(OUT):: p_best(n) 
 DOUBLE PRECISION, INTENT(OUT):: z_best

 z_best = 1.d30 
 DO k = 1 , n     
   DO j = 1 , n 
     aa(:,:) = a    
     l(:) = 0   
     f(:) = 0    
     IF (q == 1) THEN
       f(1) = n1 
       l(1) = n2 
     ELSE IF (q == 0) THEN
       f(1) = k 
       l(1) = j 
     END IF      
      p(f(1)) = l(1)      
     aa(f(1) , f(1)) = -1.d30       
     DO i = 1 , n - 1  
        IF (MAXVAL(aa(f(i) , :), DIM = 1) > MAXVAL(aa(: , f(i)), DIM = 1)) THEN                       
          f(i + 1) = MAXLOC(aa(f(i),:) , DIM = 1)  
        ELSE
          f(i + 1) = MAXLOC(aa(:,f(i)) , DIM = 1) 
        END IF        
        v1 = 1.d30
        DO s = 1,n       
 !        check if s is an index of l (to avoid assign locations already assigned)     
          IF ( MINVAL(ABS(l-s)) /= 0) THEN       
            v0 = DOT_PRODUCT(a(f(i+1),f(1:i)) , b(s,l(1:i)))&
              + DOT_PRODUCT(a(f(1:i),f(i+1)) , b(l(1:i),s))            
             IF (v0 < v1) THEN
              v1 = v0 
              r = s    
            END IF
          END IF  
        END DO      
       l(i+1) = r 
       aa(f(i+1), f(1:i+1)) = -1.d30 
       aa(f(1:i+1), f(i+1)) = - 1.d30
       p(f(i + 1)) = l(i + 1)
     END DO  
     CALL FUNOBJ(p , z)
     IF (q==1) THEN
       p_best = p
       z_best = z
       RETURN
     END IF
     IF (z < z_best)  THEN
       z_best = z
       p_best = p
     END IF
       
   END DO
 END DO

 END SUBROUTINE greedy1 
  
 
! =====================================================================  
! Greedy2
 SUBROUTINE greedy2(q,n1,n2,p_best,z_best)    
 USE datos    
!.. Scalar Arguments ..           
 INTEGER, INTENT(IN) :: q, n1, n2
!.. Local Scalars .. 
 INTEGER :: i , j, k, r, s  , p(n)
 DOUBLE PRECISION :: v0, v1, z_p 
!.. Local Arrays ..                  
 INTEGER :: f(n), l(n), aux(n), vec(n)        
 DOUBLE PRECISION :: bb(n,n)
!.. Outputs .. 
 INTEGER, INTENT(OUT) :: p_best(n) ! Final permutation: vec(a)=b means assigning location b to facility a 
 DOUBLE PRECISION, INTENT(OUT) :: z_best
 
 z_best = 1.d30
 DO k = 1,n
  DO j = 1,n
    bb(:,:) = b     
    f(:) = 0            
    f(1) = j 
    l(1) = k 
    IF (q == 1) THEN 
      f(1) = n1
      l(1) = n2 
    ENDIF 
    p(f(1)) = l(1)
    bb(l(1),l(1)) = 1.d30       
    DO i = 1,n-1  
       IF (MINVAL(bb(l(i),:), DIM=1) < MINVAL(bb(:,l(i)),DIM=1)) THEN                       
        l(i+1) = MINLOC(bb(l(i),:),DIM=1)  
       ELSE
        l(i+1) = MINLOC(bb(:,l(i)),DIM=1) 
       ENDIF
      v1 = -1.d30
      DO s = 1,n       
         IF (MINVAL(ABS(f-s)) /=  0) THEN       
           v0 = DOT_PRODUCT(b(l(i+1),l(1:i)) , a(s,f(1:i))) &
                + DOT_PRODUCT(b(l(1:i),l(i+1)) , a(f(1:i),s))           
           IF (v0 > v1) THEN
             v1 = v0 
             r = s    
           ENDIF
         ENDIF  
      ENDDO      
      f(i+1) = r               
      p(f(i + 1)) = l(i + 1)
      bb(l(i+1),l(1:i+1))=1.d30    ! actualizar a para no asignar negocios ya asignados
      bb(l(1:i+1),l(i+1))=1.d30
    ENDDO 
    CALL funobj(p , z_p)
    IF (q==1)THEN
      p_best = p
      z_best = z_p
      RETURN
    END IF
    IF (z_p < z_best) THEN
      z_best = z_p
      p_best = p
    END IF
      
  ENDDO
 ENDDO
 ENDSUBROUTINE greedy2 
                          
! =====================================================================  
! Greedy3      
 
SUBROUTINE greedy3(q,n1,n2,p_best,z_best)  
USE datos  
!.. Scalar Arguments ..           
 INTEGER, INTENT(IN):: q,n1,n2
!.. Local Scalars .. 
 INTEGER i,j,k,r,s , p(n)
 DOUBLE PRECISION v0,v1 , z_p
!.. Local Arrays ..                  
 INTEGER f(n),l(n), aux(n),vec(n)        
 DOUBLE PRECISION aa(n,n), bb(n,n)
!.. Outputs .. 
 INTEGER, INTENT(OUT):: p_best(n) ! Final permutation: vec(a)=b means assigning location b to facility a 
 DOUBLE PRECISION, INTENT(OUT) :: z_best
 

z_best = 1.d30 
p = 0
DO i = 1,n 
 DO j = 1,n 
   f(:) = 0
   l(:) = 0
   IF (q == 1)THEN
       f(1) = n1 
       l(1) = n2    
   ELSE IF (q == 0) THEN
       f(1) = i 
       l(1) = j
   ENDIF                 
   p(f(1)) = l(1)
   aa = a 
   aa(f(1),f(1)) = -1.d30
   bb = b           
   bb(l(1),l(1)) = 1.d30 
   DO s = 1,n-1     
     IF (sim == 1) THEN
        f(s+1) = MAXLOC( aa(f(s),:) , DIM = 1) 
        l(s+1) = MINLOC( bb(l(s),:) , DIM = 1)
     ELSE            
        IF (MAXVAL( aa(f(s),:) , DIM = 1) > MAXVAL( aa(:,f(s)) , DIM = 1)) THEN                       
           f(s+1) = MAXLOC( aa(f(s),:),DIM=1)  
        ELSE
           f(s+1) = MAXLOC( aa(:,f(s)),DIM=1) 
        ENDIF
        IF (MINVAL(bb(l(s),:), DIM=1) < MINVAL(bb(:,l(s)),DIM=1)) THEN                       
           l(s+1) = MINLOC( bb(l(s),:),DIM=1)  
        ELSE
           l(s+1) = MINLOC( bb(:,l(s)),DIM=1) 
        ENDIF
     ENDIF

     p(f(s + 1)) = l(s + 1)
     aa(f(1:s+1), f(s+1)) = -1.d30    
     aa(f(s+1), f(1:s+1)) = -1.d30         
     bb(l(s+1), l(1:s+1)) =  1.d30    
     bb(l(1:s+1), l(s+1)) =  1.d30    
     ENDDO    
     CALL funobj(p , z_p)

     IF (q == 1) THEN
       p_best = p
       z_best = z_p
       RETURN
     END IF
     IF (z_p < z_best) THEN
       z_best = z_p
       p_best = p
     END IF
          
 ENDDO 
ENDDO   

 
ENDSUBROUTINE greedy3


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
INTEGER i1, i2 ,aux, ne, nr
DOUBLE PRECISION temp
!.. Outputs ..
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm   
!.. External functions ..
DOUBLE PRECISION, EXTERNAL:: delta
  
m = p 
zm = zp  
ne = 0
nr = 0
50 DO i1 = 1, n-1 
     DO i2 = i1+1, n        
       ne = ne + 1     
       temp = delta(m,i1,i2)
       IF (temp < 0.d0 ) THEN    
         nr    = nr + 1
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
! Delta^1
  
DOUBLE PRECISION FUNCTION delta1(p,i1,i2,i3)  ! for {i1,i2,i3} --> {i2,i3,i1}
USE datos 
!.. Scalar arguments ..
INTEGER, INTENT(IN) :: i1,i2,i3
!.. Array arguments ..
INTEGER, INTENT(IN) :: p(n)
!.. Local scalars ..  
INTEGER :: j 

 
IF (sim == 0) THEN ! si el problema no es simetrico 
   delta1 = A(i1,i1)*(B(p(i2),p(i2))-B(p(i1),p(i1)))+A(i1,i2)*(B(p(i2),p(i3))-B(p(i1),p(i2)))+&
            A(i1,i3)*(B(p(i2),p(i1))-B(p(i1),p(i3)))+A(i2,i1)*(B(p(i3),p(i2))-B(p(i2),p(i1)))+&
            A(i2,i2)*(B(p(i3),p(i3))-B(p(i2),p(i2)))+A(i2,i3)*(B(p(i3),p(i1))-B(p(i2),p(i3)))+ &
            A(i3,i1)*(B(p(i1),p(i2))-B(p(i3),p(i1)))+A(i3,i2)*(B(p(i1),p(i3))-B(p(i3),p(i2)))+ &
            A(i3,i3)*(B(p(i1),p(i1))-B(p(i3),p(i3)))
          
    DO j=1,n   
        IF((j /= i1) .AND. (j /= i2). AND. (j /= i3) ) THEN
          delta1 = delta1 +  A(i1,j)*(B(p(i2),p(j))-B(p(i1),p(j)))+A(i2,j)*(B(p(i3),p(j))-B(p(i2),p(j)))+&
                             A(i3,j)*(B(p(i1),p(j))-B(p(i3),p(j)))+A(j,i1)*(B(p(j),p(i2))-B(p(j),p(i1)))+&
                             A(j,i2)*(B(p(j),p(i3))-B(p(j),p(i2)))+A(j,i3)*(B(p(j),p(i1))-B(p(j),p(i3)))
        ENDIF
    ENDDO   
ELSEIF (sim == 1) THEN ! si el problema es simetrico 
    delta1 = A(i1,i2)*(B(p(i2),p(i3))-B(p(i1),p(i2)))+A(i1,i3)*(B(p(i1),p(i2))-B(p(i1),p(i3)))+&
             A(i2,i3)*(B(p(i1),p(i3))-B(p(i2),p(i3)))
    DO j = 1,n  
          IF((j /= i1) .AND. (j /= i2). AND. (j /= i3) ) THEN
             delta1 = delta1 +  A(j,i1)*(B(p(j),p(i2))-B(p(j),p(i1)))+A(j,i2)*(B(p(j),p(i3))-B(p(j),p(i2)))+&
                                A(j,i3)*(B(p(j),p(i1))-B(p(j),p(i3)))
           ENDIF
    ENDDO  
    delta1 = 2.d0 * delta1
ENDIF      

ENDFUNCTION  DELTA1 

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
    delta2 = 2.d0 * delta2
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



! =====================================================================  
! Simulated Annealing  

SUBROUTINE TEMPLE(x,zx,xstar,z)
USE datos                         
!.. Scalar arguments ..          
DOUBLE PRECISION :: zx           
!.. Array arguments ..  
INTEGER :: x(n) 
!.. Local scalars ..
INTEGER :: aux, it , i1, i2
DOUBLE PRECISION ::  T , d, alfa, beta, t1, t2
!.. Local Arrays ..
DOUBLE PRECISION :: m(2)    
!.. Outputs ..
DOUBLE PRECISION, INTENT(OUT) :: z       
INTEGER, INTENT(OUT) :: xstar(n)  
!.. External functions..
DOUBLE PRECISION, EXTERNAL :: DELTA 



IF (n <= 49) THEN 
   T = 5.d3
   alfa = 9.d-1
   beta = 1.d7
ELSE 
   T = 5.d3
   alfa = 95.d-2
   beta = 1.d8
ENDIF    

CALL CPU_TIME(t1)
xstar = x 
z = zx   
CALL CPU_TIME(t2)        
DO WHILE (abs(t1-t2) < t_max) !criterio de parada: tiempo maximo
      it=0
    DO WHILE(it < beta/T .AND. abs(t1-t2)< t_max)
     CALL RANDOM_NUMBER(m) ! genera dos numeros reales aleatorios entre 0 y 1  
     i1= NINT((n-1)*m(1)+1)  !numeros enteros aleatorios entre 0 y n para 2-intercambios
     i2= NINT((n-1)*m(2)+1) 
     
     DO    
        IF (i2.NE.i1) EXIT  !para tener i1 y i2 distintos, si no no hay ningun intercambio
        CALL RANDOM_NUMBER(m(2))
        i2=NINT((n-1)*m(2)+1)    
     ENDDO
     
     d=DELTA(x,i1,i2)
     IF (d < 0.d0) THEN ! es decir si con el 2-intercambio mejoro x
        aux=x(i1)
        x(i1)=x(i2)
        x(i2)=aux 
        zx=zx+d
        IF (zx < z) THEN  
           xstar=x  
           z=zx  
        ENDIF 
        ELSE
        CALL RANDOM_NUMBER(m(1)) ! genero numero aleatorio entre 0 y 1
        
        IF (m(1)<EXP(-d/T)) THEN
            aux=x(i1)
            x(i1)=x(i2)
            x(i2)=aux 
            zx=zx+d 
        ENDIF
      ENDIF
      it=it+1
      CALL CPU_TIME(t2)
    ENDDO 
  T=alfa*T
  CALL CPU_TIME(t2)
ENDDO        
ENDSUBROUTINE TEMPLE            

! =====================================================================  
! Tabu Search 

SUBROUTINE TABU(p_best , z_best)
USE datos            
!.. Scalar arguments .. 
DOUBLE PRECISION  :: z_start  
!.. Array arguments ..
INTEGER:: p_start(n)
!.. Local Scalars ..   
INTEGER tabu_size, temp, aux, i1, i2, j1, j2, ii, h
DOUBLE PRECISION :: mu, k, zp, dmin, t1, t2 , r
!.. Local Arrays .. 
DOUBLE PRECISION Bmod(n,n)
INTEGER, DIMENSION(:), ALLOCATABLE :: tabu_list  ! tabu list        
INTEGER :: p(n), p_temp(n),  LTM(n,n)           
!.. Outputs ..
INTEGER :: p_best(n)
DOUBLE PRECISION, INTENT(OUT) :: z_best      
!.. External functions ..
EXTERNAL :: DELTA




IF (n < 50) THEN
 tabu_size = 20
  k = -1.8d0
ELSE 
  tabu_size = 20
  k = -1.4d0
END IF 

CALL greedy3(0 , 1 , 1 , p_start, z_start)

  mu = k * (1.0 * SUM(B)/n**2)      
      
ALLOCATE (tabu_list(tabu_size))
       
p_best = p_start   
p = p_start 
p_temp = p_start 
zp = z_start          
z_temp = zp
z_best = z_start

ii = 1            
LTM = 0
Bmod = B
tabu_list = 0    
CALL CPU_TIME(t1)
CALL CPU_TIME(t2)

DO WHILE (abs(t1-t2) < t_max)    
555 dmin = 0.d0        

! Local Search : seek the best 2-exchange starting from p         

loop1:   DO i1 = 1, n-1
 loop2 :  DO i2 = i1+1, n 
            temp = DELTA(p, i1, i2)  
            IF (temp < dmin) THEN                         
              h = search(i1, i2)    
              IF (h == 0 .OR. zp + temp < z_best)  THEN    
                p_temp = p
                j1 = i1
                j2 = i2    
                p_temp(j1) = p(j2)
                p_temp(j2) = p(j1)
                z_temp = zp + dmin   
                dmin = temp                 
              ENDIF
            ENDIF         
           ENDDO loop2  
         ENDDO loop1    
             
         ! If tabu list is full, start again
         IF (ii >= tabu_size) ii = 1      
         ! Update tabu list              
         tabu_list(ii) = p(j1)
         tabu_list(ii+1) = p(j2)
         ii = ii + 2    
         
         IF (dmin < 0.d0) THEN 
           p = p_temp 
           zp = z_temp    
           LTM(j1,j2) = LTM(j1,j2) + 1
           LTM(j2,j1) = LTM(j2,j1) + 1  
!           WRITE(12,*) NINT(zp), ','  

           IF (zp < z_best ) THEN
             z_best = zp       
             p_best = p      
           ENDIF  
           CALL CPU_TIME(t2)                          
           GOTO 555  
         ENDIF       
  Bmod = Bmod + mu * LTM
! CALL greedy3(0,1,2,p,zp)
  CALL greedymod(p,zp)
  CALL CPU_TIME(t2)               
ENDDO 

DEALLOCATE(tabu_list)
6969 FORMAT (12(i2,1x))

CONTAINS

  INTEGER FUNCTION search(i,j)
  ! this function looks for (i,j) in tabu_list;
  ! search =  1 if they are in the list, 0 if not.
  INTEGER i, j, it
  it = 1 
  search = 0       
  DO   
    IF ((tabu_list(it)==i .AND. tabu_list(it+1)==j) .OR. (tabu_list(it)==j .AND. tabu_list(it+1)==i) )  THEN
      search = 1 
      RETURN
    ENDIF
    it = it + 2   
    IF (it >= tabu_size) RETURN
  ENDDO      
  ENDFUNCTION search
  


 
 
  SUBROUTINE greedymod(p_best,z_best)   
!SUBROUTINE greedy3(q,n1,n2,p,z_p)  
USE datos  
!.. Scalar Arguments ..         
!.. Local Scalars .. 
 INTEGER i,j,k,r,s,p(n)
 DOUBLE PRECISION v0,v1
!.. Local Arrays ..                  
 INTEGER f(n),l(n), aux(n),vec(n)        
 DOUBLE PRECISION aa(n,n), bb(n,n),z
!.. Outputs .. 
 INTEGER, INTENT(OUT):: p_best(n) ! Final permutation: vec(a)=b means assigning location b to facility a 
 DOUBLE PRECISION, INTENT(OUT) :: z_best
 

z_best = 1.d30  
p = 0
DO i = 1,n 
 DO j = 1,n 
   f(:) = 0
   l(:) = 0
   f(1) = i 
   l(1) = j             
   p(f(1)) = l(1)
   aa = a 
   aa(f(1),f(1)) = -1.d30
   bb = b           
   bb(l(1),l(1)) = 1.d30 
   DO s = 1,n-1     
     IF (sim == 1) THEN
        f(s+1) = MAXLOC( aa(f(s),:) , DIM = 1) 
        l(s+1) = MINLOC( bb(l(s),:) , DIM = 1)
     ELSE            
        IF (MAXVAL( aa(f(s),:) , DIM = 1) > MAXVAL( aa(:,f(s)) , DIM = 1)) THEN                       
           f(s+1) = MAXLOC( aa(f(s),:),DIM=1)  
        ELSE
           f(s+1) = MAXLOC( aa(:,f(s)),DIM=1) 
        ENDIF
        IF (MINVAL(bb(l(s),:), DIM=1) < MINVAL(bb(:,l(s)),DIM=1)) THEN                       
           l(s+1) = MINLOC( bb(l(s),:),DIM=1)  
        ELSE
           l(s+1) = MINLOC( bb(:,l(s)),DIM=1) 
        ENDIF
     ENDIF
     p(f(s + 1)) = l(s + 1)
     aa(f(1:s+1), f(s+1)) = -1.d30    
     aa(f(s+1), f(1:s+1)) = -1.d30         
     bb(l(s+1), l(1:s+1)) =  1.d30    
     bb(l(1:s+1), l(s+1)) =  1.d30    
     ENDDO    
     CALL funobj(p , z)   
     IF (z < z_best) THEN
       z_best = z
       p_best = p
     END IF
 ENDDO 
ENDDO     
!ENDSUBROUTINE greedy3
  ENDSUBROUTINE greedymod      
  
  
ENDSUBROUTINE TABU                                                                                          
  
! =====================================================================  
! Ant Colony


SUBROUTINE ACO(p_best , z_best)                                 
USE datos      
!.. Parameters ..
INTEGER, PARAMETER :: m = 5 , R = 2                
INTEGER ::  S
DOUBLE PRECISION, PARAMETER :: alpha = 0.25d0 , q = 0.85d0                           
!.. Outputs ..
INTEGER, INTENT(OUT):: p_best(n)       ! final permutation
!DOUBLE PRECISION, INTENT(OUT):: z_best ! objective function value  
DOUBLE PRECISION :: z_best
!.. Local Arrays ..
INTEGER :: pi(n) , p_hat(n), p_tilde(n) , p_first(n) , p(n) , p_star(n)     
INTEGER :: Ants(n,m) , P_stock(n , m)
DOUBLE PRECISION :: prob(n)   
DOUBLE PRECISION :: T(n,n) , M_update(n,n) , M_out(n,n)
DOUBLE PRECISION :: z_stock1(m), z_stock2(m) 
DOUBLE PRECISION :: Tensor(n,n,m)
!.. Local Scalars..
INTEGER :: i , j , k ,  it , i_best
DOUBLE PRECISION :: z, z_hat, z_tilde, temp , u , t1 , t2    
LOGICAL :: intensification
!.. External..
EXTERNAL :: random_permutation, Delta

  
S = 5*n

CALL CPU_TIME(t1)
CALL CPU_TIME(t2)
                 
  M_update = 0.d0
  z_best = 1.d30
  Ants(: , :)  = 0
   


! Initialisation of Pheromone Trail 
CALL pheromone_initialisation(T)

it = 0
!DO 
DO 
  IF ( ABS(t1 - t2) > t_max) EXIT

   
     
   CALL CPU_TIME(t2)
   it = it + 1 

 DO k = 1 , m              
  ! Initialisation of solutions   
    CALL random_permutation(p)
    CALL funobj(p, z)    
    CALL opt2first(p , z , p , z)  
    IF (z < z_best) THEN
      z_best = z
      i_best = k 
      p_best(:) = p(:)
    ENDIF                       
    
  ! Manipulation of solutions
  DO i = 1 , R                             
    CALL apply_2exchange(p_hat , z_hat) 
  END DO    
    CALL pheromone_update(p_tilde , z_tilde , M_out) 
    ants(: , k) = p_tilde(:) 
    IF (z_tilde < z_best) THEN
      z_best = z_tilde
      p_best(:) = p_tilde(:)
      intensification = .TRUE.
      M_update(: , :) = M_out(: , :)
      i_best = k
    END IF   
    z_stock1(k) = z
    z_stock2(k) = z_tilde
    P_stock(: , k) = p(:)
  END DO 

  ! Intensification
   IF (intensification = .TRUE.) THEN 
     DO k = 1 , m
       IF (z_stock1(k) < z_stock2(k))    ants(: , k) = P_stock(: , k)
     END DO
   END IF
   intensification = .FALSE.
 
 ! Pheromone trail update
   T(: , :) = (1 - alpha) * T(: , :)
   T(: , :) = T(: , :) + M_update(: , :)

!print*, p_best 

 ! Diversification
   IF (it > S) THEN
     T(: , :) = 0d0
     CALL pheromone_initialisation(T)
     Ants(: , :) = 0
   END IF
   CALL CPU_TIME(t2) 
END DO  
      


 
   9000 FORMAT(4(I1,1x))  
   9001 FORMAT(4(f3.0,1x))  
   
   
CONTAINS


SUBROUTINE pheromone_initialisation(T)
INTEGER :: i, j , i1, i2, aux
INTEGER :: pi(n)
DOUBLE PRECISION :: z, temp
DOUBLE PRECISION :: M1(n,n) , M2(n,n)
DOUBLE PRECISION, INTENT(OUT) :: T(n,n)
 M1 = 0.d0
 M2 = 0.d0  
 DO j = 1 , m 
 ! Setting pi = (i, i+1, ..., n, 1, ...., i-1)
  DO i = 1,n
    pi(i) = MODULO(i + j - 1, n)     
  END DO
  pi(n - j + 1) = n
 ! Local search on pi 
  CALL funobj(pi,z)
  56 DO i1 = 1 , n - 1 
         DO i2 = i1 + 1 , n             
           temp = DELTA(pi , i1 , i2)
           IF (temp < 0.d0 ) THEN                 
             M2(i1 , pi(i1)) = M2(i1 , pi(i1)) - 1d0
             M2(i2 , pi(i2)) = M2(i2 , pi(i2)) - 1d0
             M2(i1 , pi(i2)) = M2(i1 , pi(i2)) + 1d0   
             M2(i2 , pi(i1)) = M2(i2 , pi(i1)) + 1d0
             aux = pi(i1)
             pi(i1) = pi(i2)
             pi(i2) = aux     
             z = z + temp
             GO TO 56 
           END IF
         END DO      
       END DO     
     DO i = 1 , n
      M1(i , pi(i)) = M1(i , pi(i)) + 1
     END DO                           
END DO    
 M2 = M2 + ABS(MINVAL(M2))
  T = M1 + M2  
END SUBROUTINE pheromone_initialisation


SUBROUTINE apply_2exchange(p_hat , z_hat) 
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_hat(n)
DOUBLE PRECISION, INTENT(OUT) :: z_hat
!.. Local scalars ..
INTEGER :: it,  r , s
DOUBLE PRECISION :: u , aux , sum_prob, sum_cumulative
    ! Choosing r
    CALL RANDOM_NUMBER(u)  
    r = NINT(u * (n-1) + 1 )       
    ! Choosing s     
    CALL RANDOM_NUMBER(u)
    IF ( u < q) THEN   
      aux = 0d0
      DO it = 1 , n
        IF ( T(r , p(it)) + T(it , p(r)) > aux .AND. it /= r) THEN
          aux = T(r , p(it)) + T(it , p(r))      
          s = it 
        END IF
      END DO    
    ELSE       
    ! Exploration case                                                     
      sum_prob = 0.d0
      DO i = 1 , n            
        IF (i /= r) sum_prob = sum_prob + T(r , p(i)) + T(i, p(r))
      END DO
      DO s = 1 , n        
        prob(s) = T(r , p(s)) + T(s , p(r))
        prob(s) = prob(s) / sum_prob
      END DO                               
      prob(r) = 0.d0
      CALL RANDOM_NUMBER(u)  
      i = 1     
      sum_cumulative  = 0d0

      DO 
        sum_cumulative = sum_cumulative + prob(i)
        IF (u < sum_cumulative ) THEN
          s = i
          EXIT
        END IF  
        i = i + 1
      END DO
   END IF
   p_hat(:) = p(:)
   p_hat(r) = p(s)  
   p_hat(s) = p(r)  
   z_hat = z + DELTA(p , r , s)

END SUBROUTINE apply_2exchange

SUBROUTINE pheromone_update(p_tilde , z_tilde , M_out)
USE datos
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_tilde(n)
DOUBLE PRECISION, INTENT(OUT) :: z_tilde , M_out(n,n)
!.. Local Scalars ..
INTEGER ::   i1, i2 , ii 
DOUBLE PRECISION :: z_tilde , temp 
M_out = 0d0
p_tilde = p_hat   
z_tilde = z_hat     
551 DO i1 = 1 , n - 1 
     DO i2 = i1 + 1 , n             
       temp = DELTA(p_tilde , i1 , i2)
         IF (temp < 0.d0 ) THEN      
            M_out(i1 , p_tilde(i1)) =  M_out(i1 , p_tilde(i1)) - 1d0
            M_out(i2 , p_tilde(i2)) =  M_out(i2 , p_tilde(i2)) - 1d0
            M_out(i1 , p_tilde(i2)) =  M_out(i1 , p_tilde(i2)) + 1d0   
            M_out(i2 , p_tilde(i1)) =  M_out(i2 , p_tilde(i1)) + 1d0                
           ii = p_tilde(i1)
           p_tilde(i1) = p_tilde(i2)
           p_tilde(i2) = ii     
           z_tilde  = z_tilde + temp
           GO TO 551 
         END IF
       END DO      
     END DO   
     M_out = M_out + ABS(MINVAL(M_out))
     DO i = 1 , n
       M_out(i , p_tilde(i)) = M_out(i , p_tilde(i)) + 1d0
     END DO 
END SUBROUTINE pheromone_update   

 
END SUBROUTINE ACO    

 

! =====================================================================  
! VND first improvement
SUBROUTINE VNDfirst(p_start,z_start,p_best,z_best) ! VND con algoritmos 2 y 3 optimos de primera mejora
USE datos

!.. Scalar Arguments ..
DOUBLE PRECISION, INTENT(IN) :: z_start 
!.. Array Arguments .. 
INTEGER,INTENT(IN) :: p_start(n)   
!.. Local Scalars ..
INTEGER :: j
DOUBLE PRECISION :: zx
INTEGER, PARAMETER :: ja = 2
!.. Local Arrays ..
INTEGER :: x(n)  
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_best(n)
DOUBLE PRECISION, INTENT(OUT):: z_best 
!.. External functions ..
EXTERNAL :: opt2first, opt3first 

p_best = p_start
z_best = z_start
j = 1
DO
 IF (j > ja) THEN 
   EXIT
 ELSE IF (j == 1) THEN 
   CALL opt2first(p_best, z_best, x, zx)
 ELSE IF (j == 2) THEN
   CALL opt3first(p_best, z_best, x, zx)
 END IF 
 IF (zx < z_best) THEN
     p_best = x 
     z_best = zx
     j = 1
 ELSE
     j = j+1
 END IF
END DO
END SUBROUTINE VNDfirst
                              
! =====================================================================  
! GVNS first improvement
SUBROUTINE GVNSfirst(p,zp,xstar,z) ! GVNS con algoritmos  de primera mejora
USE datos

!.. Scalar Arguments ..
DOUBLE PRECISION, INTENT(IN):: zp ! valor funcion objectivo inicial
!  Array Arguments .. 
INTEGER,INTENT(IN) :: p(n)    ! permutacion inicial
!.. Local Scalars ..
INTEGER, PARAMETER :: jb=3
INTEGER it, itmax ,j  , n1, n2, n3 
DOUBLE PRECISION zx , t1,t2  
!.. Local Arrays ..
INTEGER x(n)
DOUBLE PRECISION :: ran(3)
!.. Outputs ..
INTEGER, INTENT(OUT):: xstar(n)
DOUBLE PRECISION, INTENT(OUT):: z 
!.. External Functions ..
DOUBLE PRECISION, EXTERNAL:: DELTA, DELTA1
EXTERNAL ::  VNDPRIMERA,FUNOBJ


xstar=p
z=zp 
CALL CPU_TIME(t1) 
it=1 
itmax=10000 !numero di iteraciones maximo sin mejorar 
 
CALL CPU_TIME(t2)   
DO WHILE (it <= itmax .AND. ABS(t1 - t2)< t_max)
70 j = 1
 DO      
!   print*,'inizio ciclo' 
  IF (j > jb) THEN   
    EXIT           
  ELSEIF (j == 1) THEN  ! 2-intercambio en comparicion con la mejor solucion xstar         
!    print*,'caso 1'
    CALL RANDOM_NUMBER(ran)
    n1= NINT((n-1)*ran(1)+1)
    n2= NINT((n-1)*ran(2)+1)
        
    DO    
      IF (n2 /= n1) EXIT  !para tener n1 y n2 distintos, si no no hay ningun intercambio
      CALL RANDOM_NUMBER(ran(2))
      n2=NINT((n-1)*ran(2)+1)    
    ENDDO
    
 
   x = xstar
    x(n1) = xstar(n2)
    x(n2) = xstar(n1)
    zx = z + DELTA(xstar,n1,n2)
!    print*,'fine del IF'
  ELSEIF (j == 2) THEN  ! 3-intercambio en comparicion con la mejor solucion xstar    
!   print*,'caso 2'
    CALL RANDOM_NUMBER(ran)
    n1 = NINT((n-1)*ran(1)+1)
    n2 = NINT((n-1)*ran(2)+1)
       
    DO    
      IF (n2.NE.n1) EXIT  !para n1 y n2 distintos
      CALL RANDOM_NUMBER(ran(2))
      n2 = NINT((n-1)*ran(2)+1)    
    ENDDO
    
    n3 = NINT((n-1)*ran(3)+1) 
    
    DO
      IF (n3.NE.n1 .AND. n3.NE.n2) EXIT ! para tener n3 distintos de n1 y n2
       
     CALL RANDOM_NUMBER(ran(3))
     n3=NINT((n-1)*ran(3)+1)
     
    ENDDO    
    x = xstar
    x(n1) = xstar(n2)
    x(n2) = xstar(n3)
    x(n3) = xstar(n1) 
    zx = z+DELTA1(xstar,n1,n2,n3)
    
  ELSEIF (j==3) THEN  ! primera mmitad se pone segunda y al reves     
!   print*,'caso 3'
    
      x(1:CEILING(1.*n/2)) = xstar(n-CEILING(1.*n/2)+1:n)
      x(CEILING(1.*n/2)+1:n) = xstar(1 : n-CEILING(1.*n/2)) 
       
      CALL FUNOBJ(x,zx)
     
  ENDIF                     

    CALL VNDfirst(x,zx,x,zx)  
    IF (zx < z) THEN
       xstar = x
       z = zx 
       it = 0 ! aqui esta mejorando --> pongo it=0
       GOTO 70
    ENDIF
    j = j+1
 ENDDO
 it = it+1 
 CALL CPU_TIME(t2)     
ENDDO
ENDSUBROUTINE GVNSfirst



! =====================================================================  
! VNS best improvement
SUBROUTINE VNDbest(p_start,z_start,p_best,z_best) ! VND con algoritmos 2 y 3 optimos de primera mejora
USE datos

!.. Scalar Arguments ..
DOUBLE PRECISION, INTENT(IN) :: z_start 
!.. Array Arguments .. 
INTEGER,INTENT(IN) :: p_start(n)   
!.. Local Scalars ..
INTEGER :: j
DOUBLE PRECISION :: zx
INTEGER, PARAMETER :: ja = 2
!.. Local Arrays ..
INTEGER :: x(n)  
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_best(n)
DOUBLE PRECISION, INTENT(OUT):: z_best 
!.. External functions ..
EXTERNAL :: opt2best, opt3best

p_best = p_start
z_best = z_start
j = 1
   
DO
 IF (j > ja) THEN
   EXIT
 ELSE IF (j == 1) THEN 
   CALL opt2best(p_best, z_best, x, zx)
 ELSE IF (j == 2) THEN
   CALL opt3best(p_best, z_best, x, zx)
 END IF 
 IF (zx < z_best) THEN
     p_best = x 
     z_best = zx
     j = 1
 ELSE
     j = j+1
 END IF
END DO
END SUBROUTINE VNDbest                                              

! =====================================================================  
! GVNS best improvement
SUBROUTINE GVNSbest(p,zp,xstar,z) ! GVNS con algoritmos de mayor mejora
USE datos
 
INTEGER,INTENT(IN):: p(n)    ! permutacion inicial
DOUBLE PRECISION, INTENT(IN):: zp ! valor funcion objectivo inicial
INTEGER, INTENT(OUT):: xstar(n)
DOUBLE PRECISION, INTENT(OUT):: z 
INTEGER, PARAMETER:: jb=3
INTEGER it, itmax ,j  , x(n), n1, n2, n3 
DOUBLE PRECISION zx , ran(3),t1,t2
EXTERNAL VNDbest,FUNOBJ 
DOUBLE PRECISION, EXTERNAL:: DELTA, DELTA1


xstar = p
z=zp 
CALL CPU_TIME(t1)
it=1 
itmax=10000 !numero di iteraciones maximo sin mejorar 
 
CALL CPU_TIME(t2)
DO WHILE (it <= itmax .AND. ABS(t1-t2)<t_max) !tiempo maximo : tmax
75 j=1
 DO
  IF (j > jb) THEN 
   EXIT           
   
  ELSEIF (j == 1) THEN  ! 2-intercambio en comparicion con la mejor solucion xstar  
    CALL RANDOM_NUMBER(ran)
    n1 = NINT((n-1) * ran(1)+1)
    n2 = NINT((n-1) * ran(2)+1)
       
    DO    
      IF (n2 /= n1) EXIT  !!para tener n1 y n2 distintos, si no no hay ningun intercambio
      CALL RANDOM_NUMBER(ran(2))
      n2 = NINT((n-1) * ran(2) + 1)    
    ENDDO    
    
    zx = z + DELTA(xstar,n1,n2)
    x = xstar
    x(n1) = xstar(n2)
    x(n2) = xstar(n1)
    
  ELSEIF (j==2) THEN  ! 3-intercambio en comparicion con la mejor solucion xstar
    CALL RANDOM_NUMBER(ran)
    n1= NINT((n-1)*ran(1)+1)
    n2= NINT((n-1)*ran(2)+1)
    
    DO    
      IF (n2.NE.n1) EXIT  !!para tener n1 y n2 distintos, si no no hay ningun intercambio
      CALL RANDOM_NUMBER(ran(2))
      n2=NINT((n-1)*ran(2)+1)    
    ENDDO
    
    n3=NINT((n-1)*ran(3)+1) 
    
    DO
      IF (n3.NE.n1 .AND. n3.NE.n2) EXIT! para tener n3 distintos de n1 y n2
       
     CALL RANDOM_NUMBER(ran(3))
     n3=NINT((n-1)*ran(3)+1)
     
    ENDDO    
    x=xstar
    zx=z+DELTA1(xstar,n1,n2,n3)   
    x(n1)=xstar(n2)
    x(n2)=xstar(n3)
    x(n3)=xstar(n1)
    
  ELSEIF (j==3) THEN  ! primera mitad se pone segunda y al reves    
  
      x(1:CEILING(1.*n/2))=xstar(n-CEILING(1.*n/2)+1:n)
      x(CEILING(1.*n/2)+1:n)=xstar(1 : n-CEILING(1.*n/2))
      CALL FUNOBJ(x,zx)
  ENDIF  
    CALL VNDbest(x,zx,x,zx) 
    IF (zx < z) THEN
       xstar=x
       z=zx 
       it=0 ! aqui esta mejorando --> pongo it=0
       GOTO 75
    ENDIF
    j=j+1
 ENDDO
 it=it+1 
 CALL CPU_TIME(t2)   
ENDDO
ENDSUBROUTINE GVNSbest 


SUBROUTINE random_permutation(p)
USE datos                             
INTEGER, INTENT(OUT) :: p(n)
INTEGER :: i
  p = (/ (i, i = 1,n)/)
  CALL G05EHF(p,n,0) 
END SUBROUTINE  

