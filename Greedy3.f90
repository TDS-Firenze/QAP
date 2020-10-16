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
