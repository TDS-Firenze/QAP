! =====================================================================  
 SUBROUTINE greedy(q,n1,n2,p,z_best)        
 USE datos       
!.. Scalar Arguments ..           
 INTEGER, INTENT(IN):: q,n1,n2
!.. Local Scalars .. 
 INTEGER i,j,k,r,s
 DOUBLE PRECISION v0,v1,zp 
!.. Local Arrays ..                  
 INTEGER f(n),l(n), aux(n),vec(n)        
 DOUBLE PRECISION aa(n,n)
!.. Outputs .. 
 INTEGER, INTENT(OUT):: p(n) ! Final permutation: p(a)=b means assigning location b to facility a 
 DOUBLE PRECISION, INTENT(OUT):: z_best

 z_best=1.d30 
 DO k = 1,n     
  DO j = 1,n      
    tmpl = 0            
    aa = a        
    f(1) = k 
    l(1) = j  
    IF (q==1) THEN
      f(1) = n1 
      l(1) = n2 
    ENDIF
    aa(f(1),f(1)) = 1.d30       
    DO i=1,n-1  
       IF (MINVAL(aa(f(i),:), DIM=1) < MINVAL(aa(:,f(i)),DIM=1)) THEN                       
         f(i+1) = MINLOC( aa( f(i),:),DIM=1)  
       ELSE
         f(i+1) = MINLOC( aa(:,f(i)),DIM=1) 
       ENDIF        
      v1 = 1.d30
      DO s=1,n       
!       check if s is an index of l (to avoid assign locations already assigned)     
        IF (MINVAL(ABS(l-s)).NE. 0) THEN       
          v0 = DOT_PRODUCT(a(f(i+1),f(1:i)) , b(s,l(1:i)))+&
               DOT_PRODUCT(a(f(1:i),f(i+1)) , b(l(1:i),s))            
          IF (v0 < v1) THEN
            v1 = v0 
            r = s    
          ENDIF
        ENDIF  
      ENDDO      
      l(i+1)=r 
      aa(f(i+1),f(1:i+1))=1.d30 
      aa(f(1:i+1),f(i+1))=1.d30
    ENDDO  
    CALL ISORT@(aux,f,n) 
    vec = l(aux)        
    CALL FUNOBJ(vec,zp) 
    IF (zp < z_best) THEN
       z_best = zp
       p = vec
    ENDIF
    IF (q==1) EXIT
  ENDDO
 ENDDO
 ENDSUBROUTINE greedy 
  
