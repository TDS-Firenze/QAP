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
