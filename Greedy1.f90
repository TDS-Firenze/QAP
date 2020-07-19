!!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO GREEDY 1  ! 
!!!!!!!!!!!!!!!!!!!!!!!!
  
 SUBROUTINE greedy(q,n1,n2,vec,z_best) ! vec: mejor permutacion encontrada por el greedy; z: su valor de la funcion objectivo       
 USE datos                   
 INTEGER,INTENT(IN):: q     ! q=1 si quiero eligir n1 e n2; q=0 para hacer todas las combinaciones de n1 y n2 
 INTEGER, INTENT(IN)::n1,n2 ! para q=1 n1 negocio (planta) inicial, n2 lugar inicial 
 
 INTEGER i,j,k , r,s                     
 INTEGER tmpn(n),tmpl(n), aux(n),p(n)        
 DOUBLE PRECISION aa(n,n),v0,v1
 DOUBLE PRECISION zp    
 
 INTEGER, INTENT(OUT):: vec(n) !Permutacion final: por ejemplo vec(a)=b significa que al negocio (planta) a asocio lugar b 
 DOUBLE PRECISION, INTENT(OUT):: z_best
   
 z_best=1.d30
 
 DO k = 1,n     
  DO j = 1,n     
 
    tmpl = 0            
    aa = a      !aa = matriz plantas   
    tmpn(1) = k !vector temporal para negocios (plantas) asignados
    tmpl(1) = j ! vector temporal para lugares asignados   
    IF (q==1) THEN
      tmpn(1) = n1 ! vector temporal para negocios (plantas) asignados
      tmpl(1) = n2 ! vector temporal para luugares asigandos
    ENDIF
    aa(tmpn(1),tmpn(1)) = 1.d30       
    DO i=1,n-1  
       IF (MINVAL(aa(tmpn(i),:), DIM=1) < MINVAL(aa(:,tmpn(i)),DIM=1)) THEN                       
         tmpn(i+1) = MINLOC( aa( tmpn(i),:),DIM=1)  
       ELSE
         tmpn(i+1) = MINLOC( aa(:, tmpn(i)),DIM=1) 
       ENDIF
        
      v1 = 1.d30
      DO s=1,n       
         IF (MINVAL(ABS(tmpl-s)).NE. 0) THEN   ! es decir, si s no es un index de tmpl (para no asignar lugares ya asignados) 
         
         v0 = DOT_PRODUCT(a(tmpn(i+1),tmpn(1:i)) , b(s,tmpl(1:i)))+&
              DOT_PRODUCT(a(tmpn(1:i),tmpn(i+1)) , b(tmpl(1:i),s))  
          
         IF (v0 < v1) THEN
             v1=v0 
             r=s    
         ENDIF
        ENDIF  
      ENDDO      
      tmpl(i+1)=r  !nuevo lugar
      aa(tmpn(i+1),tmpn(1:i+1))=1.d30 ! actualizar a para no asignar negocios ya asignados
      aa(tmpn(1:i+1),tmpn(i+1))=1.d30
    ENDDO 
    CALL ISORT@(aux, tmpn,n) 
    p=tmpl(aux) !ordino      
   
    CALL FUNOBJ(p,zp) 
    IF (zp < z_best) THEN
       z_best = zp
       vec = p
    ENDIF
    IF (q==1) EXIT
  ENDDO
 ENDDO
 ENDSUBROUTINE greedy 
