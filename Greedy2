 !!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO GREEDY  ! !version 2
!!!!!!!!!!!!!!!!!!!!!!
  
 SUBROUTINE greedy2(q,n1,n2,vec,z) ! vec: mejor permutacion encontrada por el greedy; z: su valor de la funcion objectivo        
 USE datos     
 INTEGER,INTENT(IN):: q ! q=1 si quiero eligir n1 e n2; q=0 para hacer todas las combinaciones de n1 y n2  
 INTEGER, INTENT(IN)::n1,n2 ! para q=1 n1 negocio (planta) inicial, n2 lugar inicial 
   
 INTEGER i,j,k, r,s                      
 INTEGER tmpn(n),tmpl(n), aux(n),p(n)              
 DOUBLE PRECISION bb(n,n),v0,v1
 DOUBLE PRECISION zp 
      
  INTEGER, INTENT(OUT):: vec(n)     !Permutacion final: por ejemplo vec(a)=b significa que al negocio (planta) a asocio lugar b 
 DOUBLE PRECISION, INTENT(OUT):: z
   
  
 z=1.d30
 
 DO k=1,n
  DO j=1,n
    tmpn=0            
    bb = b      !bb = matriz temporal lugares 
    tmpn(1) = j !vector temporal para negocios (plantas) asignados
    tmpl(1) = k ! vector temporal para lugares asignados 
    IF (q==1) THEN 
      tmpn(1) = n1 !vector temporal para negocios (plantas) asignados
      tmpl(1) = n2 ! vector temporal para lugares asignados 
    ENDIF 
  
    bb(tmpl(1),tmpl(1)) = 1.d30     
  

    DO i=1,n-1  
       IF (MINVAL(bb(tmpl(i),:), DIM=1) < MINVAL(bb(:,tmpl(i)),DIM=1)) THEN                       
        tmpl(i+1) = MINLOC( bb( tmpl(i),:),DIM=1)  
       ELSE
        tmpl(i+1) = MINLOC( bb(:, tmpl(i)),DIM=1) 
       ENDIF
        
      v1 = 1.d30
      DO s=1,n       
         IF (MINVAL(ABS(tmpn-s)).NE. 0) THEN    ! es decir, si s no es un index de tmpn (para no asignar negocios ya asignados)    
   
         v0 = DOT_PRODUCT(b(tmpl(i+1),tmpl(1:i)) , a(s,tmpn(1:i)))+&
              DOT_PRODUCT(b(tmpl(1:i),tmpl(i+1)) , a(tmpn(1:i),s)) 
          
         IF (v0 < v1) THEN
             v1=v0 
             r=s    
         ENDIF
        ENDIF  
      ENDDO      
      tmpn(i+1)=r        !nuevo negocio
      bb(tmpl(i+1),tmpl(1:i+1))=1.d30    ! actualizar a para no asignar negocios ya asignados
      bb(tmpl(1:i+1),tmpl(i+1))=1.d30
    ENDDO 
    CALL ISORT@(aux, tmpn,n) 
    p=tmpl(aux) !ordino      
   
    CALL FUNOBJ(p,zp) 
    IF (zp < z) THEN
       z=zp
       vec=p
    ENDIF   
    IF (q==1) GOTO 40
  ENDDO
 ENDDO
 40 ENDSUBROUTINE greedy2 
                         
