!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO GREEDY  ! !version 3
!!!!!!!!!!!!!!!!!!!!!!  
      
 
SUBROUTINE greedy3(q,n1,n2,vec,z_best)  
USE datos  
INTEGER,INTENT(IN):: q ! q=1 si quiero eligir n1 e n2; q=0 para hacer todas las combinaciones de n1 y n2   
INTEGER, INTENT(IN):: n1,n2 ! para q=1 n1 negocio (planta) inicial, n2 lugar inicial 
INTEGER i,j,s,k  
INTEGER tmpn(n),tmpl(n), aux(n),p(n) 
INTEGER, INTENT(OUT):: vec(n)
DOUBLE PRECISION aa(n,n), bb(n,n)
DOUBLE PRECISION zp
DOUBLE PRECISION, INTENT(OUT)::z_best   
z_best=1.d30 
tmpn = 0
tmpl = 0
DO i=1,n 
 DO j=1,n 
   tmpn(1)=i 
   tmpl(1)=j
   IF (q==1)THEN
       tmpn(1)=n1 
       tmpl(1)=n2
   ENDIF  
   aa=a  !aa = matriz temporal plantas 
   aa(tmpn(1),tmpn(1)) = -1.d30
   bb=b   !bb = matriz temporal lugares         
   bb(tmpl(1),tmpl(1))=1.d30 
   DO s=1,n-1     
     IF (sim==1) THEN
        tmpn(s+1) = MAXLOC( aa( tmpn(s),:),DIM=1) 
        tmpl(s+1) = MINLOC( bb( tmpl(s),:),DIM=1)
     ELSE            
        IF (MAXVAL(aa(tmpn(s),:),DIM=1) > MAXVAL(aa(:,tmpn(s)),DIM=1)) THEN                       
           tmpn(s+1) = MAXLOC( aa( tmpn(s),:),DIM=1)  
        ELSE
           tmpn(s+1) = MAXLOC( aa(:, tmpn(s)),DIM=1) 
        ENDIF
        IF (MINVAL(bb(tmpl(s),:), DIM=1) < MINVAL(bb(:,tmpl(s)),DIM=1)) THEN                       
           tmpl(s+1) = MINLOC( bb( tmpl(s),:),DIM=1)  
        ELSE
           tmpl(s+1) = MINLOC( bb(:, tmpl(s)),DIM=1) 
        ENDIF
     ENDIF
     aa(tmpn(1:s+1),tmpn(1:s+1))=-1.d30     ! actualizo aa para no asignar negocios (plantas) ya asignados
     bb(tmpl(1:s+1),tmpl(1:s+1))= 1.d30     ! actualizo bb para no asignar lugares ya asignados  
   ENDDO    
   123 FORMAT(12(I2,1x))
   WRITE(*,123) tmpn
   CALL ISORT@(aux, tmpn,n) 
   WRITE(*,123) aux
   p = tmpl(aux) !ordino 
   CALL FUNOBJ(p,zp)
   IF (zp < z_best) THEN
      z_best=zp
      vec=p
   ENDIF 
   IF (q==1) RETURN 
 ENDDO 
ENDDO     
ENDSUBROUTINE greedy3
