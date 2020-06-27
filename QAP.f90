MODULE datos
  INTEGER  n , sim ! sim = 1 para un problema simetrico, sim=0 si no
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: A, B 
ENDMODULE datos 


PROGRAM QAP   
USE datos
                       
IMPLICIT NONE  
INTEGER m,i,j,d
INTEGER r,s,t, itOpt 
INTEGER, EXTERNAL:: MODDISTANCIA 
DOUBLE PRECISION, EXTERNAL:: DELTA, DELTA1, DELTA2   
EXTERNAL greedy, greedy2, greedy3  
EXTERNAL OPT2PRIMERA, OPT2MAYOR, OPT3PRIMERA, OPT3MAYOR 
EXTERNAL FUNOBJ, TEMPLE

 
OPEN(12,file='Solucion.sol')  
OPEN(11,file='nug12.dat')
READ(11,*) n                       
WRITE(*,*)
WRITE(*,9998) 'Dimension = ', n
WRITE(*,*)
    
  ALLOCATE(A(n,n),B(n,n))
            
  PRINT*, '1 if the problem is symmetric, 0 if it is not'
  
  READ(*,*) sim
  
  DO i=1,n
    READ(11,*) A(i,:)
  ENDDO    
  
  DO i=1,n
    READ(11,*) B(i,:)
  ENDDO           
      
  CALL resolucion     
  DEALLOCATE(A,B)         
  
  9998 FORMAT(T5,(A),T20,i5)
  
CONTAINS 


SUBROUTINE resolucion  
USE datos             
DOUBLE PRECISION s(2)
INTEGER itOpt 
INTEGER x(n), p(n),  i,j ,ii, tabu_size, xstar(n)
INTEGER p_greedy(n), p_random(n)
DOUBLE PRECISION z,ztemp, ztemp2(4)    
INTEGER itmax


itOpt=1

WRITE(12,998)
WRITE(12,999)
WRITE(12,998)






print*,'=============='
print*,'|   RANDOM   |'
print*,'=============='

p_random = (/ (i, i=1, n) /)

CALL G05EHF(p_random,n,0)

CALL FUNOBJ(p_random,z)
WRITE(*,9999) 'Random', z


print*,'================'
print*,'|   GREEDY     |'
print*,'================'

CALL greedy (0,0,0,p_greedy,ztemp)
WRITE(*,9999) 'Greedy 1', ztemp
CALL greedy2 (0,0,0,p_greedy,ztemp)
WRITE(*,9999) 'Greedy 2', ztemp
CALL greedy3 (0,0,0,p_greedy,ztemp)
WRITE(*,9999) 'Greedy 3', ztemp

    
                    
ztemp2 = 1.d30  
                                                                                
!Busqueda Local 

print*,'======================'
print*,'|   LOCAL SEARCH     |'
print*,'======================'

!DO ii=1,itOpt
!  WRITE(*,9998)'iteracion : ', ii
  CALL RANDOM_NUMBER(s)
  i= FLOOR((n-1)*s(1)+1)  !numeri interi aleatori fra 0 e n 
  j= FLOOR((n-1)*s(2)+1)            
  CALL greedy(1,i,j,p,ztemp)
  
  CALL OPT2PRIMERA(p,ztemp,x,z)
  IF(z < ztemp2(1)) ztemp2(1) = z 
  CALL OPT2MAYOR(p,ztemp,x,z)   
  IF(z < ztemp2(2)) ztemp2(2) = z
  CALL OPT3PRIMERA(p,ztemp,x,z)
  IF(z < ztemp2(3)) ztemp2(3) = z
  CALL OPT3MAYOR(p,ztemp,x,z)
  IF(z < ztemp2(4)) ztemp2(4) = z 
!ENDDO
  
  WRITE(*,9999) 'Opt2 First Improvement', ztemp2(1)
  WRITE(*,9999) 'Opt2 Best  Improvement', ztemp2(2)
  WRITE(*,9999) 'Opt3 First Improvement', ztemp2(3)
  WRITE(*,9999) 'Opt3 Best  Improvement', ztemp2(4)  

! Temple Simulado

!print*,'=============================='
!print*,'|   INICIO TEMPLE SIMULADO   |'
!print*,'=============================='

ztemp2= 1.d30
!DO ii = 1,itOpt
!      WRITE(*,9998)'iteracion : ', ii
!      CALL RANDOM_NUMBER(s)
!      i= FLOOR((n-1)*s(1)+1)  !numeri interi aleatori fra 0 e n 
!      j= FLOOR((n-1)*s(2)+1)
!      CALL greedy(1,i,j,p,ztemp)
!      CALL TEMPLE(p,ztemp,x,z)
!      IF(z < ztemp2(1)) ztemp2(1) = z
!   ENDDO

!PRINT*, ' '
!WRITE(12,9999) 'Simulated annealing ',  ztemp2(1)

!!!!!!!!!!!!!!!!!
!  TABU SEARCH  ! 
!!!!!!!!!!!!!!!!!


print*,'====================='
print*,'|   TABU SEARCH     |'
print*,'====================='

     
  tabu_size = 50
  itmax = 60
 CALL TABU(p_greedy,z, tabu_size, itmax)    
  WRITE(*,9999) 'Tabu Search' , z







! Colonia de hormigas
!  WRITE(12,9999) 'Colonia di formiche' , ztemp2(1)
! Grasp 
!  WRITE(12,9999) 'GRASP' , ztemp2(1)
!  CALL FORMICHE(x,z) 



998  FORMAT(t4,44('-'))  
999  FORMAT(t6, 'Algoritmo',t34,'Risultato')
9998 FORMAT(T6,(A),T20,i5 \)
9999 FORMAT(T6,(A),T30,F10.2)

ENDSUBROUTINE resolucion    

ENDPROGRAM QAP 

  
 



SUBROUTINE  FUNOBJ(p,zp) ! calculo funcion objectivo
USE datos 
INTEGER, INTENT(IN)::p(n)
DOUBLE PRECISION, INTENT(OUT)::zp
INTEGER i,j
 
zp=0.d0 
DO i=1,n
   DO j=1,n  
      zp = zp + A(i,j)*B(p(i),p(j)) 
   ENDDO
ENDDO
ENDSUBROUTINE FUNOBJ 
  
!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO GREEDY  ! !version 1
!!!!!!!!!!!!!!!!!!!!!!
  
 SUBROUTINE greedy(q,n1,n2,vec,z) ! vec: mejor permutacion encontrada por el greedy; z: su valor de la funcion objectivo       
 USE datos                   
 INTEGER,INTENT(IN):: q     ! q=1 si quiero eligir n1 e n2; q=0 para hacer todas las combinaciones de n1 y n2 
 INTEGER, INTENT(IN)::n1,n2 ! para q=1 n1 negocio (planta) inicial, n2 lugar inicial 
 
 INTEGER i,j,k , r,s                     
 INTEGER tmpn(n),tmpl(n), aux(n),p(n)        
 DOUBLE PRECISION aa(n,n),v0,v1
 DOUBLE PRECISION zp    
 
 INTEGER, INTENT(OUT):: vec(n) !Permutacion final: por ejemplo vec(a)=b significa que al negocio (planta) a asocio lugar b 
 DOUBLE PRECISION, INTENT(OUT):: z
   
 z=1.d30
 
 DO k=1,n
  DO j=1,n
    tmpl=0            
    aa = a      !aa = matriz plantas   
    tmpn(1) = k !vector temporal para negocios (plantas) asignados
    tmpl(1) = j ! vector temporal para lugares asignados   
    IF (q==1) THEN
      tmpn(1) = n1 !vector temporal para negocios (plantas) asignados
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
    IF (zp < z) THEN
       z=zp
       vec=p
    ENDIF
    IF (q==1) GOTO 35
  ENDDO
 ENDDO
35 ENDSUBROUTINE greedy 
  
 
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
                         
 
 
!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO GREEDY  ! !version 3
!!!!!!!!!!!!!!!!!!!!!!  
      
 
SUBROUTINE greedy3(q,n1,n2,vec,z)  
USE datos  
INTEGER,INTENT(IN):: q ! q=1 si quiero eligir n1 e n2; q=0 para hacer todas las combinaciones de n1 y n2   
INTEGER, INTENT(IN):: n1,n2 ! para q=1 n1 negocio (planta) inicial, n2 lugar inicial 
INTEGER i,j ,s  
INTEGER tmpn(n),tmpl(n), aux(n),p(n) 
INTEGER, INTENT(OUT):: vec(n)
DOUBLE PRECISION aa(n,n), bb(n,n)
DOUBLE PRECISION zp
DOUBLE PRECISION, INTENT(OUT)::z  
 
z=1.d30 

DO i=1,n 
 DO j=1,n 
   tmpn(1)=i 
   tmpl(1)=j
   IF (q==1)THEN
       tmpn(1)=n1 
       tmpl(1)=n2
   ENDIF
   
   aa=a  !aa = matriz temporal plantas 
   aa(tmpn(1),tmpn(1))=1.d30
   bb=b   !bb = matriz temporal lugares 
   bb(tmpl(1),tmpl(1))=1.d30 
    
   DO s=1,n-1 
    
     tmpn(s+1)=MAXLOC(aa(tmpn(s),:),DIM=1)
     tmpl(s+1)=MINLOC(bb(tmpl(s),:),DIM=1) 
    
     aa(tmpn(s+1),tmpn(1:s+1))=1.d30    ! actualizo aa para no asignar negocios (plantas) ya asignados
     bb(tmpl(s+1),tmpl(1:s+1))=1.d30    ! actualizo bb para no asignar lugares ya asignados
   ENDDO 
   
    CALL ISORT@(aux, tmpn,n) 
    p = tmpl(aux) !ordino 
    CALL FUNOBJ(p,zp)
    IF (zp < z) THEN
       z=zp
       vec=p
    ENDIF 
    IF (q==1) GOTO 45
 ENDDO
ENDDO     
45 ENDSUBROUTINE greedy3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO BUSQUEDA LOCAL  2-OPTIMO !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION DELTA(p,i1,i2) !variacion de la funcion objectivo para (i1,i2)-->(i2,i1) 
USE datos 
INTEGER, INTENT(IN):: i1,i2,p(*)  
INTEGER j

IF (sim==0)  THEN   ! si el problema no es simetrico
   DELTA= (A(i1,i1)-A(i2,i2))*(B(p(i2),p(i2))-B(p(i1),p(i1)))+ &
          (A(i2,i2)-A(i2,i1))*(B(p(i2),p(i1))-B(p(i1),p(i2))) 
   DO j=1,n   
      IF((j .NE. i1) .AND. (j .NE. i2) ) THEN
         DELTA=DELTA+ (A(i1,j)-A(i2,j))*(B(p(i2),p(j))-B(p(i1),p(j)))+& 
                      (A(j,i1)-A(j,i2))*(B(p(j),p(i2))-B(p(j),p(i1))) 
      ENDIF
   ENDDO       
  
ELSEIF (sim==1) THEN   ! si el problema es simetrico
  DELTA=0.D0

  DO j=1,n   
     IF((j .NE. i1) .AND. (j .NE. i2) ) THEN
       DELTA=DELTA+ (A(j,i1)-A(j,i2))*(B(p(j),p(i2))-B(p(j),p(i1)))  
     ENDIF
  ENDDO
      
  DELTA= DELTA*2 
ENDIF     

ENDFUNCTION  DELTA

SUBROUTINE OPT2PRIMERA(p,zp,m,zm)  ! primera mejora 2-optimo
USE datos
DOUBLE PRECISION, INTENT(IN):: zp
INTEGER, INTENT(IN):: p(n)
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm   
INTEGER i1, i2 ,aux 
DOUBLE PRECISION temp
DOUBLE PRECISION, EXTERNAL:: DELTA
  
m=p 
zm=zp
50 DO i1=1, n-1 
   DO i2=i1+1,n             
         temp=DELTA(m,i1,i2)
      IF (temp < 0.d0 ) THEN 
         zm=zm+temp
!         print*, zm
         aux=m(i1)
         m(i1)=m(i2)
         m(i2)=aux 
         GOTO 50 
      ENDIF
   ENDDO 
ENDDO
ENDSUBROUTINE OPT2PRIMERA   
 
SUBROUTINE OPT2MAYOR(p,zp,m,zm) ! mayor mejora 2-optimo 
USE datos
DOUBLE PRECISION, INTENT(IN):: zp
INTEGER, INTENT(IN):: p(n)
INTEGER, INTENT(OUT):: m(n) 
DOUBLE PRECISION, INTENT(OUT)::zm
INTEGER i1, i2 , j1,j2, aux
DOUBLE PRECISION, EXTERNAL:: DELTA
DOUBLE PRECISION:: temp, dmin

m=p 
zm=zp 

55 dmin=0.d0 

DO i1=1,n-1
   DO i2=i1+1,n 
      temp=DELTA(m,i1,i2)
      IF(temp < dmin)  THEN
        j1=i1
        j2=i2
        dmin=temp
      ENDIF
   ENDDO
ENDDO

IF (dmin <0.d0) THEN 
   aux=m(j1)
   m(j1)=m(j2)
   m(j2)=aux  
   zm=zm+dmin 
   GOTO 55
ENDIF

ENDSUBROUTINE OPT2MAYOR
                                               

    




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO BUSQUEDA LOCAL  3-OPTIMO !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
DOUBLE PRECISION FUNCTION DELTA1(p,i1,i2,i3)  ! para (i1,i2,i3) --> (i2,i3,i1)
USE datos 
INTEGER, INTENT(IN):: i1,i2,i3,p(*)  
INTEGER j  

 
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


DOUBLE PRECISION FUNCTION DELTA2(p,i1,i2,i3) ! para (i1,i2,i3) --> (i3,i1,i2) 
USE datos 
INTEGER, INTENT(IN):: i1,i2,i3,p(*)  
INTEGER j  

 
IF (sim==0) THEN ! si el problema no es simetrico 
   DELTA2=A(i1,i1)*(B(p(i3),p(i3))-B(p(i1),p(i1)))+A(i1,i2)*(B(p(i3),p(i1))-B(p(i1),p(i2)))+&
          A(i1,i3)*(B(p(i3),p(i2))-B(p(i1),p(i3)))+A(i2,i1)*(B(p(i1),p(i3))-B(p(i2),p(i1)))+&
          A(i2,i2)*(B(p(i1),p(i1))-B(p(i2),p(i2)))+A(i2,i3)*(B(p(i1),p(i2))-B(p(i2),p(i3)))+ &
          A(i3,i1)*(B(p(i2),p(i3))-B(p(i3),p(i1)))+A(i3,i2)*(B(p(i2),p(i1))-B(p(i3),p(i2)))+ &
          A(i3,i3)*(B(p(i2),p(i2))-B(p(i3),p(i3)))
          
    DO j=1,n   
        IF((j .NE. i1) .AND. (j .NE. i2). AND. (j.NE.i3) ) THEN
          DELTA2=DELTA2+  A(i1,j)*(B(p(i3),p(j))-B(p(i1),p(j)))+A(i2,j)*(B(p(i1),p(j))-B(p(i2),p(j)))+&
                          A(i3,j)*(B(p(i2),p(j))-B(p(i3),p(j)))+A(j,i1)*(B(p(j),p(i3))-B(p(j),p(i1)))+&
                          A(j,i2)*(B(p(j),p(i1))-B(p(j),p(i2)))+A(j,i3)*(B(p(j),p(i2))-B(p(j),p(i3)))
        ENDIF
    ENDDO   
ELSEIF (sim==1) THEN ! si el problema es simetrico 
    DELTA2=A(i1,i2)*(B(p(i1),p(i3))-B(p(i1),p(i2)))+A(i1,i3)*(B(p(i2),p(i3))-B(p(i1),p(i3)))+&
           A(i2,i3)*(B(p(i1),p(i2))-B(p(i2),p(i3)))
    DO j=1,n  
          IF((j .NE. i1) .AND. (j .NE. i2). AND. (j.NE.i3) ) THEN
             DELTA2=DELTA2+  A(j,i1)*(B(p(j),p(i3))-B(p(j),p(i1)))+A(j,i2)*(B(p(j),p(i1))-B(p(j),p(i2)))+&
                             A(j,i3)*(B(p(j),p(i2))-B(p(j),p(i3)))
           ENDIF
    ENDDO  
    DELTA2= DELTA2*2
ENDIF      

ENDFUNCTION  DELTA2

SUBROUTINE OPT3PRIMERA(p,zp,m,zm)  ! primera mejora 3-optimo
USE datos
DOUBLE PRECISION, INTENT(IN):: zp
INTEGER, INTENT(IN):: p(n)
INTEGER, INTENT(OUT)::  m(n)
DOUBLE PRECISION, INTENT(OUT):: zm  
INTEGER i1, i2 ,i3,aux
DOUBLE PRECISION, EXTERNAL:: DELTA1, DELTA2
DOUBLE PRECISION temp
  
m=p 
zm=zp
  
60 DO i1=1,n-2 
    DO i2=i1+1,n-1
      DO i3=i2+1,n 
           temp=DELTA1(m,i1,i2,i3)  
          IF (temp < 0.d0) THEN
            zm=zm+temp
            aux=m(i1)
            m(i1)=m(i2) 
            m(i2)=m(i3)
            m(i3)=aux
            GOTO 60 
          ENDIF 
          
           temp=DELTA2(m,i1,i2,i3)
          IF (temp <0.d0) THEN
            zm=zm+temp
            aux=m(i1)
            m(i1)=m(i3)
            m(i3)=m(i2)
            m(i2)=aux
            GOTO 60
          ENDIF
      ENDDO
   ENDDO
ENDDO 

ENDSUBROUTINE  OPT3PRIMERA   

SUBROUTINE OPT3MAYOR(p,zp,m,zm)  ! mayor mejora 3-optimo
USE datos
DOUBLE PRECISION, INTENT(IN):: zp
INTEGER, INTENT(IN):: p(n)
INTEGER, INTENT(OUT):: m(n)
DOUBLE PRECISION, INTENT(OUT)::zm
INTEGER i1, i2,i3, j1,j2,j3,l, aux 
DOUBLE PRECISION, EXTERNAL:: DELTA1, DELTA2 
DOUBLE PRECISION temp, dmin

 
m=p 
zm=zp 

65 dmin=0.d0 

DO i1=1,n-2
   DO i2=i1+1,n-1
      DO i3=i2+1,n
         temp=DELTA1(m,i1,i2,i3)
         IF(temp < dmin)  THEN
           j1=i1
           j2=i2 
           j3=i3
           l=1
           dmin=temp
         ENDIF  
         temp=DELTA2(m,i1,i2,i3)
         IF(temp < dmin)  THEN
           j1=i1
           j2=i2 
           j3=i3
           l=2
           dmin=temp
         ENDIF 
      ENDDO
   ENDDO
ENDDO

IF (dmin <0.d0) THEN
   aux=m(j1)
   IF (l==1) THEN
       m(j1)=m(j2)
       m(j2)=m(j3)
       m(j3)=aux 
       
   ELSEIF(l==2) THEN 
       m(j1)=m(j3)
       m(j3)=m(j2)
       m(j2)=aux
   ENDIF 
   zm=zm+dmin 
   GOTO 65
ENDIF

ENDSUBROUTINE OPT3MAYOR      



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ALGORITMO TEMPLE SIMULADO !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

SUBROUTINE TEMPLE(x,zx,xstar,z)
USE datos  
INTEGER, INTENT(OUT):: xstar(n)  ! permutacion final
DOUBLE PRECISION, INTENT(OUT):: z ! funcion objectivo en la solucion inicial
INTEGER x(n),aux
DOUBLE PRECISION zx, T , m(2),d, alfa, beta, tmax 
DOUBLE PRECISION t1, t2 ! scegliamo come criterio di arresto il tempo 
DOUBLE PRECISION, EXTERNAL:: DELTA 
INTEGER it , i1, i2

IF (n <= 49) THEN !tamano del problema: pequeno --> utilizo su calibracion
   T=5.d3
   alfa=9.d-1
   beta=1.d7
   tmax= 5.d0
ELSE !tamano del problema : grande
   T=5.d3
   alfa=95.d-2
   beta=1.d8
   tmax=5.d0
ENDIF
CALL CPU_TIME(t1)
xstar=x 
z=zx   
CALL CPU_TIME(t2)        
DO WHILE (abs(t1-t2) < tmax) !criterio de parada: tiempo maximo
      it=0
    DO WHILE(it < beta/T .AND. abs(t1-t2)< tmax)
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

                       
                                              
!!!!!!!!!!!!!!!!!
!  TABU SEARCH  !
!!!!!!!!!!!!!!!!!

SUBROUTINE TABU(p,z_best,tabu_size,max_iter)
USE datos                  
INTEGER tabu_size, i,j, temp, aux, k, i1, i2, j1, j2, ii, h      
INTEGER, INTENT(IN)  :: max_iter  ! maximum number of iterations
INTEGER tabu_list(tabu_size)  ! tabu list        
INTEGER p(n)       ! starting permutation              
INTEGER p_best(n)  ! best known solution 
DOUBLE PRECISION zp  
INTEGER loc_count  ! local counter
INTEGER M          ! parameter
INTEGER loc_best   ! best neighbor
DOUBLE PRECISION, INTENT(OUT):: z_best ! best known obj. function value
DOUBLE PRECISION dmin


CALL FUNOBJ(p,zp)     
! initial solution 
z_best = zp         
M = 10
loc_best  = - M
loc_count = 0
ii = 1            

tabu_list = 0          

DO k = 1, max_iter  
! print*,'k=',k      
  555 dmin = 0.d0         
!  print*,'iterazione',k
loop1:   DO i1 = 1, n-1
 loop2 :  DO i2 = i1+1, n 
            temp = DELTA(p, i1, i2)               
            h = search(i1, i2, tabu_size, tabu_list)    
            IF ( (temp < dmin .AND. h == 0 ) .OR. (temp < dmin .AND. h == 1 .AND. zp + temp < z_best) )  THEN   
              j1 = i1
              j2 = i2
              dmin = temp   
              ! adding p(i1), p(i2) in the tabu_list  
              tabu_list(ii) = p(i1)
              tabu_list(ii+1) = p(i2)
              ii = ii + 2                   
            ENDIF      
            IF (ii >= tabu_size) EXIT loop1   
            IF (dmin < 0.d0) THEN 
              aux   = p(j1)
              p(j1) = p(j2)
              p(j2) = aux  
              z_best = z_best + dmin         
              print*,'z_best = ', z_best
              GOTO 555    
            ENDIF
          ENDDO loop2  
         ENDDO loop1                                                          
 
 ! scramble p(1) and p(2) and put ( p(1), p(2) ) in the tabu_list      
! aux  = p(1)
! p(1) = p(2)
! p(2) = aux
! tabu_list(ii) = p(1)
! tabu_list(ii+1) = p(2)
   
ENDDO 



CONTAINS
INTEGER FUNCTION search(i,j,n,v)  
! this function looks for (i,j) in the list;
! search =  1 if they are in the list, 0 if not.
INTEGER i, j, it, n
INTEGER v(n)
it = 1 
search = 0       
DO   
  IF ( (v(it)==i .AND. v(it+1)==j) .OR. (v(it)==j .AND. v(it+1)==i) )  THEN
    search = 1 
    RETURN
  ENDIF
  it = it + 2   
  IF (it >= n) RETURN
ENDDO      
ENDFUNCTION search

ENDSUBROUTINE TABU 

                                                                                                            
                                                                                                                                                                                                                                                                                                                                    






!!!!!!!!!!!!!!!!!!
!! INIZIO GRASP !!
!!!!!!!!!!!!!!!!!!

SUBROUTINE GRASP(x,xstar,z)
USE datos
INTEGER, INTENT(OUT):: xstar(n)  ! permutacion final
DOUBLE PRECISION, INTENT(OUT):: z ! funcion objectivo en la solucion inicial  
DOUBLE PRECISION x
ENDSUBROUTINE GRASP
  
!!!!!!!!!!!!!!!!!!!!!!!
!! INIZIO ANT COLONY !!
!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FORMICHE(xstar,z)
USE datos
INTEGER, INTENT(OUT):: xstar(n)  ! permutacion final
DOUBLE PRECISION, INTENT(OUT):: z ! funcion objectivo en la solucion inicial  
DOUBLE PRECISION x 
INTEGER pi(n),m(n),permu(n), Mat1(n,n),Mat2(n,n),aux
INTEGER i,j,ii,zpi,zm, minimo,  medias, zmej,zmin, i1,i2
DOUBLE PRECISION temp
DOUBLE PRECISION, EXTERNAL:: DELTA               

minimo = 9999999
medias = 0
zmej   = 9999999
zmin   = 0   
Mat1 = 0 
Mat2 = 0  
permu = 0
DO i=1,n

  ! Definizione del vettore p_iniziale    
  ! pi = [i,i+1,...,n,...,i-1]
  DO j=1,n
    pi(j) = MODULO(j + i - 1, n)     
  ENDDO     
  pi(n-i+1) = n
    
  CALL FUNOBJ(pi,zpi) 
  CALL OPT2PRIMERA(pi,zpi,m,zm)
  
  DO j=1,n
     Mat1(j,m(j)) = Mat1(j,m(j)) + 1
  ENDDO
 
 
ENDDO        

WRITE(12,*)
WRITE(12,*)

  
m=pi 
50 DO i1=1, n-1 
   DO i2=i1+1,n             
         temp=DELTA(m,i1,i2)
      IF (temp < 0.d0 ) THEN     
      print*, i1,m(i1)
      print*, i2,m(i2)
         Mat2(i1,m(i1)) = Mat2(i1,m(i1)) - 1
         Mat2(i1,m(i2)) = Mat2(i1,m(i2)) + 1   
         Mat2(i2,m(i1)) = Mat2(i2,m(i1)) + 1
         Mat2(i2,m(i2)) = Mat2(i2,m(i2)) - 1
         aux=m(i1)
         m(i1)=m(i2)
         m(i2)=aux     
         GOTO 50 
      ENDIF
   ENDDO      
ENDDO      

DO i=1,n
  WRITE(12,9000) Mat1(i,:)
ENDDO                           
WRITE(12,*) 
WRITE(12,*)

DO i=1,n
  WRITE(12,9000) Mat2(i,:)             
ENDDO
 
   9000 FORMAT(30(I2,1x))  
ENDSUBROUTINE FORMICHE
 


