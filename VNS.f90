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
