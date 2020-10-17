
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
ENDSUBROUTINE greedymod      
  
  
ENDSUBROUTINE TABU       
