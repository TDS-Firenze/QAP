SUBROUTINE TABU(p,z_best,tabu_size,mu)
USE datos                  
INTEGER tabu_size, temp, aux, k, i1, i2, j1, j2, ii, h   
DOUBLE PRECISION Bmod(n,n)
INTEGER tabu_list(tabu_size)  ! tabu list        
INTEGER p(n)       ! starting permutation              
INTEGER p_best(n)  ! best known solution  
INTEGER p_temp(n)       
INTEGER LTM(n,n)   ! Long Term Memory Matrix     
DOUBLE PRECISION mu
DOUBLE PRECISION zp  
DOUBLE PRECISION z_best ! best known obj. function value
DOUBLE PRECISION dmin
DOUBLE PRECISION t1, t2



zp = z_best           
! initial solution      
p_best = p    
p_temp = p         
z_temp = zp
ii = 1            
LTM = 0
Bmod = B
tabu_list = 0     ! initializing the tabu list     
CALL CPU_TIME(t1)
CALL CPU_TIME(t2)

DO WHILE (abs(t1-t2) < t_max) 
555 dmin = 0.d0        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local Search : seek the best 2-exchange starting from p !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

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
                dmin = temp         
                z_temp = zp + dmin           
              ENDIF
            ENDIF         
           ENDDO loop2  
         ENDDO loop1    
         IF (ii >= tabu_size) ii = 1                    
         tabu_list(ii) = p(j1)
         tabu_list(ii+1) = p(j2)
         ii = ii + 2
         IF (dmin < 0.d0) THEN 
           p = p_temp 
           zp = z_temp
           LTM(j1,j2) = LTM(j1,j2) + 1
           LTM(j2,j1) = LTM(j2,j1) + 1   
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

p = p_best

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
  


 
  SUBROUTINE greedymod(vec,z)   
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
     aa=a  
     aa(tmpn(1),tmpn(1))=1.d30
     bb=bmod   
     bb(tmpl(1),tmpl(1))=1.d30 
     DO s=1,n-1     
       tmpn(s+1)=MAXLOC(aa(tmpn(s),:),DIM=1)
       tmpl(s+1)=MINLOC(bb(tmpl(s),:),DIM=1)     
       aa(tmpn(s+1),tmpn(1:s+1))=1.d30    
       bb(tmpl(s+1),tmpl(1:s+1))=1.d30    
     ENDDO 
   
    CALL ISORT@(aux, tmpn,n) 
    p = tmpl(aux) !ordino 
    CALL FUNOBJ(p,zp)
    IF (zp < z) THEN
       z=zp
       vec=p
    ENDIF 

   ENDDO
  ENDDO     
  ENDSUBROUTINE greedymod
ENDSUBROUTINE TABU                                                                                          
                                               
