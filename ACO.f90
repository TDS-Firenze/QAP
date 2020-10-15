! =====================================================================  
! Ant Colony


SUBROUTINE ACO(m , R, q ,  alpha , S , p_best , z_best)                                 
USE datos      
!.. Inputs ..
INTEGER, INTENT(IN) :: m , S , R
DOUBLE PRECISION, INTENT(IN) :: alpha , q                            
!.. Outputs ..
INTEGER, INTENT(OUT):: p_best(n)       ! final permutation
!DOUBLE PRECISION, INTENT(OUT):: z_best ! objective function value  
DOUBLE PRECISION :: z_best
!.. Local Arrays ..
INTEGER :: pi(n) , p_hat(n), p_tilde(n) , p_first(n) , p(n) , p_star(n)     
INTEGER :: Ants(n,m) , P_stock(n , m)
DOUBLE PRECISION :: prob(n)   
DOUBLE PRECISION :: T(n,n) , M_update(n,n) , M_out(n,n)
DOUBLE PRECISION :: z_stock1(m), z_stock2(m) 
DOUBLE PRECISION :: Tensor(n,n,m)
!.. Local Scalars..
INTEGER :: i , j , k ,  it , i_best
DOUBLE PRECISION :: z, z_hat, z_tilde, temp , u , t1 , t2    
LOGICAL :: intensification
!.. External..
EXTERNAL :: random_permutation, Delta


CALL CPU_TIME(t1)
CALL CPU_TIME(t2)
                 
  M_update = 0.d0
  z_best = 1.d30
  Ants(: , :)  = 0
   


! Initialisation of Pheromone Trail 
CALL pheromone_initialisation(T)

it = 0
!DO 
DO 
  IF ( ABS(t1 - t2) > t_max) EXIT

   
     
   CALL CPU_TIME(t2)
   it = it + 1 

 DO k = 1 , m              
  ! Initialisation of solutions   
    CALL random_permutation(p)
    CALL funobj(p, z)    
    CALL opt2first(p , z , p , z)  
    IF (z < z_best) THEN
      z_best = z
      i_best = k 
      p_best(:) = p(:)
    ENDIF                       
    
  ! Manipulation of solutions
  DO i = 1 , R                             
    CALL apply_2exchange(p_hat , z_hat) 
  END DO    
    CALL pheromone_update(p_tilde , z_tilde , M_out) 
    ants(: , k) = p_tilde(:) 
    IF (z_tilde < z_best) THEN
      z_best = z_tilde
      p_best(:) = p_tilde(:)
      intensification = .TRUE.
      M_update(: , :) = M_out(: , :)
      i_best = k
    END IF   
    z_stock1(k) = z
    z_stock2(k) = z_tilde
    P_stock(: , k) = p(:)
  END DO 

  ! Intensification
   IF (intensification = .TRUE.) THEN 
     DO k = 1 , m
       IF (z_stock1(k) < z_stock2(k))    ants(: , k) = P_stock(: , k)
     END DO
   END IF
   intensification = .FALSE.
 
 ! Pheromone trail update
   T(: , :) = (1 - alpha) * T(: , :)
   T(: , :) = T(: , :) + M_update(: , :)

!print*, p_best 

 ! Diversification
   IF (it > S) THEN
     T(: , :) = 0d0
     CALL pheromone_initialisation(T)
     Ants(: , :) = 0
   END IF
   CALL CPU_TIME(t2) 
END DO  
      


 
   9000 FORMAT(4(I1,1x))  
   9001 FORMAT(4(f3.0,1x))  
   
   
CONTAINS


SUBROUTINE pheromone_initialisation(T)
INTEGER :: i, j , i1, i2, aux
INTEGER :: pi(n)
DOUBLE PRECISION :: z, temp
DOUBLE PRECISION :: M1(n,n) , M2(n,n)
DOUBLE PRECISION, INTENT(OUT) :: T(n,n)
 M1 = 0.d0
 M2 = 0.d0  
 DO j = 1 , m 
 ! Setting pi = (i, i+1, ..., n, 1, ...., i-1)
  DO i = 1,n
    pi(i) = MODULO(i + j - 1, n)     
  END DO
  pi(n - j + 1) = n
 ! Local search on pi 
  CALL funobj(pi,z)
  56 DO i1 = 1 , n - 1 
         DO i2 = i1 + 1 , n             
           temp = DELTA(pi , i1 , i2)
           IF (temp < 0.d0 ) THEN                 
             M2(i1 , pi(i1)) = M2(i1 , pi(i1)) - 1d0
             M2(i2 , pi(i2)) = M2(i2 , pi(i2)) - 1d0
             M2(i1 , pi(i2)) = M2(i1 , pi(i2)) + 1d0   
             M2(i2 , pi(i1)) = M2(i2 , pi(i1)) + 1d0
             aux = pi(i1)
             pi(i1) = pi(i2)
             pi(i2) = aux     
             z = z + temp
             GO TO 56 
           END IF
         END DO      
       END DO     
     DO i = 1 , n
      M1(i , pi(i)) = M1(i , pi(i)) + 1
     END DO                           
END DO    
 M2 = M2 + ABS(MINVAL(M2))
  T = M1 + M2  
END SUBROUTINE pheromone_initialisation


SUBROUTINE apply_2exchange(p_hat , z_hat) 
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_hat(n)
DOUBLE PRECISION, INTENT(OUT) :: z_hat
!.. Local scalars ..
INTEGER :: it,  r , s
DOUBLE PRECISION :: u , aux , sum_prob, sum_cumulative
    ! Choosing r
    CALL RANDOM_NUMBER(u)  
    r = NINT(u * (n-1) + 1 )       
    ! Choosing s     
    CALL RANDOM_NUMBER(u)
    IF ( u < q) THEN   
      aux = 0d0
      DO it = 1 , n
        IF ( T(r , p(it)) + T(it , p(r)) > aux .AND. it /= r) THEN
          aux = T(r , p(it)) + T(it , p(r))      
          s = it 
        END IF
      END DO    
    ELSE       
    ! Exploration case                                                     
      sum_prob = 0.d0
      DO i = 1 , n            
        IF (i /= r) sum_prob = sum_prob + T(r , p(i)) + T(i, p(r))
      END DO
      DO s = 1 , n        
        prob(s) = T(r , p(s)) + T(s , p(r))
        prob(s) = prob(s) / sum_prob
      END DO                               
      prob(r) = 0.d0
      CALL RANDOM_NUMBER(u)  
      i = 1     
      sum_cumulative  = 0d0

      DO 
        sum_cumulative = sum_cumulative + prob(i)
        IF (u < sum_cumulative ) THEN
          s = i
          EXIT
        END IF  
        i = i + 1
      END DO
   END IF
   p_hat(:) = p(:)
   p_hat(r) = p(s)  
   p_hat(s) = p(r)  
   z_hat = z + DELTA(p , r , s)

END SUBROUTINE apply_2exchange

SUBROUTINE pheromone_update(p_tilde , z_tilde , M_out)
USE datos
!.. Outputs ..
INTEGER, INTENT(OUT) :: p_tilde(n)
DOUBLE PRECISION, INTENT(OUT) :: z_tilde , M_out(n,n)
!.. Local Scalars ..
INTEGER ::   i1, i2 , ii 
DOUBLE PRECISION :: z_tilde , temp 
M_out = 0d0
p_tilde = p_hat   
z_tilde = z_hat     
551 DO i1 = 1 , n - 1 
     DO i2 = i1 + 1 , n             
       temp = DELTA(p_tilde , i1 , i2)
         IF (temp < 0.d0 ) THEN      
            M_out(i1 , p_tilde(i1)) =  M_out(i1 , p_tilde(i1)) - 1d0
            M_out(i2 , p_tilde(i2)) =  M_out(i2 , p_tilde(i2)) - 1d0
            M_out(i1 , p_tilde(i2)) =  M_out(i1 , p_tilde(i2)) + 1d0   
            M_out(i2 , p_tilde(i1)) =  M_out(i2 , p_tilde(i1)) + 1d0                
           ii = p_tilde(i1)
           p_tilde(i1) = p_tilde(i2)
           p_tilde(i2) = ii     
           z_tilde  = z_tilde + temp
           GO TO 551 
         END IF
       END DO      
     END DO   
     M_out = M_out + ABS(MINVAL(M_out))
     DO i = 1 , n
       M_out(i , p_tilde(i)) = M_out(i , p_tilde(i)) + 1d0
     END DO 
END SUBROUTINE pheromone_update   

 
END SUBROUTINE ACO    
