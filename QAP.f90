MODULE datos
  INTEGER  n , sim ! sim = 1 if A and B are symmetric, sim = 0 otherwise
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:)     
  DOUBLE PRECISION, PARAMETER :: t_max = 1.d0 
  DOUBLE PRECISION :: BKS
ENDMODULE datos 

!=====================================================================
 PROGRAM QAP   
 USE datos   
 IMPLICIT NONE                         
!.. Local Scalars ..    
 DOUBLE PRECISION tcpu1, tcpu2!, BKS 
 CHARACTER(LEN=20) nomfile, nomfile2
 INTEGER m,i,j,useless,r,s,t,itOpt  
!.. External Functions ..
 DOUBLE PRECISION, EXTERNAL:: DELTA, DELTA1, DELTA2   
 INTEGER, EXTERNAL:: MODDISTANCIA 
!.. External Subroutines .. 
 EXTERNAL funobj
 EXTERNAL greedy, greedy2, greedy3  
 EXTERNAL opt2first, opt2best, opt3first, opt3best
 EXTERNAL TEMPLE

OPEN (10,FILE='todos.txt')  
OPEN (90,FILE='todos_sln.txt')
OPEN (12,file='solucion.sol')       

CALL CPU_TIME (tcpu1)
DO
  READ(10,'(A)',END=10) nomfile
  OPEN(11, FILE=nomfile)
  READ(11,*) sim, n                       
  WRITE(*,'(A)') nomfile        
!  WRITE(*,9998) 'Dimension = ', n        
  WRITE(12,'(/A,3X,I0)') nomfile, n   
  READ(90,'(A)',END=10) nomfile2   
  OPEN(91, FILE=nomfile2)  
  READ(91,*) useless, BKS        
!  WRITE(*,'(/A,3X,I10)')'Best known solution', INT(BKS)
  ALLOCATE(A(n,n),B(n,n))           
  DO i=1,n
    READ(11,*) A(i,:)
  ENDDO      
  DO i=1,n
    READ(11,*) B(i,:)
  ENDDO             
  CALL resolucion     
  DEALLOCATE(A,B)           
ENDDO  
10 CALL CPU_TIME (tcpu2)
WRITE (*,'(A,F10.3)') 'total time = ', tcpu2-tcpu1
9998 FORMAT(T5,(A),T20,i5)
CONTAINS                            

!=====================================================================
SUBROUTINE resolucion  
USE datos  
!.. Local Scalars .. 
 INTEGER itOpt, it  
 INTEGER  i,j ,ii, tabu_size, n1, n2, q, i1, i2 ,i3 ,i4
 DOUBLE PRECISION mu,k          
 DOUBLE PRECISION z_greedy1, z_greedy2, z_greedy3     
 DOUBLE PRECISION z,ztemp,z_random     
!.. Local Arrays ..              
DOUBLE PRECISION s(2) , q_ant, alpha_ant
INTEGER :: m_ant , S_ant , R_ant
INTEGER :: pop(n,m)
INTEGER x(n), p(n),  xstar(n)  , p_ant(n)         
INTEGER p_greedy1(n),p_greedy2(n),p_greedy3(n), p_random(n)    
DOUBLE PRECISION :: SOLUTIONS(4)
DOUBLE PRECISION  z_localsearch(4), z_ant
!Greedy
Double precision:: TEMPO(6)    

q = 0
CALL greedy3 (q , n1 , n2 , p_greedy3 , z_greedy3)  


itOpt = 1
SOLUTIONS(:) = 1.d30

!-------------!
! Tabu Search !
!-------------!  


DO it = 1 , itopt
  CALL TABU(p, z)        
  IF (z < SOLUTIONS(1)) SOLUTIONS(1) = z
ENDDO
 SOLUTIONS(1) =          (SOLUTIONS(1) - BKS) / BKS * 100.d0   


!-----!
! ACO ! 
!-----!    

DO it = 1 , itopt                                 
  CALL ACO(p , z_ant)   
  IF (z_ant < SOLUTIONS(2)) SOLUTIONS(2) = z_ant
END DO  
SOLUTIONS(2) = (SOLUTIONS(2) - BKS)*100d0/BKS                   


!------!
! GVNS ! 
!------!

DO it = 1 , itopt        
   CALL GVNSfirst(p_greedy3, z_greedy3, x, z)   
   IF (z < SOLUTIONS(3)) SOLUTIONS(3) = z
END DO
SOLUTIONS(3) = (SOLUTIONS(3) - BKS) / BKS * 100.d0   

DO it = 1 , itopt        
   CALL GVNSbest(p_greedy3, z_greedy3, x, z)   
   IF (z < SOLUTIONS(4)) SOLUTIONS(4) = z
END DO
   SOLUTIONS(4) =          (SOLUTIONS(4) - BKS) / BKS * 100.d0   



WRITE(*,1235) nomfile, SOLUTIONS



949  FORMAT(t4,A,1x,1(' & ',f6.2,1x),'\\')         
998  FORMAT(t4,44('-'))  
999  FORMAT(t6, 'Algorithm',t34,'Solution')
9999 FORMAT(T6,(A),T30,F13.2)
1235  FORMAT(A,1x,4('&'1x,f5.2,1x),'\\')
1234 FORMAT('\textt{'A'}',1x'&'1x,A,1x'&'1x,4(f4.2,1x'&'1x),'\\')

END SUBROUTINE resolucion    
END PROGRAM QAP 
