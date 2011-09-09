! linslv.f                                                              
!    from bnrchemv5.f                                                   
!                                                                       
! This contains the subroutines LINSLV and RESOLV                       
!   for inverting a matrix.                                             
!   (Dates from Prather, 1988, and Numerical Recipes)                   
!                                                                       
! It solves for X in the matrix equation AX=B.                          
!                                                                       
!     call LINSLV(A, B, X, N) for AX=B, used dimension N                
!          A is updated to reduced form  by LINSLV.                     
!          NOTE: declared dimension of A, B, and X are set at 100 (not N
!                                                                       
!     call RESOLV(A, B, X, N) for AX=B, with A reduced from LINSLV.     
!          (if error in last solution is small)                         
!                                                                       
! THE ORIGINAL CALL:                                                    
!       err=abs(xoo(i)/xr( 1 ,ic))                                      
!      if(errxo.lt.ermat) lnumat=.false.                                
!      if(lnumat) call LINSLV(fxo,xoo,nchem)                            
!      if(.not.lnumat) call RESOLV(fxo,xoo,nchem)                       
!                                                                       
! MODIFICATIONS FROM bnrchemv5.f:                                       
!   include 'commq5.f' and 'commb5.f' are removed.                      
!   matrix A is parameter of subroutine.                                
!                                                                       
!                                                                       
! ------------------------------------------------------                
                                                                        
      SUBROUTINE LINSLV(A,B,X,N) 
!  *** SUB -LINSLV- SOLVES THE FOLLOWING MATRIX EQUATION--              
!  ***     A(N,N)*X(N) = B(N)  BY REDUCING THE A-MATRIX IN PLACE.       
!  *** THE ENTRY -RESOLV- ASSUMES THAT THE A-MATRIX HAS BEEN PROPERLY   
!  ***     REDUCED AND JUST SOLVES FOR X(N).  THIS OPTION SAVES TIME    
!  ***     WHEN THE SYSTEM IS TO BE RESOLVED WITH A NEW B-VECTOR..  (Pra
                                                                        
      IMPLICIT NONE 
                                                                        
!                                                                       
      DIMENSION B(100),X(100), A(100,100) 
      DIMENSION S(100) 
      DIMENSION IPA(100) 
                                                                        
      DOUBLE PRECISION A, B, X 
      DOUBLE PRECISION S, SMAX,DIV 
      INTEGER N, KR, K, KRM1, KRMAX, KRP1, J, JP, JP1, I, IPA 
                                                                        
      DO 10 KR=1,N 
        DO    K=1,N 
          S(K) = A(K,KR) 
                             !DO    K=1,N                               
        ENDDO 
        IF(KR.EQ.1) GO TO 14 
        KRM1 = KR - 1 
        DO    J=1,KRM1 
          JP = IPA(J) 
          A(J,KR) = S(JP) 
          S(JP) = S(J) 
          JP1 = J + 1 
          DO    I=JP1,N 
            S(I) = S(I) - A(I,J)*A(J,KR) 
                                        !DO    I=JP1,N                  
          ENDDO 
                                     !DO    J=1,KRM1                    
        ENDDO 
   14   KRMAX = KR 
        SMAX =  ABS(S(KR)) 
        DO 15 I=KR,N 
          IF( ABS(S(I)).Le.SMAX) GO TO 15 
          KRMAX = I 
          SMAX =  ABS(S(I)) 
                                  !DO 15 I=KR,N                         
   15   CONTINUE 
        IPA(KR) = KRMAX 
        A(KR,KR) = S(KRMAX) 
        DIV = 1.0E0/S(KRMAX) 
        S(KRMAX) = S(KR) 
        IF(KR.EQ.N) GO TO 10 
        KRP1 = KR + 1 
        DO    I=KRP1,N 
          A(I,KR) = S(I)*DIV 
                                !DO    I=KRP1,N                         
        ENDDO 
                              !DO 10 KR=1,N                             
   10 END DO 
      CALL RESOLV(A,B,X,IPA,N) 
      RETURN 
      END                                           
!    -------------------------------------------                        
                                                                        
      SUBROUTINE RESOLV(A,B,X, IPA, N) 
! ENTRY FOR BACK SOLUTION WITH DIFFERENT B-VALUE                        
!                                                                       
                                                                        
      IMPLICIT NONE 
                                                                        
      DIMENSION B(100),X(100), A(100,100) 
      DIMENSION S(100) 
      DIMENSION IPA(100) 
                                                                        
      DOUBLE PRECISION A, B, X 
      DOUBLE PRECISION S, SUMM 
      INTEGER I, IP, IP1, II, IIP1, J, IPA, N 
                                                                        
                                                                        
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                        
!      CALL CPUTIM(-11)                                                 
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                        
      DO    I=1,N 
        S(I) = B(I) 
                              !DO    I=1,N                              
      ENDDO 
      DO    I=1,N 
        IP = IPA(I) 
        X(I) = S(IP) 
        S(IP) = S(I) 
        IF(I.EQ.N) GO TO 23 
        IP1 = I + 1 
        DO    J=IP1,N 
          S(J) = S(J) - A(J,I)*X(I) 
                                !DO    J=IP1,N                          
        ENDDO 
                              !DO    I=1,N                              
      ENDDO 
   23 CONTINUE 
                                                                        
      DO    I=1,N 
        II = N + 1 - I 
        SUMM = X(II) 
! NOTE POSSIBLE ERROR IN THIS NEXT CHANGED LINE.                        
        IF(II.LT.N) THEN 
          IIP1 = II + 1 
          DO    J=IIP1,N 
            SUMM = SUMM - A(II,J)*X(J) 
                                    !DO    J=IIP1,N                     
          ENDDO 
                                   !IF(II.LT.N) THEN                    
        ENDIF 
        X(II) = SUMM/A(II,II) 
                              !DO    I=1,N                              
      ENDDO 
                                                                        
      RETURN 
      END                                           
