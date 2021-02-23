PROGRAM LU_Fact
        IMPLICIT NONE
        INTEGER :: i, j, k
        INTEGER, PARAMETER :: n = 3 , m = 5
        REAL(8) :: S1, S2
        REAL(8) , DIMENSION(n) :: x, y , B 
        REAL(8),  DIMENSION(n,n) :: A , Id, U, L 

        OPEN(UNIT = 3, FILE = 'mat1.dat')
        OPEN(UNIT = 4, FILE = 'vect1.dat')
         DO k = 1,  n
         READ(3,*) A(k,:)
         READ(4,*) B(k)
         END DO
        CLOSE(3)
        CLOSE(4)
        Id = 0.0
        DO i = 1, n 
        Id(i,i) = 1.0
        END DO 
         U = A ; L = Id

         DO k = 1, n
         DO j = k+1, n 
         L(j,k) = U(j,k)/U(k,k)
         U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n)
         END DO 
         END DO

         PRINT*, 'The lower triangular matrix is'

         DO i = 1, n
         PRINT*,L(i,:)
         END DO
         
         PRINT*, 'The upper triangular matrix is'

         DO i = 1, n
         PRINT*,U(i,:)
         END DO

         PRINT*,' ------ Factorization A = LU Completed ---------'

         PRINT*, '---------------------------------------------------'

         ! We have now to solve two system 
         ! First Solve for y in  Ly = B by forward substitution
         y(1) = B(1)/L(1,1)

         DO i = 2,  n 
        S1 = 0
        DO j = 1, i-1
       S1 = S1 + L(i,j)*y(j)
       END DO 
      y(i) =( B(i) - S1 )/L(i,i)
      END DO 

     ! Now solve for x in Ux = y by backward substitution
     
    x(n) = y(n)/U(n,n)
    DO i = n-1, 1, -1
     S2 = 0.0 
     DO j = i+1, n 
     S2 = S2 + U(i,j)*x(j)
     END DO
     x(i) =( y(i) - S2)/U(i,i)
     END DO 

     PRINT*,'The solution of the system is',  x 
     
       END PROGRAM LU_Fact 

