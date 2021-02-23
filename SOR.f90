PROGRAM Gauss_Seidel
      IMPLICIT NONE
      REAL(8), ALLOCATABLE :: A(:,:), B(:,:), X(:), Xo(:)
      REAL(8) :: S, tol , S1, S2, S3
      INTEGER :: i, j,  k, n, m , p
      REAL(8)::w
       n = 3 ;  tol=1.0E-6 ;  m = 100

      OPEN(UNIT = 1 , FILE = 'mat1.dat')
      OPEN(UNIT = 2, FILE = 'vect1.dat')
      ALLOCATE(A(n,n))
      ALLOCATE(B(n,1))
      ALLOCATE(X(n))
      ALLOCATE(Xo(n))

      Xo = 1.0
      print*,'Enter de value of w between 1 and 2'
      read*,w

      DO i = 1, n
      READ(1,*) A(i,:)
      READ(2,*) B(i,:)
      END DO
      CLOSE(1)
      CLOSE(2)

! Implement the convergence condition 
      DO i = 1, n
      S1 = 0.0
      DO j = 1, n
      IF (j .NE. i ) THEN
             S1 = S1 + ABS( A(i,j))
     END IF
     END DO
     END DO

     IF (ABS(A(1,1)) > S1 .AND.ABS(A(2,2))>S1 .AND.ABS(A(3,3))>S1 .AND. A(4,4)> S1 ) THEN

 ! Perform Gauss-Seidel  method if the convergence condition is satisfied 
      DO k  = 1, m
       DO i = 1 , n 
       !Compute the first sum of the G-S method 
       S2 = 0.0 
       DO j = 1, i-1
       S2 = S2 + A(i,j)*X(j)
       END DO 
       ! Compute the second sum of the G-S formula
       S3 = 0.0 
       DO j = i, n 
       S3 = S3 + A(i,j)*Xo(j)
       END DO 
       S = S2 + S3 
       ! Compute the solution X
       X(i) = Xo(i)+w*( B(i,1) - S )/A(i,i)
       END DO 

       IF (NORM2( X - Xo) < tol) EXIT 
       Xo = X
       PRINT*, k, X 

       END DO

      PRINT*, 'The solution is', X 
             ELSE
                PRINT*, 'The convergence condition is not satisfied'

        END IF

        END PROGRAM Gauss_Seidel
 
              

