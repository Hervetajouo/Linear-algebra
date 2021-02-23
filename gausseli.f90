PROGRAM GAUSSELI
      IMPLICIT NONE
      INTEGER, PARAMETER:: n =3
      INTEGER:: i,j,k
      real(8):: A(n,n+1),x(n), sum_x
      open(1, file='mat1.dat', status='unknown')
      !Reading the matrix from file
      do i =1,n
     
      read(1,*) A(i,:)
      
      enddo
      !GAUSS-ELImation algorithm
      do k=1,n-1
      do i=k+1,n
      do j = k+1,n+1
      A(i,j)=A(i,j)-(A(i,k)/A(k,k))*A(k,j)
      enddo
      enddo
      enddo
      x(n) = A(n,n+1)/A(n,n)
      
      do k=n-1,1,-1
      sum_x=0d0
      do j=k+1,n
      sum_x=sum_x+A(k,j)*x(j)
      enddo
      x(k)=(A(k,n+1)-sum_x)/A(k,k)
      enddo
      print*, x
      
      end program GAUSSELI
