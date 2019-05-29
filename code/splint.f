! *************************************************

      SUBROUTINE splint(xa,ya,y2a,N,x,y,dt)
      IMPLICIT REAL*8 (A-H,K,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      INTEGER :: N
      REAL(kind=DP), dimension(N) :: xa,y2a,ya
! ------------------------------------------------
!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xai’s in order), and given the array y2a(1:n), which is the output from spline above,and given a value of x, this routine returns a cubic-spline interpolated value y.
! -------------------------------------------------
      INTEGER :: k,khi,klo
      REAL(kind=DP) :: a,b,h,x1

! We will find the right place in the table by means of bisection. This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely spaced, one would do better to store previous values of klo and khi and test if they remain appropriate on the next call.

      klo=1 
      khi=n
 1    IF (khi-klo.GT.1) THEN
         k=(khi+klo)/2
         IF(xa(k).GT.x)THEN
            khi=k
         ELSE
            klo=k
         END IF
         GO TO 1
      END IF
      
      h=xa(khi)-xa(klo)

      IF (x <= xa(1)) THEN
         y = ya(1)
         RETURN
      END IF
      IF (x >= xa(n)) THEN
         y = ya(n)
         RETURN
      END IF
      
!     The xa’s must be distinct.
! Cubic spline polynomial is now evaluated.
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      
      RETURN
      END

! *******************************************
