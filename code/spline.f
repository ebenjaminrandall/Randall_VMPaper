!*********************************************
!

      SUBROUTINE spline(x,y,N,yp1,ypn,y2)
      IMPLICIT REAL*8 (A-H,K,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      INTEGER :: N,NMAX
      REAL(kind=DP), dimension(N) :: x,y,y2
      REAL(kind=DP) :: yp1,ypn
      PARAMETER (NMAX=100000)
! ---------------------------------------------------------------------
!      Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!      x1 < x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
!      function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!      length n which contains the second derivatives of the interpolating function at the tabulated
!      points xi. If yp1 and/or ypn are equal to 1 × 10**30 or larger, the routine is signaled to set
!      the corresponding boundary condition for a natural spline, with zero second derivative on
!      that boundary.
!      Parameter: NMAX is the largest anticipated value of n.
! ----------------------------------------------------------------------      
      INTEGER :: i,k
      REAL(KIND=DP) :: p,qn,sig,un
      REAL(kind=DP), dimension(NMAX) :: u
! The lower boundary condition is set either to be "natural" or else to have a specified derivative

      IF (yp1 .GT. .99e30) THEN  
         y2(1)=0
         u(1)=0.
      ELSE
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      END IF
! This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary
! storage of the decomposed factors.
      DO i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      END DO
! The upper boundary condition is set either to be "natural" or else to have a specified first derivative.
      IF (ypn .GT. .99e30) THEN
         qn=0.
         un=0.
      ELSE
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      END IF 
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
!     This is the backsubstitution loop of the tridiagonal algorithm.
      DO k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      END DO
      
      RETURN
      END

!***********************************************
