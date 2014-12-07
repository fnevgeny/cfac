C     calculate dirac coulomb function for positive and negtive energies
C     INPUT:
C     z Nuclear charge. positive for electron-ion system.
C     e energy.
C     k kappa
C     r radius
C     OUTPUT:
C     p large component of regular solution
C     q small component of regular solution
C     p1 iregular large for e > 0. ignore for e < 0
C     q1 iregular small for e > 0, ignore for e < 0
C     ierr error code returned by coulcc
      subroutine dcoul(z, e, k, r, p, q, p1, q1, ierr)
      implicit none
      integer k, ierr, kfn
      double precision z, e, r, p, q, p1, q1, c, ki, zp, gam
      double precision lambda, y, qi, x0, b1, b2, np
      complex*16 x, eta, zlmin, omega, a, pp, qq, mu, nu, IONE
      complex*16 fc(1), gc(1), fcp(1), gcp(1), sig(1), clogam, lam0
      double precision SL, SL2, TSL2, ALPHA
      parameter (SL=137.036D0,SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      real*8 HALFPI
      parameter (HALFPI = 1.5707963268D0)

      IONE = dcmplx(0.0, 1.0)
      c = 1.0+0.5*e/SL2
      ki = sqrt(2.0*abs(e)*c)
      zp = z*ALPHA
      gam = sqrt(k*k - zp*zp)
      lambda = gam - 0.5
      qi = sqrt(c/ki)
      y = (1.0+e/SL2)*z/ki

      x0 = ki*r
      if (e .lt. 0) then
         x = dcmplx(0.0, x0)
         eta = dcmplx(1D-5*(0.5+y), 0.5+y)
         mu = dcmplx(k - z/ki, 0.0)
         nu = dcmplx(0.5+y-x0, 0.0)
      else
         x = dcmplx(x0, 0.0)
         eta = dcmplx(-y, 0.5)
         mu = dcmplx(k, -z/ki)
         nu = IONE*(x - eta)
      endif

      zlmin = dcmplx(lambda, 0.0)
      kfn = 0

      ierr = 0
      call coulcc(x, eta, zlmin, 1, fc, gc, fcp, gcp, sig,
     +     11, kfn, ierr)
      if (ierr .ne. 0) then
          return
      endif
      if (e .lt. 0) then
         omega = IONE*(HALFPI*(lambda - y - 0.5) - sig(1))
         if (z > 0) then
            np = y - gam
            b1 = sqrt(z*c*(z/ki-k))*ki/z
            lam0 = np+1.0
            b2 = dble(clogam(lam0))
            lam0 = np + 1.0 + 2.0*gam
            b2 = b2 + dble(clogam(lam0))
            b2 = -0.5*b2
            omega = omega + b2
            a = exp(omega)*b1
            a = a/(mu*sqrt(2.0*x0))
         else
            a = exp(IONE*dimag(omega))
            a = a/mu
            a = a*qi/sqrt(2.0*x0)
         endif
         pp = a*((mu + nu)*gc(1) - x*gcp(1))
         qq = (ALPHA*e/ki)*a*((mu - nu)*gc(1) + x*gcp(1))
         p = dble(pp)
         q = dble(qq)
         p1 = dble(omega)
         q1 = dimag(omega)
      else
         a = qi/sqrt(2.0*IONE*mu*x)
         pp = a*((mu + nu)*gc(1) - x*gcp(1))
         qq = (IONE*e*ALPHA/ki)*a*((mu - nu)*gc(1) + x*gcp(1))
         p1 = dble(pp)
         q1 = dble(qq)
         p = dimag(pp)
         q = dimag(qq)
      endif

      end
