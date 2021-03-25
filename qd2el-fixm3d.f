      PROGRAM diag2el
C     ********************************
C     * 2-el. quantum dot in 3D      *
C     * full (CM+relative) motion    *
C     * fixed M = m1+m2 = m3+m4      *
C     * diagonalization of 1/|r1-r2| *
C     * in the 3D HO basis           *
C     ********************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SM(50,50),EN(50)
      real*8 meff,l0,lam,intm(50,50),laguerre
      COMPLEX W(50),EIGFUN(50,50)
      OPEN(1,FILE='enlev.dat',STATUS='UNKNOWN')
      OPEN(2,FILE='en-b.dat',STATUS='UNKNOWN')
C
      write(*,*)' om0[meV], omz/om0[2.5], m*[0.067], eps_r[12]:'
      read(*,*) om0mev,omz,meff,epsr
C      
      om0au = om0mev/27211.4d0
      l0 = 1.d0/sqrt(meff*om0au)
      cc = 1.d0/(epsr*l0*om0au)
      lam = dsqrt(2.d0)*cc
      write(*,111) l0,cc,lam
  111 format(' l0 =',f8.2,' a0, k =',f8.4,' [lam =',f8.4,']')
      write(*,*)'--------------------------------------------'
C
      WRITE(*,*)' M, Nmax[1], Nz_max[1], eps[1e-2(-3)]:'
      READ(*,*) M,nmax,nzmax,eps
      WRITE(*,*)' rho_max[5], Ngr1<50>:'
      READ(*,*) xmax, ngr1
      WRITE(*,*)' z_max[5], Ngr2<50>:'
      READ(*,*) ymax, ngr2
      WRITE(*,*)' oml_min, oml_max, noml:'
      READ(*,*) oml1,oml2,noml
C
      NN = (nzmax+1)**2*(nmax+1)**2*(M+3)
      write(*,*)' N =',NN
      dx = xmax/ngr1
      dy = ymax/ngr2
      step = (oml2-oml1)/noml   
C
C     ************************************************
C     Calculation of eigenenergies as functions of OmL
C     ************************************************
C
      DO ioml=0,noml
        oml = oml1+ioml*step
	om2 = oml*oml+1.d0
	om = dsqrt(om2)
C       ---------------------------------------------------- 
C       Interaction matrix elements: 
C       <n1,m1,nz1;n2,M-m1,nz2|1/r12|n3,m3,nz3;n4,M-m3,nz4>
C       ----------------------------------------------------
        do nz4=0,nzmax
        do nz3=0,nzmax
        do nz2=0,nzmax
        do nz1=0,nzmax
        do n4=0,nmax
        do n3=0,nmax
        do n2=0,nmax
        do n1=0,nmax
        do m3=-1,m+1
        do m1=-1,m+1
	i = m1+2 + (M+3)*(n1+(nmax+1)*n2)
     *    + (M+3)*(nmax+1)**2*(nz1+(nzmax+1)*nz2)
        j = m3+2 + (M+3)*(n3+(nmax+1)*n4)
     *    + (M+3)*(nmax+1)**2*(nz3+(nzmax+1)*nz4)
        m2 = M-m1
        m4 = M-m3
        call INTMEL(om,omz,n1,m1,nz1,n2,m2,nz2,n3,m3,nz3,n4,m4,nz4,
     *  ngr1,ngr2,eps,dx,dy,intm(i,j))
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
	enddo
	enddo
	enddo
	enddo
C       ---------------
C       Diagonalization
C       ---------------
        do i=1,NN
        do j=1,NN
          SM(i,j) = cc*intm(i,j)
        enddo
        enddo
C
        do nz2=0,nzmax
        do nz1=0,nzmax
	do n2=0,nmax
	do n1=0,nmax
        do m1=-1,M+1
	  i = m1+2 + (M+3)*(n1+(nmax+1)*n2)
     *    + (M+3)*(nmax+1)**2*(nz1+(nzmax+1)*nz2)
	  m2 = M-m1
	  E1 = om*(2*n1+abs(m1)+1) + omz*(nz1+.5d0)
	  E2 = om*(2*n2+abs(m2)+1) + omz*(nz2+.5d0)
	  E0 = E1 + E2 - oml*M
          SM(i,i) = SM(i,i) + E0
        enddo
	enddo
	enddo
        enddo
        enddo
        CALL EIG(SM,W,EIGFUN)
        DO I=1,NN
          EN(I) = REAL(W(I))
        ENDDO
        CALL ORDER(NN,EN)
        B = 17.27598*meff*om0mev*oml
	WRITE(*,*) oml,B
        WRITE(1,200) oml,(EN(I),I=1,20)	
        WRITE(2,200) B,(om0mev*EN(I),I=1,25)	
      ENDDO
C     ************************************************
  200 FORMAT(26(F10.4,1X))
      STOP
      END


      SUBROUTINE INTMEL(om,omz,n1,m1,nz1,n2,m2,nz2,n3,m3,nz3,n4,m4,nz4,
     *ngr1,ngr2,eps,dx,dy,ime)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      real*8 laguerre,ime
      pi = 3.1415926535898d0
      drho = dx/dsqrt(om)
      dz = dy/dsqrt(omz)
      ime = 0.d0
      IF (m1+m2.eq.m3+m4) THEN
        m1a = abs(m1)
        m2a = abs(m2)
        m3a = abs(m3)
        m4a = abs(m4)
        c1 = dsqrt(2.d0*om*fact(n1)/fact(n1+m1a))*
     *       dsqrt(dsqrt(omz/pi)/(2**nz1*fact(nz1)))
        c2 = dsqrt(2.d0*om*fact(n2)/fact(n2+m2a))*
     *       dsqrt(dsqrt(omz/pi)/(2**nz2*fact(nz2)))
        c3 = dsqrt(2.d0*om*fact(n3)/fact(n3+m3a))*
     *       dsqrt(dsqrt(omz/pi)/(2**nz3*fact(nz3)))
        c4 = dsqrt(2.d0*om*fact(n4)/fact(n4+m4a))*
     *       dsqrt(dsqrt(omz/pi)/(2**nz4*fact(nz4)))
	do k1=0,ngr1
	  x1 = k1*dx
	  rho1 = x1/dsqrt(om)
          do k2=0,ngr1
	    x2 = k2*dx+dx/2
	    rho2 = x2/dsqrt(om)
	    do k3=-ngr2,ngr2
	      y1 = k3*dy
	      z1 = y1/dsqrt(omz)
	      do k4=-ngr2,ngr2
	        y2 = k4*dy+dy/2
	        z2 = y2/dsqrt(omz)
                x1sq = x1*x1
	        x2sq = x2*x2
	        rwf1 = c1*pow(x1,m1a)*exp(-x1sq/2)*Laguerre(n1,m1a,x1sq)
     *                *exp(-y1*y1/2)*Hermite(nz1,y1)
	        rwf2 = c2*pow(x2,m2a)*exp(-x2sq/2)*Laguerre(n2,m2a,x2sq)
     *                *exp(-y2*y2/2)*Hermite(nz2,y2)
	        rwf3 = c3*pow(x1,m3a)*exp(-x1sq/2)*Laguerre(n3,m3a,x1sq)
     *                *exp(-y1*y1/2)*Hermite(nz3,y1)
                rwf4 = c4*pow(x2,m4a)*exp(-x2sq/2)*Laguerre(n4,m4a,x2sq)
     *                *exp(-y2*y2/2)*Hermite(nz4,y2)
                r1sq = rho1*rho1 + z1*z1
                r2sq = rho2*rho2 + z2*z2
                a = (2*rho1*rho2)/(r1sq+r2sq-2*z1*z2)
                call aint(m3-m1,a,eps,ai)
                angm = ai/dsqrt(r1sq+r2sq-2*z1*z2)    
                ime = ime + angm*
     *                rwf1*rwf2*rwf3*rwf4*rho1*rho2*drho*drho*dz*dz
              enddo
	    enddo
          enddo
        enddo
      ENDIF
      RETURN
      END


      FUNCTION POW(X,N)
      real*8 x,pow
      if (n.eq.0) then
        pow = 1.d0
      else
        pow = x**n
      endif
      RETURN
      END
      

      SUBROUTINE ORDER(NN,EN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION EN(50)      
      do i=1,nn
      do j=i+1,nn
        if (en(i).gt.en(j)) then
	  temp = en(i)
	  en(i) = en(j)
	  en(j) = temp
	endif
      enddo
      enddo
      RETURN
      END


      SUBROUTINE EIG(SM,W,Z)
      REAL*8 SM(50,50)
      COMPLEX W(50),Z(50,50)
      DIMENSION A(50,50),WK(200)
      DO 10 I=1,50
      DO 10 J=1,50
   10 A(I,J) = SM(I,J)
      IA = 50
      IZ = 50
      N = 50
      IJOB = 1
      CALL EIGRF(A,N,IA,IJOB,W,Z,IZ,WK,IER)
      RETURN
      END


      FUNCTION LAGUERRE(N,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      real*8 laguerre
      laguerre = fact(n+m)/(fact(n)*fact(m)) 
      xk = 1.d0
      do k=1,n
        c = fact(n+m)/(fact(n-k)*fact(m+k)*fact(k)) 
	xk = -xk*x
	laguerre = laguerre + c*xk
      enddo
      RETURN
      END
      

      FUNCTION HERMITE(N,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      if (n.eq.0) h = 1.d0
      if (n.eq.1) h = 2*x     
      if (n.eq.2) h = 4*x*x - 2.d0
      if (n.eq.3) h = 8*x*x*x - 12*x
      if (n.eq.4) h = 16*x*x*x*x - 48*x*x + 12.d0
      if (n.eq.5) h = 32*x**5 - 160*x*x*x + 120*x
      if (n.eq.6) h = 64*x**6 - 480*x*x*x*x + 720*x*x - 120.d0
      if (n.eq.7) h = 128*x**7 - 1344*x**5 + 3360*x*x*x - 1680*x
      if (n.eq.8) h = 256*x**8 - 3584*x**6 + 13440*x**4 - 13440*x*x + 
     *                1680.d0
      if (n.eq.9) h = 512*x**9 - 9216*x**7 + 48384*x**5 - 80640*x*x*x +
     *                30240*x      
      if (n.eq.10) h = 1024*x**10 - 23040*x**8 + 161280*x**6 -
     *                 403200*x*x*x*x + 302400*x*x - 30240.d0
      hermite = h
      RETURN
      END

      
      SUBROUTINE AINT(m,x,eps,ai)
      implicit double precision(a-h,o-z)
      ai = 0.d0
      k = (abs(m)+m)/2
      if (2*k.eq.m) then
        y = 1.d0
      else
        y = (x/4)**(2*k-m)
      endif
      term = fact2(4*k-2*m-1)/(fact(k)*fact(k-m))*y
   10 continue
	ai = ai + term
	if (term.lt.eps) RETURN
	k = k+1
	term = term*(dfloat(4*k-2*m-1)/dfloat(k))*
     *              (dfloat(4*k-2*m-3)/dfloat(k-m))*x*x/16
      go to 10
      END


      FUNCTION FACT(N)
      real*8 fact
      fact = 1.d0
      if (n.le.1) return
      do i = 1,n
        fact = fact*i
      enddo
      RETURN
      END


      FUNCTION FACT2(N)
      real*8 fact2
      fact2 = 1.d0
      if (n.le.1) return
      do i = 1,n,2
        fact2 = fact2*i
      enddo
      RETURN
      END
