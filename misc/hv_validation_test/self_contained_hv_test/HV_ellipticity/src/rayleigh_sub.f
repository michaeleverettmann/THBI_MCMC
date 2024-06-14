c
c     rayleigh wave eigenvalue (r-k-g version)
c
      subroutine  raydsp(diffeq,h,rho,vp,vs,ap,ae,l,
     *                   w,cmn,cmx,dc,tol,itr,ia,c,u,ek,
     1                   y0,yij,ier)
c
c     input
c       diffeq : subroutine to integrate eq. of motion
c       h    : layer thickness.  h(l) is arbitrary
c       rho  : density
c       vp   : compressional wave velocity
c       vs   : shear wave velocity
c       ap   : anisotropy factor phi
c       ae   : anisotropy factor eta
c       l    : number of layers including the bottom
c              half-space
c       w    : angular frequency
c       cmn  : lower limit of phase velocity
c       cmx  : upper limit of phase velocity
c       dc   : increment   of phase velocity
c       tol  : relative accuracy of phase velocity
c       itr  : maximum number of iterations
c       ia   : = 0 ;   isotropic model
c              = 1 ; anisotropic model
c     output
c       c    : phase velocity
c       u    : group velocity by differentiation
c       ekd  : energy integral by differentiation
c              2*k**2*i3
c       y0   : surface values of eigenfunction
c              y0(1) ; y1 (scale factor)
c              y0(2) ; y2/abs(y1), dispersion function
c              y0(3) ; y3/y1
c       yij  : surface values of yij (compounds),
c              c*dyij/dc, and w*dyij/dw
c       ier  : return code
c              < 0 ; input error
c              = 0 ; no error
c              = 1 ; slow convergence
c              = 2 ; root not found
c
c     subroutine : diffeq (raymrx or rayrkg)
c
c     disper-80,  ver-86
c
c     m. saito  23/vi/79
c     revised   25/ix/85
c
      implicit real*8 (a-h,o-z)
      dimension  h(l),rho(l),vp(l),vs(l),ap(l),ae(l),
     *           y0(3),yij(15)
c
c     initialization
c
      if( l.le.2 .or. w.le.0 .or. dc.eq.0 )  go to  90
      one = 1
c
c     machine epsilon
c
      eps = 1
    1 if( (1+eps).le.1 )  go to  2
        eps = eps/2
        go to  1
    2 eps = eps/2
c
      tol1 = tol
      if( tol1.gt.0 )  tol1 = max(tol1,eps)
      ier1 = 0
c      write(6,3)  w
c    3   format(/7x,'w',15x,'rayleigh wave'/
c     *         1pe18.6//7x,'c',17x,'y2',16x,'y3')
c     *         1pd18.6//7x,'c',17x,'y2',16x,'y3')
      c3 = cmn
      if( c3.le.0 )  go to  92
      call  diffeq(h,rho,vp,vs,ap,ae,l,w,c3,ia,1,u,
     *             ek,y0,yij,ier)
      if( ier.ne.0 )  return
      f3 = y0(2)
c      write(6,4)  c3,f3,y0(3)
c    4   format(1p3e18.6)
c    4   format(1p3d18.6)
      if( f3.eq.0 .and. tol1.gt.0 )  go to  9
c
c     find a zero-cross
c
      kx = (cmx - cmn)/dc + 0.5
      kx = max(1,kx)
      do  5  k=1,kx
        cc = cmn + k*dc
        c1 = c3
        f1 = f3
        c3 = cc
        if( c3.le.0 )  go to  92
        call  diffeq(h,rho,vp,vs,ap,ae,l,w,c3,ia,1,u,
     *               ek,y0,yij,ier)
        if( ier.ne.0 )  return
        f3 = y0(2)
c        write(6,4)  c3,f3,y0(3)
        if(  tol1.le.0 )  go to  5
        if( f3*sign(one,f1).le.0 )  go to  6
    5 continue
c
      ier = 2
      if( tol1.le.0 )  return
      go to  94
c
c     interpolation
c
    6 if( f3.eq.0 )  go to  9
      c2 = c3
      f2 = f3
      e  = c1 - c2
      d  = e/2
      c3 = c2 + d
      kx = max(1,itr)
c
      do  7  k=1,kx
        call  diffeq(h,rho,vp,vs,ap,ae,l,w,c3,ia,1,u,
     *               ek,y0,yij,ier)
        f3 = y0(2)
c        write(6,4) c3,f3,y0(3)
        if( f3*sign(one,f2).gt.0 )  then
          ff = c1
          c1 = c2
          c2 = ff
          ff = f1
          f1 = f2
          f2 = ff
        end if
        if( abs(f3).gt.abs(f2) )  then
          ff = c2
          c2 = c3
          c3 = ff
          ff = f2
          f2 = f3
          f3 = ff
        end if
        e = c2 - c3
        if( f3.eq.0 )  go to  9
        tolc = c3*tol1
        dd = d
        f32 = f3/f2
        f31 = f3/f1
        f21 = f2/f1
        q = f32*(e*(1 - f31) + f21*(f31 - f21)*(c1 - c3))
        s = (f21 - 1)*(f32 - 1)*(f31 - 1)
c
c       test range
c
        if( q.lt.0 )  s =-s
        q = abs(q)
        if( q.ge.(e*s-abs(tolc*s)) )  then
c
c       linear interpolation
c
          d = e*f32/(f32 - 1)
c
c       inverse quadratic interpolation
c
        else
          d = q/s
        end if
c
c       test convergence
c
        c1 = c2
        f1 = f2
        c2 = c3
        f2 = f3
        c3 = c2 + d
        if( abs(e).le.tolc )  go to  9
        if( abs(d).le.tolc )  then
c
c       bisection
c
          if( abs(dd).le.tolc )  then
            d = e/2
            c3 = c2 + d
          else
            c3 = c2 + sign(tolc,d)
          end if
        end if
c
    7 continue
c
c     slow convergence
c
      write(6,8)  kx
    8   format(20x,5('?'),3x,'(raydsp)   slow conv. ',
     *         'after',i5,' iterations',3x,5('?'))
      ier = 1
c
c     root is found
c
    9 call  diffeq(h,rho,vp,vs,ap,ae,l,w,c3,ia,3,u,
     *             ek,y0,yij,ier1)
c      write(6,4)  c3,y0(2),y0(3)
      c   = c3
      return
c
c     input error
c
   90 write(6,91)  l,w,dc
   91   format(20x,5('?'),3x,'(raydsp)   input error',3x,
     *         'l  =',i5,3x,'w  =',1pe13.6,3x,
     1         'dc =',e13.6,3x,5('?'))
      ier =-1
      return
c
   92 write(6,93)  c3
   93   format(20x,5('?'),3x,'(raydsp)   input error',3x,
     *         'c  =',1pe13.6,3x,5('?'))
      ier =-1
      return
c
c     no root
c
   94 write(6,95)
   95   format(20x,5('?'),3x,'(raydsp)',3x,
     *         'root not found',3x,5('?'))
      return
      end
c
c     rayleigh wave runge-kutta-gill integration
c
      subroutine  rayrkg(h,rho,vp,vs,ap,ae,l,w,c,ia,ig,u,
     *                   ek,y0,yij,ier)
c
c     input
c       h    : grid interval.  h(l) is arbitrary
c              every two successive h must be equal
c       rho  : density
c       vp   : p wave velocity
c       vs   : s wave velocity
c       ap   : anisotropy factor phi
c       ae   : anisotropy factor eta
c       l    : number of grid points
c       w    : angular frequency
c       c    : phase velocity
c       ia   : = 0 ;   isotropic model
c              = 1 ; anisotropic model
c       ig   : = 1 ; to compute y only
c              = 2 ; to compute y and c*cy/dc
c              = 3 ; to compute y, c*dy/dc, and w*dy/dw
c     output
c       u    : group velocity by differentiation (ig=3)
c       ek   : energy integral 2*k**2*i3 by
c              differentiation  (ig>=2)
c       y0   : surface values of eigenfunctions
c              y0(1) ; y1 (scale factor)
c              y0(2) ; y2/abs(y1), dispersion function
c              y0(3) ; y3/y1
c       yij  : surface values of yij (compounds),
c              c*dyij/dc, and w*dyij/dw
c              for solid layer
c                (ij) = (12),(13),(14),(23),(24)
c       ier  : return code
c              < 0 ; invalid input
c              = 0 ; no error
c
c     disper-80,  ver-86
c
c     m. saito  30/vii/79
c     revised   10/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  h(l),rho(l),vp(l),vs(l),ap(l),ae(l),
     *           y0(3),yij(15),
     1           yy(15),dy(15),qq(15),yyy(15),qqq(15)
      data  eps/1.e-30/,big/1.e+30/
     *     ,rkgp/1.707107/,rkgm/0.2928932/
c
      rkgp = 1.d0 + dsqrt(0.5d0)
      rkgm = 0.5d0/rkgp
c
      if( l.le.2 .or. w.le.0 .or. c.le.0 )  go to  90
      i   = l
      if( rho(i).le.0. .or. c.ge.vp(i) )  go to  92
      ier = 0
      ih  = 0
      cc  = c*c
      wn  = w/c
      igg = max(1,min(3,ig))
c
c     initial value (isotropic model)
c
      ro  = rho(i)
      roc = ro*cc
      sv  = vs(i)
      el  = ro*sv**2
      cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
      do  1  j=1,15
        yy(j) = 0
        qq(j) = 0
    1 continue
c
c     liquid bottom
c
      if( sv.le.0 )  then
        jx = 2*igg
        yy(1) = eps*ra
        yy(2) =-eps*roc
        yy(3) =-eps*cp**2/ra
        yy(4) =-eps*roc*2
c
c     solid bottom
c
      else
        if( c.ge.sv )  go to  92
        jx  = 5*igg
        cs  = c/sv
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = 2*el
        yy(3) =-eps*ra
        yy(4) =-eps*rb
        yy(2) =-eps*(cp**2*rbb + cs**2)/(roc*(ra*rb + 1))
        yy(1) = rg*yy(2) + eps
        yy(5) =-rg*(yy(1) + eps) + eps*roc
c
        if( igg.ge.2 )  then
          yy( 8) = eps*cp**2/ra
          yy( 9) = eps*cs**2/rb
          yy( 7) =-(rb*yy(8) + ra*yy(9))/roc - yy(2)*2
          yy( 6) = rg*yy(7)
          yy(10) =-rg*yy(6) + eps*roc*2
        endif
      endif
c
c     integrate upward
c
    2 if( h(i-1).gt.0 )  go to  6
      i = i - 1
      if( i.le.1 )  go to  92
      if( (vs(i).gt.0 .and. el.gt.0) .or.
     *    (vs(i).le.0 .and. el.le.0) )  go to  6
      do  3  j=1,15
        yyy(j) = yy(j)
        qqq(j) = qq(j)
         yy(j) = 0
         qq(j) = 0
    3 continue
      jx1 = 0
c
c     solid  -----> liquid boundary
c
      if( el.gt.0 )  then
        y0(3) =-yyy(1)/yyy(3)
        jx = 2*igg
        do  4  j=1,jx,2
          yy(j  ) = yyy(jx1+3)
          qq(j  ) = qqq(jx1+3)
          yy(j+1) = yyy(jx1+5)
          qq(j+1) = qqq(jx1+5)
          jx1 = jx1 + 5
    4   continue
c
c     liquid -----> solid boundary
c
      else
        jx = 5*igg
        do  5  j=1,jx,5
          yy(j+1) = yyy(jx1+1)
          qq(j+1) = qqq(jx1+1)
          yy(j+3) = yyy(jx1+2)
          qq(j+3) = qqq(jx1+2)
          jx1 = jx1 + 2
    5   continue
      endif
c
    6 hh  = wn*(h(i-1) + h(i-2))
      h5  = hh/2
      jmp = 0
c
c     derivatives  (in-line subroutine)
c
    7 if( i.le.0  )  go to  92
      if( i.ne.ih )  then
        ih = i
        ro = rho(i)
        roc= ro*cc
        ea = ro*vp(i)**2
        if( ea.le.0 )  go to  92
        el = ro*vs(i)**2
        ec = ea
        ef = ea - el-el
        if( ia.gt.0 )  then
          ec = ec*ap(i)
          ef = ef*ae(i)
        endif
c
        if( el.le.0 )  then
          a12 = 1/ea - 1/roc
          if( igg.ge.2 )  c12 =-a12 + 2/roc
c
        else
          a12 = 1/ec
          a13 = ef*a12
          b31 = a13+a13
          a34 = 1/el
          a43 =-roc + ea - ef*a13
          if( igg.ge.2 )  c43 =-a43 - roc-roc
        endif
      endif
c
c     liquid layer
c
    8 if( el.le.0 )  then
        dy(1) = a12*yy(2)
        dy(2) =-roc*yy(1)
c
        if( igg.ge.2 )  then
          dy(3) = a12*yy(4) + c12*yy(2)
          dy(4) =-roc*(yy(3) + yy(1))
c
          if( igg.gt.2 )  then
            dy(5) = a12*(yy(6) + yy(2))
            dy(6) =-roc*(yy(5) + yy(1))
          endif
        endif
c
c     solid layer
c
      else
        dy( 1) =     yy(3) - a13*yy(4)
        dy( 2) = a34*yy(3) + a12*yy(4)
        dy( 3) =-b31*yy(1) + a43*yy(2) + a12*yy(5)
        dy( 4) =   2*yy(1) - roc*yy(2) + a34*yy(5)
        dy( 5) =-roc*yy(3) + a43*yy(4)
c
        if( igg.ge.2 )  then
          dy( 6) = yy(8) - yy(3) - a13*(yy(9) - yy(4))
          dy( 7) = a34*(yy(8) - yy(3))
     *           + a12*(yy(9) - yy(4))
          dy( 8) =-b31*(yy(6) - yy(1))
     *           + a43*yy(7) + c43*yy(2)
     1           + a12*(yy(10) - yy(5))
          dy( 9) =(yy(6) - yy(1))*2 - roc*(yy(7) + yy(2))
     *           + a34*(yy(10) - yy(5))
          dy(10) =-roc*(yy(8) + yy(3)) + a43*yy(9)
     *           + c43*yy(4)
c
          if( igg.gt.2 )  then
            dy(11) = yy(13) + yy(3) - a13*(yy(14) + yy(4))
            dy(12) = a34*(yy(13) + yy(3))
     *             + a12*(yy(14) + yy(4))
            dy(13) =-b31*(yy(11) + yy(1))
     *             + a43*(yy(12) + yy(2))
     1             + a12*(yy(15) + yy(5))
            dy(14) =(yy(11) + yy(1))*2
     *             - roc*(yy(12) + yy(2))
     1             + a34*(yy(15) + yy(5))
            dy(15) =-roc*(yy(13) + yy(3))
     *             + a43*(yy(14) + yy(4))
          endif
        endif
      endif
c
      jmp = jmp + 1
      go to  (9,11,13,15),jmp
c
c     runge-kutta-gill
c
    9 do  10  j=1,jx
        wk = h5*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)
        qq(j) = wq + (yy(j) - wy)*3 - wk
   10 continue
      i = i - 1
      go to  7
c
   11 do  12  j=1,jx
        wk = hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgm
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgm
   12 continue
      go to  8
c
   13 do  14  j=1,jx
        wk = hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgp
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgp
   14 continue
      i = i - 1
      go to  7
c
   15 do  16 j=1,jx
        wk = h5*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)/3
        qq(j) = wq + (yy(j) - wy)*3 - wk
   16 continue
c
c     normalization
c
      z1 = 0
      do  17  j=1,jx
   17   z1 = max(z1, abs(yy(j)))
      if( z1.gt.big )  then
        do  18  j=1,jx
          yy(j) = eps*yy(j)
          qq(j) = eps*qq(j)
   18   continue
      endif
c
      if( i.gt.2 )  go to  2
      if( i.ne.1 )  go to  92
c
c     normal exit
c
      do  19  j=1,jx
   19   yij(j) = yy(j)
c
c     liquid surface
c
      if( el.le.0 )  then
        y0(1) = yy(1)
        if( abs(yy(2))*eps.le.abs(yy(1)) )  then
          y0(2) = yy(2)/abs(yy(1))
        else
          y0(2) = sign(big,yy(2))
        endif
        if( igg.ge.2 )  then
          ek =-wn*yy(4)/yy(1)
          if( igg.gt.2 )  u = c*yy(4)/(yy(4) + yy(6))
        endif
c
c     solid surface
c
      else
        y0(1) = yy(3)
        if( abs(yy(5))*eps.le.abs(yy(3)) )  then
          y0(2) = yy(5)/abs(yy(3))
          y0(3) =-yy(1)/yy(3)
        else
          y0(2) = sign(big,yy(5))
        endif
        if( igg.ge.2 )  then
          ek =-wn*yy(10)/yy(3)
          if( igg.gt.2 )  u = c*yy(10)/(yy(10) + yy(15))
        endif
      endif
c
      return
c
c     input error
c
   90 write(6,91)  l,w,c
   91   format(20x,5('?'),3x,'(rayrkg)   input error',3x,
     *         'l  =',i5,3x,'wn =',1pe13.6,3x,
     1         'c  =',e13.6,3x,5('?'))
      ier =-1
      return
c
   92 write(6,93)  i
   93   format(20x,5('?'),3x,'(rayrkg)   input error',
     *         'at layer',i5,3x,5('?'))
      ier =-max(1,i)
      return
      end
c
c     rayleigh wave matrix method integration
c
      subroutine  raymrx(h,rho,vp,vs,ap,ae,l,w,c,
     *                   ia,ig,u,ek,y0,y,ier)
c
c     input
c       h    : layer thickness.  h(l) is arbitrary
c       rho  : density
c       vp   : p wave velocity
c       vs   : s wave velocity
c       l    : number of layers inclding the bottom
c              half-space
c       w    : angular frequency
c       c    : phase velocity
c       ig   : = 1 ; to compute y only
c              = 2 ; to compute y and c*cy/dc
c              = 3 ; to compute y, c*dy/dc, and w*dy/dw
c     output
c       u    : group velocity by differentiation (ig=3)
c       ek   : energy integral 2*k**2*i3  (ig>=2)
c       y0   : surface values of eigenfunctions
c              y0(1) ; y1 (scale factor)
c              y0(2) ; y2/abs(y1), dispersion function
c              y0(3) ; y3/y1
c       y    : surface values of yij (compounds),
c              c*dyij/dc, and w*dyij/dw
c              for solid layer
c                (ij) = (12),(13),(14),(23),(24)
c       ier  : return code
c              < 0 ; input error
c              = 0 ; no error
c
c     isotropic model only.  ap, ae, ia are dummy.
c
c     disper-80,  ver-86
c
c     m. saito  30/vii/79
c     revised   10/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  h(l),rho(l),vp(l),vs(l),y0(3),y(15),z(15)
     *          ,ap(1),ae(1)
      data  eps/1.e-30/,big/1.e+30/
c
c     define sinh(x)/x and (cosh(x)-sinh(x)/x)/x**2
c
c      sh0(x) = 0.9999997 + x*(0.1666667 + x*(0.0083361
c     *       + x*0.0001984))
c      sh1(x) = 0.3333333 + x*(0.0333333 + x*(0.0011907
c     *       + x*0.0000220))
c
c     for double precision use
c 
      sh0(x) = 1.0d0 + x*(1.6666 66666 66666 7d-1
     *       + x*(8.33 33333 33334 0d-3
     1       + x*(1.9 84126 98412 7d-4
     2       + x*(2.7557 31918 9d-6 + x*(2.50 12108 4d-8
     3       + x*(1.60596 1d-10 + x*7.64 7d-13))))))
      sh1(x) = 3.3333 33333 33333 3d-1
     *       + x*(3.333 33333 33333 3d-2
     1       + x*(1.19 04761 90476 2d-3
     2       + x*(2.20458 55379 2d-5
     3       + x*(2.505 21083 7d-7 + x*(1.9 27085 3d-9
     4       + x*(1.0706 3d-11 + x*4.50d-14))))))
c 
c     initial value
c
      if( l.le.0 .or. w.le.0 .or. c.le.0 )  go to  90
      i   = l
      ro  = rho(i)
      if( ro.le.0 .or. c.ge.vp(i) )  go to  92
      ier = 0
      cc  = c*c
      wn  = w/c
      igg = max(1,min(3,ig))
      roc = ro*cc
      sv  = vs(i)
      cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
      do  1  j=1,15
    1   y(j) = 0
c
c     liquid bottom
c
      if( sv.le.0 )  then
        y(1) = ra*eps
        y(2) =-roc*eps
        jx = 2
c
        if( igg.ge.2 )  then
          y(3) =-cp**2*eps/ra
          y(4) =-2*roc*eps
          jx = 4
c
          if( igg.gt.2 )  jx = 6
        endif
c
c     solid bottom
c
      else
        if( c.ge.sv )  go to  92
        cs  = c/sv
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = 2*ro*vs(i)**2
        y(3) =-ra*eps
        y(4) =-rb*eps
        y(2) =-eps*(cp**2*rbb + cs**2)/(roc*(ra*rb + 1))
        y(1) = rg*y(2) + eps
        y(5) =-rg*(y(1) + eps) + roc*eps
        jx = 5
c
        if( igg.ge.2 )  then
          y(8) = eps*cp**2/ra
          y(9) = eps*cs**2/rb
          y(7) =-(rb*y(8) + ra*y(9))/roc - 2*y(2)
          y(6) = rg*y(7)
          y(10)=-rg*y(6) + eps*roc*2
          jx = 10
c
          if( igg.gt.2 )  jx = 15
        endif
      endif
c
c     integrate upward
c
      if( l.le.1 )  go to  8
      do  7  ii=2,l
        i  = i - 1
        ro = rho(i)
        roc= ro*cc
        pv = vp(i)
        sv = vs(i)
        if( pv.le.0 .or. ro.le.0 .or.
     *      h(i).le.0 )  go to  92
        do  2  j=1,15
    2     z(j) = y(j)
        if( (sv.le.0 .and. vs(i+1).le.0) .or.
     *      (sv.gt.0 .and. vs(i+1).gt.0) )  go to  3
c
c       solid  -----> liquid boundary
c
        if( sv.le.0 )  then
          y0(3) =-y(1)/y(3)
          z(1) = y(3)
          z(2) = y(5)
          jx = 2
c
          if( igg.ge.2 )  then
            z(3) = y(8)
            z(4) = y(10)
            jx = 4
c
            if( igg.gt.2 )  then
              z(5) = y(13)
              z(6) = y(15)
              jx = 6
            endif
          endif
c
c       liquid -----> solid boundary
c
        else
          z( 2) = y(1)
          z( 4) = y(2)
          z( 1) = 0
          z( 3) = 0
          z( 5) = 0
          jx = 5
c
          if( igg.ge.2 )  then
            z( 7) = y(3)
            z( 9) = y(4)
            z( 6) = 0
            z( 8) = 0
            z(10) = 0
            jx = 10
c
            if( igg.gt.2 )  then
              z(12) = y(5)
              z(14) = y(6)
              z(11) = 0
              z(13) = 0
              z(15) = 0
              jx = 15
            endif
          endif
        endif
c
    3   r2  = 1/roc
        cp  = c/pv
        raa = (1 + cp)*(1 - cp)
        cs  = 0
        if( sv.gt.0 )  cs = c/sv
        rbb = (1 + cs)*(1 - cs)
        hk  = h(i)*wn
        hkk = hk**2
        xx  = raa*hkk
        one = 1
c
c       sinh(x)/x
c
        do  4  k=1,2
          cha = chb
          sha = shb
          dha = dhb
          aa  = abs(xx)
          if( aa.le.1 )  then
            shb = sh0(xx)
            chb = 1 + xx*sh0(xx/4)**2/2
            if( igg.ge.2 )  dhb = sh1(xx)*hkk
          else
            aa  = sqrt(aa)
            if( xx.le.0 )  then
              chb = cos(aa)
              shb = sin(aa)/aa
            else
              if( aa.gt.100 )  one = 0
              if( aa.le.100 )  one = one/cosh(aa)
              chb = 1
              shb = tanh(aa)/aa
            endif
            if( igg.ge.2 )  dhb = (hkk/xx)*(chb - shb)
          endif
          xx  = hkk*rbb
          shb = hk*shb
    4   continue
c
c     layer matrices
c
c       liquid layer
c
        if( sv.le.0 )  then
          b11 = cha
          b12 =-raa*sha*r2
          b21 =-roc*sha
c
          y(1) = b11*z(1) + b12*z(2)
          y(2) = b21*z(1) + b11*z(2)
c
          if( igg.ge.2 )  then
            c11 =-hk*sha
            c12 = (hk*cha + (1 + raa)*sha)*r2
            c21 = (hk*dha - sha)*roc
c
            y(3) = b11*z(3) + b12*z(4) + c11*z(1)
     *           + c12*z(2)
            y(4) = b21*z(3) + b11*z(4) + c21*z(1)
     *           + c11*z(2)
c
            if( igg.gt.2 )  then
              w11 = hk*raa*sha
              w12 =-hk*raa*cha*r2
              w21 =-hk*roc*cha
c
              y(5) = b11*z(5) + b12*z(6) + w11*z(1)
     *             + w12*z(2)
              y(6) = b21*z(5) + b11*z(6) + w21*z(1)
     *             + w11*z(2)
            endif
          endif
c
c       solid layer
c
        else
          g1  = 2/cs**2
          rg  = g1*roc
          r4  = rg - roc
          e1  = cha*chb
          e2  = e1 - one
          e3  = sha*shb
          e5  = sha*chb
          e6  = shb*cha
          f1  = e2 - e3
          f2  = r2*f1
          f3  = g1*f1 + e3
          b33 = e1
          b34 = raa*e3
          b43 = rbb*e3
          b25 =-r2*(f2 + r2*(e2 - raa*b43))
          b15 = rg*b25 + f2
          b16 =-rg*b15 - f3
          b22 = b16 + e1
          b12 = rg*b16 - r4*f3
          b52 =-rg*b12 + r4*(rg*f3 + r4*e3)
          b23 = r2*(e5 - rbb*e6)
          b13 = rg*b23 - e5
          b42 =-rg*b13 + r4*e5
          b24 = r2*(e6 - raa*e5)
          b14 = rg*b24 - e6
          b32 =-rg*b14 + r4*e6
          b11 = one - b16-b16
          b21 = b15+b15
          b31 = b14+b14
          b41 = b13+b13
          b51 = b12+b12
c
          y(1) = b11*z(1) + b12*z(2) + b13*z(3)
     *         + b14*z(4) + b15*z(5)
          y(2) = b21*z(1) + b22*z(2) + b23*z(3)
     *         + b24*z(4) + b25*z(5)
          y(3) = b31*z(1) + b32*z(2) + b33*z(3)
     *         + b34*z(4) + b24*z(5)
          y(4) = b41*z(1) + b42*z(2) + b43*z(3)
     *         + b33*z(4) + b23*z(5)
          y(5) = b51*z(1) + b52*z(2) + b42*z(3)
     *         + b32*z(4) + b22*z(5)
c
          if( igg.ge.2 )  then
            raac =-2*cp*cp
            rbbc =-2*cs*cs
            r1c = roc+roc
            e1c =-hk*(e5 + e6)
            e3c =-e3-e3 - hk*(dha*shb + dhb*sha)
            e5c =-e5 - hk*(dha*chb + e3)
            e6c =-e6 - hk*(dhb*cha + e3)
            f1c = e1c - e3c
            f2c = r2*(f1c - f1-f1)
            f3c = g1*(f1c - f1-f1) + e3c
            c33 = e1c
            c34 = raa*e3c + raac*e3
            c43 = rbb*e3c + rbbc*e3
            c25 =-r2*(f2c + r2*(e1c - raa*c43 - raac*b43))
     *          - 2*(b25+b25 + r2*f2)
            c15 = rg*c25 + f2c
            c16 =-rg*c15 - f3c
            c22 = c16 + e1c
            c12 = rg*c16 + r1c*f3 - r4*f3c
            c52 =-rg*c12 + r4*(rg*f3c + r4*e3c)
     *          - r1c*(rg*f3 + 2*r4*e3)
            c23 = r2*(e5c - rbb*e6c - rbbc*e6) - b23-b23
            c13 = rg*c23 - e5c
            c42 =-rg*c13 + r4*e5c - r1c*e5
            c24 = r2*(e6c - raa*e5c - raac*e5) - b24-b24
            c14 = rg*c24 - e6c
            c32 =-rg*c14 + r4*e6c - r1c*e6
            c11 =-c16-c16
            c21 = c15+c15
            c31 = c14+c14
            c41 = c13+c13
            c51 = c12+c12
c
            y( 6) = b11*z(6) + b12*z(7) + b13*z(8)
     *            + b14*z(9) + b15*z(10)
     *            + c11*z(1) + c12*z(2) + c13*z(3)
     *            + c14*z(4) + c15*z(5)
            y( 7) = b21*z(6) + b22*z(7) + b23*z(8)
     *            + b24*z(9) + b25*z(10)
     *            + c21*z(1) + c22*z(2) + c23*z(3)
     *            + c24*z(4) + c25*z(5)
            y( 8) = b31*z(6) + b32*z(7) + b33*z(8)
     *            + b34*z(9) + b24*z(10)
     *            + c31*z(1) + c32*z(2) + c33*z(3)
     *            + c34*z(4) + c24*z(5)
            y( 9) = b41*z(6) + b42*z(7) + b43*z(8)
     *            + b33*z(9) + b23*z(10)
     *            + c41*z(1) + c42*z(2) + c43*z(3)
     *            + c33*z(4) + c23*z(5)
            y(10) = b51*z(6) + b52*z(7) + b42*z(8)
     *            + b32*z(9) + b22*z(10)
     *            + c51*z(1) + c52*z(2) + c42*z(3)
     *            + c32*z(4) + c22*z(5)
c
            if( igg.gt.2 )  then
              e1w = hk*(raa*e5 + rbb*e6)
              e3w = hk*(e5 + e6)
              e5w = hk*(e1 + b43)
              e6w = hk*(e1 + b34)
              f1w = e1w - e3w
              f2w = r2*f1w
              f3w = g1*f1w + e3w
              w33 = e1w
              w34 = raa*e3w
              w43 = rbb*e3w
              w25 =-r2*(f2w + r2*(e1w - raa*w43))
              w15 = rg*w25 + f2w
              w16 =-rg*w15 - f3w
              w22 = w16 + e1w
              w12 = rg*w16 - r4*f3w
              w52 =-rg*w12 + r4*(rg*f3w + r4*e3w)
              w23 = r2*(e5w - rbb*e6w)
              w13 = rg*w23 - e5w
              w42 =-rg*w13 + r4*e5w
              w24 = r2*(e6w - raa*e5w)
              w14 = rg*w24 - e6w
              w32 =-rg*w14 + r4*e6w
              w11 =-w16-w16
              w21 = w15+w15
              w31 = w14+w14
              w41 = w13+w13
              w51 = w12+w12
c
              y(11) = b11*z(11) + b12*z(12) + b13*z(13)
     *              + b14*z(14) + b15*z(15)
     1              + w11*z( 1) + w12*z( 2) + w13*z( 3)
     2              + w14*z( 4) + w15*z( 5)
              y(12) = b21*z(11) + b22*z(12) + b23*z(13)
     *              + b24*z(14) + b25*z(15)
     1              + w21*z( 1) + w22*z( 2) + w23*z( 3)
     2              + w24*z( 4) + w25*z( 5)
              y(13) = b31*z(11) + b32*z(12) + b33*z(13)
     *              + b34*z(14) + b24*z(15)
     1              + w31*z( 1) + w32*z( 2) + w33*z( 3)
     2              + w34*z( 4) + w24*z( 5)
              y(14) = b41*z(11) + b42*z(12) + b43*z(13)
     *              + b33*z(14) + b23*z(15)
     1              + w41*z( 1) + w42*z( 2) + w43*z( 3)
     2              + w33*z( 4) + w23*z( 5)
              y(15) = b51*z(11) + b52*z(12) + b42*z(13)
     *              + b32*z(14) + b22*z(15)
     1              + w51*z( 1) + w52*z( 2) + w42*z( 3)
     2              + w32*z( 4) + w22*z( 5)
            endif
          endif
        endif
c
c     normalization
c
        z1 = 0
        do  5  j=1,jx
    5     z1 = max(z1,abs(y(j)))
        if( z1.gt.big )  then
          do  6  j=1,jx
    6       y(j) = eps*y(j)
        endif
c
    7 continue
c
c     normal exit
c
c     liquid surface
c
    8 if( sv.le.0 )  then
        y0(1) = y(1)
        if( abs(y(2))*eps.le.abs(y(1)) )  then
          y0(2) = y(2)/abs(y(1))
        else
          y0(2) = sign(big,y(2))
        endif
c
        if( igg.ge.2 )  then
          ek =-wn*y(4)/y(1)
          if( igg.gt.2 )  u = c*y(4)/(y(4) + y(6))
        endif
c
c     solid surface
c
      else
        y0(1) = y(3)
        if( abs(y(5))*eps.le.abs(y(3)) )  then
          y0(2) = y(5)/abs(y(3))
          y0(3) =-y(1)/y(3)
        else
          y0(2) = sign(big,y(5))
        endif
c
        if( igg.ge.2 )  then
          ek =-wn*y(10)/y(3)
          if( igg.gt.2 )  u = c*y(10)/(y(10) + y(15))
        endif
      endif
      return
c
c     input error
c
   90 write(6,91)  l,w,c
   91   format(20x,5('?'),3x,'(raymrx)   input error',3x,
     *         'l  =',i5,3x,'w  =',1pe13.6,3x,
     1         'c  =',e13.6,3x,5('?'))
      ier =-1
c
c     dummy statements
c
      ap(1) = ap(1)
      ae(1) = ae(1)
      ia = ia
      return
c
   92 write(6,93)  i
   93   format(20x,5('?'),3x,'(raymrx)   input error ',
     *         'at layer',i5,3x,5('?'))
      ier =-max(1,i)
      return
      end
c
c     rayleigh wave eigenfunction and partials by
c        runge-kutta-gill integration method
c
      subroutine  rayefr(h,rho,vp,vs,ap,ae,qp,qs,l,ly,
     *                   w,c,u,ia,iq,q,ekd,eki,er,
     1                   iy,y,ip,p,ny,ws,ier)
c
c     input
c       h    : grid interval.  h(l) is arbitrary
c              every 4 h's must be equal
c       rho  : density
c       vp   : p-wave velocity
c       vs   : s-wave velocity
c       ap   : anisotropy factor phi
c       ae   : anisotropy factor eta
c       qp   : inverse of p-wave q
c       qs   : inverse of s-wave q
c       l    : number of grid points
c       ly   : deepest coordinate point where
c              eigenfunction is to be computed
c       w    : angular frequency
c       c    : phase velocity
c       ekd  : energy integral 2*k**2*i3 computed by
c              rayrkg.  for normalization only
c       ia   : = 0 ;   isotropic model
c              = 1 ; anisotropic model
c       iy   : = 0 ; not to compute eigenfunction
c              = 1 ;     to compute eigenfunction
c       iq   : = 0 ; not to compute attenuation
c              = 1 ;     to compute attenuation
c       ip   : = 0 ; not to compute partials
c              = 1 ;     to compute partials
c     output
c       u    : group velocity by energy integrals
c       q    : inverse of q
c       eki  : energy integral 2*k**2*i3 by integration
c       er   : relative error in energy integrals
c       y    : eigenfunction
c       p    : partial derivatives
c              p(1,j) ; dc/drho
c              p(2,j) ; dc/dvp
c              p(3,j) ; dc/dvs
c              p(4,j) ; dc/dphi
c              p(5,j) ; dc/deta
c       ny   : length of y, p
c       ier  : return code
c              < 0 ; input error
c              = 0 ; no error
c     working spaces
c       ws   : 5*(l - no. of discon)/2 words
c
c     subprogram : rayeng
c
c     disper-80,  ver-86
c
c     m. saito  20/xi/80
c     revised   10/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  h(l),rho(l),vp(l),vs(l),ap(l),ae(l),
     *           qp(l),qs(l),y(4,l),p(5,l),ws(5,l),
     1           yy(5),qq(5),dy(5),yyy(5),qqq(5),
     2           f(4,3),s(4)
      data  eps/1.e-30/,big/1.e+30/
     *     ,rkgp/1.707107/,rkgm/0.2928932/
c
      rkgp = 1.d0 + dsqrt(0.5d0)
      rkgm = 0.5d0/rkgp
c
      if( l.le.2 .or. w.le.0 .or. c.le.0 )  go to  90
      i   = l
      if( rho(i).le.0 .or. c.ge.vp(i) )  go to  92
      ih  = 0
      lly = max(1,min(l,ly))
      ier = 0
      cc  = c*c
      wn  = w/c
      ro  = rho(i)
      roc = ro*cc
      sv  = vs(i)
      el  = rho(i)*sv**2
      do  1  j=1,5
        yy(j) = 0
        qq(j) = 0
    1 continue
      cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
c
      if( sv.le.0 )  then
        jx = 2
        yy(1) = eps*ra
        yy(2) =-eps*roc
c
      else
        if( c.ge.sv )  go to  92
        jx  = 5
        cs  = c/sv
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = 2*el
        yy(3) =-eps*ra
        yy(4) =-eps*rb
        yy(2) =-eps*(cp**2*rbb + cs**2)/(roc*(ra*rb + 1))
        yy(1) = rg*yy(2) + eps
        yy(5) =-rg*(yy(1) + eps) + eps*roc
      endif
c
      do  2  j=1,5
    2   ws(j,1) = yy(j)
      jy = 1
c
c     integrate upward
c
    3 if( h(i-1).gt.0 )  go to  6
      i = i - 1
      if( i.le.1 )  go to  92
      if( (vs(i).le.0 .and. el.le.0) .or.
     *    (vs(i).gt.0 .and. el.gt.0) )  go to  6
      do  4  j=1,5
        yyy(j) = yy(j)
        qqq(j) = qq(j)
         yy(j) = 0
         qq(j) = 0
    4 continue
c
c     solid  -----> liquid boundary
c
      if( el.gt.0 )  then
        jx = 2
        yy(1) = yyy(3)
        qq(1) = qqq(3)
        yy(2) = yyy(5)
        qq(2) = qqq(5)
c
c     liquid -----> solid boundary
c
      else
        jx = 5
        yy(2) = yyy(1)
        qq(2) = qqq(1)
        yy(4) = yyy(2)
        qq(4) = qqq(2)
      endif
c
      jy = jy + 1
      do  5  j=1,5
    5   ws(j,jy) = yy(j)
c
    6 hh  = wn*(h(i-1) + h(i-2))
      h5  = hh/2
      jmp = 0
c
    7 if( i.le.0  )  go to  92
      if( i.ne.ih )  then
        ih = i
        ro = rho(i)
        roc= ro*cc
        ea = ro*vp(i)**2
        if( ea.le.0 )  go to  92
        el = ro*vs(i)**2
        ec = ea
        ef = ea - el-el
        if( ia.gt.0 )  then
          ec = ec*ap(i)
          ef = ef*ae(i)
        endif
c
        if( el.le.0 )  then
          a12 = 1/ea - 1/roc
c
        else
          a12 = 1/ec
          a13 = ef*a12
          b31 = a13+a13
          a34 = 1/el
          a43 =-roc + ea - ef*a13
        endif
      endif
c
    8 if( el.le.0 )  then
        dy(1) = a12*yy(2)
        dy(2) =-roc*yy(1)
c
      else
        dy(1) =     yy(3) - a13*yy(4)
        dy(2) = a34*yy(3) + a12*yy(4)
        dy(3) =-b31*yy(1) + a43*yy(2) + a12*yy(5)
        dy(4) =   2*yy(1) - roc*yy(2) + a34*yy(5)
        dy(5) =-roc*yy(3) + a43*yy(4)
      endif
c
      jmp = jmp + 1
      go to  (9,11,13,15),jmp
c
    9 do  10  j=1,jx
        wk = h5*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)
        qq(j) = wq + (yy(j) - wy)*3 - wk
   10 continue
      i = i - 1
      go to  7
c
   11 do  12  j=1,jx
        wk = hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgm
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgm
   12 continue
      go to  8
c
   13 do  14  j=1,jx
        wk = hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgp
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgp
   14 continue
      i = i - 1
      go to  7
c
   15 do  16 j=1,jx
        wk = h5*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)/3
        qq(j) = wq + (yy(j) - wy)*3 - wk
   16 continue
c
      jy = jy + 1
      do  17  j=1,5
   17   ws(j,jy) = yy(j)
c
c     normalization
c
      if( el.gt.0 )  then
        z1 = 0
        do  18  j=1,5
   18     z1 = max(z1, abs(yy(j)))
        if( z1.gt.big )  then
          do  19  j=1,5
            yy(j) = eps*yy(j)
            qq(j) = eps*qq(j)
   19     continue
        endif
      endif
c
      if( i.gt.2 )  go to  3
c
      ny = 1
      do  20  j=1,4
          s(j) = 0
        f(j,1) = 0
        f(j,2) = 0
        f(j,3) = 0
   20 continue
c
c     liquid surface
c
      if( el.le.0 )  then
        a = 1/yy(1)
        yy(1) = 1
        yy(2) = 0
        ec = ea
        ef = ea - 2*el
c
c     solid surface
c
      else
        yy(3) =-yy(1)/yy(3)
        yy(1) = 1
        yy(2) = 0
        yy(4) = 0
      endif
c
      call  rayeng(ro,ea,ec,el,ef,qp(i),qs(i),c,yy,f,
     *             iy,y,ip,p,iq)
c
c     integrate downward
c
   21 do  33  k=2,3
c
      hh  = wn*(h(i) + h(i+1))
      h5  = hh/2
      jmp = 0
c
      if( el.le.0 )  then
        jy = jy - 1
        yy(1) = a*ws(1,jy)
        yy(2) = a*ws(2,jy)
        i  = i + 2
      endif
c
   22 if( i.ne.ih )  then
        ih  = i
        ro  = rho(i)
        ea  = ro*vp(i)**2
        el  = ro*vs(i)**2
        roc = ro*cc
        ec  = ea
        ef  = ea - el-el
        if( el.le.0 )  go to  32
        if( ia.gt.0 )  then
          ec = ec*ap(i)
          ef = ef*ae(i)
        endif
        a12 = 1/ec
        a13 = ef*a12
        a34 = 1/el
        a43 =-roc + ea - ef*a13
      endif
c
   23 dy(1) = a12*yy(2) + a13*yy(3)
      dy(2) =-roc*yy(1) +     yy(4)
      dy(3) =    -yy(1) + a34*yy(4)
      dy(4) =-a13*yy(2) + a43*yy(3)
c
      jmp = jmp + 1
      go to  (24,26,28,30),jmp
c
   24 do  25  j=1,4
        wk =-h5*dy(j)
        wy = yy(j)
        yy(j)= wy + wk
        qq(j) = (yy(j) - wy)*3 - wk
   25 continue
      i = i + 1
      go to  22
c
   26 do  27  j=1,4
        wk =-hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgm
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgm
   27 continue
      go to  23
c
   28 do  29  j=1,4
        wk =-hh*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)*rkgp
        qq(j) = wq + (yy(j) - wy)*3 - wk*rkgp
   29 continue
      i = i + 1
      go to  22
c
   30 do  31 j=1,4
        wk =-h5*dy(j)
        wy = yy(j)
        wq = qq(j)
        yy(j) = wy + (wk - wq)/3
        qq(j) = wq + (yy(j) - wy)*3 - wk
   31 continue
c
      jy = jy - 1
      yy(2) = (yy(1)*ws(4,jy) + yy(3)*ws(1,jy))/ws(2,jy)
      yy(4) = (yy(1)*ws(1,jy) + yy(3)*ws(3,jy))/ws(2,jy)
c
   32 ny = ny + 1
      call  rayeng(ro,ea,ec,el,ef,qp(i),qs(i),c,
     *             yy,f(1,k),iy,y(1,ny),ip,p(1,ny),iq)
c
   33 continue
c
      h1 = h(i-3) + h(i-4)
      h2 = h(i-2) + h(i-1)
      hh = (h1 + h2)/6
      s1 = hh*(2 - h2/h1)
      s2 = hh*(4 + (1 - h1/h2)*(h2/h1 - 1))
      s3 = hh*(2 - h1/h2)
      do  34  j=1,4
          s(j) = s(j) + s1*f(j,1) + s2*f(j,2) + s3*f(j,3)
        f(j,1) = f(j,3)
   34 continue
c
      if(  i.ge.lly )  go to  36
      if( h(i).gt.0 )  go to  21
      i = i + 1
      if( (vs(i).le.0 .and. el.le.0) .or.
     *    (vs(i).gt.0 .and. el.gt.0) )  go to  35
c
c     liquid -----> solid
c
      jy = jy - 1
      if( el.le.0 )  then
        yy(3) =-yy(1)*(ws(1,jy)/ws(3,jy))
        yy(4) = 0
c
c     solid  -----> liquid
c
      else
        a  = yy(1)/ws(1,jy)
      endif
c
   35 ny = ny + 1
      ro = rho(i)
      roc= ro*cc
      ea = ro*vp(i)**2
      el = ro*vs(i)**2
      ec = ea
      ef = ea - el-el
      if( ia.gt.0 .and. el.gt.0 )  then
        ec = ec*ap(i)
        ef = ef*ae(i)
      endif
      call  rayeng(ro,ea,ec,el,ef,qp(i),qs(i),
     *             c,yy,f,iy,y(1,ny),ip,p(1,ny),iq)
      if( i.le.lly-4 )  go to  21
c
c     normal exit
c
c     half space integral (isotropic model)
c
   36 cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
c
c     liquid half space
c
      if( el.le.0 )  then
        aa = (yy(1)/ra)**2/(2*wn*ra)
        s(1) = s(1) + ro*(1 + raa)*aa
        s(2) = s(2) + roc*cp**2*aa
        s(3) = s(3) + roc*aa
        if( iq.ne.0 )  s(4) = s(4) + qp(i)*cp**2*aa
c
c     solid half sapce
c
      else
        cs  = c/vs(i)
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = 2*el
        cb  = (ra*rb + 1)/(cp**2 + cs**2*raa)
        ca  = (yy(3) - rb*yy(1))*cb
        cb  = (yy(1) - ra*yy(3))*cb
        aa  = ca**2/(2*wn*ra)
        ab  = ca*cb/(wn*(ra + rb))
        bb  = cb**2/(2*wn*rb)
        ya2 = rg - roc
        ya4 = rg*ra
        yb2 = rg*rb
        s11 = raa*aa + 2*ra*ab + bb
        s22 = ya2**2*aa + 2*ya2*yb2*ab + yb2**2*bb
        s33 = aa + 2*rb*ab + rbb*bb
        s44 = ya4**2*aa + 2*ya4*ya2*ab + ya2**2*bb
        s14 = ra*ya4*aa + (ra*ya2 + ya4)*ab + ya2*bb
        s23 = ya2*aa + (rb*ya2 + yb2)*ab + rb*yb2*bb
        ef  = ea - el*2
        af  = (ea - ef*ef/ea)*s33
        s(1) = s(1) + ro*(s11 + s33)
        s(2) = s(2) + s22/ea + s44/el + af
        s(3) = s(3) + af + s14 - ef*s23/ea
        if( iq.ne.0 )  then
          dp  = (s22 - 4*el*(s23 - el*s33))/ea
          ds  = s44/el + 4*el*(s23 + ef*s33)/ea
          s(4) = s(4) + qp(i)*dp + qs(i)*ds
        endif
      endif
c
      u   = s(3)/(c*s(1))
      er  = 1 - s(2)/(cc*s(1))
      eki = 2*wn*wn*s(3)
      if( iq.ne.0 )  q = s(4)/s(3)
c
c     normalization
c
      if( ip.ne.0 )  then
        if( ekd.gt.0 )  then
          z1 = wn**2/ekd
        else
          z1 = 1/(2*s(3))
        endif
        do  37  i=1,5
        do  37  j=1,ny
          p(i,j) = z1*p(i,j)
   37   continue
      endif
c
c     conversion to the conventional notation
c
      if( iy.ne.0 )  then
        z2 = wn
        do  38  j=1,ny
          y(2,j) = z2*y(2,j)
          y(4,j) = z2*y(4,j)
   38   continue
      endif
c
      return
c
c     input error
c
   90 write(6,91)  l,w,c
   91   format(20x,5('?'),3x,'(rayefr)   input error',3x,
     *         'l  =',i5,3x,'w =',1pe13.6,3x,
     1         'c  =',e13.6,3x,5('?'))
      ier =-1
      return
c
   92 write(6,93)  i
   93   format(20x,5('?'),3x,'(rayefr)   input error ',
     *         'at layer',i5,3x,5('?'))
      ier =-max(1,i)
      return
      end
c
c     rayleigh wave eigenfunction and partial derivatives
c                  by matrix method integration
c
      subroutine  rayefx(h,rho,vp,vs,qp,qs,l,ly,ld,w,c,
     *                   u,iq,q,ekd,eki,er,iy,y,
     1                   ip,p,ny,ws,ier)
c
c     input
c       h    : layer thickness.  h(l) is arbitrary
c       rho  : density
c       vp   : p wave velocity
c       vs   : s wave velocity
c       qp   : inverse of p-wave q
c       qs   : inverse of s-wave q
c       l    : number of layers inclding the bottom
c              half-space
c       ly   : deepest layer where eigenfunction has
c              to be computed
c       ld   : number of sublayers/layer
c              actually, 2*ld sublayers are introduced
c       w    : angular frequency
c       c    : phase velocity
c       ekd  : energy integral 2*k**2*i3
c              computed by raymrx
c       iy   : = 0 ; not to compute eigenfunction
c              = 1 ;     to compute eigenfunction
c       ip   : = 0 ; not to compute partials
c              = 1 ;     to compute partials
c       iq   : = 0 ; not to compute attenuation
c              = 1 ;     to compute attenuation
c     output
c       u    : group velocity by energy integrals
c       q    : inverse of rayleigh wave q
c       eki  : energy integral 2*k**2*i3 by integration
c       er   : relative error in energy integrals
c       y    : eigenfunction
c       p    : partials
c              p(1,j) ; dc/drho
c              p(2,j) ; dc/dvp
c              p(3,j) ; dc/dvs
c              p(4,j) ; dc/dphi
c              p(5,j) ; dc/deta
c       ny   : length of y and/or p
c       ier  : return code
c              < 0 ; input error
c              = 0 ; no error
c              > 0 ; overflow
c     working space
c       ws   : 5*l words
c
c     isotropic model only
c
c     subprogram : rayeng
c
c     disper-80,  ver-86
c
c     m. saito  30/vii/79
c     revised   10/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  h(l),rho(l),vp(l),vs(l),qp(l),qs(l),
     *           y(4,*),p(5,*),ws(5,*),yy(5),z(5),
     1           f(4,3),s(4)
      data  eps/1.e-30/,big/1.e+30/
c
c     define sinh(x)/x
c
c      sh0(x) = 0.9999997 + x*(0.1666667 + x*(0.0083361
c     *       + x*0.0001984))
c
c     for double precision use
c
      sh0(x) = 1.0d0 + x*(1.6666 66666 66666 7d-1
     *       + x*(8.33 33333 33334 0d-3
     1       + x*(1.9 84126 96412 7d-4
     2       + x*(2.7557 31918 9d-6 + x*(2.50 12108 4d-8
     3       + x*(1.60596 1d-10 + x*7.64 7d-13))))))
c 
c     initial value
c
      if( l.le.1 .or. c.le.0 .or. w.le.0 )  go to  90
      i   = l
      ro  = rho(i)
      if( ro.le.0 .or. c.ge.vp(i) )  go to  92
      ier = 0
      ny  = 0
      lly = max(1,min0(l,ly))
      cc  = c*c
      wn  = w/c
      roc = ro*cc
      sv  = vs(i)
      cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
      do  1  j=1,5
    1   yy(j) = 0
c
      if( sv.le.0 )  then
        yy(1) = ra
        yy(2) =-roc
        jx = 2
c
      else
        if( c.ge.sv )  go to  92
        cs  = c/sv
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = 2*roc/cs**2
        yy(3) =-eps*ra
        yy(4) =-eps*rb
        yy(2) =-eps*(cp**2*rbb + cs**2)/(roc*(ra*rb + 1))
        yy(1) = rg*yy(2) + eps
        yy(5) =-rg*(yy(1) + eps) + eps*roc
        jx = 5
      endif
c
      do  2  j=1,5
    2   ws(j,i) = yy(j)
c
c     integrate upward
c
      if( l.le.1 )  go to  17
      do  9  ii=2,l
        i  = i - 1
        ro = rho(i)
        roc= ro*cc
        pv = vp(i)
        sv = vs(i)
        if( pv.le.0 .or. ro.le.0 .or.
     *      h(i).le.0 )  go to  92
        do  3  j=1,5
    3     z(j) = yy(j)
        if( (sv.le.0 .and. vs(i+1).le.0) .or.
     *      (sv.gt.0 .and. vs(i+1).gt.0) )  go to  4
c
c     solid  -----> liquid
c
        if( sv.le.0 )  then
          z(1) = yy(3)
          z(2) = yy(5)
          jx = 2
c
c     liquid -----> solid
c
        else
          z(2) = yy(1)
          z(4) = yy(2)
          z(1) = 0
          z(3) = 0
          z(5) = 0
          jx = 5
        endif
c
    4   r2  = 1/roc
        cp  = c/pv
        raa = (1 + cp)*(1 - cp)
        cs  = 0
        if( sv.gt.0 )  cs = c/sv
        rbb = (1 + cs)*(1 - cs)
        hk  = h(i)*wn
        hkk = hk**2
        xx  = raa*hkk
        one = 1
c
c       sinh(x)/x
c
        do  5  k=1,2
          cha = chb
          sha = shb
          ax  = abs(xx)
          if( ax.le.1 )  then
            shb = sh0(xx)
            chb = 1 + xx*sh0(xx/4)**2/2
          else
            ax  = sqrt(ax)
            if( xx.le.0 )  then
              chb = cos(ax)
              shb = sin(ax)/ax
            else
              if( ax.gt.100 )  then
                one = 0
              else
                one = one/cosh(ax)
              endif
              chb = 1
              shb = tanh(ax)/ax
            endif
          endif
          xx  = hkk*rbb
          shb = hk*shb
    5   continue
c
c     liquid layer
c
        if( sv.le.0 )  then
          b11 = cha
          b12 =-raa*sha*r2
          b21 =-roc*sha
c
          yy(1) = b11*z(1) + b12*z(2)
          yy(2) = b21*z(1) + b11*z(2)
c
c     solid layer
c
        else
          g1  = 2/cs**2
          rg  = g1*roc
          r4  = rg - roc
          e1  = cha*chb
          e2  = e1 - one
          e3  = sha*shb
          e5  = sha*chb
          e6  = shb*cha
          f1  = e2 - e3
          f2  = r2*f1
          f3  = g1*f1 + e3
          b33 = e1
          b34 = raa*e3
          b43 = rbb*e3
          b25 =-r2*(f2 + r2*(e2 - raa*b43))
          b15 = rg*b25 + f2
          b16 =-rg*b15 - f3
          b22 = b16 + e1
          b12 = rg*b16 - r4*f3
          b52 =-rg*b12 + r4*(rg*f3 + r4*e3)
          b23 = r2*(e5 - rbb*e6)
          b13 = rg*b23 - e5
          b42 =-rg*b13 + r4*e5
          b24 = r2*(e6 - raa*e5)
          b14 = rg*b24 - e6
          b32 =-rg*b14 + r4*e6
          b11 = one - b16-b16
          b21 = b15+b15
          b31 = b14+b14
          b41 = b13+b13
          b51 = b12+b12
c
          yy(1) = b11*z(1) + b12*z(2) + b13*z(3)
     *          + b14*z(4) + b15*z(5)
          yy(2) = b21*z(1) + b22*z(2) + b23*z(3)
     *          + b24*z(4) + b25*z(5)
          yy(3) = b31*z(1) + b32*z(2) + b33*z(3)
     *          + b34*z(4) + b24*z(5)
          yy(4) = b41*z(1) + b42*z(2) + b43*z(3)
     *          + b33*z(4) + b23*z(5)
          yy(5) = b51*z(1) + b52*z(2) + b42*z(3)
     *          + b32*z(4) + b22*z(5)
        endif
c
        do  6  j=1,5
    6     ws(j,i) = yy(j)
c
c     normalization
c
        if( el.gt.0 )  then
          z1 = 0
          do  7  j=1,jx
    7       z1 = max(z1,abs(yy(j)))
          if( z1.gt.big )  then
            do  8  j=1,jx
    8         yy(j) = eps*yy(j)
          endif
        endif
c
    9 continue
c
c     integrate downward
c
      do  10  j=1,4
          s(j) = 0
        f(j,1) = 0
        f(j,2) = 0
        f(j,3) = 0
   10 continue
      lx  = max(1,ld)
      div = 2*lx
      if( sv.gt.0 )  yy(3) =-yy(1)/yy(3)
      yy(1) = 1
      yy(2) = 0
      yy(4) = 0
c
      if( lly.le.1 )  go to  17
      do  16  ii=2,lly
        ro  = rho(i)
        roc = ro*cc
        pv  = vp(i)
        sv  = vs(i)
        hh  = h(i)/div
        ea  = ro*pv**2
        el  = ro*sv**2
        ef  = ea - el-el
        r2  = 1/roc
        cp  = c/pv
        raa = (1 + cp)*(1 - cp)
        cs  = 0
        if( sv.gt.0 )  cs = c/sv
        rbb = (1 + cs)*(1 - cs)
        hk  = hh*wn
        hkk = hk**2
        xx  = raa*hkk
c
        do  11  k=1,2
          cha = chb
          sha = shb
          ax  = abs(xx)
          if( ax.le.1 )  then
            shb = sh0(xx)
            chb = 1 + xx*sh0(xx/4)**2/2
          else
            ax  = sqrt(ax)
            if( xx.le.0 )  then
              chb = cos(ax)
              shb = sin(ax)/ax
            else
              if( ax.gt.100 )  go to  94
              ex  = exp(ax)
              exi = 1/ex
              chb = (ex + exi)/2
              shb = (ex - exi)/(2*ax)
            endif
          endif
          xx = hkk*rbb
          shb =-hk*shb
   11   continue
c
c     layer matrix
c
c       liquid layer
c
        if( sv.le.0 )  then
          b11 = cha
          b12 =-raa*sha*r2
          b21 =-roc*sha
c
c       solid layer
c
        else
          rg  = el+el
          r4  = rg - roc
          b14 = (cha - chb)*r2
          b11 =-rg*b14 + cha
          b33 = rg*b14 + chb
          b23 = rg*(b33 - cha)
          b34 = (sha - rbb*shb)*r2
          b24 = rg*b34 - sha
          b21 =-rg*b24 + r4*sha
          b12 = (shb - raa*sha)*r2
          b13 =-rg*b12 + shb
          b43 = rg*b13 + r4*shb
        endif
c
        ny = ny + 1
        call  rayeng(ro,ea,ea,el,ef,qp(i),qs(i),c,
     *               yy,f,iy,y(1,ny),ip,p(1,ny),iq)
c
        do  15  ll=1,lx
c
          do  13  k=2,3
            do  12  j=1,4
   12         z(j) = yy(j)
c
            if( sv.le.0 )  then
              yy(1) = b11*z(1) + b12*z(2)
              yy(2) = b21*z(1) + b11*z(2)
c
            else
              yy(1) = b11*z(1) + b12*z(2) + b13*z(3)
     *              + b14*z(4)
              yy(2) = b21*z(1) + b11*z(2) + b23*z(3)
     *              + b24*z(4)
              yy(3) =-b24*z(1) - b14*z(2) + b33*z(3)
     *              + b34*z(4)
              yy(4) =-b23*z(1) - b13*z(2) + b43*z(3)
     *              + b33*z(4)
            endif
c
            ny = ny + 1
            call  rayeng(ro,ea,ea,el,ef,qp(i),qs(i),c,yy,
     *                   f(1,k),iy,y(1,ny),ip,p(1,ny),iq)
c
   13     continue
c
          h6 = hh/3
          do  14  j=1,4
              s(j) = s(j)
     *             + h6*(f(j,1) + 4*f(j,2) + f(j,3))
            f(j,1) = f(j,3)
   14     continue
c
   15   continue
c
        i = i + 1
c
        if( vs(i).le.0 )  then
          yy(2) = (ws(2,i)/ws(1,i))*yy(1)
c
        else
          if( sv.le.0 )  yy(3) =-(ws(1,i)/ws(3,i))*yy(1)
          yy(2) = (ws(4,i)*yy(1) + ws(1,i)*yy(3))/ws(2,i)
          yy(4) = (ws(1,i)*yy(1) + ws(3,i)*yy(3))/ws(2,i)
        endif
c
   16 continue
c
c     half space integral
c
   17 ro  = rho(i)
      roc = ro*cc
      ea  = ro*vp(i)**2
      el  = ro*vs(i)**2
      ef  = ea - el-el
      ny  = ny + 1
      call  rayeng(ro,ea,ea,el,ef,qp(i),qs(i),c,
     *             yy,f,iy,y(1,ny),ip,p(1,ny),iq)
      cp  = c/vp(i)
      raa = (1 + cp)*(1 - cp)
      ra  = sqrt(raa)
c
c     liquid half space
c
      if( vs(i).le.0 )  then
        aa = (yy(1)/ra)**2/(2*wn*ra)
        s(1) = s(1) + ro*(1 + raa)*aa
        s(2) = s(2) + roc*cp**2*aa
        s(3) = s(3) + roc*aa
        if( iq.ne.0 )  s(4) = s(4) + qp(i)*cp**2*aa
c
c     solid half space
c
      else
        cs  = c/vs(i)
        rbb = (1 + cs)*(1 - cs)
        rb  = sqrt(rbb)
        rg  = el+el
        cb  = (ra*rb + 1)/(cp**2 + cs**2*raa)
        ca  = (yy(3) - rb*yy(1))*cb
        cb  = (yy(1) - ra*yy(3))*cb
        aa  = ca**2/(2*wn*ra)
        ab  = ca*cb/(wn*(ra + rb))
        bb  = cb**2/(2*wn*rb)
        ya2 = rg - roc
        ya4 = rg*ra
        yb2 = rg*rb
        s11 = raa*aa + 2*ra*ab + bb
        s22 = ya2**2*aa + 2*ya2*yb2*ab + yb2**2*bb
        s33 = aa + 2*rb*ab + rbb*bb
        s44 = ya4**2*aa + 2*ya4*ya2*ab + ya2**2*bb
        s14 = ra*ya4*aa + (ra*ya2 + ya4)*ab + ya2*bb
        s23 = ya2*aa + (rb*ya2 + yb2)*ab + rb*yb2*bb
        af  = (ea - ef*ef/ea)*s33
        s(1) = s(1) + ro*(s11 + s33)
        s(2) = s(2) + s22/ea + s44/el + af
        s(3) = s(3) + af + s14 - ef*s23/ea
        if( iq.ne.0 )  then
          dp = (s22 - 4*el*(s23 - el*s33))/ea
          ds = s44/el + 4*el*(s23 + ef*s33)/ea
          s(4) = s(4) + qp(i)*dp + qs(i)*ds
        endif
      endif
c
      u   = s(3)/(c*s(1))
      er  = 1 - s(2)/(cc*s(1))
      eki = 2*wn*wn*s(3)
      if( iq.ne.0 )  q = s(4)/s(3)
c
      if( ip.ne.0 )  then
        if( ekd.gt.0 )  then
          z1 = wn**2/ekd
        else
          z1 = 1/(2*s(3))
        endif
        do  18  k=1,5
        do  18  j=1,ny
          p(k,j) = z1*p(k,j)
   18   continue
      endif
c
      if( iy.ne.0 )  then
        do  19  j=1,ny
          y(2,j) = wn*y(2,j)
          y(4,j) = wn*y(4,j)
   19   continue
      endif
      return
c
c     input error
c
   90 write(6,91)  l,w,c
   91   format(20x,5('?'),3x,'(rayefx)   input error',3x,
     *         'l  =',i5,3x,'w  =',1pe13.6,3x,
     1         'c  =',e13.6,3x,5('?'))
      ier =-1
      return
c
   92 write(6,93)  i
   93   format(20x,5('?'),3x,'(rayefx)   input error ',
     *         'at layer',i5,3x,5('?'))
      ier =-max(1,i)
      return
c
c     overflow
c
   94 write(6,95)  i
   95   format(20x,5('?'),3x,'(rayefx)   overflow  ',
     *         'at layer',i5,3x,5('?'))
      ier = i
      return
      end
c------------------------------------------------------------
c
c     rayleigh wave energy integral
c
      subroutine  rayeng(rho,ea,ec,el,ef,qp,qs,c,yy,f,
     *                   iy,y,ip,p,iq)
c
c     disper-80,  ver-86
c
c     m. saito  30/vii/79
c     revised   27/ix/85
c
      implicit real*8 (a-h,o-z)
      dimension  yy(4),f(4),y(4),p(5)
c
      cc = c*c
      y1 = yy(1)
      y2 = yy(2)
      y3 = yy(3)
      y4 = yy(4)
c
      if( el.le.0 )  then
        y3 =-y2/(rho*cc)
        y4 = 0
      endif
c
      if( iy.ne.0 )  then
        y(1) = y1
        y(2) = y2
        y(3) = y3
        y(4) = y4
      endif
c
      y42  = 0
      if( el.gt.0 )  y42 = y4**2/el
      fc = ef/ec
      f3 = (ea - ef*fc)*y3**2
      f(1) = rho*(y1**2 + y3**2)
      f(2) = y2**2/ec + y42 + f3
      f(3) = f3 + y1*y4 - fc*y2*y3
      if( ip.eq.0 .and. iq.eq.0 )  return
c
      et = ef/(ea - el-el)
      dph= y2 + ef*y3
      dp = (y2 - 2*et*el*y3)**2/ec
     *   + ea*(1 - (ea/ec)*et**2)*y3**2
      ds = y42 + 4*et*(el/ec)*y3*dph
      if( iq.ne.0 )  f(4) = qp*dp + qs*ds
c
      if( ip.ne.0 )  then
        p(1) = f(2) - cc*f(1)
        p(2) = dp+dp
        p(3) = ds+ds
        p(4) = dph**2/ec
        p(5) =-2.0d0*fc*y3*dph
      endif
c
      return
      end
c
c     find grid number for a given depth
c
      subroutine  dspdep(h,n,dep,k,l)
c
c     input
c       h    : grid interval
c       n    : number of grid points
c       dep  : depth
c       k    : search every k grid point
c              = 2 ; to find grid point consistent
c                    to r-k-g scheme
c              = 4 ; to find grid point consistent
c                    to simpson scheme
c                    (to compute energy integrals)
c     output
c       l    : grid number
c
c     disper-80,  ver-86
c
c     m. saito   28/i/78
c     revised     2/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  h(n)
c
      if( k.le.1 )  then
        kk = 1
      else if( k.le.3 )  then
        kk = 2
      else
        kk = 4
      endif
      dd = 0
      i  = 1
c
    1 if( i.ge.n )  go to  2
      if( h(i).eq.0 )  i = i + 1
      if( dep-dd.le.h(i)*0.1 )  go to  3
      dd = dd + h(i)
      if( kk.gt.1 )  then
        dd = dd + h(i+1)
        if( kk.eq.4 )  dd = dd + h(i+2) + h(i+3)
      endif
      i  = i + kk
      go to  1
c
    2 if( i.gt.n )  i = i - kk
    3 l  = i
      return
      end
c
c     prints eigenfunction or partials
c
      subroutine  dsprnt(y,l,m,h,d0)
c
c     disper-80,  ver-86
c
c     m. saito   13/xi/81
c     revised    18/xi/86
c
      implicit real*8 (a-h,o-z)
      dimension  y(l,m),h(m)
c
      dd = d0
      i  =-1
      j  = 0
c
    1 i  = i + 2
      j  = j + 1
      write(6,2)  i,dd,(y(k,j),k=1,l)
    2   format(i5,1p7e15.6)
c   2   format(i5,1p7d15.6)
      if( j.ge.m )  return
      if( h(i).gt.0 )  then
        dd = dd + h(i) + h(i+1)
      else
        i  = i - 1
      endif
      go to  1
c
      end
c
c     print eigenfunction or partial
c
      subroutine  dsprnx(y,k,ny,ld,h,d0)
c
c     disper-80,  ver-86
c
c     m. saito  15/xi/81
c     revised    3/xii/86
c
      implicit real*8 (a-h,o-z)
      dimension  y(k,ny),h(1)
c
      j   = 0
      i   = 0
      d   = d0
      ld2 = ld*2
      div = ld2
c
    1 j = j + 1
      i = i + 1
      write(6,2)  i,d,(y(m,j),m=1,k)
    2   format(i5,1p5e15.6)
c   2   format(i5,1p5d15.6)
      if( j.ge.ny )  return
      hh = h(i)/div
      do  3  l=1,ld2
        d = d + hh
        j = j + 1
        write(6,2)  i,d,(y(m,j),m=1,k)
    3 continue
      if( j.lt.ny )  go to  1
c
      return
      end
