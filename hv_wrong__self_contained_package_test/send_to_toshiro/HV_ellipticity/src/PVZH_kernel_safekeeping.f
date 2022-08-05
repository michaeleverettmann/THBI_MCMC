c
c program PVZH_kernel.f
c
c	rayleigh wave flat-layer version
c  	calls rayleigh_sub.f
c
c
      implicit real*8(a-h,o-z)
      character modfile*80 
      character ofile*80      
      common /struc/ r(50000),h(50000),rho(50000),vp(50000),
     2 vs(50000),qp(50000), qs(50000), nlay
      common /struc2/ r2(50000),h2(50000),rho2(50000),vp2(50000),
     2 vs2(50000),qp2(50000), qs2(50000), nlay2
      common /raysol/ sol1(4,50000), sol2(4,50000)
      dimension y0(3),dep(50000),
     *          y(4,50000),p(5,50000),ap(50000),ae(50000),
     *          x(4,50000),yyy2(50000),dum(50000),
     *          yyy4(50000),yyy1(50000),yyy3(50000),yyy5(50000),
     *          yyy6(50000),yij(15),ws(5,50000)
      dimension regn(50000)
      dimension knot(50000)
      dimension dummy(50000),rnew(50000),rold(50000)
      dimension rho_org(50000),vp_org(50000),vs_org(50000),r_org(50000)
      dimension qp_org(50000), qs_org(50000)
      dimension depth1(50000),depth2(50000)
      dimension vpave(50000),vsave(50000),rhoave(50000)
      dimension akrh0(50000),akvp(50000),akvs(50000)
      dimension akrh(50000),akk(50000),akmu(50000)
      dimension pv_rh(50000), pv_vp(50000), pv_vs(50000)
      dimension zhrho(50000),zhvp(50000),zhvs(50000)
      dimension zh_rh(50000), zh_vp(50000), zh_vs(50000)
      data  pi/3.141592653589793d0/
      external  rayrkg
c
      pi2=2.0d0*pi
c
c     Input model : The Original Model (referenced at 1 Hz)
c
 101  format(a80)

      print *,' input model file : '
      read(*,101) modfile
      open(3,file=modfile)

      print *,' output model file : '
      read(*,101) ofile
      open(6,file=ofile)

      rewind 3
C     First line of infile is frequency range to output
      read(3,*)fstart,fstop,dfq
C     Second line of infile is knot description 
      read(3,*)dh,ncase0
C     Third line of infile is number of layers in input model 
      read(3,*)nlay0
      do i=1,nlay0
      read(3,*)r_org(i),rho_org(i),vp_org(i),vs_org(i),
     2 qp_org(i),qs_org(i)
      if(r_org(i).gt.5000.) r_org(i)=6371.-r_org(i)
      qp_org(i)=1./300.
      qs_org(i)=1./200.
      enddo
      close(3)
c
c knot determination - now defined in infile
c	perturbation of layers with thickness dh
c	ncase is the number of layers
c
C     ncase0=10
C     dh=1.0
      do i=1,ncase0
	depth1(i)=dh*(i-1)
        depth2(i)=depth1(i)+dh
      enddo
c
c average density and velocities
c	They are required only for derivation of akk and amu
c
      do j=1,ncase0
      do ipar=1,3
      xd=0.0
      kount=0
      do i=1,nlay0
      yd=r_org(i)
      if(yd.ge.depth1(j).and.yd.lt.depth2(j)) then
	kount=kount+1
	if(ipar.eq.1) xd=xd+vs_org(i)
	if(ipar.eq.2) xd=xd+vp_org(i)
	if(ipar.eq.3) xd=xd+rho_org(i)
      endif
      enddo
      if(ipar.eq.1) vsave(j)=xd/kount
      if(ipar.eq.2) vpave(j)=xd/kount
      if(ipar.eq.3) rhoave(j)=xd/kount
      enddo
      enddo
c
      write(6,9000) ncase0
 9000 format(i10)
      write(6,9001)(i,depth1(i),depth2(i),rhoave(i),vpave(i),vsave(i),
     2 i=1,ncase0)
 9001 format(i5,5e15.7)
c
c Frequency range and increment - now get these from infile
c
C     fstart=0.13d0
C     fstop=0.37d0
C     dfq=0.01d0
      nfq=(fstop-fstart)/dfq+1.1
c
      write(6,780)fstart, fstop, dfq
  780 format(3e15.7)
c
c tol : integration tolerance
c itr : iteration for mode search
      tol = 0.00001
      itr = 40
      ia  = 0
c

c big loop over frq
      do ifq=1,nfq
      frq= fstart+dfq*(ifq-1)
      om = pi2*frq
c
c determination of starting depth for integration
      mlay=nlay0
c
c determine the eigenfrequency and eigenfunctions of the reference state
c
      do i=1,mlay
      r(i)=r_org(i)
      rho(i)=rho_org(i)
      vp(i)=vp_org(i)*(1.0d0+qp_org(i)/pi*qcoe)
      vs(i)=vs_org(i)*(1.0d0+qs_org(i)/pi*qcoe)
      rold(i)=r(i)
      qp(i)=qp_org(i)
      qs(i)=qs_org(i)
      enddo
c
c interpolation to increase the numebr of layers four-fold
c
      call interp(mlay,r,rnew,rho,dummy,kount)
      do i=1,kount
	rho(i)=dummy(i)
      enddo
      call interp(mlay,r,rnew,vp,dummy,kount)
      do i=1,kount
        vp(i)=dummy(i)
      enddo
      call interp(mlay,r,rnew,qp,dummy,kount)
      do i=1,kount
        qp(i)=dummy(i)
      enddo
      call interp(mlay,r,rnew,qs,dummy,kount)
      do i=1,kount
        qs(i)=dummy(i)
      enddo
      call interp(mlay,r,rnew,vs,dummy,kount)
      do i=1,kount
        vs(i)=dummy(i)
	r(i)=rnew(i)
      enddo
c
c nlay is defined here
      nlay=kount
c
      do i=1,nlay-1
	h(i)=r(i+1)-r(i)
	if(dabs(h(i)).lt.1.0e-10) h(i)=0.0d0
      enddo
      h(nlay)=0.0
c
      cmn = 0.001
      cmx = 5.00
      dc  = 0.005
c
      call  raydsp(rayrkg,h,rho,vp,vs,ap,ae,nlay,om,cmn,cmx,
     *             dc,tol,itr,ia,pvel,u,ekd,y0,yij,ier)
c
      iy=1
      ip=1
      iq=1
      ly=nlay
c
      call  rayefr(h,rho,vp,vs,ap,ae,qp,qs,nlay,ly,om,pvel,
     *             u,ia,iq,q,ekd,eki,er,
     1             iy,y,ip,p,ny,ws,ier)
c
      pvel0=pvel
      gvel0=u
      ak=om/pvel
      zh0=y(1,1)/y(3,1)
c
      call sol_rayl(om,ak)
c
      call egn_knots(ny,h,nlay2,knot)
      do i=1,nlay2
        j=knot(i)
        regn(i)=r(j)
      enddo
      do i=1,nlay2
	xd=regn(i)
	call search(xd,r,nlay,ians)
	r2(i)=r(ians)
	rho2(i)=rho(ians)
	vp2(i)=vp(ians)
	vs2(i)=vs(ians)
	do j=1,4
	   x(j,i)=sol1(j,ians)
	enddo
      enddo
c
      do i=1,nlay2-1
	h2(i)=r2(i+1)-r(i)
      enddo
      h2(nlay2)=0.0
c
      call ZHker(om,ak,x,y,zhrho,zhvp,zhvs)
c
      write(6,790)frq,er,gvel0,pvel0,zh0
  790 format('frq=',5e15.7)
c
c averaging analytical kernels
c reduce the number of knots to the original
c input file
c
      do j=1,ncase0
      do ipar=1,6
      xd=0.0d0
      kount=0
      do i=1,nlay2
      yd=r2(i)
      if(yd.ge.depth1(j).and.yd.lt.depth2(j)) then
        kount=kount+1
        if(ipar.le.3) xd=xd+p(ipar,i)
	if(ipar.eq.4) xd=xd+zhrho(i)
	if(ipar.eq.5) xd=xd+zhvp(i)
	if(ipar.eq.6) xd=xd+zhvs(i)
      endif
      enddo
      if(ipar.eq.1) pv_rh(j)=xd/kount*dh
      if(ipar.eq.2) pv_vp(j)=xd/kount*dh
      if(ipar.eq.3) pv_vs(j)=xd/kount*dh
      if(ipar.eq.4) zh_rh(j)=xd/kount*dh
      if(ipar.eq.5) zh_vp(j)=xd/kount*dh
      if(ipar.eq.6) zh_vs(j)=xd/kount*dh
      enddo
      enddo
c
      write(6,9000) ncase0
      do i=1,ncase0
      write(6,860)depth1(i),depth2(i),pv_rh(i),
     2 pv_vp(i),pv_vs(i),zh_rh(i),zh_vp(i),zh_vs(i)
  860 format(2f11.6,6e15.7)
      enddo
c
      enddo
c
  999 continue
      stop
      end
c-------------------------------------------------------------------------
      subroutine ZHker(om,ak,x,y,zhrho,zhvp,zhvs)
c
c ZH kernels for a mode at om (omega) and ak (wavenumber)
c The number of points for depth is nlay2
c
      implicit real*8(a-h,o-z)
      common /struc2/ r2(50000),h2(50000),rho2(50000),vp2(50000),
     2 vs2(50000),qp2(50000), qs2(50000), nlay2
      dimension zhrho(1), zhvp(1), zhvs(1)
      dimension x(4,1), y(4,1), dxdx(10), dydx(10), dum(10)
c
c normalization check
c
c at the surface, we must have y3=1 and x2=1
c
      if(dabs(y(3,1)-1.0d0).gt.1.0e-5) then
	anorm=y(3,1)
	do i=1,nlay2
	do j=1,4
		y(j,i)=y(j,i)/anorm
	enddo
	enddo
      endif
      if(dabs(x(2,1)-1.0d0).gt.1.0e-5) then
        anorm=x(2,1)
        do i=1,nlay2
        do j=1,4
                x(j,i)=x(j,i)/anorm
        enddo
        enddo
      endif
c
      call intJK(om,ak,4,x,y,ansxy)
      call intJK(om,ak,4,y,y,ansyy)
      ratJK=ansxy/ansyy
c
c ZH ratio
      zh0=y(1,1)/y(3,1)
c
      do i=1,nlay2
	ala2mu=rho2(i)*vp2(i)*vp2(i)
        amu=rho2(i)*vs2(i)*vs2(i)
        ala=ala2mu-2.0d0*amu
        rom2=rho2(i)*om*om
	coevp=2.0d0*ala2mu/zh0
	coevs=2.0d0*amu/zh0
	coerho=rho2(i)/zh0
        do j=1,4
           dum(j)=x(j,i)
        enddo
        call dy_rayl(ak,rom2,ala,amu,dum,dxdx)
        do j=1,4
           dum(j)=y(j,i)
        enddo
        call dy_rayl(ak,rom2,ala,amu,dum,dydx)
c
c zhvp
	xd1=ratJK*(dydx(1)-ak*y(3,i))*(dydx(1)-ak*y(3,i))
	xd2=(dxdx(1)-ak*x(3,i))*(dydx(1)-ak*y(3,i))
	zhvp(i)=coevp*(xd1-xd2)
c
c zhvs
	xd1=ratJK*((dydx(3)+ak*y(1,i))**2+4.0d0*ak*dydx(1)*y(3,i))
 	xdum1=(dxdx(3)+ak*x(1,i))*(dydx(3)+ak*y(1,i))
	xdum2=2.0d0*ak*(dxdx(1)*y(3,i)+dydx(1)*x(3,i))
	zhvs(i)=coevs*(xd1-xdum1-xdum2)
c
c zhrho
	xdum1=vp2(i)*vp2(i)*(dydx(1)-ak*y(3,i))*(dydx(1)-ak*y(3,i))
	xdum2=vs2(i)*vs2(i)*(dydx(3)+ak*y(1,i))*(dydx(3)+ak*y(1,i))
	xdum3=4.0d0*ak*vs2(i)*vs2(i)*dydx(1)*y(3,i)
	xdum4= -om*om*(y(1,i)*y(1,i)+y(3,i)*y(3,i))
	xd1=ratJK*(xdum1+xdum2+xdum3+xdum4)
	xdum1=vp2(i)*vp2(i)*(dxdx(1)-ak*x(3,i))*(dydx(1)-ak*y(3,i))
	xdum2=vs2(i)*vs2(i)*(dxdx(3)+ak*x(1,i))*(dydx(3)+ak*y(1,i))
	xdum3=2.0d0*ak*vs2(i)*vs2(i)*(dxdx(1)*y(3,i)+dydx(1)*x(3,i))
	xdum4= -om*om*(x(1,i)*y(1,i)+x(3,i)*y(3,i))
	xd2= xdum1+xdum2+xdum3+xdum4
	zhrho(i)=coerho*(xd1-xd2)
c
      enddo
      return
      end
c-------------------------------------------------------------------------
      subroutine intJK(om,ak,n,x,y,ans)
      implicit real*8(a-h,o-z)
      common /struc2/ r2(50000),h2(50000),rho2(50000),vp2(50000),
     2 vs2(50000),qp2(50000), qs2(50000), nlay2
      dimension x(4,nlay), y(4,nlay), dum(10), dxdx(10), dydx(10)
c
      ans=0.0d0
      do i=1,nlay2-1
	ala2mu=rho2(i)*vp2(i)*vp2(i)
	amu=rho2(i)*vs2(i)*vs2(i)
	ala=ala2mu-2.0d0*amu
	rom2=rho2(i)*om*om
	do j=1,4
	   dum(j)=x(j,i)
	enddo
	call dy_rayl(ak,rom2,ala,amu,dum,dxdx)
	do j=1,4
           dum(j)=y(j,i)
        enddo
        call dy_rayl(ak,rom2,ala,amu,dum,dydx)
	xd1=-x(3,i)*(dydx(1)-ak*y(3,i))-y(3,i)*(dxdx(1)-ak*x(3,i))
	xd2=x(1,i)*(dydx(3)+ak*y(1,i))+y(1,i)*(dxdx(3)+ak*x(1,i))
	xd3=dxdx(1)*y(3,i)+x(3,i)*dydx(1)
	x1=ala2mu*xd1+amu*(xd2+2.0d0*xd3)
c
	ala2mu=rho2(i+1)*vp2(i+1)*vp2(i+1)
        amu=rho2(i+1)*vs2(i+1)*vs2(i+1)
        ala=ala2mu-2.0d0*amu
        rom2=rho2(i+1)*om*om
        do j=1,4
           dum(j)=x(j,i+1)
        enddo
        call dy_rayl(ak,rom2,ala,amu,dum,dxdx)
        do j=1,4
           dum(j)=y(j,i+1)
        enddo
        call dy_rayl(ak,rom2,ala,amu,dum,dydx)
        xd1=-x(3,i+1)*(dydx(1)-ak*y(3,i+1))
     2		-y(3,i+1)*(dxdx(1)-ak*x(3,i+1))
        xd2=x(1,i+1)*(dydx(3)+ak*y(1,i+1))
     2		+y(1,i+1)*(dxdx(3)+ak*x(1,i+1))
        xd3=dxdx(1)*y(3,i+1)+x(3,i+1)*dydx(1)
        x2=ala2mu*xd1+amu*(xd2+2.0d0*xd3)
c
	ans=ans+0.5d0*(x1+x2)*h2(i)
      enddo
      return
      end
c--------------------------------------------------------------------------
      subroutine  egn_knots(ny,h,k,knot)
c
c 	from input ny and h()
c	knot numbers for eigenfunctions
c	for saito's algorithm are returned
c	
c	kk : total number of knots
c	knot() : contains layer numbers for eigenfunction output
c	referred to the original depth knots
c
c	i=knot(j) for j=1:k
c	give the depths (radius) that eigenfunctions are computed
c
      implicit real*8(a-h,o-z)
      dimension  h(1),knot(1)
c
      k=0
      i  =-1
      j  = 0
c
    1 i  = i + 2
      j  = j + 1
      k = k + 1
      knot(k)=i
      if( j.ge.ny )  return
      if( h(i).le.0.0d0 )  then
        i  = i - 1
      endif
      go to  1
c
      end
c-------------------------------------------------------------------------
      subroutine search(r,dum,ny,ians)
      implicit real*8(a-h,o-z)
      dimension dum(1)
c
      do i=1,ny
      dd=dabs(r-dum(i))
      if(dd.le.1.0e-7) go to 100
      enddo
      write(6,600)r
  600 format('Search failed', '   r=',e15.7)
      stop
  100 ians=i
      return
      end
c-----------------------------------------------------------------
c interpolation at three internal points
c  r(i) is increasing with depth --- this is actually depth
c  for rayleigh wave version
c
c interpolation to increase four-fold
c from v(i) to dum(i)
c
      subroutine interp(nlay,r,rnew,v,dum,kount)
      implicit real*8(a-h,o-z)
      dimension r(nlay),rnew(nlay),v(1),dum(nlay)
c
      kount=0
      do 10 i=1,nlay-1
      hh=r(i+1)-r(i)
      kount=kount+1
      dum(kount)=v(i)
      rnew(kount)=r(i)
      if(hh.le.1.0e-10) go to 10
      kount=kount+1
      rnew(kount)=r(i)+0.25*hh
      dum(kount)=0.75*v(i)+0.25*v(i+1)
      kount=kount+1
      rnew(kount)=r(i)+0.5*hh
      dum(kount)=0.5*v(i)+0.5*v(i+1)
      kount=kount+1
      rnew(kount)=r(i)+0.75*hh
      dum(kount)=0.25*v(i)+0.75*v(i+1)
   10 continue
c
      kount=kount+1
      rnew(kount)=r(nlay)
      dum(kount)=v(nlay)
      return
      end
c--------------------------------------------------------------------------
      subroutine sol_rayl(om,ak)
c
c for a given set of omega (om) and wavenumber (ak),
c get two independent solutions for a Rayleigh 4x4 system
c
c two solutions sol1 and sol2 are returned
c at depths specified by dep(i).
c dep(nlay) is the deepest point below which there is halfspace.
c dep(1) is the surface (zero)
c
      implicit real*8(a-h,o-z)
      common /struc/dep(50000),thk(50000),rho(50000),vp(50000),
     2		    vs(50000),qp(50000),qs(50000),nlay
      common /raysol/ sol1(4,50000),sol2(4,50000)
      dimension y(10),yout(10),dydx(10)
      dimension dyk1(10),dyk2(10),dyk3(10),dyk4(10)
      dimension yk1(10), yk2(10), yk3(10)
c
      neq=4
      pi=3.141592653589793d0
      pi2=2.0*pi
      eps=1.0d-1
c
      am0=rho(nlay)*vs(nlay)*vs(nlay)
c
      xds=om/vs(nlay)
      xdp=om/vp(nlay)
      beta=dsqrt(ak*ak-xds*xds)
      alpha=dsqrt(ak*ak-xdp*xdp)
c
      do icase=1,2
      if(icase.eq.1) then
	y(1)= alpha+ak
      	y(2)= am0*(beta+ak)*(beta+ak)
	y(3)= beta+ak
	y(4)= am0*(2.0*ak*ak+2.*ak*alpha-xds*xds)
	call normalize(4,y)
	do j=1,4
	    sol1(j,nlay)=y(j)
	enddo
      endif
      if(icase.eq.2) then
	yk1(1)= alpha+ak
        yk1(2)= am0*(beta+ak)*(beta+ak)
        yk1(3)= beta+ak
        yk1(4)= am0*(2.0*ak*ak+2.*ak*alpha-xds*xds)
	y(1)= alpha-ak
	y(2)= am0*(beta-ak)*(beta-ak)
	y(3)= ak-beta
	y(4)= am0*(xds*xds-xdp*xdp-(xdp*xdp/(alpha+ak)**2))
	call orthogonal(4,yk1,y)
	do j=1,4
            sol2(j,nlay)=y(j)
        enddo
      endif
c
      do ii=1,nlay-1
      ilay=nlay-ii
      nint=thk(ilay)/0.5
      if(nint.lt.1) nint=1
      h=thk(ilay)/nint
      hh=0.5d0*h
      h6=h/6.0d0
c
      rho1=rho(ilay+1)
      rho2=rho(ilay)
      vp1=vp(ilay+1)
      vp2=vp(ilay)
      vs1=vs(ilay+1)
      vs2=vs(ilay)
c
      anint=nint
      do i=1,nint
	rat=(i-1)/anint
	arat=1.0d0-rat
	rrho=arat*rho1+rat*rho2
	vvp=arat*vp1+rat*vp2
	vvs=arat*vs1+rat*vs2
	rom2=rrho*om*om
	amu=rrho*vvs*vvs
	ala2mu=rrho*vvp*vvp
	ala=ala2mu-2.0d0*amu
	call dy_rayl(ak,rom2,ala,amu,y,dyk1)
	do j=1,neq
		yk1(j)=y(j)+hh*dyk1(j)
	enddo
	rat=(i-0.5d0)/anint
	arat=1.0d0-rat
	rrho=arat*rho1+rat*rho2
        vvp=arat*vp1+rat*vp2
        vvs=arat*vs1+rat*vs2
        rom2=rrho*om*om
        amu=rrho*vvs*vvs
        ala2mu=rrho*vvp*vvp
        ala=ala2mu-2.0d0*amu
	call dy_rayl(ak,rom2,ala,amu,yk1,dyk2)
	do j=1,neq
                yk2(j)=y(j)+hh*dyk2(j)
        enddo
	call dy_rayl(ak,rom2,ala,amu,yk2,dyk3)
        do j=1,neq
                yk3(j)=y(j)+h*dyk3(j)
        enddo
	rat=i/anint
        arat=1.0d0-rat
        rrho=arat*rho1+rat*rho2
        vvp=arat*vp1+rat*vp2
        vvs=arat*vs1+rat*vs2
        rom2=rrho*om*om
        amu=rrho*vvs*vvs
        ala2mu=rrho*vvp*vvp
        ala=ala2mu-2.0d0*amu
	call dy_rayl(ak,rom2,ala,amu,yk3,dyk4)
	do j=1,neq
          y(j)=y(j)+h6*(dyk1(j)+dyk4(j)+2.0d0*(dyk2(j)+dyk3(j)))
        enddo
      enddo
c
      do i=1,neq
	if(icase.eq.1) sol1(i,ilay)=y(i)
	if(icase.eq.2) sol2(i,ilay)=y(i)
      enddo
c
      enddo
      enddo
c
      anorm1=sol1(1,1)
      anorm2=sol2(1,1)
      do i=1,nlay
	do j=1,4
	   sol1(j,i)=sol1(j,i)/anorm1
	   sol2(j,i)=sol2(j,i)/anorm2
	enddo
      enddo
c
      return
      end
c---------------------------------------------------------------------
      subroutine normalize(n,y)
c
c normalize a vector in y(i) with length n
c
      implicit real*8(a-h,o-z)
      dimension y(1)
      xd=0.0d0
      do i=1,n
	xd=xd+y(i)*y(i)
      enddo
      anorm=dsqrt(xd)
      do i=1,n
	y(i)=y(i)/anorm
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine orthogonal(n,x,y)
c
c make y orthogonal to x and normalize it
c both vectors (x and y) have length n
c
      implicit real*8(a-h,o-z)
      dimension x(1), y(1)
      call normalize(n,x)
      call normalize(n,y)
      xd=0.0
      do i=1,n
	xd=xd+x(i)*y(i)
      enddo
      do i=1,n
	y(i)=y(i)-xd*x(i)
      enddo
      call normalize(n,y)
      xd=0.0
      do i=1,n
	xd=xd+y(i)*x(i)
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine dy_rayl(ak,rom2,ala,amu,y,dydx)
c
c calculation of dy/dx for a 4x4 Rayleigh wave system
c called in the Runge-Kutta integration scheme
c
      implicit real*8(a-h,o-z)
      dimension y(10),dydx(10)
c
      xd=ala+2.0d0*amu
      dydx(1)=y(2)/xd+ak*ala*y(3)/xd
      dydx(2)= -rom2*y(1)+ak*y(4)
      dydx(3)= -ak*y(1)+y(4)/amu
      xdum=xd-ala*ala/xd
      xd2=ak*ak*xdum-rom2
      dydx(4)= -ak*ala*y(2)/xd+xd2*y(3)
c
      return
      end
