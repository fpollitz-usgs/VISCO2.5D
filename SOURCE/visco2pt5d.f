c	Test matrel subroutines
	parameter (ncmax=47)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (N2=N*N)
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
	character*1 gincl
	character*80 recfil,disfil,progfil
c***
	real*8 phi,dphi
	real*8 phi2,dphi2x,dphi2z,wt2
	parameter (lmax2=(lmax+2)*(lmax+2))
	common/basisfns/phi(lmax+2,lmax+2),dphi(lmax+2,lmax+2),
     &	phi2(lmax2,lmax2),dphi2x(lmax2,lmax2),dphi2z(lmax2,lmax2),
     &	wt2(lmax2)
c	phi(i,j) has phi_i (x_j)
c	dphi(i,j) has [d(phi_i)/dx]|x=x_j
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	common/grdpar/NX,NZ
	real*8 ang
	real*8 vp,vs,rho,eta,mupr,eta1
	real*8 mu0,mu2,s0,tau1,tau2,tes,flj,xj,yj
	common/struc/dx(ncmax),dz(ncmaz),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
	real*8 xgdble,phdble,cdelt,sdelt,sfhi,cfhi,delta,sphi,cphi,phirec,fhi,cpsi,spsi,psi
	real*8 dlat,plat,plon,phiref
c	real*8 pi,rad
	common/sphpar/dlat,plat,plon,phiref
	common/minfo/minc,mmax
	common/ginfo/gincl
c***
	complex*16 d,dsrec,dsrecp,ds,dsverp,dsvrpp,dsvert,dsvrtp
	common/datavec/d(3*nptmax,4)
	parameter (maxrec=12000)

c	gravz has the free air gravity field
c	along a swath of surface points surrounding phi=0.
	real*8 drho,rlp2,y7g,bigg
c	Assume a maximum spherical harmonic degree / azimuthal order = 2000
c	when considering free air gravity computations.
	real*8 p(2001*2002/2)
	real*8 theta,ctheta,stheta
	complex*16 gravz(ncmax*N,61,14),pb(ncmax*N,2001)
	complex*16 sumg1(2001),sumg2(2001)
	complex*16 dilap,dilam,bmplup,bmminp,bmplux,bmminx,bmpluz,bmminz
	parameter (bigg=6.673d-11)

	real*8 slat,slon
	real*8 reclat(maxrec),reclon(maxrec),recz(maxrec)
	real*8 zthres
	dimension ncxs(maxrec),nczs(maxrec)
	dimension xrecd(maxrec),xrec(maxrec),zrec(maxrec)
	dimension phir(maxrec),deltar(maxrec)
	dimension dsverp(3*nptmax),dsvrpp(3*nptmax)
	common/vertpr/dsvert(3*nptmax,14),dsvrtp(3*nptmax,14)
	common/msumr/dsrec(3*N*N*maxrec,14),dsrecp(3*N*N*maxrec,14)
	common/msum/ds(3*nptmax,14)
c	complex*16 g
c	common/gmat/g(icompr,3*nptmax)
c***
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
	integer*8 Ap
	integer(8), ALLOCATABLE :: Ai(:)
	integer(8), ALLOCATABLE :: jAi(:)
	complex*16, ALLOCATABLE :: AA0(:) , AA1(:) , AA2(:)
	dimension Ap(nmax)
C***
	real*8 phiv,dphivx,dphivz
	complex*16 s,sval,ui,facp,facm,facg
	real*8 lapl
	real*8 tm0,tm1,tm2,tm1o,tm2o,tm1p,tm2p
	real*8 corrf0(5),corrf(10)
	common/tmvals/tm1,tm2
	complex*16 tmpx(14),tmpy(14),tmpz(14),tmpg(14)
	complex*16 dstrxx(14),dstrxy(14),dstrxz(14),dstryy(14),dstryz(14),dstrzz(14)
	real*8 dspx,dspy,dspz,dspg
	real*8 xtxx,xtxy,xtxz
	real*8 xtyy,xtyz,xtzz
	common/flapl/sval(14)
	real*8 ky
	real*8 twopi,rtwopi,efac
	parameter (efac=dlog(100.d0)/9.d0)
	parameter (bigr=6371.)
c
	pi=3.1415926535
	twopi=2.d0*3.14159265358979d0
	rtwopi=dsqrt(twopi)
	rad=360./real(twopi)
	ui=dcmplx(0.,1.)
	call init
	nmat=(NX+1)*(NZ+1)*N2*N2*9
	allocate ( Ai(nmat) )
	allocate ( jAi(nmat) )
	allocate ( AA0(nmat) )
	allocate ( AA1(nmat) )
	allocate ( AA2(nmat) )
		write(6,*)'nmat=',nmat
		write(6,*)'AA size of Ai=',size(Ai)
		write(6,*)'AA size of jAi=',size(jAi)
		write(6,*)'AA Ai(20)=',Ai(20)
		write(6,*)'AA jAi(20)=',jAi(20)
cUMF
	write(6,*)'Determine matrix ordering? (yes=1)'
	write(6,*)'(no, and do all inverse LT (vertical slice; receivers; earths surface)=2)'
	write(6,*)'(no, and do inverse LT for receivers only=3)'
	write(6,*)'(no, and do not do LT=any other value)'
	read(5,*) imatr
	if(imatr.eq.1) then
	write(6,*)'calling sporder'
	call sporder(nmat,Ap,Ai,jAi)
	stop
	else
	write(6,*)'Reading in previously calculated matrix ordering'
	open(2,file='sporder.out',form='unformatted')
	read(2) Ap
	read(2) Ai
	read(2) jAi
	close(2)
	endif
c-- 
c	Read in lat,lon,dep of receivers
	open(2,file='receivers-latlondep.txt',status='old',err=25)
	read(2,*) nrec
c	Receivers cannot be located at depth 0 or above Earth's surface.  
c	Introduce a small distance below the surface as a limiting depth.
	zthres=dble(dz(NZ))*1.d-4
	do nr=1,nrec
	read(2,*) alatin,alonin,depin
	reclat(nr)=(pi/2.d0-alatin/rad)
	reclon(nr)=alonin/rad
	recz(nr)=depin
	if(recz(nr).gt.-zthres) recz(nr)=-zthres
	enddo
	close(2)
	write(6,*)'line C'
c
	  do nr=1,nrec
	slat=reclat(nr)
	slon=reclon(nr)
		write(6,*)'receiver colat,lon=',slat,slon
		write(6,*)'pole colat,lon=',plat,plon
	zs=recz(nr)
c*	Determine angular distance (rad) and azimuth (rad CCW from due N)
c	of the receiver from the pole of the spherical coordinate system.
c	the azimuth is referenced to the azimuth of the origin (theta=dlat,z=0)
c	of the 2D grid.
        cdelt=dcos(slat)*dcos(plat)+dsin(slat)*dsin(plat)*dcos(plon-slon)
        delta=dacos(cdelt)
	deltar(nr)=delta
	sdelt=dsin(delta)
        sphi=dsin(slat)*dsin(plon-slon)/sdelt
        cphi=(dcos(slat)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	phirec=datan2(sphi,cphi)-phiref
	phir(nr)=phirec
c--
		write(6,*)'RECEIVER: phi_receiver=',rad*phirec
c*	The theta-coord of the receiver is the angular distance from (plat,plon)
	xs=real(delta)
		write(6,*)'RECEIVER: theta-coordinate x[receiver]=',xs
		write(6,*)'RECEIVER: z-coordinate z[receiver]=',zs
	xrecd(nr)=xs
c	Determine which cell number the receiver is located in
	ncxs(nr)=0
	nczs(nr)=0
	do ncz=1,NZ
	do ncx=1,NX
c	lower left corner
	ilx=1
	ilz=1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x1=xg(ig)
	z1=zg(ig)
c		write(6,*)'RECEIVER: cell #',ncx,ncz,'dimensionless lower left corner=',x1,z1
c	upper right corner
	ilx=N
	ilz=N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x2=xg(ig)
	z2=zg(ig)
	tx=(xs-x1)*(xs-x2)
	tz=(zs-z1)*(zs-z2)
	  if(tx.le.0.0.and.tz.le.0.0) then
	  ncxs(nr)=ncx
	  nczs(nr)=ncz
	  xrec(nr)=(2./dx(ncx))*(xs-x1)-1.
	  zrec(nr)=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xrec,zrec=',xrec(nr),zrec(nr)
c	  ilc=N*((N/2)-1)+(N/2)
c	  igc=igrd(ncx,ncz,ilc)
c		write(6,*)'RECEIVER #',nr,': gridpoint of cell center=',igc
	  endif
	enddo
	enddo
	write(6,*)'RECEIVER #',nr,': Cell #s of receiver =(',ncxs(nr),nczs(nr),')'
	if(ncxs(nr).eq.0.or.nczs(nr).eq.0) then
	write(6,*)'RECEIVER #',nr,': Receiver location lies outside model domain'
	stop
	endif
	write(6,*)'Local dimensionless receiver coords=(',xrec(nr),zrec(nr),')'
	  enddo
c	if(i.ne.9999) stop
c---
cOLD	write(6,*)'entering diaglim'
cOLD	call diaglim(KU,KL)
cOLD	write(6,*)'out of diaglim'
cOLD	write(6,*)'KU=',KU,'KL=',KL
	ky=0.d0
c	Find smallest inverse relaxation time
	s0=0.
	do ncz=1,NZ
	do ncx=1,NX
	do il=1,N*N
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	mu0=0.1*rho(ncx,ncz,il)*vs(ncx,ncz,il)**2
	mu2=mupr(ncx,ncz,il)*mu0/(mu0-mupr(ncx,ncz,il))
	  tau1=mu0/eta1(ncx,ncz,il)
          tau2=mu2/eta(ncx,ncz,il)
cOLD	tes=mu0/eta1(ncx,ncz,il)
	tes=0.5*(tau1+tau2+mu0/eta(ncx,ncz,il) + sqrt((tau1+tau2+mu0/eta(ncx,ncz,il))**2-4.d0*tau1*tau2))
		write(6,*)'il=',il
		write(6,*)'eta=',eta(ncx,ncz,il),'tau2=',tau2
		write(6,*)'vs=',vs(ncx,ncz,il),'mupr=',mupr(ncx,ncz,il)
		write(6,*)'mu0,eta1=',mu0,eta1(ncx,ncz,il),'tau1=',tau1
		write(6,*)'tes=',tes
	if(tes.gt.s0) s0=tes
c	tes=mu2/eta(ncx,ncz,il)
c		write(6,*)'mu2,eta=',mu2,eta(ncx,ncz,il)
c		write(6,*)'tes=',tes
c	if(tes.gt.s0) s0=tes
c	tes=mu0/eta(ncx,ncz,il)
c		write(6,*)'mu2,eta=',mu2,eta(ncx,ncz,il)
c		write(6,*)'tes=',tes
c	if(tes.gt.s0) s0=tes
        enddo
	enddo
	enddo
c       Determine s-values to be used in sampling Laplace transformed 
c	displacement field.
        do j=1,14
        flj=dble(j)
	xj=s0/5.6-0.055*(flj-1.0)*s0
        yj=0.11*(flj-1.0)*s0
        sval(j)=xj+ui*yj
c--
	write(6,*)'sval(',j,')=',sval(j)
	enddo
c	sval(13)=sval(1)/3.
c        sval(14)=sval(1)/9.
c		pause
c
	write(6,*)'minc,mmax=',minc,mmax
cTE	Reset mmax to do just ky=0 and ky=minc
c	minc=2000
c	mmax=minc
c--
	write(6,*)'Number of wavenumber components=',1+2*(mmax/minc)
	write(6,*)'File for storing displacement spectra at receivers?'
	read(5,10) recfil
	write(6,*)'File for storing displacement spectra at surface?'
	read(5,10) disfil
	write(6,*)'Progress file?'
	read(5,10) progfil
	write(6,*)'Beginning and ending fraction of s-value indices?'
	read(5,*) f1,f2
	is1=int(f1*14.)+1
	is2=int(f2*14.)
	write(6,*)'Doing s-indices ',is1,'to ',is2
cTE
	if(imatr.eq.2.or.imatr.eq.3) then
	write(6,*)'Will calculate inverse LTs after one call to SOURCE'
	s=sval(1)
	ky=0.d0
	call source(s,ky)
	go to 20
	endif
c--
c	Zero out gravity array.
	    do indh=1,NX*N
	    do iphi=1,61
	    do j=1,14
	    gravz(indh,iphi,j)=0.d0
	    enddo
	    enddo
	    enddo
c--
	open(4,file=disfil,form='unformatted')
	open(8,file=recfil,form='unformatted')
c
	do j=is1,is2
	s=sval(j)
c	Zero out the dsverp-array
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	dsverp(i1)=0.d0
	dsverp(i2)=0.d0
	dsverp(i3)=0.d0
	dsvrpp(i1)=0.d0
        dsvrpp(i2)=0.d0
        dsvrpp(i3)=0.d0
	enddo
c--
c	Start loop over azimuthal order number
	  do m=minc,mmax+minc,minc
cTE	should be +
	ky=dble(m-minc)
c--
	write(6,*)'MLOOP: m=',m,'out of',mmax+minc,'with increment',minc
	write(6,*)'MLOOP: ky=',ky
c - - -
	open(2,file=progfil,status='old',access='append')
	write(2,*)'Doing j=',j,'out of',14
	write(2,*) 'MLOOP: m=',m,'out of',mmax+minc,'with increment',minc
	close(2)
c - - -
	call source(s,ky)
	write(6,*)'After source: tm1,tm2=',tm1,tm2
	write(6,*)'entering matrel, j=',j,'out of',14
	write(6,*)'s=',s,'ky=',ky
cOLD	write(6,*)'KU,KL=',KU,KL
	call matrel(s,ky,nmat,Ap,Ai,jAi,AA0,AA1,AA2)
	write(6,*)'out of matrel'
c**	Store values of Laplace-transformed and phi-transformed displacement
c	at the cells containing receivers.
	do nr=1,nrec
	ncz=nczs(nr)
	do ilz=1,N
	ncx=ncxs(nr)
c		write(6,*)'line A: ncx,ncz=',ncx,ncz
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
c	x- and z-displacements are symmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	x- and z-displacements are antisymmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	y-displacement is antisymmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	y-displacement is symmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	The situation is reversed for the components with
c	-ui*sin(ky*phisrc) dependence (i.e., those with d(_,3) and d(_,4) output)
	write(8) d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4),d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4),
     &	d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4)
	if(m.gt.minc) write(8) d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4),-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4),
     &	d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4)
	enddo
	enddo
	enddo
c	Skip the gravity computations if gravity anomaly is not desired.
	if(gincl.ne.'y') go to 45
c--
c	Update the gravz array for azimuthal order m-minc and degrees l=m-minc through mmax-minc.
c	Compute the integrals over theta-r
c	Do spherical harmonic degrees up to min(2000,(8*mmax/5))
	  mval=m-minc
	  lmaxi=min(2000,mmax)
c - - - - - -
	  do l=1,lmaxi+1
	  sumg1(l)=0.d0
	  sumg2(l)=0.d0
	  enddo
	 
	indh=0
	do ncx=1,NX
	do ilx=1,N
	indh=indh+1
	  igtmp=igrd(ncx,NZ,ilx)
c	  Next line: angular distance from pole in radians
	  theta=dble(xg(igtmp))
	  stheta=dsin(theta)
	  ctheta=dcos(theta)
	  call plmon(p, mval, lmaxi, ctheta, 1, 1)
	  l=mval-1
50	  l=l+1
	ind=l*(l+1)/2 + mval + 1
	pb(indh,l+1)=p(ind)*rtwopi
	do ncz=1,NZ

	do ilz=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1

c	First the contribution to free air gravity from vertical uplift.
c	Next line: dimensionless radius
	y7g=dble(l+2)*dble(log(1.0+zg(ig)/bigr))
	rlp2=0.d0
	if(y7g.gt.-20.d0) rlp2=dexp(y7g)
c	rlp2 is r**(l+2)
	if(ilz.gt.1) then
c	Use the same cell and the z-point just below the current one to evaluate drho/dr.
	il2=N*(ilz-2)+ilx
	drho=(rho(ncx,ncz,il)-rho(ncx,ncz,il2))/((dz(ncz)/2.d0)*(x(ilz)-x(ilz-1)))
	endif
	if(ilz.eq.1.and.ncz.eq.1) drho=0.d0
	if(ilz.eq.1.and.ncz.gt.1) then
c	Use the z-point just above the current one.
	il2=N*(ilz)+ilx
	drho=rho(ncx,ncz,il2)-rho(ncx,ncz,il)
	endif
c	Factor of 1.e+4 converts the displacement coefficient to units of meters.
	bmpluz=(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*1.d+4
	bmminz=(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*1.d+4
	sumg1(l+1)=sumg1(l+1) + rlp2 * pb(indh,l+1) * bmpluz * drho * wt(ilx) * dble(dx(ncx)/2.)
     &	* stheta
	if(m.gt.minc) sumg2(l+1)=sumg2(l+1) + rlp2 * pb(indh,l+1) * bmminz * drho * wt(ilx) * dble(dx(ncx)/2.)
     &	* stheta

c	Next the contribution to free air gravity from dilatation.
        fac1=1.d0/dble(bigr+zg(ig))
c       For strains at gridpoint #il in cell (ncx,ncz), use Lagrangian interpolation.
c	Start with the eyy contribution (which doesn't need this interpolation)
	bmplup=(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*(1.d+4)*ui*dble(mval)
	bmminp=-(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*(1.d+4)*ui*dble(mval)
	dilap=bmplup/stheta * fac1
	dilam=bmminp/stheta * fac1
c		dilap=0.
c		dilam=0.
        do ikz=1,N
        do ikx=1,N
        k=N*(ikz-1)+ikx
        igk=igrd(ncx,ncz,k)
        k1=3*igk-2
        k2=k1+1
        k3=k2+1
c	Factor of 1.e+4 converts the displacement coefficient to units of meters.
	bmplux=(d(k1,1)+d(k1,2)+d(k1,3)+d(k1,4))*1.d+4
	bmminx=(d(k1,1)-d(k1,2)-d(k1,3)+d(k1,4))*1.d+4
	bmpluz=(d(k3,1)+d(k3,2)+d(k3,3)+d(k3,4))*1.d+4
	bmminz=(d(k3,1)-d(k3,2)-d(k3,3)+d(k3,4))*1.d+4
c	Next the exx contribution to dilatation
        dilap=dilap + bmplux*dphi2x(k,il)*dble(2./dx(ncx)) * fac1 + bmplux*phi2(k,il)*(ctheta/stheta) * fac1
        dilam=dilam + bmminx*dphi2x(k,il)*dble(2./dx(ncx)) * fac1 + bmminx*phi2(k,il)*(ctheta/stheta) * fac1
c	Next the ezz contribution to dilatation
	dilap=dilap + bmpluz*dphi2z(k,il)*dble(2./dz(ncz))
	dilam=dilam + bmminz*dphi2z(k,il)*dble(2./dz(ncz))
        enddo
        enddo
c	Next the 2*u_r/r part of the exx contribution to dilatation
	bmpluz=(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*1.d+4
	bmminz=(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*1.d+4
	dilap=dilap + 2.d0*bmpluz * fac1
	dilam=dilam + 2.d0*bmminz * fac1
cTE
	sumg1(l+1)=sumg1(l+1) + rlp2 * pb(indh,l+1) * dilap * rho(ncx,ncz,il) * wt2(il) * dble(dx(ncx)/2.)
     &	* dble(dz(ncz)/2.) * stheta
	if(m.gt.minc) sumg2(l+1)=sumg2(l+1) + rlp2 * pb(indh,l+1) * dilam * rho(ncx,ncz,il) * wt2(il)
     &	* dble(dx(ncx)/2.) * dble(dz(ncz)/2.) * stheta

	enddo

c	At a horizontal boundary (ilz=N), add a term due to the jump in density at the boundary,
c	which is part of the contribution to free air gravity from vertical uplift
	il=N*(N-1)+ilx
c	Restore values of bmpluz,bmminz that may have been overridden in do-loops above for strains
	ig=igrd(ncx,ncz,il)
	i3=3*ig
	bmpluz=(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*1.d+4
	bmminz=(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*1.d+4
	if(ncz.lt.NZ) then
c	Use the z-point in the cell just above the current one.
	il2=ilx
	drho=rho(ncx,ncz+1,il2)-rho(ncx,ncz,il)
	elseif(dabs(rho(ncx,ncz,il)-1.000d0).lt.1.d-5) then 
c	We're at the free surface -- implement a density jump unless water is in the shallowest layer.
	drho=0.0
	else
	drho=-rho(ncx,ncz,il)
	endif
	sumg1(l+1)=sumg1(l+1) + rlp2 * pb(indh,l+1) * bmpluz * drho * wt(ilx) * dble(dx(ncx)/2.) * stheta
	if(m.gt.minc) sumg2(l+1)=sumg2(l+1) + rlp2 * pb(indh,l+1) * bmminz * drho * wt(ilx) * dble(dx(ncx)/2.) * stheta
	enddo
	  if(l.lt.lmaxi) go to 50

	enddo
	enddo
c - - - - - -
c	Loop over the swath of surface observation points
	ncz=NZ
	ilz=N
	indh=0
	do ncx=1,NX
	do ilx=1,N
	indh=indh+1
c	il=N*(ilz-1)+ilx
c	ig=igrd(ncx,ncz,il)
c	theta=dble(xg(ig))
c	  ctheta=dcos(theta)
c	  call plmon(p, lmaxi, ctheta, 1, 1)
c	Determine the associated normalized Legendre function P_l^(m-minc) (cos theta).
c	(They are already stored in the pb-array.)
	  l=mval-1
60	  l=l+1
c* * * *
c	With facf apply tapering over last half of frequency range.
	fl=real(l)
	flmax=real(lmaxi)
c		write(6,*)'fl,flmax=',fl,flmax
	iw0m=lmaxi/2
	iw0x=lmaxi-iw0m+1
	fiw0m=real(iw0m)
c		write(6,*)'iw0m,iw0x,fiw0m=',iw0m,iw0x,fiw0m
	facf=1.0
c		write(6,*)'l,iw0x=',i,iw0x
	if(l.ge.iw0x) facf=1.0-cos((flmax-fl)/(fiw0m-1.) * real(pi)/2.)
c		write(6,*)'facf=',facf
c* * * *
	  facg = -(dble(l)+1.d0)/(2.d0*dble(l)+1.d0) * (2.d0*twopi) * bigg
	  ind=l*(l+1)/2 + mval + 1
c	  pb(indh,l+1)=p(ind)*rtwopi
	    do iphi=1,61
	    phival=0.20*real(twopi)/(real(minc)*real(dsin(dlat)))*real(iphi-31)/30.
	    ang=dble(phival)*dble(m-minc)
	    facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	    facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	    gravz(indh,iphi,j)=gravz(indh,iphi,j) + facg * pb(indh,l+1) * (sumg1(l+1)*facp + sumg2(l+1)*facm) * dble(facf)
c		write(6,*)'gravz(',indh,iphi,j,')=',gravz(indh,iphi,j),'facf=',facf
	    enddo 
	  if(l.lt.lmaxi) go to 60
	enddo
	enddo
c--
45	continue
c	Store values of Laplace-transformed and phi-transformed displacement
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
c	x- and z-displacements are symmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	x- and z-displacements are antisymmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	y-displacement is antisymmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	y-displacement is symmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	The situation is reversed for the components with
c	-ui*sin(ky*phisrc) dependence (i.e., those with d(_,3) and d(_,4) output)
	write(4) d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4),d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4),
     &	d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4)
	if(m.gt.minc) write(4) d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4),-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4),
     &	d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4)
cOLD	if(m.gt.minc) write(4) d(i1,1)-d(i1,2),-d(i2,1)+d(i2,2),d(i3,1)-d(i3,2)
	if(ig.eq.1.or.ig.eq.2308) then
	write(6,*)'After matrel ig=',ig,'s=',s,'m=',m,'ky=',ky
	write(6,*)'After matrel d(',i1,1,')=',d(i1,1)
	write(6,*)'After matrel d(',i2,1,')=',d(i2,1)
	write(6,*)'After matrel d(',i3,1,')=',d(i3,1)
	write(6,*)'After matrel d(',i1,2,')=',d(i1,2)
	write(6,*)'After matrel d(',i2,2,')=',d(i2,2)
	write(6,*)'After matrel d(',i3,2,')=',d(i3,2)
	write(6,*)'After matrel d(',i1,3,')=',d(i1,3)
	write(6,*)'After matrel d(',i2,3,')=',d(i2,3)
	write(6,*)'After matrel d(',i3,3,')=',d(i3,3)
	write(6,*)'After matrel d(',i1,4,')=',d(i1,4)
	write(6,*)'After matrel d(',i2,4,')=',d(i2,4)
	write(6,*)'After matrel d(',i3,4,')=',d(i3,4)
	endif
	enddo	
c	Update dsverp and dsverm arrays with sum over wavenumber at
c	presumed phi=0.
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	  facp=dble(minc)/twopi
	  facm=dble(minc)/twopi
	dsverp(i1)=dsverp(i1)+(d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4))*dble(minc)/twopi
	dsverp(i2)=dsverp(i2)+(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*dble(minc)/twopi
	dsverp(i3)=dsverp(i3)+(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*dble(minc)/twopi
	dsvrpp(i1)=dsvrpp(i1)+(d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4))*facp*ui*dble(m-minc)
        dsvrpp(i2)=dsvrpp(i2)+(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*facp*ui*dble(m-minc)
        dsvrpp(i3)=dsvrpp(i3)+(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*facp*ui*dble(m-minc)
	if(m.gt.minc) then
	dsverp(i1)=dsverp(i1)+(d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4))*dble(minc)/twopi
	dsverp(i2)=dsverp(i2)+(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*dble(minc)/twopi
	dsverp(i3)=dsverp(i3)+(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*dble(minc)/twopi
	dsvrpp(i1)=dsvrpp(i1)-(d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4))*facm*ui*dble(m-minc)
        dsvrpp(i2)=dsvrpp(i2)-(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*facm*ui*dble(m-minc)
        dsvrpp(i3)=dsvrpp(i3)-(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*facm*ui*dble(m-minc)
	endif
	enddo
	  enddo
c	Write out latest dsverp-array and dsvrpp-array
	if(j.eq.is1) then
	open(2,file='vertp-j',form='unformatted')
	else
	open(2,file='vertp-j',form='unformatted',status='old',access='append')
	endif
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	write(2) dsverp(i1),dsverp(i2),dsverp(i3)
	write(2) dsvrpp(i1),dsvrpp(i2),dsvrpp(i3)
	enddo
	close(2)

	enddo
	close(4)
	close(8)

c	Write out gravity coefficients
	open(8,file='grav-coeff')
	do j=is1,is2
	do indh=1,NX*N
	do i=1,61
	write(8,*) gravz(indh,i,j)
	enddo
	enddo
	enddo
	close(8)

cOLD	if(imatr.ne.2) stop
	stop
c****************
20	continue
	write(6,*)'after 20 continue'
	deallocate (Ai)
	deallocate (jAi)
	deallocate (AA0)
	deallocate (AA1)
	deallocate (AA2)
c*****
c       There is coupling between static and postseismic (exponentially decaying) terms
c       when LAPL is called to evaluate the inverse Laplace transform.  Evaluate the
c       effect of the static term on the postseismic terms.
        do iom=1,14
        tmpx(iom)=1.d0/sval(iom)
        enddo
c       Do postseismic displacements for 10 time intervals
        do np=1,10
        tm1p=tm1*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
        tm2p=tm2*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
        corrf(np)=lapl(tmpx,tm2p)-lapl(tmpx,tm1p)
        write(6,*)'corrf(',np,')=',corrf(np)
        enddo
	tm1o=tm1
	tm2o=tm2
	do np=1,5
	tm1p=dble(np)*tm1o/3.16881d0
	tm2p=dble(np)*tm2o/3.16881d0
        corrf0(np)=lapl(tmpx,tm2p)-lapl(tmpx,tm1p)
        write(6,*)'corrf0(',np,')=',corrf0(np)
        enddo
c*****
c*-*-*-*
c	Evaluate static and postseismic displacements at receiver points
	open(2,file='visco2pt5d-stat-rec.gmt')
	open(8,file='visco2pt5d-post-rec.gmt')
c*
	write(6,*)'entering addmr'
	call addmr(nrec,phir)
	write(6,*)'after addmr'
c       Evaluate inverse Laplace transform. 
		do nr=1,nrec
		write(6,*)'doing nr=',nr,'out of',nrec
	xs=xrec(nr)
	zs=zrec(nr)
c	Do x-, y-, and z-displacements
	do j=1,14
	tmpx(j)=0.d0
	tmpy(j)=0.d0
	tmpz(j)=0.d0
	enddo
c	For displacements at (xrec,zrec) in cell (ncxz,nczs), use Lagrangian interpolation.
c	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do j=1,14
	tmpx(j)=tmpx(j)+dsrec(i1,j)*phiv(ilx,ilz,xs,zs)
	tmpy(j)=tmpy(j)+dsrec(i2,j)*phiv(ilx,ilz,xs,zs)
	tmpz(j)=tmpz(j)+dsrec(i3,j)*phiv(ilx,ilz,xs,zs)
c		write(6,*)'ilz=',ilz,'ilx=',ilx
c		write(6,*)'dsrec(',i1,j,')=',dsrec(i1,j)
c		write(6,*)'phiv(',ilx,ilz,xs,zs,')=',phiv(ilx,ilz,xs,zs)
c		write(6,*)'------------------'
	enddo
	enddo
	enddo
c       Do for times tm0 only (STATIC deformation).
	tm0=0.d0
	dspx=lapl(tmpx,tm0)
	dspy=lapl(tmpy,tm0)
	dspz=lapl(tmpz,tm0)

c	Write out in GMT format using geographic longitude and latitude.
	xgdble=deltar(nr)
	phirec=phir(nr)
	phdble=phiref+dble(phirec)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c		write(2,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon
	write(2,*) tm0,tm0,real(plon-psi)*rad,90.-rad*real(delta),recz(nr),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
c       Do for postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	dspx=lapl(tmpx,tm2p)-lapl(tmpx,tm1p) - corrf(np)*lapl(tmpx,tm0)
        dspy=lapl(tmpy,tm2p)-lapl(tmpy,tm1p) - corrf(np)*lapl(tmpy,tm0)
        dspz=lapl(tmpz,tm2p)-lapl(tmpz,tm1p) - corrf(np)*lapl(tmpz,tm0)
c	Write out in GMT format using geographic longitude and latitude.
	write(8,*) tm1p*3.16881,tm2p*3.16881,real(plon-psi)*rad,90.-rad*real(delta),recz(nr),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
	enddo
		enddo
	close(2)
	close(8)
c*-*-*-*
c	Evaluate static and postseismic strains at receiver points
	open(2,file='visco2pt5d-statstrains-rec.gmt')
	open(8,file='visco2pt5d-poststrains-rec.gmt')
c*
c	write(6,*)'entering addmr'
c	call addmr(nrec,phir)
c	write(6,*)'after addmr'
c       Evaluate inverse Laplace transform. 
		do nr=1,nrec
		write(6,*)'doing nr=',nr,'out of',nrec
	xs=xrec(nr)
	zs=zrec(nr)
	fac1=1.d0/(dble(bigr)+recz(nr))
c       Next tensor strain components in the order 
c	exx=d(x-displacement)/dx
c	exy=0.5*[d(x-displacement)/dy + d(y-displacement)/dx]
c	exz=0.5*[d(x-displacement)/dz + d(z-displacement)/dx]
c	eyy=d(y-displacement)/dy
c	eyz=0.5*[d(y-displacement)/dz + d(z-displacement)/dy]
c	ezz=d(z-displacement)/dz
	do j=1,14
	dstrxx(j)=0.d0
	dstrxy(j)=0.d0
	dstrxz(j)=0.d0
	dstryy(j)=0.d0
	dstryz(j)=0.d0
	dstrzz(j)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
        ncx=ncxs(nr)
        ncz=nczs(nr)
        ic=N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do j=1,14
	dstrxx(j)=dstrxx(j) + dsrec(i1,j)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	dstrxy(j)=dstrxy(j) + 0.5d0*dsrecp(i1,j)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1
     &  + 0.5d0*dsrec(i2,j)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	dstrxz(j)=dstrxz(j) + 0.5d0*dsrec(i1,j)*dphivz(ilx,ilz,xs,zs)*dble(2./dz(ncz))
     &  + 0.5d0*dsrec(i3,j)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	dstryy(j)=dstryy(j) + dsrecp(i2,j)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1
	dstryz(j)=dstryz(j) + 0.5d0*dsrec(i2,j)*dphivz(ilx,ilz,xs,zs)*dble(2./dz(ncz))
     &  + 0.5d0*dsrecp(i3,j)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1
	dstrzz(j)=dstrzz(j) + dsrec(i3,j)*dphivz(ilx,ilz,xs,zs)*dble(2./dz(ncz))
c		if(j.eq.1.and.ilz.eq.2) then
c		write(6,*)'dsrecp(',i1,j,')=',dsrec(i1,j)
c		write(6,*)'phiv=',phiv(ilx,ilz,xs,zs)
c		write(6,*)'dphivx=',dphivx(ilx,ilz,xs,zs)
c		write(6,*)'dx(',ncx,')=',dx(ncx)
c		write(6,*)'dsrec(',i2,j,')=',dsrec(i2,j),'fac1=',fac1
c		write(6,*)'latest dstrxx(',j,')=',dstrxx(j)
c		write(6,*)'----------'
c		endif
	enddo
	enddo
	enddo
c       Do for times tm0 only (STATIC deformation).
	tm0=0.d0
	xtxx=lapl(dstrxx,tm0)
c		write(6,*)'out of lapl: dstrxx=',dstrxx,'xtxx=',xtxx
        xtxy=lapl(dstrxy,tm0)
        xtxz=lapl(dstrxz,tm0)
        xtyy=lapl(dstryy,tm0)
        xtyz=lapl(dstryz,tm0)
        xtzz=lapl(dstrzz,tm0)
c	Write out in GMT format using geographic longitude and latitude.
	xgdble=deltar(nr)
	phirec=phir(nr)
	phdble=phiref+dble(phirec)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c       Rotate strain components into primed coord system:
c       xhat' = -sfhi*xhat + cfhi*yhat
c       yhat' = -cfhi*xhat - sfhi*yhat
c       zhat = zhat
        expxp = xtxx*sfhi**2 - 2.*sfhi*cfhi*xtxy + xtyy*cfhi**2
        eypyp = xtxx*cfhi**2 + 2.*sfhi*cfhi*xtxy + xtyy*sfhi**2
        expyp = (xtxx-xtyy)*sfhi*cfhi + xtxy*(sfhi**2-cfhi**2)
        expzp = -xtxz*sfhi + xtyz*cfhi
        eypzp = -xtxz*cfhi - xtyz*sfhi
        ezpzp = xtzz
c* * *  Units: Original U=F/K displacement has units of 10^4 m; division by 10^3 meters
c       yields a dimensionless factor of 10.  Hence multiplication by 1.e+7
c* * *  to yield strain in units of microstrain.
	write(2,*) tm0,tm0,real(plon-psi)*rad,90.-rad*real(delta),recz(nr),(1.e+7)*expxp,
     &  (1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
c       Do for postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	xtxx=lapl(dstrxx,tm2p)-lapl(dstrxx,tm1p) - corrf(np)*lapl(dstrxx,tm0)
        xtxy=lapl(dstrxy,tm2p)-lapl(dstrxy,tm1p) - corrf(np)*lapl(dstrxy,tm0)
        xtxz=lapl(dstrxz,tm2p)-lapl(dstrxz,tm1p) - corrf(np)*lapl(dstrxz,tm0)
        xtyy=lapl(dstryy,tm2p)-lapl(dstryy,tm1p) - corrf(np)*lapl(dstryy,tm0)
        xtyz=lapl(dstryz,tm2p)-lapl(dstryz,tm1p) - corrf(np)*lapl(dstryz,tm0)
        xtzz=lapl(dstrzz,tm2p)-lapl(dstrzz,tm1p) - corrf(np)*lapl(dstrzz,tm0)
        expxp = xtxx*sfhi**2 - 2.*sfhi*cfhi*xtxy + xtyy*cfhi**2
        eypyp = xtxx*cfhi**2 + 2.*sfhi*cfhi*xtxy + xtyy*sfhi**2
        expyp = (xtxx-xtyy)*sfhi*cfhi + xtxy*(sfhi**2-cfhi**2)
        expzp = -xtxz*sfhi + xtyz*cfhi
        eypzp = -xtxz*cfhi - xtyz*sfhi
        ezpzp = xtzz

c	Write out in GMT format using geographic longitude and latitude.
	write(8,*) tm1p*3.16881,tm2p*3.16881,real(plon-psi)*rad,90.-rad*real(delta),recz(nr),(1.e+7)*expxp,
     &  (1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
	enddo
		enddo
	close(2)
	close(8)
	if(imatr.eq.3) stop
c* * *
c	Evaluate static [and postseismic] displacements on a vertical profile
		write(6,*)'calling readp'
	  call readp
		write(6,*)'out of readp'
	open(2,file='visco2pt5d-stat_vertp.gmt')
	open(8,file='visco2pt5d-post_vertp.gmt')
c       Evaluate inverse Laplace transform. 
	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	do j=1,14
	tmpx(j)=dsvert(i1,j)
	tmpy(j)=dsvert(i2,j)
	tmpz(j)=dsvert(i3,j)
	enddo
c       Do for times tm0 only (STATIC deformation).
	tm0=0.d0
	dspx=lapl(tmpx,tm0)
        dspy=lapl(tmpy,tm0)
        dspz=lapl(tmpz,tm0)
c	(xa[km],za[km])=Cartesian coordinates along phi=0 plane.
c	stheta=sin(xg(ig)-real(dlat))
c	ctheta=cos(xg(ig)-real(dlat))
c	xa=(bigr+zg(ig))*stheta
c	ya=(bigr+zg(ig))*phival
c	za=(bigr+zg(ig))*ctheta
c	Write out in GMT format using geographic longitude and latitude.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
	write(2,*) tm0,tm0,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
c       Do for postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	dspx=lapl(tmpx,tm2p)-lapl(tmpx,tm1p) - corrf(np)*lapl(tmpx,tm0)
        dspy=lapl(tmpy,tm2p)-lapl(tmpy,tm1p) - corrf(np)*lapl(tmpy,tm0)
        dspz=lapl(tmpz,tm2p)-lapl(tmpz,tm1p) - corrf(np)*lapl(tmpz,tm0)

c	Write out in GMT format using geographic longitude and latitude.
	write(8,*) tm1p*3.16881,tm2p*3.16881,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
	enddo

	enddo
	enddo
	enddo
	enddo
	close(2)
	close(8)
c* * *
c	Evaluate static [and postseismic] strains on a vertical profile
c		write(6,*)'calling readp'
c	  call readp
c		write(6,*)'out of readp'
	open(2,file='visco2pt5d-statstrains_vertp.gmt')
	open(8,file='visco2pt5d-poststrains_vertp.gmt')
c       Evaluate inverse Laplace transform. 
	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	xgdble=dble(xg(ig))
        fac1=1.d0/dble(bigr+zg(ig))
c       Next tensor strain components in the order 
c	exx=d(x-displacement)/dx
c	exy=0.5*[d(x-displacement)/dy + d(y-displacement)/dx]
c	exz=0.5*[d(x-displacement)/dz + d(z-displacement)/dx]
c	eyy=d(y-displacement)/dy
c	eyz=0.5*[d(y-displacement)/dz + d(z-displacement)/dy]
c	ezz=d(z-displacement)/dz
        do j=1,14
        dstrxx(j)=0.d0
        dstrxy(j)=0.5d0*dsvrtp(i1,j)/dsin(xgdble) * fac1
        dstrxz(j)=0.d0
        dstryy(j)=dsvrtp(i2,j)/dsin(xgdble) * fac1
	dstryz(j)=0.5d0*dsvrtp(i3,j)/dsin(xgdble) * fac1
        dstrzz(j)=0.d0
        enddo
c       For strains at gridpoint #il in cell (ncx,ncz), use Lagrangian interpolation.
        do ikz=1,N
        do ikx=1,N
        k=N*(ikz-1)+ikx
        igk=igrd(ncx,ncz,k)
        k1=3*igk-2
        k2=k1+1
        k3=k2+1
        do j=1,14
        dstrxx(j)=dstrxx(j) + dsvert(k1,j)*dphi2x(k,il)*dble(2./dx(ncx)) * fac1
	dstrxy(j)=dstrxy(j) + 0.5d0*dsvert(k2,j)*dphi2x(k,il)*dble(2./dx(ncx)) * fac1
	dstrxz(j)=dstrxz(j) + 0.5d0*dsvert(k1,j)*dphi2z(k,il)*dble(2./dz(ncz))
     &  + 0.5d0*dsvert(k3,j)*dphi2x(k,il)*dble(2./dx(ncx)) * fac1
c	No interpolation for yy strain component.
	dstryz(j)=dstryz(j) + 0.5d0*dsvert(k2,j)*dphi2z(k,il)*dble(2./dz(ncz))
	dstrzz(j)=dstrzz(j) + dsvert(k3,j)*dphi2z(k,il)*dble(2./dz(ncz))
        enddo
        enddo
        enddo
c       Do for times tm0 only (STATIC deformation).
	tm0=0.d0
	xtxx=lapl(dstrxx,tm0)
        xtxy=lapl(dstrxy,tm0)
        xtxz=lapl(dstrxz,tm0)
        xtyy=lapl(dstryy,tm0)
        xtyz=lapl(dstryz,tm0)
        xtzz=lapl(dstrzz,tm0)
c	(xa[km],za[km])=Cartesian coordinates along phi=0 plane.
c	stheta=sin(xg(ig)-real(dlat))
c	ctheta=cos(xg(ig)-real(dlat))
c	xa=(bigr+zg(ig))*stheta
c	ya=(bigr+zg(ig))*phival
c	za=(bigr+zg(ig))*ctheta
c	Write out in GMT format using geographic longitude and latitude.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c       Rotate strain components into primed coord system:
c       xhat' = -sfhi*xhat + cfhi*yhat
c       yhat' = -cfhi*xhat - sfhi*yhat
c       zhat = zhat
        expxp = xtxx*sfhi**2 - 2.*sfhi*cfhi*xtxy + xtyy*cfhi**2
        eypyp = xtxx*cfhi**2 + 2.*sfhi*cfhi*xtxy + xtyy*sfhi**2
        expyp = (xtxx-xtyy)*sfhi*cfhi + xtxy*(sfhi**2-cfhi**2)
        expzp = -xtxz*sfhi + xtyz*cfhi
        eypzp = -xtxz*cfhi - xtyz*sfhi
        ezpzp = xtzz
c* * *  Units: Original U=F/K displacement has units of 10^4 m; division by 10^3 meters
c       yields a dimensionless factor of 10.  Hence multiplication by 1.e+7
c* * *  to yield strain in units of microstrain.
	write(2,*) tm0,tm0,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),(1.e+7)*expxp,
     &  (1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
c       Do for postseismic displacements for 10 time intervals
	do np=1,10
	tm1p=tm1*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	tm2p=tm2*(0.15d0*dexp(efac*dble(np-1)))/3.16881d0
	xtxx=lapl(dstrxx,tm2p)-lapl(dstrxx,tm1p) - corrf(np)*lapl(dstrxx,tm0)
        xtxy=lapl(dstrxy,tm2p)-lapl(dstrxy,tm1p) - corrf(np)*lapl(dstrxy,tm0)
        xtxz=lapl(dstrxz,tm2p)-lapl(dstrxz,tm1p) - corrf(np)*lapl(dstrxz,tm0)
        xtyy=lapl(dstryy,tm2p)-lapl(dstryy,tm1p) - corrf(np)*lapl(dstryy,tm0)
        xtyz=lapl(dstryz,tm2p)-lapl(dstryz,tm1p) - corrf(np)*lapl(dstryz,tm0)
        xtzz=lapl(dstrzz,tm2p)-lapl(dstrzz,tm1p) - corrf(np)*lapl(dstrzz,tm0)
        expxp = xtxx*sfhi**2 - 2.*sfhi*cfhi*xtxy + xtyy*cfhi**2
        eypyp = xtxx*cfhi**2 + 2.*sfhi*cfhi*xtxy + xtyy*sfhi**2
        expyp = (xtxx-xtyy)*sfhi*cfhi + xtxy*(sfhi**2-cfhi**2)
        expzp = -xtxz*sfhi + xtyz*cfhi
        eypzp = -xtxz*cfhi - xtyz*sfhi
        ezpzp = xtzz

c	Write out in GMT format using geographic longitude and latitude.
	write(8,*) tm1p*3.16881,tm2p*3.16881,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),(1.e+7)*expxp,
     &  (1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
	enddo

	enddo
	enddo
	enddo
	enddo
	close(2)
	close(8)

c* * *
25	continue
c***
c	Read in gravity coefficients
	open(8,file='grav-coeff')
	rewind(8)
	do j=1,14
	do indh=1,NX*N
	do i=1,61
	read(8,*) gravz(indh,i,j)
	enddo
	enddo
	enddo
	close(8)
c***
c	Evaluate static and postseismic displacements at Earth's surface at longitudes ranging from
c	-0.20 to +0.20 x the spatial periodicity, i.e.
c	-0.20*real(twopi)/(real(minc)*real(dsin(dlat))) to
c	+0.20*real(twopi)/(real(minc)*real(dsin(dlat)))
	open(2,file='visco2pt5d-stat-phival.gmt')
	    do iphi=1,61
	phival=0.20*real(twopi)/(real(minc)*real(dsin(dlat)))*real(iphi-31)/30.
	call addm(phival)
c       Evaluate inverse Laplace transform.  Do coseismic
	tm0=0.d0
c       Do for times tm0 only (STATIC deformation).
	ncz=NZ
	indh=0
	do ncx=1,NX
	ilz=N
	do ilx=1,N
	indh=indh+1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	do j=1,14
	tmpx(j)=ds(i1,j)
	tmpy(j)=ds(i2,j)
	tmpz(j)=ds(i3,j)
	tmpg(j)=gravz(indh,iphi,j)
	enddo
	dspx=lapl(tmpx,tm0)
        dspy=lapl(tmpy,tm0)
        dspz=lapl(tmpz,tm0)
        dspg=lapl(tmpg,tm0)
c	(xa[km],za[km])=Cartesian coordinates along phi=0 plane.
c	stheta=sin(xg(ig)-real(dlat))
c	ctheta=cos(xg(ig)-real(dlat))
c	xa=(bigr+zg(ig))*stheta
c	ya=(bigr+zg(ig))*phival
c	za=(bigr+zg(ig))*ctheta
c	Write out in GMT format using geographic longitude and latitude.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
	write(2,*) tm0,tm0,real(plon-psi)*rad,90.-rad*real(delta),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz),(1.e+8)*real(dspg)
c	(1.e+8)*real(dspg) is the gravity anomaly in mGal.
	enddo
	enddo
	    enddo
	close(2)
c* * *
	tm1o=tm1
	tm2o=tm2
	open(2,file='visco2pt5d-post-phival.gmt')
	    do iphi=1,61
	phival=0.20*real(twopi)/(real(minc)*real(dsin(dlat)))*real(iphi-31)/30.
	call addm(phival)
c       Evaluate inverse Laplace transform.  
c       Do for times tm1 and tm2 (POSTSEISMIC deformation).
	ncz=NZ
	indh=0
	do ncx=1,NX
	ilz=N
	do ilx=1,N
	indh=indh+1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	do j=1,14
	tmpx(j)=ds(i1,j)
	tmpy(j)=ds(i2,j)
	tmpz(j)=ds(i3,j)
	tmpg(j)=gravz(indh,iphi,j)
	enddo
	do np=1,5
	tm1p=dble(np)*tm1o/3.16881d0
	tm2p=dble(np)*tm2o/3.16881d0
	dspx=lapl(tmpx,tm2p)-lapl(tmpx,tm1p) - corrf0(np)*lapl(tmpx,tm0)
        dspy=lapl(tmpy,tm2p)-lapl(tmpy,tm1p) - corrf0(np)*lapl(tmpy,tm0)
        dspz=lapl(tmpz,tm2p)-lapl(tmpz,tm1p) - corrf0(np)*lapl(tmpz,tm0)
        dspg=lapl(tmpg,tm2p)-lapl(tmpg,tm1p) - corrf0(np)*lapl(tmpg,tm0)
c	(xa[km],za[km])=Cartesian coordinates along phi=0 plane.
c	stheta=sin(xg(ig)-real(dlat))
c	ctheta=cos(xg(ig)-real(dlat))
c	xa=(bigr+zg(ig))*stheta
c	ya=(bigr+zg(ig))*phival
c	za=(bigr+zg(ig))*ctheta
c	Write out in GMT format using geographic longitude and latitude.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
	write(2,*) tm1p*3.16881,tm2p*3.16881,real(plon-psi)*rad,90.-rad*real(delta),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz),(1.e+8)*real(dspg)
c	(1.e+8)*real(dspg) is the gravity anomaly in mGal.
	enddo
	enddo
	enddo
	    enddo
	close(2)
c* * *
10	format(a80)
	end
