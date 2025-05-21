c	Progeram visco3d
	parameter (ncmax=13)
	parameter (ncmay=13)
	parameter (ncmaz=13)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=3)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmay-(ncmay-1))*(N*ncmaz-(ncmaz-1)))
	parameter (N3=N*N*N)
	parameter (kmax=ncmax*ncmay*ncmaz*N3*N3*9)
	character*1 gincl
	character*80 recfil,disfil,progfil
c***
	real*8 phi,dphi
	real*8 phi2,dphi2x,dphi2y,dphi2z,wt2
	parameter (lmax3=(lmax+2)*(lmax+2)*(lmax+2))
	common/basisfns/phi(lmax+2,lmax+2),dphi(lmax+2,lmax+2),
     &	phi2(lmax3,lmax3),dphi2x(lmax3,lmax3),dphi2y(lmax3,lmax3),dphi2z(lmax3,lmax3),
     &	wt2(lmax3)
c	phi(i,j) has phi_i (x_j)
c	dphi(i,j) has [d(phi_i)/dx]|x=x_j
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	common/grdpar/NX,NY,NZ
	real*8 ang
	real*8 vp,vs,rho,eta,mupr,eta1
	real*8 mu0,mu2,tau1,tau2,tes,flj,xj,yj
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
c	xg, yg and zg have the (theta[radians],z[km]) coordinates at each of the global gridpoints
c	assuming dimensions of (dx)x(dy)x(dz) for each cell 
c	The origin of the 3D grid is (theta=dlat,y=0,z=0)
c	(dimensionless x, y and z each run from -1 to +1).
	common/glbxy/xg(nptmax),yg(nptmax),zg(nptmax)
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
	real*8 xgdble,phdble,cdelt,sdelt,sfhi,cfhi,delta,sphi,cphi,phirec,fhi,cpsi,spsi,psi
	real*8 ygdble,sphivl,cphivl,phival
	real*8 dlat,plat,plon,phiref
c	real*8 pi,rad
	common/sphpar/dlat,plat,plon,phiref
c	common/minfo/minc,mmax
	common/ginfo/gincl
c***
	complex*16 d,dsrec,dsrecp,ds,dsvert,dsvrtp
	complex*16 ds1,ds2,ds3
	common/datavec/d(3*nptmax,4)
	parameter (maxrec=3000)

c	gravz has the free air gravity field
c	along a swath of surface points surrounding phi=0.
c	real*8 drho,rlp2,y7g,bigg
c	Assume a maximum spherical harmonic degree / azimuthal order = 2000
c	when considering free air gravity computations.
c	real*8 p(2001*2002/2)
	real*8 theta,ctheta,stheta
c	complex*16 gravz(ncmax*N,41,28),pb(ncmax*N,2001)
c	complex*16 sumg1(2001),sumg2(2001)
c	complex*16 dilap,dilam,bmplup,bmminp,bmplux,bmminx,bmpluz,bmminz
c	parameter (bigg=6.673d-11)

	real*8 slat,slon
	real*8 reclat(maxrec),reclon(maxrec),recz(maxrec)
	real*8 zthres
	dimension ncxs(maxrec),ncys(maxrec),nczs(maxrec)
	dimension xrecd(maxrec),xrec(maxrec),yrec(maxrec),zrec(maxrec)
	dimension phir(maxrec),deltar(maxrec)
	common/vertpr/dsvert(3*nptmax,28),dsvrtp(3*nptmax,28)
	common/msumr/dsrec(3*N*N*N*maxrec,28),dsrecp(3*N*N*N*maxrec,28)
c	common/msum/ds(3*nptmax,28)
c	complex*16 g
c	common/gmat/g(icompr,3*nptmax)
c***
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmay-(ncmay-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	integer*8 Ap
	integer(8), ALLOCATABLE :: Ai(:)
	integer(8), ALLOCATABLE :: jAi(:)
	dimension Ap(nmax)
C***
	real*8 phiv,dphivx,dphivz
	complex*16 s,ui,facp,facm,facg
	real*8 lapl
	real*8 tm0,tm1,tm2,tm1o,tm2o,tm1p,tm2p
	real*8 corrf0(5),corrf(10)
	common/tmvals/tm1,tm2
	complex*16 tmpx(28),tmpy(28),tmpz(28),tmpg(28)
	complex*16 dstrxx(28),dstrxy(28),dstrxz(28),dstryy(28),dstryz(28),dstrzz(28)
	real*8 dspx,dspy,dspz,dspg
	real*8 xtxx,xtxy,xtxz
	real*8 xtyy,xtyz,xtzz
c - - - - - - - - -
	complex*16 sval
	common/flapl/sval(28)
	real*8 s0
	common/s0val/s0
	real*8 eig,vm,bigb
	common/flapl2/eig(12),vm(12,12),bigb(12,12)
c - - - - - - - - -
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
	nmat=(NX+1)*(NY+1)*(NZ+1)*N3*N3*9
	allocate ( Ai(nmat) )
	allocate ( jAi(nmat) )
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
c	of any horizontal slice of the 3D grid.
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
	ys=phirec
		write(6,*)'RECEIVER: phi-coordinate y[receiver]=',ys
		write(6,*)'RECEIVER: z-coordinate z[receiver]=',zs
	xrecd(nr)=xs
c	Determine which cell number the receiver is located in
	ncxs(nr)=0
	ncys(nr)=0
	nczs(nr)=0
	do ncz=1,NZ
	do ncy=1,NY
	do ncx=1,NX
c	lower left corner
	ilx=1
	ily=1
	ilz=1
	il=N*(N*(ilz-1)+ilx-1)+ily
	ig=igrd(ncx,ncy,ncz,il)
	x1=xg(ig)
	y1=yg(ig)
	z1=zg(ig)
c		write(6,*)'SOURCE: cell #',ncx,ncy,ncz,'dimensionless lower left corner=',x1,y1,z1
c	upper right corner
	ilx=N
	ily=N
	ilz=N
	il=N*(N*(ilz-1)+ilx-1)+ily
	ig=igrd(ncx,ncy,ncz,il)
	x2=xg(ig)
	y2=yg(ig)
	z2=zg(ig)
	tx=(xs-x1)*(xs-x2)
	ty=(ys-y1)*(ys-y2)
	tz=(zs-z1)*(zs-z2)
	  if(tx.le.0.0.and.ty.le.0.0.and.tz.le.0.0) then
	  ncxs(nr)=ncx
	  ncys(nr)=ncy
	  nczs(nr)=ncz
	  xrec(nr)=(2./dx(ncx))*(xs-x1)-1.
	  yrec(nr)=(2./dy(ncy))*(ys-y1)-1.
	  zrec(nr)=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'ys,y1=',ys,y1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xrec,zrec=',xrec(nr),zrec(nr)
c	  ilc=N*((N/2)-1)+(N/2)
c	  igc=igrd(ncx,ncz,ilc)
c		write(6,*)'RECEIVER #',nr,': gridpoint of cell center=',igc
	  endif
	enddo
	enddo
	enddo
	write(6,*)'RECEIVER #',nr,': Cell #s of receiver =(',ncxs(nr),ncys(nr),nczs(nr),')'
	if(ncxs(nr).eq.0.or.ncys(nr).eq.0.or.nczs(nr).eq.0) then
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
c - - - - - - - - -
c	Read in various arrays needed for the inverse Laplace transform
	open(2,file='INVLAPL.PARAM')
	rewind(2)
	read(2,*) (sval(j), j=1,28)
	read(2,*) bigb
	read(2,*) vm
	read(2,*) eig
	close(2)
c - - - - - - - - -
c	Find smallest inverse relaxation time
	s0=0.
	do ncz=1,NZ
	do ncy=1,NY
	do ncx=1,NX
	do il=1,N*N*N
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	mu0=0.1*rho(ncx,ncy,ncz,il)*vs(ncx,ncy,ncz,il)**2
	mu2=mupr(ncx,ncy,ncz,il)*mu0/(mu0-mupr(ncx,ncy,ncz,il))
	  tau1=mu0/eta1(ncx,ncy,ncz,il)
          tau2=mu2/eta(ncx,ncy,ncz,il)
cOLD	tes=mu0/eta1(ncx,ncy,ncz,il)
	tes=0.5*(tau1+tau2+mu0/eta(ncx,ncy,ncz,il) + sqrt((tau1+tau2+mu0/eta(ncx,ncy,ncz,il))**2-4.d0*tau1*tau2))
		write(6,*)'il=',il
		write(6,*)'eta=',eta(ncx,ncy,ncz,il),'tau2=',tau2
		write(6,*)'vs=',vs(ncx,ncy,ncz,il),'mupr=',mupr(ncx,ncy,ncz,il)
		write(6,*)'mu0,eta1=',mu0,eta1(ncx,ncy,ncz,il),'tau1=',tau1
		write(6,*)'tes=',tes
	if(tes.gt.s0) s0=tes
c	tes=mu2/eta(ncx,ncy,ncz,il)
c		write(6,*)'mu2,eta=',mu2,eta(ncx,ncy,ncz,il)
c		write(6,*)'tes=',tes
c	if(tes.gt.s0) s0=tes
c	tes=mu0/eta(ncx,ncz,il)
c		write(6,*)'mu2,eta=',mu2,eta(ncx,ncz,il)
c		write(6,*)'tes=',tes
c	if(tes.gt.s0) s0=tes
        enddo
	enddo
	enddo
	enddo
	open(4,file='minimum_relaxation_time.txt')
	write(4,*)'The minimum material relaxation time mu/eta is '
	write(4,*)'tau_0 =',3.16881/s0,' years'
	write(4,*)'Numerical results are accurate out to ~1000*tau_0'
	close(4)	
c - - - - - - - - -
c	Re-scale the inverse Laplace transform arrays that depend on s0
	do j=1,28
	sval(j)=sval(j)*s0
	enddo
	do i=1,12
	eig(i)=eig(i)*(0.5/s0)**4
	enddo
c - - - - - - - - -
	write(6,*)'File for storing displacement spectra at receivers?'
	read(5,10) recfil
	write(6,*)'File for storing displacement spectra at surface?'
	read(5,10) disfil
	write(6,*)'Progress file?'
	read(5,10) progfil
	write(6,*)'Beginning and ending fraction of s-value indices?'
	read(5,*) f1,f2
	is1=int(f1*28.)+1
	is2=int(f2*28.)
	write(6,*)'Doing s-indices ',is1,'to ',is2
cTE
	if(imatr.eq.2.or.imatr.eq.3) then
	write(6,*)'Will calculate inverse LTs after one call to SOURCE'
	s=sval(1)
	call source(s)
	go to 20
	endif
c--
	open(4,file=disfil,form='unformatted')
	open(8,file=recfil,form='unformatted')
c
	do j=is1,is2
	s=sval(j)
	open(2,file=progfil,status='old',access='append')
	write(2,*)'Doing j=',j,'out of',28
	close(2)
c - - -
	call source(s)
	write(6,*)'After source: tm1,tm2=',tm1,tm2
	write(6,*)'entering matrel, j=',j,'out of',28
	write(6,*)'s=',s
cOLD	write(6,*)'KU,KL=',KU,KL
	call matrel(s,nmat,Ap,Ai,jAi)
	write(6,*)'out of matrel'
c**	Store values of Laplace-transformed displacement
c	at the cells containing receivers.
          write(6,*)'Store values of LT displacement at the rec cells'
	do nr=1,nrec
          write(6,*)'Receiver #',nr
	ncx=ncxs(nr)
	ncy=ncys(nr)
	ncz=nczs(nr)
	do ilz=1,N
	do ily=1,N
c		write(6,*)'line A: ncx,ncz=',ncx,ncz
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	ig=igrd(ncx,ncy,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	write(8) d(i1,1),d(i2,1),d(i3,1)
          write(6,*)'ilx,ily,ilz=',ilx,ily,ilz,'di2,di2,di3=',d(i1,1),d(i2,1),d(i3,1)
	enddo
	enddo
	enddo
	enddo
c	Store values of Laplace-transformed displacement
	do ig=1,(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	write(4) d(i1,1),d(i2,1),d(i3,1)
	if(ig.eq.1.or.ig.eq.2308) then
	write(6,*)'After matrel ig=',ig,'s=',s
	write(6,*)'After matrel d(',i1,1,')=',d(i1,1)
	write(6,*)'After matrel d(',i2,1,')=',d(i2,1)
	write(6,*)'After matrel d(',i3,1,')=',d(i3,1)
	endif
	enddo	

	enddo
	close(4)
	close(8)

c	Write out gravity coefficients
c	open(8,file='grav-coeff')
c	do j=is1,is2
c	do indh=1,NX*N
c	do i=1,41
c	write(8,*) gravz(indh,i,j)
c	enddo
c	enddo
c	enddo
c	close(8)

cOLD	if(imatr.ne.2) stop
	stop
c****************
20	continue
	write(6,*)'after 20 continue'
	deallocate (Ai)
	deallocate (jAi)
c*****
c       There is coupling between static and postseismic (exponentially decaying) terms
c       when LAPL is called to evaluate the inverse Laplace transform.  Evaluate the
c       effect of the static term on the postseismic terms.
        do iom=1,28
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
	open(2,file='visco3d-stat-rec.gmt')
	open(8,file='visco3d-post-rec.gmt')
c*
	open(4,file='rec-j-m',form='unformatted')
	rewind(4)
	do j=1,28
	ic=0
	do nr=1,nrec
	do idum=1,N*N*N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	dsrec(i1,j)=ds1
	dsrec(i2,j)=ds2
	dsrec(i3,j)=ds3
c                write(6,*)'READIN: dsrec(',i1,j,')=',dsrec(i1,j)
c                write(6,*)'READIN: dsrec(',i2,j,')=',dsrec(i2,j)
c                write(6,*)'READIN: dsrec(',i3,j,')=',dsrec(i3,j)
	enddo
	enddo
	enddo
	close(4)
c       Evaluate inverse Laplace transform. 
		do nr=1,nrec
		write(6,*)'doing nr=',nr,'out of',nrec
	xs=xrec(nr)
	ys=yrec(nr)
	zs=zrec(nr)
c	Do x-, y-, and z-displacements
	do j=1,28
	tmpx(j)=0.d0
	tmpy(j)=0.d0
	tmpz(j)=0.d0
	enddo
c	For displacements at (xrec,yrec,zrec) in cell (ncxs,ncys,nczs), use Lagrangian interpolation.
c	ncz=nczs(nr)
	ic=N*N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ily=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do j=1,28
	tmpx(j)=tmpx(j)+dsrec(i1,j)*phiv(ilx,ily,ilz,xs,ys,zs)
	tmpy(j)=tmpy(j)+dsrec(i2,j)*phiv(ilx,ily,ilz,xs,ys,zs)
	tmpz(j)=tmpz(j)+dsrec(i3,j)*phiv(ilx,ily,ilz,xs,ys,zs)
c		write(6,*)'ilz=',ilz,'ilx=',ilx
c		write(6,*)'dsrec(',i1,j,')=',dsrec(i1,j)
c		write(6,*)'dsrec(',i2,j,')=',dsrec(i2,j)
c		write(6,*)'dsrec(',i3,j,')=',dsrec(i3,j)
c		write(6,*)'phiv(',ilx,ily,ilz,xs,ys,zs,')=',phiv(ilx,ily,ilz,xs,ys,zs)
c               write(6,*)'latest tmpx(',j,')=',tmpx(j)
c               write(6,*)'latest tmpy(',j,')=',tmpy(j)
c               write(6,*)'latest tmpz(',j,')=',tmpz(j)
c		write(6,*)'------------------'
	enddo
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
	if(imatr.eq.3) stop
c* * *
c	Evaluate static [and postseismic] displacements at all nodes
	open(4,file='displ-j-m',form='unformatted')
	rewind(4)
	igmax=(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
	do j=1,28
	do ig=1,igmax
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	dsvert(i1,j)=ds1
	dsvert(i2,j)=ds2
	dsvert(i3,j)=ds3
c		write(6,*)'readp: dsvert(',i1,j,')=',dsvert(i1,j)
c		write(6,*)'readp: dsvert(',i2,j,')=',dsvert(i2,j)
	enddo
	enddo
	close(4)
	open(2,file='visco3d-stat_all_nodes.gmt')
	open(8,file='visco3d-post_all_nodes.gmt')
c       Evaluate inverse Laplace transform. 
	do ncx=1,NX
	do ncy=1,NY
	do ncz=1,NZ
	do ilz=1,N
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	ig=igrd(ncx,ncy,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	do j=1,28
	tmpx(j)=dsvert(i1,j)
	tmpy(j)=dsvert(i2,j)
	tmpz(j)=dsvert(i3,j)
	enddo
                if(ig.eq.6922) then
                write(6,*)'all_nodes lines: ig=',ig,'tmpx=',tmpx
                write(6,*)'all_nodes lines: ig=',ig,'tmpy=',tmpy
                write(6,*)'all_nodes lines: ig=',ig,'tmpz=',tmpz
                endif
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
c	phival=0.
	xgdble=dble(xg(ig))
	ygdble=dble(yg(ig))
	sphivl=dsin(ygdble)/dsin(xgdble)
	cphivl=(dcos(ygdble)-dcos(xgdble)**2)/(dsin(xgdble)**2)
	phival=datan2(sphivl,cphivl)
	phdble=phiref+phival
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
	enddo
	enddo
	close(2)
	close(8)
c* * *
25	continue
c* * *
10	format(a80)
	end
