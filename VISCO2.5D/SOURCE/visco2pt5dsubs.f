	subroutine sporder(Ap,Ai,jAi)
c	For the matrix G to be constructed by subroutine MATREL,
c	characterize it in terms of row-counter matrices Ap and Ai
c	used in UMFPACK.  
c	Use the same loops as in subroutine MATREL to determine this
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (N2=N*N)
	common/grdpar/NX,NZ
	common/glbgrd/igrd(ncmax,ncmax,N*N)
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
	dimension irow(kmax),icol(kmax)
	integer*8 nzmax,nmax
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
	integer*8 Ap,Ai,jAi
	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
c
	k=0
c--	Following lines follow structure of subroutine MATREL
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncx,ncz,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N2) go to 20
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncz)
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
c--
c	Do the elements (i1,j1),(i1,j2),(i1,j3),(i2,j1),(i2,j2),(i2,j3),
c	(i3,j1),(i3,j2),(i3,j3) in that order.
	k=k+1
	irow(k)=i1
	icol(k)=j1
	k=k+1
	irow(k)=i1
	icol(k)=j2
	k=k+1
	irow(k)=i1
	icol(k)=j3
	k=k+1
	irow(k)=i2
	icol(k)=j1
	k=k+1
	irow(k)=i2
	icol(k)=j2
	k=k+1
	irow(k)=i2
	icol(k)=j3
	k=k+1
	irow(k)=i3
	icol(k)=j1
	k=k+1
	irow(k)=i3
	icol(k)=j2
	k=k+1
	irow(k)=i3
	icol(k)=j3
	go to 10
25	continue
	if(i1.eq.j1) then
	k=k+1
	irow(k)=i1
	icol(k)=j1
	endif
	if(i2.eq.j2) then 
	k=k+1
	irow(k)=i2
	icol(k)=j2
	endif
	if(i3.eq.j3) then
	k=k+1
	irow(k)=i3
	icol(k)=j3
	endif
	go to 10
30	continue
c	Do the elements (i1,j1),(i2,j2),(i3,j3) in that order
	k=k+1
	irow(k)=i1
	icol(k)=j1
	k=k+1
	irow(k)=i2
	icol(k)=j2
	k=k+1
	irow(k)=i3
	icol(k)=j3
	go to 10
20	continue
	enddo
	enddo
	enddo
	enddo	
	ktot=k

c	Assume 0-based Ap and Ai as in UMFPACK conventions.
	ncol=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
		write(6,*)'ncol=',ncol
	do ic=1,ncol+1
	Ap(ic)=0
	enddo
c
c	do k=1,kmax
c	Ai(k)=-1
c	enddo
c
	do k=1,ktot
	ir=irow(k)
	ic=icol(k)
		k5000=k-50000*(k/50000)
		if(k5000.eq.0) then
		write(6,*)'SPORDER: k=',k,'out of',ktot,'; kmax=',kmax
		write(6,*)'row #, column #=',ir,ic
		endif
c	See if (ir,ic) is already represented in column j.  If not,
c	then add it to Ap / Ai

	if(Ap(ic+1).eq.Ap(ic)) then
c	Initiate a new column of non-zero entries
	j0=Ap(ic)+1
	Ai(j0)=ir-1
c		write(6,*)'new column: Ai(',j0,')=',Ai(j0)
	do ip=ic+1,ncol+1
	Ap(ip)=Ap(ip)+1
c		write(6,*)'new Ap(',ip,')=',Ap(ip)
	enddo
	go to 15
	endif	

c	Find smallest j s.t. Ai(j) >= (ir-1) 
c		write(6,*)'Find smallest j s.t. Ai(j) >=',ir-1 
	j0=Ap(ic+1)
	do j=Ap(ic+1),Ap(ic)+1,-1
	if(Ai(j).ge.(ir-1)) j0=j
	enddo
c		write(6,*)'k=',k,'smallest j=',j0

	if(Ai(j0).eq.(ir-1)) then
c		write(6,*)'matches pre-existing entry'
	go to 15
	endif

	if(Ai(j0).gt.(ir-1)) then
c		write(6,*)'Ai(j0) >',ir-1
c	Here Ai(j0) > (ir-1)
c	Add new entry to Ai at index j0; shift everything to its right.
	  do j1=Ap(ncol+1),j0,-1
	  Ai(j1+1)=Ai(j1)
	  enddo
	Ai(j0)=ir-1
	do ip=ic+1,ncol+1
	Ap(ip)=Ap(ip)+1
	enddo
	go to 15
	endif	

c	Here Ai(j0) < (ir-1)
c		write(6,*)'Ai(j0) <',ir-1
c	Add new entry to Ai at index j0+1; shift everything to its right.
c		write(6,*)'Ap(',ic+1,')=',Ap(ic+1),'; j0=',j0
c		write(6,*)'Expect above two numbers to be equal'
	  if(j0.lt.Ap(ncol+1)) then
	  do j1=Ap(ncol+1),j0+1,-1
	  Ai(j1+1)=Ai(j1)
	  enddo
	  endif
	Ai(j0+1)=ir-1
	do ip=ic+1,ncol+1
	Ap(ip)=Ap(ip)+1
	enddo

15	continue
c		if(k.eq.1200) then
c		write(6,*)'After doing k=1200:'
c		write(6,*)'Ap(1 thru 300)='
c		do kdum=1,300
c		write(6,*) Ap(kdum)
c		enddo
c		write(6,*)'Ai(1 thru 2400)='
c		do kdum=1,2400
c		write(6,*) Ai(kdum)
c		enddo
c		endif
c		write(6,*)'-------'
	enddo
c		write(6,*)'ncol=',ncol
c		write(6,*)'Final Ap(1 thru 1200)='
c		do kdum=1,ncol
c		write(6,*) Ap(kdum)
c		enddo
c		write(6,*)'Final Ai(1 thru 2400)='
c		do kdum=1,Ap(ncol+1)
c		write(6,*) Ai(kdum)
c		enddo

c	Assign an element of Ai to every G-matrix update #k
c	The jAi-index will be used to assign values to the AA array in MATREL
	do k=1,ktot
	ir=irow(k)
	ic=icol(k)
	do j=Ap(ic)+1,Ap(ic+1)
	if(Ai(j).eq.(ir-1)) jAi(k)=j
	enddo
c		write(6,*)'Ai element # associated with k=',k,'=',jAi(k)
	enddo
c	Write out Ap,Ai,jAi
	open(2,file='sporder.out',form='unformatted')
	write(2) Ap
	write(2) Ai
	write(2) jAi
	close(2)

	return
	end

	subroutine globgrd(dlat)
c	Determine global grid numbers as a function of
c	(ncx,ncz), cell number indices
c	and local grid number il = N*(ilz-1) + ilx,
c	where (ilx,ilz) are local grid number indices
c
c	OUTPUT
c	igrd(ncx,ncz,il) array of global grid numbers
c
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dz(ncmax),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
c	xg and zg have the (theta[radians],z[km]) coordinates at each of the global gridpoints
c	assuming dimensions of (dx)x(dz) for each cell 
c	The origin of the 2D grid is (theta=dlat,z=0)
c	(dimensionless x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
c
	common/grdpar/NX,NZ
	common/glbgrd/igrd(ncmax,ncmax,N*N)
	real*8 dlat
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c*
	sumzn=0.
	do iz=1,NZ
	sumzn=sumzn+dz(iz)/2.
	enddo
c*
	ig=0
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	inew=1
c
c	When ilx=1 and ncx>1, then reset ilx-->N and ncx-->ncx-1
	if(ilx.eq.1.and.ncx.gt.1) then
	il1=N
	nc1=ncx-1
	il0=N*(ilz-1)+il1
	igrd(ncx,ncz,il)=igrd(nc1,ncz,il0)
c		write(6,*)'ilx=1: igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	inew=0
	endif
c
c	When ilz=1 and ncz>1, then reset ilz-->N and ncz-->ncz-1
	if(ilz.eq.1.and.ncz.gt.1) then
	il2=N
	nc2=ncz-1
	il0=N*(il2-1)+ilx
	igrd(ncx,ncz,il)=igrd(ncx,nc2,il0)
c		write(6,*)'ilz=1: igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	inew=0
	endif
c
	if(inew.eq.1) then
	ig=ig+1
	igrd(ncx,ncz,il)=ig
		write(6,*)'igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	sumx=real(dlat)/2.
	do ix=1,ncx
	sumx=sumx+dx(ix)/2.
	enddo
	xg(ig)=2.*sumx + (dx(ncx)/2.)*(real(x(ilx))-1.)
	sumz=0.
	do iz=1,ncz
	sumz=sumz+dz(iz)/2.
	enddo
	zg(ig)=2.*(sumz-sumzn) + (dz(ncz)/2.)*(real(x(ilz))-1.)
	endif
		write(6,*)'----------------'
c
	enddo
	enddo
	enddo
	enddo	

	return
	end

	subroutine init
	character*1 gincl
	character*80 b80,sfile
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c***
	common/grdpar/NX,NZ
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dz(ncmax),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
	common/glbgrd/igrd(ncmax,ncmax,N*N)
	real*8 elat,elon,plat,plon,phiref,pi,rad
	real*8 dlat,dlon,cdelt,sdelt,sphi,cphi,delta,phi
	common/sphpar/dlat,plat,plon,phiref
	common/minfo/minc,mmax
	common/ginfo/gincl
	real*8 g0
	common/gval/g0
	parameter (m0max=100)
	parameter (bigr=6371.)
c***
	pi=3.14159265358979d0
	rad=180.d0/pi
	write(6,*)'INIT: Reading in simulation info'
	open(2,file='simulation-spherical.info')
	rewind(2)
	read(2,5) b80
	read(2,*) elat,elon
	read(2,5) b80
	read(2,*) dlat,dlon
c*	Determine (colat,lon) of pole of the spherical coordinate system.
	elat=pi/2.d0-elat/rad
	elon=elon/rad
	dlat=dlat/rad
	dlon=dlon/rad
        cdelt=dcos(elat)*dcos(dlat)+dsin(elat)*dsin(dlat)*dcos(dlon)
		write(6,*)'elat=',elat,'dlon=',dlon
		write(6,*)'cdelt=',cdelt
        delta=dacos(cdelt)
		write(6,*)'delta=',delta
	sdelt=dsin(delta)
        sphi=dsin(dlat)*dsin(dlon)/sdelt
        cphi=(dcos(dlat)-dcos(elat)*cdelt)/(dsin(elat)*sdelt)
	phi=datan2(sphi,cphi)
	plat=delta
	plon=elon+phi
		write(6,*)'INIT: geographic plat,plon=',90.-rad*plat,rad*plon
cTE
c	Determine azimuth from the pole of the spherical coordinate system
c	to the reference point (theta=dlat,z=0) of the 2D grid in rad CCW from due N.
        sphi=dsin(elat)*dsin(dlon)/sdelt
        cphi=(dcos(elat)-dcos(dlat)*cdelt)/(dsin(dlat)*sdelt)
	phiref=datan2(sphi,cphi)
		write(6,*)'INIT: phi_ref=',rad*phiref
c*	
	read(2,5) b80
	read(2,*) NX
	read(2,5) b80
	read(2,*) (dx(i), i=1,NX)
	read(2,5) b80
	read(2,*) NZ
	read(2,5) b80
	read(2,*) (dz(i), i=1,NZ)
	read(2,5) b80
	read(2,*) ydist,wavmin
c	Structure and deformation is periodic with spatial periodicity ydist km, hence
	ydist=ydist*real(dsin(dlat))
	wavmin=wavmin*real(dsin(dlat))
	minc=int((2.*3.1415926535*bigr)/ydist)
	mmax=int((2.*3.1415926535*bigr)/wavmin)
	read(2,5) b80
	read(2,5) sfile
	read(2,5) b80
	read(2,*) g0	
	g0=g0*1.d-4
	  write(6,*)'g0=',g0
	  write(6,*)'Number of cells in theta-direction=',NX
	  write(6,*)'Number of cells in z-direction=',NZ
	  write(6,*)'Physical length of cell along theta-direction=',dx
	  write(6,*)'Physical length of cell along z-direction=',dz
	  write(6,*)'Structural parameter file is'
	  write(6,5) sfile
	gincl='n'
	read(2,5,end=45) b80
	read(2,10) gincl
45	continue
	close(2)
	write(6,*)'INIT: Finished reading in simulation info'
	write(6,*)'INIT: calling basis'
	call basis
	write(6,*)'INIT: out of basis'
	write(6,*)'INIT: calling globgrid'
	call globgrd(dlat)
	write(6,*)'INIT: out of globgrid'
	write(6,*)'INIT: calling gridgeom'
	call gridgeom(sfile)
	write(6,*)'INIT: out of gridgeom'
5	format(a80)
10	format(a1)
	return
	end

	subroutine matrel(s,ky,Ap,Ai,jAi)
cccc	subroutine matrel(ome,ky,KU,KL)
c	Generate elements of G (and d) for given 
c	Laplace-transform parameter s
cccc	angular frequency ome
c	and y-wavenumber ky. KU and KL are the number of diagonals
c	above and below the central diagonal, respectively, with non-zero 
c	matrix elements.
c		real*8 xtmp
	character*80 dum80
	complex*16 s
	real*8 ky
cccc	real*8 ome,ky
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
	parameter (N2=N*N)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	There are a maximum of nptmax GLL points in the 2D global grid,
c	each with three displacement components.
	common/glbxy/xg(nptmax),zg(nptmax)
c	igc stored in common will have the grid # of the last point
c	determined on the source plane -- representative of the source location.
	common/isrc/igc
c****
c	Every column of G will be expressed in a sparse matrix format
c	3*(N*nptmax+N*N) is a slightly high estimate of KL and/or KU,
c	so icompr below is slightly higher than 2*KL+KU+1,
c	the intended leading dimension of G.
cOLD	parameter (icompr=3*(N*nptmax+N*N) * 3 + 1)
	complex*16 g
cOLD	dimension g(icompr,3*nptmax*nptmax)
	complex*16 d
	common/datavec/d(3*nptmax,4)
cOLD	integer KL,KU,M,NRHS,LDAB,INFO
cOLD	integer IPIV(3*nptmax*nptmax)
c****
	real*8 fac,fac1,fac2,cott
	complex*16 ui
	parameter (bigr=6371.)
c***
	common/grdpar/NX,NZ
	real*8 vpval,vsval,rhoval,mupval,etaval,et1val,mu2val,tau1,tau2
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dz(ncmax),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
	real*8 g0
	common/gval/g0
	common/glbgrd/igrd(ncmax,ncmax,N*N)
c*****
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
	integer*8 nzmax,nmax
cOLD        parameter (nzmax = 9500000, nmax = 160000)
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
	integer*8 Ap,Ai,jAi
	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	complex*16 AA (nzmax) , XX (nmax), BB (nmax), r (nmax)
	double precision Ax (nzmax), xumf (nmax), b (nmax)
	double precision control(20),info(90)
	double precision Az (nzmax), xz (nmax), bz (nmax)
	integer*8 numf , numfz , numeric, symbolic, status, sys, filenum
c*****
	real*8 kappa0,mu0
	complex*16 lam2mu,lam,mu
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
	ui=dcmplx(0.d0,1.d0)

c	Zero out the elements of AA
	ncol=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
		write(6,*)'MATREL: ncol=',ncol
	write(6,*)'Zero-ing out elements of AA'
	write(6,*)'Ap(',ncol+1,')=',Ap(ncol+1)
	do j=1,Ap(ncol+1)
	AA(j)=dcmplx(0.,0.)
	enddo
	write(6,*)'Done zero-ing out elements of AA'

	write(6,*)'Constructing G'
	k0=0
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncx,ncz,il)-2
	i2=i1+1
	i3=i2+1
c		write(6,*)'Working on rows',i1,i2,i3
c		write(6,*)'out of',3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
c		write(6,*)'      '
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N2) go to 20
	jlz=(jl-1)/N+1
	jlx=jl-N*(jlz-1)
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncz)
cNOTE	If the test function corresponds to
c	the bottom of the grid, then enforce zero-displacement BC
c	This assumes that the i1,i2,i3 elements of the data vector d are zero,
c	i.e. the source is not in one of the bottom cells.
c	if(ilz.eq.1.and.ncz.eq.1) then
c	write(6,*)'ilz=',ilz,'ncz=',ncz
c	write(6,*)'ncx,ncz=',ncx,ncz
c	write(6,*)'will go to 25'
c	endif
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
c--
	l=0
15	l=l+1
	if(l.gt.N2) then
	k0=k
c		write(6,*)'Check k0: go to 10 w/k0=',k
	go to 10
	endif
	lz=(l-1)/N+1
	lx=l-N*(lz-1)
		if(ilx.ne.lx.and.ilz.ne.lz) go to 15
c**
c	write(6,*)'ilx,lx=(',ilx,lx,')','ilz,lz=(',ilz,lz,')'
c	write(6,*)'dphi2x(',il,l,')=',dphi2x(il,l)
c	write(6,*)'dphi2z(',il,l,')=',dphi2z(il,l)
c	write(6,*)'-----------------------------'

	ig=igrd(ncx,ncz,l)
	rhoval=rho(ncx,ncz,l)
	vsval=vs(ncx,ncz,l)
	vpval=vp(ncx,ncz,l)
	mupval=mupr(ncx,ncz,l)
	etaval=eta(ncx,ncz,l)
	et1val=eta1(ncx,ncz,l)
c		write(6,*)'MATREL: ig=',ig
c		write(6,*)'rhoval,vsval,vpval,mupval,etaval,et1val='
c		write(6,*) rhoval,vsval,vpval,mupval,etaval,et1val
c		write(6,*)'--------------'

	mu0=0.1*rhoval*vsval**2
	kappa0=0.1*rhoval*vpval**2 - (4.d0/3.d0)*mu0
c	Apply correspondence principle for Burgers body.
	mu2val=mupval*mu0/(mu0-mupval)
	  tau1=mu0/et1val
          tau2=mu2val/etaval
	  mu=mu0*s*(s+tau2)/((s+tau2)*(s+tau1)+mu0*s/etaval)
c		if(ncz.eq.17) then
c		write(6,*)'ncz,lz=',ncz,lz
c		write(6,*)'mu0=',mu0,'mupval=',mupval,'mu2val=',mu2val
c		write(6,*)'et1val,etaval=',et1val,etaval
c		write(6,*)'tau1=',tau1,'tau2=',tau2
c		write(6,*)'mu=',mu
c		endif
cOLD	mu=mu0*s/(s+mu0/etaval)
	lam=kappa0-(2.d0/3.d0)*mu
	lam2mu=lam+2.d0*mu

c	Eqn 7 of notes
c	darea=dble((dx/2.))*dble((dz/2.)) * (bigr+zg(ig))
	fac=bigr+zg(ig)
	fac1=1.d0/(bigr+zg(ig))
	fac2=fac1 * dble(1.0/sin(xg(ig)))
c		write(6,*)'fac1,fac2=',fac1,fac2,'xg(',ig,')=',xg(ig)
	cott=dble(cos(xg(ig))/sin(xg(ig)))
c		write(6,*)'cott=',cott
c	Terms arising from the theta- and r-parts of the gradient of the stress tensor
c	after integration by parts in theta and r.
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam2mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1
     &					    - wt2(l)*dphi2x(il,l)*lam*phi2(jl,l)*cott * dble(dz(ncz)/2.) * fac1
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac
     &		                            + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dble(dx(ncx)/2.0)

c		write(6,*)'line 1: AA(',j,')=',AA(j)
c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam*ui*ky*phi2(jl,l) * dble(dz(ncz)/2.) * fac2 
c		write(6,*)'line 2: AA(',j,')=',AA(j)

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam2mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac1
     &					    - wt2(l)*dphi2x(il,l)*lam*dphi2z(jl,l) * 1.d0
     &					    - wt2(l)*dphi2x(il,l)*lam*phi2(jl,l) * dble(dz(ncz)/2.) * fac1
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2x(jl,l) * 1.d0
	AA(j)=AA(j)
     &			    - wt2(l)*phi2(il,l)*rhoval*g0*dphi2x(jl,l) * dble(dz(ncz)/2.)
c		write(6,*)'line 3: AA(',j,')=',AA(j)
c	Above line: gravitational term

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*ui*ky*phi2(jl,l) * dble(dz(ncz)/2.) * fac2

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1
     &					    + wt2(l)*dphi2x(il,l)*mu*phi2(jl,l)*cott * dble(dz(ncz)/2.) * fac1
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac
     &					    + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dble(dx(ncx)/2.)

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2z(il,l)*mu*ui*ky*phi2(jl,l) * dble(dx(ncx)/2.) * fac2 * fac
	AA(j)=AA(j)
     &			    - wt2(l)*phi2(il,l)*rhoval*g0*ui*ky*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
c	Above line: gravitational term

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2z(jl,l) * 1.d0
     &					    + wt2(l)*dphi2x(il,l)*mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac1
     &		                            - wt2(l)*dphi2z(il,l)*lam*dphi2x(jl,l) * 1.d0
     &					    - wt2(l)*dphi2z(il,l)*lam*phi2(jl,l)*cott * dble(dx(ncx)/2.)
	AA(j)=AA(j)
     &			    + wt2(l)*phi2(il,l)*rhoval*g0*dphi2x(jl,l) * dble(dz(ncz)/2.)
c	Above line: gravitational term

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2z(il,l)*lam*ui*ky*phi2(jl,l) * dble(dx(ncx)/2.) * fac2 * fac
	AA(j)=AA(j)
     &			    + wt2(l)*phi2(il,l)*rhoval*g0*ui*ky*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
c	Above line: gravitational term

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1
     &		                            - wt2(l)*dphi2z(il,l)*lam2mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac
     &		                            - wt2(l)*dphi2z(il,l)*2.d0*lam*phi2(jl,l) * dble(dx(ncx)/2.)
	AA(j)=AA(j)
     &			    + wt2(l)*phi2(il,l)*rhoval*g0*2.d0*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
c	Above line: gravitational term

c	Terms arising from the phi-part of the gradient of the stress tensor.
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*dphi2x(jl,l) * dble(dz(ncz)/2.)
     &	+ wt2(l)*phi2(il,l)*mu*dphi2z(jl,l) * dble(dx(ncx)/2.)
     &	- wt2(l)*phi2(il,l)*mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac1
c	Above 2 lines: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *mu*dphi2x(jl,l) * dble(dz(ncz)/2.)
     &	- wt2(l)*phi2(il,l)*fac2*ui*ky *mu*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*dphi2z(jl,l) * dble(dx(ncx)/2.) * fac 
     &	+ wt2(l)*phi2(il,l)*mu*dphi2x(jl,l) * dble(dz(ncz)/2.) * fac1
c	Above 1 line: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam2mu*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam*dphi2x(jl,l) * dble(dz(ncz)/2.)
     &	+ wt2(l)*phi2(il,l)*fac1*cott *mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam2mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
     &	+ wt2(l)*phi2(il,l)*fac1*cott *mu*dphi2x(jl,l) * dble(dz(ncz)/2.)
     &	- wt2(l)*phi2(il,l)*fac1*cott *mu*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.)
     &	+ wt2(l)*phi2(il,l)*fac1 *mu*dphi2z(jl,l) * dble(dx(ncx)/2.) * fac
     &	- wt2(l)*phi2(il,l)*fac1 *mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam2mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *lam*dphi2z(jl,l) * dble(dx(ncx)/2.) * fac
     &	+ wt2(l)*phi2(il,l)*fac1 *mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1 *lam*dphi2x(jl,l) * dble(dz(ncz)/2.)
     &	- wt2(l)*phi2(il,l)*lam2mu*dphi2x(jl,l) * dble(dz(ncz)/2.) * fac1
     &	- wt2(l)*phi2(il,l)*lam*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.) * fac1
c	Above 2 lines: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *mu*dphi2z(jl,l) * dble(dx(ncx)/2.) * fac 
     &	- wt2(l)*phi2(il,l)*fac2*ui*ky *mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
     &	- wt2(l)*phi2(il,l)*lam*ui*ky*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac2 
c	Above 1 line: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*ky *mu*ui*ky*phi2(jl,l) 
     &	* dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac 
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1 *lam*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.)
     &	- wt2(l)*phi2(il,l)*fac1 *lam*dphi2z(jl,l) * dble(dx(ncx)/2.) * fac
     &	- wt2(l)*phi2(il,l)*lam2mu*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac1
     &	- wt2(l)*phi2(il,l)*lam*dphi2z(jl,l) * dble(dx(ncx)/2.)
     &	- wt2(l)*phi2(il,l)*lam*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * fac1
c	Above 3 lines: Additional terms arising from r- and theta-derivatives of unit vectors.

	go to 15
cNOTE	Enforce zero-displacement BC at the bottom of the grid.
25	continue
c	NOTE: if i1=j1, then it must also be true that i2=j2 and i3=j3
	if(i1.eq.j1) then
	k=k0+1
	j=jAi(k)
	AA(j)=1.d0
	k=k+1
	j=jAi(k)
	AA(j)=1.d0
	k=k+1
	j=jAi(k)
	AA(j)=1.d0
	endif
c	if(i1.eq.j1) g(KU+KL+1+i1-j1,j1)=1.d0
c	if(i2.eq.j2) g(KU+KL+1+i2-j2,j2)=1.d0
c	if(i3.eq.j3) g(KU+KL+1+i3-j3,j3)=1.d0
c		write(6,*)'g(',KU+KL+1+i1-j1,j1,')=',g(KU+KL+1+i1-j1,j1)
c		write(6,*)'g(',KU+KL+1+i2-j2,j2,')=',g(KU+KL+1+i2-j2,j2)
c		write(6,*)'g(',KU+KL+1+i3-j3,j3,')=',g(KU+KL+1+i3-j3,j3)
	k0=k
	go to 10
30	continue
	igb=igrd(ncx,ncz,il)
	roff=1.d0/(xg(igb)-xg(igc))
	k=k0+1
	j=jAi(k)
	AA(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))
	k=k+1
	j=jAi(k)
	AA(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))
	k=k+1
	j=jAi(k)
	AA(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))
	k0=k
	go to 10
20	continue

	enddo
	enddo
	enddo
	enddo

	write(6,*)'Done constructing G'
		write(6,*)'MATREL: end value of k=',k0
	numfz=Ap(ncol+1)
	do j=1,numfz
	Ax (j) = dreal(AA(j))
c		write(6,*)'Ax(',j,')=',Ax(j)
        Az (j) = dimag(AA(j))
c		write(6,*)'Az(',j,')=',Az(j)
	enddo

	numf=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))

c       ----------------------------------------------------------------
c       factor the matrix and save to a file
c       ----------------------------------------------------------------

c       set default parameters
	write(6,*)'entering umf4zdef'
        call umf4zdef (control)

c       print control parameters.  set control (1) to 1 to print
c       error messages only
        control (1) = 2
	write(6,*)'entering umf4zpcon'
        call umf4zpcon (control)

c       pre-order and symbolic analysis
	write(6,*)'entering umf4zsym'
        call umf4zsym (numf, numf, Ap, Ai, Ax, Az, symbolic, control, info)

c       print statistics computed so far
c       call umf4zpinf (control, info) could also be done.
        print 80, info (1), info (16),
     $      (info (21) * info (4)) / 2**20,
     $      (info (22) * info (4)) / 2**20,
     $      info (23), info (24), info (25)
80      format ('symbolic analysis:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, ' (sec)'/,
     $      '   estimates (upper bound) for numeric LU:', /,
     $      '   size of LU:    ', f10.2, ' (MB)', /,
     $      '   memory needed: ', f10.2, ' (MB)', /,
     $      '   flop count:    ', e10.2, /
     $      '   nnz (L):       ', f10.0, /
     $      '   nnz (U):       ', f10.0)

c       check umf4zsym error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsym: ', info (1)
            stop
        endif

c       numeric factorization
        call umf4znum (Ap, Ai, Ax, Az, symbolic, numeric, control, info)

c       print statistics for the numeric factorization
c       call umf4zpinf (control, info) could also be done.
        print 90, info (1), info (66),
     $      (info (41) * info (4)) / 2**20,
     $      (info (42) * info (4)) / 2**20,
     $      info (43), info (44), info (45)
90      format ('numeric factorization:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, /,
     $      '   actual numeric LU statistics:', /,
     $      '   size of LU:    ', f10.2, ' (MB)', /,
     $      '   memory needed: ', f10.2, ' (MB)', /,
     $      '   flop count:    ', e10.2, /
     $      '   nnz (L):       ', f10.0, /
     $      '   nnz (U):       ', f10.0)

c       check umf4znum error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4znum: ', info (1)
            stop
        endif

c       save the symbolic analysis to the file s42.umf
c       note that this is not needed until another matrix is
c       factorized, below.
cU	filenum = 42
cU        call umf4zssym (symbolic, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zssym: ', status
cU            stop
cU        endif

c       save the LU factors to the file n0.umf
cU        call umf4zsnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zsnum: ', status
cU            stop
cU        endif

c       free the symbolic analysis
        call umf4zfsym (symbolic)

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c	***** FIRST RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 1st rhs
	do i=1,numf
	BB(i)=d(i,1)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,1)=XX(i)
c	    write(6,*)'First rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** SECOND RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU       if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 2nd rhs
	do i=1,numf
	BB(i)=d(i,2)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,2)=XX(i)
c	    write(6,*)'Second rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** THIRD RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 3rd rhs
	do i=1,numf
	BB(i)=d(i,3)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,3)=XX(i)
c	    write(6,*)'Third rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** FOURTH RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 4th rhs
	do i=1,numf
	BB(i)=d(i,4)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,4)=XX(i)
c	    write(6,*)'Fourth rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c		stop

c	open(2,file='matrel.solution1')
c	write(6,*)'solution vector 1='
c	do i=1,numf
c	write(6,*) i,d(i,1)
c	write(2,*) i,real(d(i,1))
c	enddo
c	close(2)

c	open(2,file='matrel.solution2')
c	write(6,*)'solution vector 2='
c	do i=1,numf
c	write(6,*) i,d(i,2)
c	write(2,*) i,real(d(i,2))
c	enddo
c	close(2)

c	open(2,file='matrel.solution3')
c	write(6,*)'solution vector 3='
c	do i=1,numf
c	write(6,*) i,d(i,3)
c	write(2,*) i,real(d(i,3))
c	enddo
c	close(2)

c	open(2,file='matrel.solution4')
c	write(6,*)'solution vector 4='
c	do i=1,numf
c	write(6,*) i,d(i,4)
c	write(2,*) i,real(d(i,4))
c	enddo
c	close(2)

	return
	end

c=======================================================================
c== resid ==============================================================
c=======================================================================

c Compute the residual, r = Ax-b, its max-norm, and print the max-norm
C Note that A is zero-based.

        subroutine resid (n, nz, Ap, Ai, A, x, b, r)
        integer*8
     $      n, nz, Ap (n+1), Ai (n), j, i, p
        complex*16 A (nz), x (n), b (n), r (n), aij
	double precision rmax

        do 10 i = 1, n
            r (i) = -b (i)
10      continue

        do 30 j = 1,n
            do 20 p = Ap (j) + 1, Ap (j+1)
                i = Ai (p) + 1
                aij = A (p)
                r (i) = r (i) + aij * x (j)
20          continue
30      continue

        rmax = 0
        do 40 i = 1, n
            rmax = max (rmax, abs (r (i)))
40      continue

        print *, 'norm (A*x-b): ', rmax
        return
        end

	subroutine gridgeom(sfile)
	character*80 sfile
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
	parameter (iptmax=100000)
	real*4 muprin
	real*8 vp0(iptmax),vs0(iptmax),rho0(iptmax),
     &  mupr0(iptmax),eta0(iptmax),eta01(iptmax)
	dimension x0(iptmax),z0(iptmax)
	character*80 filin
	character*92 aread
	common/grdpar/NX,NZ
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dz(ncmax),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
	common/glbgrd/igrd(ncmax,ncmax,N*N)
c
	pi=3.1415926535
	rad=180./pi
	open(2,file=sfile)
	rewind(2)
		write(6,*)'GRIDGEOM: sfile=',sfile
	write(6,*)'Reading in structural parameters'
	etala=1.d+13
	i=0
5	read(2,14,end=10) aread 
c		write(6,14) aread
	read(aread,*,end=17) xin,zin,vpin,vsin,rhoin,muprin,etain,eta1in
c		write(6,*)'Used first read'
	i=i+1
	eta01(i)=eta1in
	mupr0(i)=muprin
	go to 12
17	i=i+1
	eta01(i)=1.d+13
	mupr0(i)=0.d0
	read(aread,*) xin,zin,vpin,vsin,rhoin,etain
c		write(6,*)'Used second read'
12	continue
cOLD	i=i+1
c
	if(i.gt.iptmax) then
	write(6,*)'Current number of parameter lines =',i,'for structural parameter input'
	write(6,*)'exceeds limit of iptmax=',iptmax
	stop
	endif
c
	vp0(i)=dble(vpin)
	vs0(i)=dble(vsin)
	rho0(i)=dble(rhoin)
	eta0(i)=etain
c	Convert (xin[deg.],zin[km]) to (xin[km],zin[km])
c	and reference to the initial x0-value
	x0(i)=(xin/rad)*6371.
	x0(i)=x0(i)-x0(1)
	z0(i)=zin
c		write(6,*)'x0 and z0(',i,')=',x0(i),z0(i)
	go to 5
10	continue
	close(2)
14	format(a92)
	imax=i
	write(6,*)'Finished reading structural parameters at',imax,'points'
c	Fill up model grid with values based on closest distance to just read-in values
c	igmax=(N*NX-(NX-1))*(N*NZ-(NZ-1))
c	do ig=1,igmax

	  do ncz=1,NZ
	  do ncx=1,NX
	  do ilz=1,N
	  do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)

	xval=(xg(ig)-xg(1))*6371.
	if(ilx.eq.1) xval=xval+0.01*dx(ncx)*6371.
	if(ilx.eq.N) xval=xval-0.01*dx(ncx)*6371.
	zval=zg(ig)
	if(ilz.eq.1) zval=zval+0.01*dz(ncz)
	if(ilz.eq.N) zval=zval-0.01*dz(ncz)
c	find closest read-in gridpoint to (xg(ig),yg(ig))
	do i=1,imax
	dist=(xval-x0(i))**2 + (zval-z0(i))**2
c		if(ig.eq.94) write(6,*)'A latest dist=',dist
c		if(ig.eq.94) write(6,*)'A xg,zg=',xg(ig),zg(ig)
c		if(ig.eq.94) write(6,*)'A x0,z0=',x0(i),z0(i)
	if(i.eq.1) then
	distm=dist
	im=i
	endif
	if(i.gt.1.and.dist.lt.distm) then
	distm=dist
	im=i
c		if(ig.eq.94) write(6,*)'latest dist=',dist
c		if(ig.eq.94) write(6,*)'xg,zg=',xg(ig),zg(ig)
c		if(ig.eq.94) write(6,*)'x0,z0=',x0(i),z0(i)
c		if(ig.eq.94) write (6,*)' ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^'
	endif
	enddo
c		if(ig.eq.94) write(6,*)'closest i=',im
		write(6,*)'NX=',NX,'NZ=',NZ
		write(6,*)'ig=',ig,'im=',im,'nptmax=',nptmax
	vp(ncx,ncz,il)=vp0(im)
	vs(ncx,ncz,il)=vs0(im)
	rho(ncx,ncz,il)=rho0(im)
	mupr(ncx,ncz,il)=mupr0(im)
	eta(ncx,ncz,il)=eta0(im)
		write(6,*)'GRIDGEOM: eta(',ncx,ncz,il,')=',eta(ncx,ncz,il)
	eta1(ncx,ncz,il)=eta01(im)
		write(6,*)'GRIDGEOM: x,z=',xg(ig),'km',zg(ig),'km'
		write(6,*)'vp and vs(',ncx,ncz,il,')=',vp(ncx,ncz,il),vs(ncx,ncz,il)
		write(6,*)'------------------------'
	  enddo
	  enddo
	  enddo
	  enddo

c	enddo

	return
	end


	subroutine source(s,ky)
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
c	igc stored in common will have the grid # of the last point
c	determined on the source plane -- representative of the source location.
	common/isrc/igc
c
	common/grdpar/NX,NZ
	real*8 mupval,etaval,et1val,mu2val,tau1,tau2
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dz(ncmax),vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),mupr(ncmax,ncmaz,N*N),eta(ncmax,ncmaz,N*N),eta1(ncmax,ncmaz,N*N)
	common/glbgrd/igrd(ncmax,ncmax,N*N)
	real*8 cdip,sdip,c2dip
	real*8 sstr,cstr,s2str,c2str,frak,srak,crak
	real*8 deltl,phil,deltd,delt0
	real*8 mu0
	real*8 shrm1,p1,p2,p3,p4,p5
	real*8 dmax(600),dmin(600),dip(600)
	real*8 flat(600),flon(600),fbigl(600),fstr(600),frake(600),fwt(600)
	real*8 tm1,tm2,vmult
	real*8 plat,plon,phiref,pi,rad
	real*8 dlat,cdelt,sdelt,delta,sphi,cphi,phisrc,phi1
	common/sphpar/dlat,plat,plon,phiref
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	real*8 mxx,myy,mzz,mxy,mxz,myz
	real*8 mtt,mpp,mtp,mrt,mrp,mrr
	complex*16 d
	common/datavec/d(3*nptmax,4)
	common/tmvals/tm1,tm2
c***
	real*8 slat,slon
	real*8 phis,dphisx,dphisz
	real*8 phiv,dphivx,dphivz
	real*8 ky,fac1,fac2,cott
	complex*16 s,ui
	complex*16 facc,facs
	parameter (bigr=6371.)
	pi=3.14159265358979d0
	rad=180.d0/pi
	ui=dcmplx(0.d0,1.d0)

c	Read in source parameters.
	open(2,file='source-spherical.param')
	rewind(2)
	read(2,5) b80
	write(6,5) b80
	read(2,*) tm0,tm1,tm2,vmult
	write(6,*) tm0,tm1,tm2,vmult
	tm1=(tm1-tm0)/vmult
	tm2=(tm2-tm0)/vmult
	read(2,5) b80
	write(6,5) b80
	read(2,*) iseg
	write(6,*)'iseg=',iseg
	read(2,5) b80
	read(2,5) b80
	do i=1,iseg
	read(2,*) dmax(i),dmin(i),dip(i)
	write(6,*) dmax(i),dmin(i),dip(i)
	read(2,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	write(6,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	flat(i)=(pi/2.d0-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
	enddo
	read(2,5) b80
	read(2,*) idstr
	read(2,5) b80
	read(2,*) iddip
	close(2)

c	read(2,*) mtt,mpp,mzz,mtp,mtz,mpz
c	read(2,5) b80
c	read(2,*) slat,slon,zs
c	slat=pi/2.d0-slat/rad
c	slon=slon/rad	

c	Determine data vector of inverse problem
	write(6,*)'SOURCE: Determine data vector of inverse problem'
	write(6,*)'arising from the source term'

	write(6,*)'Max size of d=',3*nptmax
	write(6,*)'N=',N,'NX,NZ=',NX,NZ
	write(6,*)'Will zero out elements up to index',3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
	do i=1,3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
	d(i,1)=0.d0
	d(i,2)=0.d0
	d(i,3)=0.d0
	d(i,4)=0.d0
	enddo

	i=0
10	i=i+1
	if(i.gt.iseg) go to 30
	cdip=dcos(dip(i)/rad)
	sdip=dsin(dip(i)/rad)
	c2dip=cdip*cdip-sdip*sdip
	sstr=dsin(fstr(i))
	cstr=dcos(fstr(i))
	s2str=2.d0*sstr*cstr
	c2str=cstr*cstr-sstr*sstr 
	frak=frake(i)
	srak=dsin(frake(i))
	crak=dcos(frake(i))
c	Divide fault plane into [idstr] x [iddip] patches
	ilen=0
15	ilen=ilen+1
	if(ilen.gt.idstr) go to 10
c	Determine position (delta,phil) of running point on lower edge
c	This involves starting at (flat(i),flon(i)) and moving along
c	a direction striking fstr(i)-pi
	deltl=(fbigl(i)/dble(bigr))*(dble(ilen)-0.5d0)/dble(idstr)
        cdelt=dcos(flat(i))*dcos(deltl)+dsin(flat(i))*dsin(deltl)*(-cstr)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sphi=dsin(deltl)*(-sstr)/sdelt
        cphi=(dcos(deltl)-dcos(flat(i))*cdelt)/(dsin(flat(i))*sdelt)
	phil=datan2(sphi,cphi)
c		write(6,*)'SOURCE: ilen=',ilen,'delta,phi=',90.-delta*rad,(phil+flon(i))*rad
	delt0=delta
	idip=0
20	idip=idip+1
	if(idip.gt.iddip) go to 15
c	Determine position (slat,slon) of running point along updip direction
c	This involves starting at (delt0,phil+flon(i)) and moving along
c	a direction striking fstr(i)-pi/2
	deltd=((dmax(i)-dmin(i))*cdip/dble(bigr))*(dble(idip)-0.5d0)/dble(iddip)
	cdelt=dcos(delt0)*dcos(deltd)+dsin(delt0)*dsin(deltd)*(sstr)
        slat=dacos(cdelt)
	sdelt=sin(slat)
        sphi=dsin(deltd)*(-cstr)/sdelt
        cphi=(dcos(deltd)-dcos(delt0)*cdelt)/(dsin(delt0)*sdelt)
	slon=datan2(sphi,cphi)+phil+flon(i)
	zs=dmax(i)+(dmin(i)-dmax(i))*(dble(idip)-0.5d0)/dble(iddip)
	zs=-zs
c		write(6,*)'SOURCE: ilen,idip=',ilen,idip
c		write(6,*)'SOURCE: slat,slon,zs=',90.-slat*rad,slon*rad,zs
c		write(6,*)'------------'

c*	Determine angular distance (rad) and azimuth (rad CCW from due N)
c	of the source from the pole of the spherical coordinate system.
c	the azimuth is referenced to the azimuth of the origin (theta=dlat,z=0)
c	of the 2D grid.
        cdelt=dcos(slat)*dcos(plat)+dsin(slat)*dsin(plat)*dcos(plon-slon)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sphi=dsin(slat)*dsin(plon-slon)/sdelt
        cphi=(dcos(slat)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	phisrc=datan2(sphi,cphi)-phiref
cTE
c	phisrc=0.
c--
c		write(6,*)'SOURCE: phi_src=',rad*phisrc
c*	The theta-coord of the source is the angular distance from (plat,plon)
	xs=real(delta)
c		write(6,*)'SOURCE: xs=',xs
c	Determine which cell number the source is located in
	ncxs=0
	nczs=0
	do ncz=1,NZ
	do ncx=1,NX
c	lower left corner
	ilx=1
	ilz=1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x1=xg(ig)
	z1=zg(ig)
c		write(6,*)'SOURCE: cell #',ncx,ncz,'dimensionless lower left corner=',x1,z1
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
	  ncxs=ncx
	  nczs=ncz
	  xsource=(2./dx(ncx))*(xs-x1)-1.
	  zsource=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xsource,zsource=',xsource,zsource
	  ilc=N*((N/2)-1)+(N/2)
	  igc=igrd(ncx,ncz,ilc)
c		write(6,*)'SOURCE: gridpoint of cell center=',igc
	  endif
	enddo
	enddo
c	write(6,*)'SOURCE: Cell #s of source =(',ncxs,nczs,'). ilen,idip=',ilen,idip
	if(ncxs.eq.0.or.nczs.eq.0) then
	write(6,*)'SOURCE: Source location lies outside model domain'
	stop
	endif
c	write(6,*)'Local dimensionless source coords=(',xsource,zsource,')'
c*
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
c		write(6,*)'strainA-momten: abs(frak)=',abs(frak)
	mu0=0.1*rho(ncxs,nczs,ilc)*vs(ncxs,nczs,ilc)**2
c		write(6,*)'1st mu0=',mu0
	mupval=mupr(ncxs,nczs,ilc)
	etaval=eta(ncxs,nczs,ilc)
	et1val=eta1(ncxs,nczs,ilc)
c	Apply correspondence principle for Burgers body.
	mu2val=mupval*mu0/(mu0-mupval)
	  tau1=mu0/et1val
          tau2=mu2val/etaval
	  mu0=mu0*s*(s+tau2)/((s+tau2)*(s+tau1)+mu0*s/etaval)
c		write(6,*)'SOURCE:'
c		write(6,*)'eta(',ncxs,nczs,ilc,')=',eta(ncxs,nczs,ilc)
c		write(6,*)'mupval=',mupval
c		write(6,*)'etaval,et1val=',etaval,et1val
c		write(6,*)'mu2val=',mu2val
c		write(6,*)'tau1,tau2=',tau1,tau2
c		write(6,*)'mu0=',mu0
c	Next line is shear moment.
	shrm1=fwt(i)*fbigl(i)*((dmax(i)-dmin(i))/sdip)*mu0*1.d-6 / dble(idstr*iddip)
c		write(6,*)'fwt=',fwt(i),'mu0=',mu0,'shrm1=',shrm1
	p1=srak*sdip*cdip*s2str+crak*sdip*c2str
	p2=crak*sdip*s2str-srak*sdip*cdip*c2str
	p3=-crak*cdip*sstr+srak*c2dip*cstr
	p4=srak*c2dip*sstr+crak*cdip*cstr
	p5=srak*sdip*cdip
c		write(6,*)'p1-5=',p1,p2,p3,p4,p5
	mrr=shrm1*2.d0*p5
	mtt=-shrm1*(p2+p5)
	mpp=shrm1*(p2-p5)
	mrt=-shrm1*p4
	mrp=-shrm1*p3
	mtp=-shrm1*p1
cTE
c	mrr=0.
c	mtt=0.5d0*(mtt-mpp)
c	mpp=-mtt
	
cTE
c	mrr=0.
c	mtt=0.
c	mpp=0.
c	mrt=0.
c	mrp=0.
c	mtp=0.
cNOTE	At this point, the moment tensor, which should be read in wrt local
c	geographic r,theta,phi coordinates, should be rotated by an angle phi1
c	about the vertical axis at the source point. 
        sphi=dsin(plat)*dsin(plon-slon)/sdelt
        cphi=(dcos(plat)-dcos(slat)*cdelt)/(dsin(slat)*sdelt)
	phi1=datan2(sphi,cphi)
c		write(6,*)'SOURCE: phi1=',rad*phi1
c	Rotate moment tensor components an amount phi1 about a vertical axis
c	at the source.
	mxx=mtt*cphi*cphi + mpp*sphi*sphi - 2.d0*mtp*cphi*sphi
	myy=mpp*cphi*cphi + mtt*sphi*sphi + 2.d0*mtp*cphi*sphi
	mxy=(mtt-mpp)*cphi*sphi + mtp*(cphi*cphi-sphi*sphi)
	mxz=mrt*cphi - mrp*sphi
	myz=mrp*cphi + mrt*sphi
	mzz=mrr
cTE
c		mxx=0.
c		myy=0.
c		mxy=0.
c		mxz=2.d-2
c		myz=0.
c		mzz=0.
c--
c		write(6,*)'SOURCE: original mtt,mpp,mrr,mtp,mrt,mrp=',mtt,mpp,mrr,mtp,mrt,mrp
c		write(6,*)'SOURCE: rotated  mxx,myy,mzz,mxy,mxz,myz=',mxx,myy,mzz,mxy,mxz,myz 
c*
cTE
c
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncxs,nczs,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncxs,nczs) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	Next interpolate the basis function phi2_[il] at x=xsource,z=zsource
	phis=phiv(ilx,ilz,xsource,zsource)
c	Next interpolate the basis function d(phi2_[il])/dx at x=xsource,z=zsource
	dphisx=dphivx(ilx,ilz,xsource,zsource) * dble(2./dx(ncxs))
c	Next interpolate the basis function d(phi2_[il])/dz at x=xsource,z=zsource
	dphisz=dphivz(ilx,ilz,xsource,zsource) * dble(2./dz(nczs))
c	Eqn 10 of notes
c		write(6,*)'phis=',phis
c		write(6,*)'dphisx=',dphisx
c		write(6,*)'dphisz=',dphisz
c		write(6,*)'i1,i2,i3=',i1,i2,i3
c		write(6,*)'3*nptmax=',3*nptmax
c		write(6,*)'term1=',(-mxx*dphisx+ui*ky*mxy*phis-mxz*dphisz)
c		write(6,*)'term2=',(-mxy*dphisx+ui*ky*myy*phis-myz*dphisz)
c		write(6,*)'term3=',(-mxz*dphisx+ui*ky*myz*phis-mzz*dphisz)
	facc=dcos(ky*phisrc)
	facs=-ui*dsin(ky*phisrc)
	fac1=1.d0/(bigr+zs)
	fac2=1.d0/((bigr+zs)*sdelt)
	cott=cdelt/sdelt
c	d(i1,1)=(-mxx*(dphisx-phis*cott)*fac1+ui*ky*mxy*phis*fac2-mxz*(dphisz-phis*fac1))/s * fac * fac2 
c	d(i2,1)=(-mxy*(dphisx-phis*cott)*fac1+ui*ky*myy*phis*fac2-myz*(dphisz-phis*fac1))/s * fac * fac2
c	d(i3,1)=(-mxz*(dphisx-phis*cott)*fac1+ui*ky*myz*phis*fac2-mzz*(dphisz-phis*fac1))/s * fac * fac2
	d(i1,1)=d(i1,1) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facc * fac2
	d(i2,1)=d(i2,1) + (ui*ky*myy*phis*fac2)/s * facc * fac2
	d(i3,1)=d(i3,1) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facc * fac2
	d(i1,2)=d(i1,2) + (ui*ky*mxy*phis*fac2)/s * facc * fac2 
	d(i2,2)=d(i2,2) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facc * fac2
	d(i3,2)=d(i3,2) + (ui*ky*myz*phis*fac2)/s * facc * fac2
	d(i1,3)=d(i1,3) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facs * fac2
	d(i2,3)=d(i2,3) + (ui*ky*myy*phis*fac2)/s * facs * fac2
	d(i3,3)=d(i3,3) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facs * fac2
	d(i1,4)=d(i1,4) + (ui*ky*mxy*phis*fac2)/s * facs * fac2 
	d(i2,4)=d(i2,4) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facs * fac2
	d(i3,4)=d(i3,4) + (ui*ky*myz*phis*fac2)/s * facs * fac2
c		write(6,*)'SOURCE: d(',i1,1,')=',d(i1,1)
c		write(6,*)'SOURCE: d(',i2,1,')=',d(i2,1)
c		write(6,*)'SOURCE: d(',i3,1,')=',d(i3,1)
c		write(6,*)'SOURCE: d(',i1,2,')=',d(i1,2)
c		write(6,*)'SOURCE: d(',i2,2,')=',d(i2,2)
c		write(6,*)'SOURCE: d(',i3,2,')=',d(i3,2)
	enddo
	enddo
cTE
c		if(ilen.eq.idstr.and.idip.eq.iddip) stop
	go to 20
30	continue
	
5	format(a80)
	return
	end

	function lapl(f,t)
c	Calculate inverse Laplace transform of function f at time t.
c	The array of complex s-values where f(s) was computed
c	is already stored in common.
c	Assume that only the imaginary part of s-values vary,
c	not the real part.
c	Further assume that the imaginary part of s varies quadratically.
	real*8 lapl
	real*8 t,v
	real*8 s0
	complex*16 ui,s,f(12),ds,fn,gn,val
	real*8 sarr(7)
	real*8 bigb(8,8),c(8),a(8)
	real*8 emax,eig(8),vm(8,8)
	complex*16 bcal,ccal,acal
	complex*16 sval
	common/flapl/sval(12)
	ui=dcmplx(0.,1.)
	s0=5.6d0*dble(sval(1))
c		write(6,*)'LAPL: s0=',s0
	sarr(1)=s0*2.d0
	sarr(2)=s0
	sarr(3)=s0/2.d0
	sarr(4)=s0/10.d0
	sarr(5)=s0/50.d0
	sarr(6)=s0/100.d0
	sarr(7)=s0/500.d0
c		write(6,*)'LAPL: sarr=',sarr
		do l=1,8
	ccal=0.d0
	do k=1,8
	bcal=0.
	do n=1,12
	if(l.le.7) fn=1.d0/(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1.d0/sval(n) 
	if(k.le.7) gn=1.d0/(sval(n)*(sval(n)+sarr(k)))
	if(k.eq.8) gn=1.d0/sval(n) 
	bcal=bcal+fn*gn
	enddo
	bigb(k,l)=real(bcal)
	enddo
	do n=1,12
	if(l.le.7) fn=1./(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1./sval(n) 
	ccal=ccal+f(n)*fn
	enddo
	c(l)=real(ccal)
		enddo
c			bigb(8,8)=bigb(8,8)+1.e+6
cNOTE
c		write(6,*)'entering ainver: bigb=',bigb
c		write(6,*)'c=',c
c	call ainver(a,bigb,c,8)
c		write(6,*)'a=',a
c		write(6,*)'bigb=',bigb
	m=8
	mf1=m
c		write(6,*) 'entering svdcmp: m=',m
	call svdcmp(bigb,m,mf1,8,8,eig,vm)
c		write(6,*)'initial eigenvalues are'
c		write(6,*) (eig(j), j=1,m)
	emax=0.d0
	do i=1,m
	if(eig(i).gt.emax) emax=eig(i)
	enddo
	do i=1,m
	eig(i)=eig(i)+(2.d-6)*emax
	enddo
c		write(6,*)'revised eigenvalues are'
c		write(6,*) (eig(j), j=1,m)
c		pause
	call svdksb(bigb,eig,vm,m,mf1,8,8,c,a)
c		write(6,*)'a=',a
c		pause
	lapl=0.
	do k=1,7
	v=-sarr(k)*t
	if(v.lt.-20.d0) v=-20.d0
	lapl=lapl+a(k)*(1.d0-dexp(v))/sarr(k)
	enddo
	lapl=lapl+a(8)
c--------------
cNOTE
c	Evaluate misfit
	do n=1,12
	val=0.
	do l=1,8
	if(l.le.7) fn=1./(sval(n)*(sval(n)+sarr(l)))
	if(l.eq.8) fn=1./sval(n) 
	val=val+fn*a(l)
c		write(6,*)'LAPL: fn=',fn,'sval(',n,')=',sval(n)
	enddo
c	write(6,*)'LAPL: f(',n,')=',f(n),'calc=',val
	enddo
c-------
c		pause
	return
	end	

	subroutine addmr(nrec,phir)
c	Do the summation of transformed displacement
c	over azimuthal order number using
c	the response functions which are stored in 'rec-j-m'
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	common/grdpar/NX,NZ
	real*8 ang
	complex*16 dsrec,dsrecp
	complex*16 ds1,ds2,ds3
	complex*16 ui,facp,facm
	common/minfo/minc,mmax
	parameter (maxrec=3000)
	dimension phir(maxrec)
	common/msumr/dsrec(3*N*N*maxrec,12),dsrecp(3*N*N*maxrec,12)
	real*8 twopi
c	m0tot=1+2*(mmax/minc)
c
	twopi=2.d0*3.14159265358979d0
	ui=dcmplx(0.,1.)
	open(4,file='rec-j-m',form='unformatted')
	rewind(4)
	do j=1,12
		write(6,*)'Doing j=',j,'out of',12
c	Zero out the dsrec-array
	do ic=1,N*N*nrec
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	dsrec(i1,j)=0.d0
	dsrec(i2,j)=0.d0
	dsrec(i3,j)=0.d0
	dsrecp(i1,j)=0.d0
        dsrecp(i2,j)=0.d0
        dsrecp(i3,j)=0.d0
	enddo
	  do m=minc,mmax+minc,minc
	ic=0
	do nr=1,nrec
	phival=phir(nr)
	  ang=dble(phival)*dble(m-minc)
	  facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	  facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	do idum=1,N*N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
c		write(6,*)'m=',m
c		write(6,*)'ADDMR1: ds1=',ds1,'ds2=',ds2,'ds3=',ds3,'facp=',facp,'facs=',facs
	dsrec(i1,j)=dsrec(i1,j)+ds1*facp
	dsrec(i2,j)=dsrec(i2,j)+ds2*facp
	dsrec(i3,j)=dsrec(i3,j)+ds3*facp
	dsrecp(i1,j)=dsrecp(i1,j)+ds1*facp*ui*dble(m-minc)
        dsrecp(i2,j)=dsrecp(i2,j)+ds2*facp*ui*dble(m-minc)
        dsrecp(i3,j)=dsrecp(i3,j)+ds3*facp*ui*dble(m-minc)
	if(m.gt.minc) then
	read(4) ds1,ds2,ds3
c		write(6,*)'m=',m
c		write(6,*)'ADDMR2: ds1=',ds1,'ds2=',ds2,'ds3=',ds3,'facp=',facp,'facs=',facs
	dsrec(i1,j)=dsrec(i1,j)+ds1*facm
	dsrec(i2,j)=dsrec(i2,j)+ds2*facm
	dsrec(i3,j)=dsrec(i3,j)+ds3*facm
	dsrecp(i1,j)=dsrecp(i1,j)-ds1*facm*ui*dble(m-minc)
        dsrecp(i2,j)=dsrecp(i2,j)-ds2*facm*ui*dble(m-minc)
        dsrecp(i3,j)=dsrecp(i3,j)-ds3*facm*ui*dble(m-minc)
	endif
	enddo
	enddo
	  enddo
	enddo
	close(4)
	return
	end

	subroutine addm(phival)
c	Do the summation of transformed displacement
c	over azimuthal order number using
c	the response functions which are stored in 'displ-j-m'
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	common/grdpar/NX,NZ
	real*8 ang
	complex*16 ds
	complex*16 ds1,ds2,ds3
	complex*16 ui,facp,facm
	common/minfo/minc,mmax
	common/msum/ds(3*nptmax,12)
	real*8 twopi
c	m0tot=1+2*(mmax/minc)
c
	twopi=2.d0*3.14159265358979d0
	ui=dcmplx(0.,1.)
	open(4,file='displ-j-m',form='unformatted')
	rewind(4)
	do j=1,12
c	Zero out the ds-array
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	ds(i1,j)=0.d0
	ds(i2,j)=0.d0
	ds(i3,j)=0.d0
	enddo
	  do m=minc,mmax+minc,minc
	  ang=dble(phival)*dble(m-minc)
	  facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	  facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	ds(i1,j)=ds(i1,j)+ds1*facp
	ds(i2,j)=ds(i2,j)+ds2*facp
	ds(i3,j)=ds(i3,j)+ds3*facp
	if(m.gt.minc) then
	read(4) ds1,ds2,ds3
	ds(i1,j)=ds(i1,j)+ds1*facm
	ds(i2,j)=ds(i2,j)+ds2*facm
	ds(i3,j)=ds(i3,j)+ds3*facm
	endif
	enddo
	  enddo
	enddo
	close(4)
	return
	end

	subroutine readp
c	Read in the response functions which are stored in 'vertp-j'
	parameter (ncmax=37)
	parameter (ncmaz=37)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	common/grdpar/NX,NZ
	complex*16 ds1,ds2,ds3
	complex*16 dsvert,dsvrtp
	common/vertpr/dsvert(3*nptmax,12),dsvrtp(3*nptmax,12)
	open(4,file='vertp-j',form='unformatted')
	rewind(4)
	igmax=(N*NX-(NX-1))*(N*NZ-(NZ-1))
	do j=1,12
	do ig=1,igmax
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	dsvert(i1,j)=ds1
	dsvert(i2,j)=ds2
	dsvert(i3,j)=ds3
	read(4) ds1,ds2,ds3
        dsvrtp(i1,j)=ds1
        dsvrtp(i2,j)=ds2
        dsvrtp(i3,j)=ds3
c		write(6,*)'readp: dsvert(',i1,j,')=',dsvert(i1,j)
c		write(6,*)'readp: dsvert(',i2,j,')=',dsvert(i2,j)
	enddo
	enddo
	close(4)
	return
	end

	subroutine PlmON(p, mval, lmax, z, csphase, cnorm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the normalized associated Legendre
!	functions up to degree lmax. The functions are initially scaled by 
!	10^280 sin^m in order to minimize the effects of underflow at large m 
!	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). 
!	On a Mac OSX system with a maximum allowable double precision value of 
!	2.225073858507203E-308 the scaled portion of the algorithm will not overflow 
!	for degrees less than or equal to 2800.
!
!	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta), 
!	with the intial value of rescalem being equal to 1/scalef (which is here equal 
!	to 10^280). This will gradually reduce this huge number to a tiny number, and will 
!	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
!	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
!	number will be multipled by the old value of rescalem.
!
!	Temporary variables in saved in an allocated array. In order to explicitly deallocate this
!	memory, call this routine with a spherical harmonic degree of -1.
!
!	Calling Parameters:
!		OUT
!			p:		A vector of all associated Legendgre polynomials evaluated at 
!					z up to lmax. The length must by greater or equal to (lmax+1)*(lmax+2)/2.
!		OPTIONAL (IN)
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			cnorm:		0: Use real normalization.
!					1: Use complex normalization.
!		IN
!			lmax:		Maximum spherical harmonic degree to compute.
!			z:		cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The employed normalization is the "orthonormalized convention." The integral of
!		(plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is 1.
!	2.	The integral of plm**2 over (-1,1) is (2 - delta(0,m))/2pi. If CNORM=1, then
!		this is equal to 1/2pi.
!	3.	The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!		the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the 
!		calling routine.
!	4. 	The default is to exclude the Condon-Shortley phase of (-1)^m.
!
!
!	Dependencies:	CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek September 25, 2005.
!
!	April 19, 2008: Added CNORM optional parameter compute complex normalized functions.
!
!!	February 12, 2017 -- modified by Fred Pollitz 
!!	I added the azimuthal order number mval to the argument list because we only want results for m=mval.
!!	In particular, the Plm for l<=m+1<=lmax will be computed only for m=mval, 
!!	otherwise the P(m+1,m) lines and do l subloop which does this task will be skipped.
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

c	use SHTOOLS, only: CSPHASE_DEFAULT

	implicit none
	integer, intent(in) ::	mval, lmax
cOLD	real*8, intent(out) ::	p(:)
!!	Permit a maximum spherical harmonic degree of 2000, and this subroutine will not be called with lmax larger than this.
	real*8 p(2001*2002/2)
       	real*8, intent(in) ::	z
       	integer, intent(in), optional :: csphase, cnorm
       	real*8 ::	pm2, pm1, pmm, plm, rescalem, phase, pi, u, scalef
      	real*8, save, allocatable ::	f1(:), f2(:), sqr(:)
      	integer ::	k, kstart, m, l, astat(3)
      	integer, save ::	lmax_old = 0

	if (lmax == -1) then
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		lmax_old = 0
		return
	endif

	pi = acos(-1.0d0)
	
	if (size(p) < (lmax+1)*(lmax+2)/2) then 
		print*, "Error --- PlmBar"
     		print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlmBar"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlmBar"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif     	
     	
     	if (present(csphase)) then
     		if (csphase == -1) then
     			phase = -1.0d0
     		elseif (csphase == 1) then
     			phase = 1.0d0
     		else
     			print*, "PlmBar --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		endif
     	else
     		phase = 1.d0
     	endif
     		
	scalef = 1.0d-280
	
	
	if (lmax /= lmax_old) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		
		allocate(sqr(2*lmax+1), stat=astat(1))
		allocate(f1((lmax+1)*(lmax+2)/2), stat=astat(2))
		allocate(f2((lmax+1)*(lmax+2)/2), stat=astat(3))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
			print*, "PlmON --- Error"
			print*, "Problem allocating arrays SQR, F1 and F2", astat(1), astat(2), astat(3)
			stop
		endif
	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax+1
			sqr(l) = sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l and m=l-1,
		!	as a different recursion is used for these two values.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		k = 3
	
		do l=2, lmax, 1
			k = k + 1
			f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
			f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
			do m=1, l-2
				k = k+1
				f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                		f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) 
     &                  			 / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
			enddo
			k = k + 2
		enddo
		
		lmax_old = lmax
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	!	Calculate P(l,0). These are not scaled.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	u = sqrt((1.0d0-z)*(1.0d0+z)) ! sin(theta)

      	pm2  = 1.0d0/sqrt(4.0d0*pi)
      	p(1) = pm2
      	
      	if (lmax == 0) return
      	
      	pm1  = sqr(3)*z/sqrt(4.0d0*pi)
      	p(2) = pm1
      		
	k = 2

      	do l = 2, lmax, 1
         	k = k+l
         	plm = f1(k)*z*pm1-f2(k)*pm2
         	p(k) = plm
         	pm2  = pm1
         	pm1  = plm
      	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate P(m,m), P(m+1,m), and P(l,m)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
 	if (present(cnorm)) then
		if (cnorm == 1) then
			pmm  = scalef/sqrt(4.0d0*pi)
		else
			pmm  = sqr(2)*scalef/sqrt(4.0d0*pi)
		endif
	else
      		pmm  = sqr(2)*scalef/sqrt(4.0d0*pi)
      	endif

      	rescalem = 1.0d0/scalef
      	kstart = 1

      	do m = 1, lmax - 1, 1
      		
      		rescalem = rescalem*u

		! Calculate P(m,m)
        	kstart = kstart+m+1
         
         	pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        	p(kstart) = pmm*rescalem
        	pm2 = pmm

!!		Here go through P(m+1,m) lines and the l-loop only for m=mval, otherwise just update the k-index.
		if(m.eq.mval) then
		! Calculate P(m+1,m)
		k = kstart+m+1
	   	pm1 = z * sqr(2*m+3) * pmm
	    	p(k) = pm1*rescalem

		! Calculate P(l,m)
               	do l = m+2, lmax, 1
               		k = k+l
                  	plm  = z*f1(k)*pm1-f2(k)*pm2
                  	p(k) = plm*rescalem
                  	pm2  = pm1
                  	pm1  = plm
               	enddo
		else
		k = kstart+m+1
		k=k+(lmax*(lmax+1)-(m+1)*(m+2))/2
		endif
              
      	enddo
      	
      	! Calculate P(lmax,lmax)
      	
      	rescalem = rescalem*u
      	
        kstart = kstart+m+1
        pmm = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax)
        p(kstart) = pmm*rescalem
      		
	end subroutine PlmON

