	subroutine sporder(nmat,Ap,Ai,jAi)
c	For the matrix G to be constructed by subroutine MATREL,
c	characterize it in terms of row-counter matrices Ap and Ai
c	used in UMFPACK.  
c	Use the same loops as in subroutine MATREL to determine this
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
	common/grdpar/NX,NY,NZ
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
	parameter (kmax=ncmax*ncmay*ncmaz*N3*N3*9)
cUSED	dimension irow(kmax),icol(kmax)
	integer(4), ALLOCATABLE :: irow(:)
	integer(4), ALLOCATABLE :: icol(:)
	integer*8, ALLOCATABLE :: Ai0(: , :)
	integer, ALLOCATABLE :: X(:)
	integer, ALLOCATABLE :: Y(:)
	integer*8 nzmax,nmax
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmay-(ncmay-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
	integer*8 Ap
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	dimension Ap(nmax)
	integer*8 Ai(nmat),jAi(nmat)
	dimension nrow(nmax-1)
c
	allocate ( irow(nmat) )
	allocate ( icol(nmat) )
	k=0
c--	Following lines follow structure of subroutine MATREL
	do ncz=1,NZ
	do ncy=1,NY
	do ncx=1,NX
	do ilz=1,N
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	i1=3*igrd(ncx,ncy,ncz,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncy,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N3) go to 20
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncy,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncy,ncz)
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
	if(ily.eq.1.and.ncy.eq.1) go to 30
	if(ily.eq.N.and.ncy.eq.NY) go to 30
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
	enddo	
	enddo	
	ktot=k

c	Assume 0-based Ap and Ai as in UMFPACK conventions.
	ncol=3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
		write(6,*)'ncol=',ncol
c	do a first pass to determine how many rows per column there are 
c	(including duplicate irow's for a given icol) and the maximum #irow's for any column
	do j=1,ncol
	nrow(j)=0
	enddo
	maxr=0
	do k=1,ktot
	ic=icol(k)
	nrow(ic)=nrow(ic)+1
	if(nrow(ic).gt.maxr) maxr=nrow(ic)
	enddo
	write(6,*)'maxinum number of row entries for any column=',maxr
	write(6,*)'Allocating Ai0'
c-----------------------------
	allocate ( Ai0(ncol,maxr) )
	allocate ( X(maxr) )
	allocate ( Y(maxr) )
c-----------------------------
	write(6,*)'Construct a preliminary Ai0 array that has the row values arranged in column indices'
	write(6,*)'including redundant entries'
c	Construct a preliminary Ai0 array that has the row values arranged in column indices
c	including redundant entries
	do j=1,ncol
	nrow(j)=0
	enddo
	do k=1,ktot
	ic=icol(k)
	nrow(ic)=nrow(ic)+1
	l=nrow(ic)
	Ai0(ic,l)=irow(k)-1
	enddo
	write(6,*)'Construct the Ai array from Ai0 while also removing the redundant entries'
	write(6,*)'Also construct Ap on the fly'
c	Construct the Ai array from Ai0 while also removing the redundanct entries
c	Also construct Ap on the fly
	l1=0
	Ap(1)=0
	do ic=1,ncol
		k=ic-10000*(ic/10000)
		if(k.eq.0) write(6,*)'Doing ic=',ic,'out of',ncol
	do l=1,nrow(ic)
c	See if row ir (ir-1 in zero-based approach) is already represented in the column ic tally
	  ival=0
	  if(l.gt.1) then
	  do l0=1,l-1
	  if(Ai0(ic,l0).eq.Ai0(ic,l)) ival=1
	  enddo
	  endif
	if(ival.eq.0) then
	l1=l1+1
c	Ai(l1)=Ai0(ic,l)
	lind=l1-Ap(ic)
	X(lind)=Ai0(ic,l)
c		write(6,*)'ic=',ic,'l1=',l1,'Ap(',ic,')=',Ap(ic),'X(',lind,')=',X(lind)
	endif
	enddo
	Ap(ic+1)=l1
c	Sort the Ai-values currently in X so they are in increasing order
	  N0=Ap(ic+1)-Ap(ic)
c	  write(6,*)'ic=',ic
c	  write(6,*)'entering SORT: X=',(X(lind), lind=1,N0)
	  call SORT(X,maxr,N0,Y)
c	  write(6,*)'out of SORT: Y=',(Y(lind), lind=1,N0)
	  do l=Ap(ic)+1,Ap(ic+1)
	  lind=l-Ap(ic)
	  Ai(l)=Y(lind)
	  enddo
	enddo

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
	deallocate (irow)
	deallocate (icol)
	deallocate (Ai0)
	deallocate (X)
	deallocate (Y)

	return
	end

	subroutine globgrd(dlat)
c	Determine global grid numbers as a function of
c	(ncx,ncy,ncz), cell number indices
c	and local grid number il = N*(N*(ilz-1) + ilx) + ily,
c	where (ilx,ily,ilz) are local grid number indices
c
c	OUTPUT
c	igrd(ncx,ncy,ncz,il) array of global grid numbers
c
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
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
c	xg, yg and zg have the (theta[radians],z[km]) coordinates at each of the global gridpoints
c	assuming dimensions of (dx)x(dy)x(dz) for each cell 
c	The origin of the 3D grid is (theta=dlat,y=0,z=0)
c	(dimensionless x, y and z each run from -1 to +1).
	common/glbxy/xg(nptmax),yg(nptmax),zg(nptmax)
c
	common/grdpar/NX,NY,NZ
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
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
	do ncy=1,NY
	do ncx=1,NX
	do ilz=1,N
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	inew=1
c
c	When ilx=1 and ncx>1, then reset ilx-->N and ncx-->ncx-1
	if(ilx.eq.1.and.ncx.gt.1) then
	il1=N
	nc1=ncx-1
	il0=N*(N*(ilz-1)+il1-1)+ily
	igrd(ncx,ncy,ncz,il)=igrd(nc1,ncy,ncz,il0)
c		write(6,*)'ilx=1: igrd(',ncx,ncy,ncz,il,')=',igrd(ncx,ncy,ncz,il)
	inew=0
	endif
c
c	When ily=1 and ncy>1, then reset ily-->N and ncy-->ncy-1
	if(ily.eq.1.and.ncy.gt.1) then
	il2=N
	nc2=ncy-1
	il0=N*(N*(ilz-1)+ilx-1)+il2
	igrd(ncx,ncy,ncz,il)=igrd(ncx,nc2,ncz,il0)
c		write(6,*)'ily=1: igrd(',ncx,ncy,ncz,il,')=',igrd(ncx,ncy,ncz,il)
	inew=0
	endif
c
c	When ilz=1 and ncz>1, then reset ilz-->N and ncz-->ncz-1
	if(ilz.eq.1.and.ncz.gt.1) then
	il3=N
	nc3=ncz-1
	il0=N*(N*(il3-1)+ilx-1)+ily
	igrd(ncx,ncy,ncz,il)=igrd(ncx,ncy,nc3,il0)
c		write(6,*)'ilz=1: igrd(',ncx,ncy,ncz,il,')=',igrd(ncx,ncy,ncz,il)
	inew=0
	endif
c
	if(inew.eq.1) then
	ig=ig+1
	igrd(ncx,ncy,ncz,il)=ig
		write(6,*)'igrd(',ncx,ncy,ncz,il,')=',igrd(ncx,ncy,ncz,il)
	sumx=real(dlat)/2.
	do ix=1,ncx
	sumx=sumx+dx(ix)/2.
	enddo
	xg(ig)=2.*sumx + (dx(ncx)/2.)*(real(x(ilx))-1.)
	sumy=0.
	do iy=1,ncy
	sumy=sumy+dy(iy)/2.
	enddo
	yg(ig)=2.*sumy + (dy(ncy)/2.)*(real(x(ily))-1.)
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
	enddo
	enddo

	return
	end

	subroutine init
	character*1 gincl
	character*80 b80,sfile
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
c***
	common/grdpar/NX,NY,NZ
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
	real*8 elat,elon,plat,plon,phiref,pi,rad
	real*8 dlat,dlon,cdelt,sdelt,sphi,cphi,delta,phi
	common/sphpar/dlat,plat,plon,phiref
c	common/minfo/minc,mmax
	common/ginfo/gincl
	real*8 g0
	common/gval/g0
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
	read(2,*) NY
	read(2,5) b80
	read(2,*) (dy(i), i=1,NY)
	read(2,5) b80
	read(2,*) NZ
	read(2,5) b80
	read(2,*) (dz(i), i=1,NZ)
c	read(2,5) b80
c	read(2,*) ydist,wavmin
cc	Structure and deformation is periodic with spatial periodicity ydist km, hence
c	ydist=ydist*real(dsin(dlat))
c	wavmin=wavmin*real(dsin(dlat))
c	minc=int((2.*3.1415926535*bigr)/ydist)
c	mmax=int((2.*3.1415926535*bigr)/wavmin)
	read(2,5) b80
	read(2,5) sfile
	read(2,5) b80
	read(2,*) g0	
	g0=g0*1.d-4
	  write(6,*)'g0=',g0
	  write(6,*)'Number of cells in theta-direction=',NX
	  write(6,*)'Number of cells in phi-direction=',NY
	  write(6,*)'Number of cells in z-direction=',NZ
	  write(6,*)'Physical length of cell along theta-direction=',dx,' radians'
	  write(6,*)'Physical length of cell along phi-direction=',dy,' radians'
	  write(6,*)'Physical length of cell along z-direction=',dz,' km'
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

	subroutine matrel(s,nmat,Ap,Ai,jAi)
cccc	subroutine matrel(ome,ky,KU,KL)
c	Generate elements of G (and d) for given 
c	Laplace-transform parameter s
cccc	angular frequency ome
c	KU and KL are the number of diagonals
c	above and below the central diagonal, respectively, with non-zero 
c	matrix elements.
c		real*8 xtmp
	character*80 dum80
	complex*16 s
c	real*8 ky
cccc	real*8 ome,ky
	parameter (ncmax=13)
	parameter (ncmay=13)
	parameter (ncmaz=13)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=3)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
	parameter (N3=N*N*N)
c	Next line is the maximum number of GLL points on a side of
c	the global grid.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmay-(ncmay-1))*(N*ncmaz-(ncmaz-1)))
c	There are a maximum of nptmax GLL points in the 2D global grid,
c	each with three displacement components.
	common/glbxy/xg(nptmax),yg(nptmax),zg(nptmax)
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
	common/grdpar/NX,NY,NZ
	real*8 vpval,vsval,rhoval,mupval,etaval,et1val,mu2val,tau1,tau2
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
	real*8 g0
	common/gval/g0
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
c*****
	parameter (kmax=ncmax*ncmay*ncmaz*N3*N3*9)
	integer*8 nzmax,nmax
cOLD        parameter (nzmax = 9500000, nmax = 160000)
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmay-(ncmay-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	integer*8 Ap
	dimension Ap(nmax)
	integer*8 Ai(nmat),jAi(nmat)
cUSED	complex*16 AA (nzmax) , XX (nmax), BB (nmax), r (nmax)
	complex(8), ALLOCATABLE :: AA(:)
	double precision, ALLOCATABLE :: Ax(:)
	double precision, ALLOCATABLE :: Az(:)
	complex*16 XX (nmax), BB (nmax), r (nmax)
cUSED	double precision Ax (nzmax), xumf (nmax), b (nmax)
	double precision xumf (nmax), b (nmax)
	double precision control(20),info(90)
cUSED	double precision Az (nzmax), xz (nmax), bz (nmax)
	double precision xz (nmax), bz (nmax)
	integer*8 numf , numfz , numeric, symbolic, status, sys, filenum
c*****
	real*8 kappa0,mu0
	complex*16 lam2mu,lam,mu
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
	allocate ( AA(nmat) )
	allocate ( Ax(nmat) )
	allocate ( Az(nmat) )
	ui=dcmplx(0.d0,1.d0)

c	Zero out the elements of AA
	ncol=3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
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
	do ncy=1,NY
	do ncx=1,NX
	do ilz=1,N
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	i1=3*igrd(ncx,ncy,ncz,il)-2
	i2=i1+1
	i3=i2+1
c		write(6,*)'Working on rows',i1,i2,i3
c		write(6,*)'out of',3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
c		write(6,*)'      '
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncy,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N3) go to 20
c-----	Next three lines not needed...-----
	jlz=(jl-1)/(N*N)+1
	jlx=(jl-1-N*N*(jlz-1))/N+1
	jly=jl-N*(N*(jlz-1)+jlx-1)
c-----
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncy,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncy,ncz)
cNOTE	If the test function corresponds to
c	the bottom of the grid, then enforce zero-displacement BC
c	This assumes that the i1,i2,i3 elements of the data vector d are zero,
c	i.e. the source is not in one of the bottom cells.
c	if(ilz.eq.1.and.ncz.eq.1) then
c	write(6,*)'ilz=',ilz,'ncz=',ncz
c	write(6,*)'ncx,ncy,ncz=',ncx,ncy,ncz
c	write(6,*)'will go to 25'
c	endif
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ily.eq.1.and.ncy.eq.1) go to 31
	if(ily.eq.N.and.ncy.eq.NY) go to 31
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
c--
	l=0
15	l=l+1
	if(l.gt.N3) then
	k0=k
c		write(6,*)'Check k0: go to 10 w/k0=',k
	go to 10
	endif
c	lz=(l-1)/N+1
c	lx=l-N*(lz-1)

	lz=(l-1)/(N*N)+1
	lx=(l-1-N*N*(lz-1))/N+1
	ly=l-N*(N*(lz-1)+lx-1)

c		if(ilx.ne.lx.and.ilz.ne.lz) go to 15
		if(ilx.ne.lx.and.ily.ne.ly.and.ilz.ne.lz) go to 15
c**
c	write(6,*)'ilx,lx=(',ilx,lx,')','ilz,lz=(',ilz,lz,')'
c	write(6,*)'dphi2x(',il,l,')=',dphi2x(il,l)
c	write(6,*)'dphi2z(',il,l,')=',dphi2z(il,l)
c	write(6,*)'-----------------------------'

	ig=igrd(ncx,ncy,ncz,l)
	rhoval=rho(ncx,ncy,ncz,l)
	vsval=vs(ncx,ncy,ncz,l)
	vpval=vp(ncx,ncy,ncz,l)
	mupval=mupr(ncx,ncy,ncz,l)
	etaval=eta(ncx,ncy,ncz,l)
	et1val=eta1(ncx,ncy,ncz,l)
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
	fac3=fac2 * dble(1.0/sin(xg(ig)))
c		write(6,*)'fac1,fac2=',fac1,fac2,'xg(',ig,')=',xg(ig)
	cott=dble(cos(xg(ig))/sin(xg(ig)))
c		write(6,*)'cott=',cott
c	Terms arising from the theta- and r-parts of the gradient of the stress tensor
c	after integration by parts in theta and r.
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam2mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1 * dble(dy(ncy)/2.)
     &			    - wt2(l)*dphi2x(il,l)*lam*phi2(jl,l)*cott * dble(dz(ncz)/2.) * fac1 * dble(dy(ncy)/2.)
     &		            - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac * dble(dy(ncy)/2.)
     &		            + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dble(dx(ncx)/2.0) * dble(dy(ncy)/2.)

c		write(6,*)'line 1: AA(',j,')=',AA(j)
c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam*(dphi2y(jl,l)) * dble(dz(ncz)/2.) * fac2 
c		write(6,*)'line 2: AA(',j,')=',AA(j)

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*lam2mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac1 * dble(dy(ncy)/2.)
     &			    - wt2(l)*dphi2x(il,l)*lam*dphi2z(jl,l) * 1.d0 * dble(dy(ncy)/2.)
     &			    - wt2(l)*dphi2x(il,l)*lam*phi2(jl,l) * dble(dz(ncz)/2.) * fac1 * dble(dy(ncy)/2.)
     &	                    - wt2(l)*dphi2z(il,l)*mu*dphi2x(jl,l) * 1.d0 * dble(dy(ncy)/2.)
	AA(j)=AA(j)
     &		    - wt2(l)*phi2(il,l)*rhoval*g0*dphi2x(jl,l) * dble(dz(ncz)/2.) * dble(dy(ncy)/2.)
c		write(6,*)'line 3: AA(',j,')=',AA(j)
c	Above line: gravitational term

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*(dphi2y(jl,l)) * dble(dz(ncz)/2.) * fac2

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1 * dble(dy(ncy)/2.)
     &			    + wt2(l)*dphi2x(il,l)*mu*phi2(jl,l)*cott * dble(dz(ncz)/2.) * fac1 * dble(dy(ncy)/2.)
     &                      - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac * dble(dy(ncy)/2.)
     &			    + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dble(dx(ncx)/2.) * dble(dy(ncy)/2.)

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2z(il,l)*mu*(dphi2y(jl,l)) * dble(dx(ncx)/2.) * fac2 * fac
	AA(j)=AA(j)
     &		    - wt2(l)*phi2(il,l)*rhoval*g0*(dphi2y(jl,l)) * dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
c	Above line: gravitational term

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2z(jl,l) * 1.d0 * dble(dy(ncy)/2.)
     &			    + wt2(l)*dphi2x(il,l)*mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac1 * dble(dy(ncy)/2.)
     &                      - wt2(l)*dphi2z(il,l)*lam*dphi2x(jl,l) * 1.d0 * dble(dy(ncy)/2.)
     &			    - wt2(l)*dphi2z(il,l)*lam*phi2(jl,l)*cott * dble(dx(ncx)/2.) * dble(dy(ncy)/2.)
	AA(j)=AA(j)
     &		    + wt2(l)*phi2(il,l)*rhoval*g0*dphi2x(jl,l) * dble(dz(ncz)/2.) * dble(dy(ncy)/2.)
     &		    + wt2(l)*phi2(il,l)*rhoval*g0*phi2(jl,l)*cott * dble(dz(ncz)*dx(ncx)/4.) * dble(dy(ncy)/2.)
c	Above line: gravitational term

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2z(il,l)*lam*(dphi2y(jl,l)) * dble(dx(ncx)/2.) * fac2 * fac
	AA(j)=AA(j)
     &		    + wt2(l)*phi2(il,l)*rhoval*g0*(dphi2y(jl,l)) * dble(dz(ncz)*dx(ncx)/4.) * fac2 * fac
c	Above line: gravitational term

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2x(il,l)*mu*dphi2x(jl,l) * dble(dz(ncz)/dx(ncx)) * fac1 * dble(dy(ncy)/2.)
     &                        - wt2(l)*dphi2z(il,l)*lam2mu*dphi2z(jl,l) * dble(dx(ncx)/dz(ncz)) * fac * dble(dy(ncy)/2.)
     &                        - wt2(l)*dphi2z(il,l)*2.d0*lam*phi2(jl,l) * dble(dx(ncx)/2.) * dble(dy(ncy)/2.)
	AA(j)=AA(j)
     &		    + wt2(l)*phi2(il,l)*rhoval*g0*2.d0*phi2(jl,l) * dble(dz(ncz)*dx(ncx)/4.) * dble(dy(ncy)/2.)
c	Above line: gravitational term

c**********************

c	Terms arising from the phi- parts of the gradient of the stress tensor
c	after integration by parts in phi.
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*mu*dphi2y(jl,l) * dble(dz(ncz)/dy(ncy)) * fac3 * dble(dx(ncx)/2.)

c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*mu*(dphi2x(jl,l)) * dble(dz(ncz)/2.) * fac2
     &	+ wt2(l)*dphi2y(il,l)*mu*phi2(jl,l)*cott * dble(dx(ncx)/2.) * dble(dz(ncz)/2.) * fac2

c	i1,j3 terms
	k=k+1
c	No contribution to theta (i1)-component from c_jl terms

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*lam2mu*phi2(jl,l)*cott *fac2 * dble(dx(ncx)/2.) * dble(dz(ncz)/2.)
     &	- wt2(l)*dphi2y(il,l)*lam*(dphi2x(jl,l)) * dble(dz(ncz)/2.) * fac2 

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*lam2mu*dphi2y(jl,l) * dble(dz(ncz)/dy(ncy)) * fac3 * dble(dx(ncx)/2.)

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*lam2mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac2 * dble(dx(ncx)/2.)
     &			    - wt2(l)*dphi2y(il,l)*lam*dphi2z(jl,l) * (fac2/fac1) * dble(dx(ncx)/2.)
     &			    - wt2(l)*dphi2y(il,l)*lam*phi2(jl,l) * dble(dz(ncz)/2.) * fac2 * dble(dx(ncx)/2.)

c	i3,j1 terms
	k=k+1
c	No contribution to z (i3)-component from a_jl terms

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*mu*dphi2z(jl,l) * (fac2/fac1) * dble(dx(ncx)/2.)
     &			    + wt2(l)*dphi2y(il,l)*mu*phi2(jl,l) * dble(dz(ncz)/2.) * fac2 * dble(dx(ncx)/2.)

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA(j)=AA(j)
     &	- wt2(l)*dphi2y(il,l)*mu*dphi2y(jl,l) * dble(dz(ncz)/dy(ncy)) * fac3 * dble(dx(ncx)/2.)

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
30      continue
        igb=igrd(ncx,ncy,ncz,il)
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
31      continue
        igb=igrd(ncx,ncy,ncz,il)
        roff=1.d0/(yg(igb)-yg(igc))
        k=k0+1
        j=jAi(k)
        AA(j)=phi2(jl,il)*roff+dphi2y(jl,il)*(2.d0/dy(ncy))
        k=k+1
        j=jAi(k)
        AA(j)=phi2(jl,il)*roff+dphi2y(jl,il)*(2.d0/dy(ncy))
        k=k+1
        j=jAi(k)
        AA(j)=phi2(jl,il)*roff+dphi2y(jl,il)*(2.d0/dy(ncy))
        k0=k
        go to 10
20      continue

	enddo
	enddo
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
	deallocate (AA)

	numf=3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))

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
        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

cNOTE	Do only the first RHS in the 3D version

c		stop

c	open(2,file='matrel.solution1')
c	write(6,*)'solution vector 1='
c	do i=1,numf
c	write(6,*) i,d(i,1)
c	write(2,*) i,real(d(i,1))
c	enddo
c	close(2)

c	deallocate (AA)
	deallocate (Ax)
	deallocate (Az)

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
c	xg, yg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),yg(nptmax),zg(nptmax)
	parameter (iptmax=1000000)
	real*4 muprin
	real(4), ALLOCATABLE :: vp0(:)
	real(4), ALLOCATABLE :: vs0(:)
	real(4), ALLOCATABLE :: rho0(:)
	real(4), ALLOCATABLE :: mupr0(:)
	real(4), ALLOCATABLE :: eta0(:)
	real(4), ALLOCATABLE :: eta01(:)
	real(4), ALLOCATABLE :: x0(:)
	real(4), ALLOCATABLE :: y0(:)
	real(4), ALLOCATABLE :: z0(:)
c	real*8 vp0(iptmax),vs0(iptmax),rho0(iptmax),
c     &  mupr0(iptmax),eta0(iptmax),eta01(iptmax)
c	dimension x0(iptmax),y0(iptmax),z0(iptmax)
	character*80 filin
	character*120 aread
	common/grdpar/NX,NY,NZ
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
c
c	allocate ( vp0(200000) )
c	allocate ( vs0(200000) )
c	allocate ( rho0(200000) )
c	allocate ( mupr0(200000) )
c	allocate ( eta0(200000) )
c	allocate ( eta01(200000) )
c	allocate ( x0(200000) )
c	allocate ( y0(200000) )
c	allocate ( z0(200000) )
	allocate ( vp0(70*NX*NY*NZ*N*N*N) )
	allocate ( vs0(70*NX*NY*NZ*N*N*N) )
	allocate ( rho0(70*NX*NY*NZ*N*N*N) )
	allocate ( mupr0(70*NX*NY*NZ*N*N*N) )
	allocate ( eta0(70*NX*NY*NZ*N*N*N) )
	allocate ( eta01(70*NX*NY*NZ*N*N*N) )
	allocate ( x0(70*NX*NY*NZ*N*N*N) )
	allocate ( y0(70*NX*NY*NZ*N*N*N) )
	allocate ( z0(70*NX*NY*NZ*N*N*N) )
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
	read(aread,*,end=17) xin,yin,zin,vpin,vsin,rhoin,muprin,etain,eta1in
c		write(6,*)'Used first read'
	i=i+1
	eta01(i)=eta1in
	mupr0(i)=muprin
	go to 12
17	i=i+1
	eta01(i)=1.d+13
	mupr0(i)=0.d0
	read(aread,*) xin,yin,zin,vpin,vsin,rhoin,etain
c		write(6,*)'Used second read'
12	continue
cOLD	i=i+1
c
c	if(i.gt.iptmax) then
c	write(6,*)'Current number of parameter lines =',i,'for structural parameter input'
c	write(6,*)'exceeds limit of iptmax=',iptmax
c	stop
c	endif

c	if(i.eq.(200000)) then
	if(i.eq.(70*NX*NY*NZ*N*N*N)) then
	write(6,*)'Current number of parameter lines =',i,'for structural parameter input'
	write(6,*)'exceeds limit of ',20*NX*NY*NZ*N*N*N
	write(6,*)'Increase the allocated space of arrays vp0, etc. in visco3dsubs.f'
	stop
	endif
c
	vp0(i)=dble(vpin)
	vs0(i)=dble(vsin)
	rho0(i)=dble(rhoin)
	eta0(i)=etain
c	Convert (xin[deg.],yin[deg.],zin[km]) to (xin[km],yin[km],zin[km])
c	and reference to the initial x0-value
	x0(i)=(xin/rad)*6371.
	x0(i)=x0(i)-x0(1)
	y0(i)=(yin/rad)*6371.
	y0(i)=y0(i)-y0(1)
	z0(i)=zin
c		write(6,*)'x0 and z0(',i,')=',x0(i),z0(i)
	go to 5
10	continue
	close(2)
14	format(a120)
	imax=i
	write(6,*)'Finished reading structural parameters at',imax,'points'
c	Fill up model grid with values based on closest distance to just read-in values
c	igmax=(N*NX-(NX-1))*(N*NZ-(NZ-1))
c	do ig=1,igmax

	  do ncz=1,NZ
		write(6,*)'Assigning structural parameters at depth level',ncz,' out of',NZ
	  do ncy=1,NY
		write(6,*)'ncy=',ncy,' out of',NY
	  do ncx=1,NX
	  do ilz=1,N
	  do ily=1,N
	  do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	ig=igrd(ncx,ncy,ncz,il)

	xval=(xg(ig)-xg(1))*6371.
	if(ilx.eq.1) xval=xval+0.01*dx(ncx)*6371.
	if(ilx.eq.N) xval=xval-0.01*dx(ncx)*6371.
	yval=(yg(ig)-yg(1))*6371.
	if(ily.eq.1) yval=yval+0.01*dy(ncy)*6371.
	if(ily.eq.N) yval=yval-0.01*dy(ncy)*6371.
	zval=zg(ig)
	if(ilz.eq.1) zval=zval+0.01*dz(ncz)
	if(ilz.eq.N) zval=zval-0.01*dz(ncz)
c	find closest read-in gridpoint to (xg(ig),yg(ig),zg(ig))
	do i=1,imax
	dist=(xval-x0(i))**2 + (yval-y0(i))**2 + (zval-z0(i))**2
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
c		if(ig.eq.94) then
c		write(6,*)'latest dist=',dist
c		write(6,*)'xg,yg,zg=',xg(ig),yg(ig),zg(ig)
c		write(6,*)'xval,yval,zval=',xval,yval,zval
c		write(6,*)'x0,z0=',x0(i),z0(i)
c		write (6,*)' ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^'
c		endif
	endif
	enddo
		if(ig.eq.94) then
		write(6,*)'xg,yg,zg=',xg(ig),yg(ig),zg(ig)
		write(6,*)'xval,yval,zval=',xval,yval,zval
		write(6,*)'closest i=',im,' with distance=',sqrt(distm),' km'
		write(6,*)'NX=',NX,'NY=',NY,'NZ=',NZ
		write(6,*)'ig=',ig,'im=',im,'nptmax=',nptmax
		write(6,*)'- - - - - - - - - - - - -'
		endif
	vp(ncx,ncy,ncz,il)=vp0(im)
	vs(ncx,ncy,ncz,il)=vs0(im)
	rho(ncx,ncy,ncz,il)=rho0(im)
	mupr(ncx,ncy,ncz,il)=mupr0(im)
	eta(ncx,ncy,ncz,il)=eta0(im)
c		write(6,*)'GRIDGEOM: eta(',ncx,ncy,ncz,il,')=',eta(ncx,ncy,ncz,il)
	eta1(ncx,ncy,ncz,il)=eta01(im)
c		write(6,*)'GRIDGEOM: x,y,z=',xg(ig),'radians',yg(ig),'radians',zg(ig),'km'
c		write(6,*)'vp and vs(',ncx,ncy,ncz,il,')=',vp(ncx,ncy,ncz,il),vs(ncx,ncy,ncz,il)
c		write(6,*)'eta and eta1(',ncx,ncy,ncz,il,')=',eta(ncx,ncy,ncz,il),eta1(ncx,ncy,ncz,il)
c		write(6,*)'------------------------'
	  enddo
	  enddo
	  enddo
	  enddo
	  enddo
	  enddo

	deallocate ( vp0 )
	deallocate ( vs0 )
	deallocate ( rho0 )
	deallocate ( mupr0 )
	deallocate ( eta0 )
	deallocate ( eta01 )
	deallocate ( x0 )
	deallocate ( y0 )
	deallocate ( z0 )
c	enddo

	return
	end

	subroutine source(s)
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
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),yg(nptmax),zg(nptmax)
c	igc stored in common will have the grid # of the last point
c	determined on the source plane -- representative of the source location.
	common/isrc/igc
c
	common/grdpar/NX,NY,NZ
	real*8 mupval,etaval,et1val,mu2val,tau1,tau2
	real*8 vp,vs,rho,eta,mupr,eta1
	common/struc/dx(ncmax),dy(ncmay),dz(ncmaz),vp(ncmax,ncmay,ncmaz,N*N*N),
     &	vs(ncmax,ncmay,ncmaz,N*N*N),rho(ncmax,ncmay,ncmaz,N*N*N),
     &	mupr(ncmax,ncmay,ncmaz,N*N*N),eta(ncmax,ncmay,ncmaz,N*N*N),eta1(ncmax,ncmay,ncmaz,N*N*N)
	common/glbgrd/igrd(ncmax,ncmay,ncmax,N*N*N)
	real*8 cdip,sdip,c2dip
	real*8 sstr,cstr,s2str,c2str,frak,srak,crak
	real*8 deltl,phil,deltd,delt0
	real*8 mu0
	complex*16 muval
	complex*16 shrm1
	real*8 p1,p2,p3,p4,p5
	real*8 dmax(600),dmin(600),dip(600)
	real*8 flat(600),flon(600),fbigl(600),fstr(600),frake(600),fwt(600)
	real*8 tm1,tm2,vmult
	real*8 plat,plon,phiref,pi,rad
	real*8 dlat,cdelt,sdelt,delta,sphi,cphi,phisrc,phi1
	common/sphpar/dlat,plat,plon,phiref
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	complex*16 mxx,myy,mzz,mxy,mxz,myz
	complex*16 mtt,mpp,mtp,mrt,mrp,mrr
	complex*16 d
	common/datavec/d(3*nptmax,4)
	common/tmvals/tm1,tm2
c***
	real*8 slat,slon
	real*8 phis,dphisx,dphisy,dphisz
	real*8 phiv,dphivx,dphivy,dphivz
	real*8 fac1,fac2,cott
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
	write(6,*)'N=',N,'NX,NY,NZ=',NX,NY,NZ
	write(6,*)'Will zero out elements up to index',3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
	do i=1,3*(N*NX-(NX-1))*(N*NY-(NY-1))*(N*NZ-(NZ-1))
	d(i,1)=0.d0
c	d(i,2)=0.d0
c	d(i,3)=0.d0
c	d(i,4)=0.d0
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
	deltd=((dmax(i)-dmin(i))*(cdip/sdip)/dble(bigr))*(dble(idip)-0.5d0)/dble(iddip)
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
c	the azimuth is referenced to the azimuth of the origin (theta=dlat,phi=0,z=0)
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
cNOTE --check that this is giving phi>=0 starting from the edge of the domain in the y-direction--
c*	The phi-coord of the source is the azimuth from (plat,plon)
	ys=phisrc
c		write(6,*)'SOURCE: xs=',xs,'ys=',ys
c	Determine which cell number the source is located in
	ncxs=0
	ncys=0
	nczs=0
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
	  ncxs=ncx
	  ncys=ncy
	  nczs=ncz
	  xsource=(2./dx(ncx))*(xs-x1)-1.
	  ysource=(2./dy(ncy))*(ys-y1)-1.
	  zsource=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xsource,ysource,zsource=',xsource,ysource,zsource
	  ilc=N*(N*((N/2)-1)+(N/2)-1)+(N/2)
	  igc=igrd(ncx,ncy,ncz,ilc)
c		write(6,*)'SOURCE: gridpoint of cell center=',igc
	  endif
	enddo
	enddo
	enddo
c	write(6,*)'SOURCE: Cell #s of source =(',ncxs,ncys,nczs,'). ilen,idip=',ilen,idip
	if(ncxs.eq.0.or.ncys.eq.0.or.nczs.eq.0) then
	write(6,*)'SOURCE: Source location lies outside model domain'
	stop
	endif
c	write(6,*)'Local dimensionless source coords=(',xsource,ysource,zsource,')'
c*
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
c		write(6,*)'strainA-momten: abs(frak)=',abs(frak)
c*-*-*-*-*-*-*-*
c	Use Lagrangian interpolation to determine shear modulus at the source point
	muval=0.
	do ilz=1,N
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	mu0=0.1*rho(ncxs,ncys,nczs,il)*vs(ncxs,ncys,nczs,il)**2
c		write(6,*)'1st mu0=',mu0
	mupval=mupr(ncxs,ncys,nczs,il)
	etaval=eta(ncxs,ncys,nczs,il)
	et1val=eta1(ncxs,ncys,nczs,il)
c	Apply correspondence principle for Burgers body.
	mu2val=mupval*mu0/(mu0-mupval)
	  tau1=mu0/et1val
          tau2=mu2val/etaval
	  muval=muval + mu0*s*(s+tau2)/((s+tau2)*(s+tau1)+mu0*s/etaval)
     &	   * phiv(ilx,ily,ilz,xsource,ysource,zsource)
	enddo
	enddo
	enddo
c*-*-*-*-*-*-*-*
c		write(6,*)'SOURCE:'
c		write(6,*)'eta(',ncxs,nczs,ilc,')=',eta(ncxs,nczs,ilc)
c		write(6,*)'mupval=',mupval
c		write(6,*)'etaval,et1val=',etaval,et1val
c		write(6,*)'mu2val=',mu2val
c		write(6,*)'tau1,tau2=',tau1,tau2
c		write(6,*)'muval=',muval
c	Next line is shear moment.
	shrm1=fwt(i)*fbigl(i)*((dmax(i)-dmin(i))/sdip)*muval*1.d-6 / dble(idstr*iddip)
c		write(6,*)'fwt=',fwt(i),'muval=',muval,'shrm1=',shrm1
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
	do ily=1,N
	do ilx=1,N
	il=N*(N*(ilz-1)+ilx-1)+ily
	i1=3*igrd(ncxs,ncys,nczs,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncxs,ncys,nczs) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	Next interpolate the basis function d(phi2_[il])/dy at x=xsource,y=ysource,z=zsource
	dphisy=dphivy(ilx,ily,ilz,xsource,ysource,zsource) * dble(2./dy(ncys))
c	Next interpolate the basis function d(phi2_[il])/dx at x=xsource,y=ysource,z=zsource
	dphisx=dphivx(ilx,ily,ilz,xsource,ysource,zsource) * dble(2./dx(ncxs))
c	Next interpolate the basis function d(phi2_[il])/dz at x=xsource,y=ysource,z=zsource
	dphisz=dphivz(ilx,ily,ilz,xsource,ysource,zsource) * dble(2./dz(nczs))
c	Eqn 10 of notes
c		write(6,*)'phis=',phis
c		write(6,*)'dphisx=',dphisx
c		write(6,*)'dphisz=',dphisz
c		write(6,*)'i1,i2,i3=',i1,i2,i3
c		write(6,*)'3*nptmax=',3*nptmax
c		write(6,*)'term1=',(-mxx*dphisx+ui*ky*mxy*phis-mxz*dphisz)
c		write(6,*)'term2=',(-mxy*dphisx+ui*ky*myy*phis-myz*dphisz)
c		write(6,*)'term3=',(-mxz*dphisx+ui*ky*myz*phis-mzz*dphisz)
cOLD	facc=dcos(ky*phisrc)
cOLD	facs=-ui*dsin(ky*phisrc)
	fac1=1.d0/(bigr+zs)
	fac2=1.d0/((bigr+zs)*sdelt)
	cott=cdelt/sdelt
c	d(i1,1)=(-mxx*(dphisx-phis*cott)*fac1+ui*ky*mxy*phis*fac2-mxz*(dphisz-phis*fac1))/s * fac * fac2 
c	d(i2,1)=(-mxy*(dphisx-phis*cott)*fac1+ui*ky*myy*phis*fac2-myz*(dphisz-phis*fac1))/s * fac * fac2
c	d(i3,1)=(-mxz*(dphisx-phis*cott)*fac1+ui*ky*myz*phis*fac2-mzz*(dphisz-phis*fac1))/s * fac * fac2
	d(i1,1)=d(i1,1) + (-mxx*(dphisx)*fac1-mxy*(dphisy)*fac2-mxz*(dphisz))/s * fac1
	d(i2,1)=d(i2,1) + (-mxy*(dphisx)*fac1-myy*(dphisy)*fac2-myz*(dphisz))/s * fac1
	d(i3,1)=d(i3,1) + (-mxz*(dphisx)*fac1-myz*(dphisy)*fac2-mzz*(dphisz))/s * fac1
c	d(i1,1)=d(i1,1) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facc * fac2
c	d(i2,1)=d(i2,1) + (ui*ky*myy*phis*fac2)/s * facc * fac2
c	d(i3,1)=d(i3,1) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facc * fac2
c	d(i2,1)=d(i2,1) + (ui*ky*myy*phis*fac2)/s * facc * fac2
c	d(i3,1)=d(i3,1) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facc * fac2
c	d(i1,2)=d(i1,2) + (ui*ky*mxy*phis*fac2)/s * facc * fac2 
c	d(i2,2)=d(i2,2) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facc * fac2
c	d(i3,2)=d(i3,2) + (ui*ky*myz*phis*fac2)/s * facc * fac2
c	d(i1,3)=d(i1,3) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facs * fac2
c	d(i2,3)=d(i2,3) + (ui*ky*myy*phis*fac2)/s * facs * fac2
c	d(i3,3)=d(i3,3) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facs * fac2
c	d(i1,4)=d(i1,4) + (ui*ky*mxy*phis*fac2)/s * facs * fac2 
c	d(i2,4)=d(i2,4) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facs * fac2
c	d(i3,4)=d(i3,4) + (ui*ky*myz*phis*fac2)/s * facs * fac2
c		write(6,*)'SOURCE: d(',i1,1,')=',d(i1,1)
c		write(6,*)'SOURCE: d(',i2,1,')=',d(i2,1)
c		write(6,*)'SOURCE: d(',i3,1,')=',d(i3,1)
c		write(6,*)'SOURCE: d(',i1,2,')=',d(i1,2)
c		write(6,*)'SOURCE: d(',i2,2,')=',d(i2,2)
c		write(6,*)'SOURCE: d(',i3,2,')=',d(i3,2)
	enddo
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
c	Determine inverse Laplace transform using vm and eig
c	previously determined by function lapl. 
	real*8 t,v,lapl
	complex*16 ui
	complex*16 f(28)
	real*8 sarr(11)
	complex*16 ccal,fn
	real*8 c(12),a(12)
c - - - - - - - - -
	complex*16 sval
	common/flapl/sval(28)
	real*8 s0
	common/s0val/s0
	real*8 eig,vm,bigb
	common/flapl2/eig(12),vm(12,12),bigb(12,12)
c - - - - - - - - -
	ui=dcmplx(0.,1.)
	sarr(1)=s0*2.d0
	sarr(2)=s0
	sarr(3)=s0/2.d0
	sarr(4)=s0/5.d0
	sarr(5)=s0/10.d0
	sarr(6)=s0/30.d0
	sarr(7)=s0/100.d0
	sarr(8)=s0/300.d0
	sarr(9)=s0/1000.d0
	sarr(10)=s0/3000.d0
	sarr(11)=s0/10000.d0


c		write(6,*)'f=',f
		do l=1,12
	ccal=0.d0
	do n=1,28
	if(l.le.11) fn=1./(sval(n)*(sval(n)+sarr(l))) 
	if(l.eq.12) fn=1./(2.d0*s0*sval(n)) 
	ccal=ccal+f(n)*fn
	enddo
	c(l)=real(ccal)
		enddo
c
	m=12
	mf1=m
c		write(6,*)'entering svdksb: c=',c
c		write(6,*)'bigb=',bigb
c		write(6,*)'eig=',eig
	call svdksb(bigb,eig,vm,m,mf1,12,12,c,a)
c		write(6,*)'out of svdksb: a=',a
c		write(6,*)'---------------'
c
	lapl=0.
	do k=1,11
	v=-sarr(k)*t
	if(v.lt.-20.d0) v=-20.d0
	lapl=lapl+a(k)*(1.d0-dexp(v))/sarr(k)
	enddo
	lapl=lapl+a(12)/(2.d0*s0)
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

      SUBROUTINE SORT(X,Mx,N,Y)
C
C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              SINGLE PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED. 
C                     --Mx     = THE SIZE OF X AND Y
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE SINGLE PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE SINGLE PRECISION VECTOR X.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU 
C                   (DEFINED AND USED INTERNALLY WITHIN
C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
C                   IF IL AND IU EACH HAVE DIMENSION K,
C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
C                   OF IL AND IU HAVE BEEN SET TO 36,
C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
C                   APPROXIMATELY 137 BILLION.
C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
C                   THEN THERE IS NO PRACTICAL RESTRICTION
C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE 
C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
C                   INTO THIS SUBROUTINE.)
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
C              WILL BE PLACED IN THE FIRST POSITION
C              OF THE VECTOR Y,
C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
C              WILL BE PLACED IN THE SECOND POSITION
C              OF THE VECTOR Y, ETC.
C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
C              THIS IS DONE BY HAVING THE SAME
C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE. 
C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
C              CALL SORT(X,N,X)
C              IS ALLOWABLE AND WILL RESULT IN
C              THE DESIRED 'IN-PLACE' SORT.
C     COMMENT--THE SORTING ALGORTHM USED HEREIN
C              IS THE BINARY SORT.
C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
C              FOLLOWING TIME TRIALS INDICATE.
C              THESE TIME TRIALS WERE CARRIED OUT ON THE
C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
C              IN AUGUST OF 1974.
C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
C              ALSO BEEN INCLUDED--
C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
C               NUMBERS SORTED
C                N = 10                 .002 SEC          .002 SEC
C                N = 100                .011 SEC          .045 SEC
C                N = 1000               .141 SEC         4.332 SEC
C                N = 3000               .476 SEC        37.683 SEC
C                N = 10000             1.887 SEC      NOT COMPUTED
C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
C                 BY RICHARD C. SINGLETON).
C               --CACM JANUARY 1970, PAGE 54.
C               --CACM OCTOBER 1970, PAGE 624.
C               --JACM JANUARY 1961, PAGE 41.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE--301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --NOVEMBER  1975. 
C
C---------------------------------------------------------------------
C
      INTEGER X(Mx),Y(Mx)
	INTEGER HOLD,AMED
      DIMENSION IU(36),IL(36) 
C
      IPR=6
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
      IF(N.LT.1)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
      IF(X(I).NE.HOLD)GOTO90
60 	CONTINUE
cc      WRITE(IPR, 9)HOLD
      DO61I=1,N
      Y(I)=X(I)
61 	CONTINUE
      RETURN
50 	CONTINUE
cc	WRITE(IPR,15) 
cc      WRITE(IPR,47)N
      RETURN
55 	CONTINUE
cc	WRITE(IPR,18) 
      Y(1)=X(1)
      RETURN
90 	CONTINUE
cc9	FORMAT(1H ,108H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUME)
cc     1NT (A VECTOR) TO THE SORT   SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6
cc     1H *****)
cc15 	FORMAT(1H , 91H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE
cc     1 SORT   SUBROUTINE IS NON-POSITIVE *****)
cc18 	FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUME
cc     1NT TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****)
cc47 	FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
C     COPY THE VECTOR X INTO THE VECTOR Y
      DO100I=1,N
      Y(I)=X(I)
100 	CONTINUE
C
C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
C
      NM1=N-1
      DO200I=1,NM1
      IP1=I+1
      IF(Y(I).LE.Y(IP1))GOTO200
      GOTO250
200 	CONTINUE
      RETURN
250 	M=1 
      I=1 
      J=N 
305 	IF(I.GE.J)GOTO370
310 	K=I 
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO320 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
320 	L=J 
      IF(Y(J).GE.AMED)GOTO340 
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO340 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO340
330 	Y(L)=Y(K)
      Y(K)=TT
340 	L=L-1
      IF(Y(L).GT.AMED)GOTO340 
      TT=Y(L)
350 	K=K+1
      IF(Y(K).LT.AMED)GOTO350 
      IF(K.LE.L)GOTO330
      LMI=L-I
      JMK=J-K
      IF(LMI.LE.JMK)GOTO360
      IL(M)=I
      IU(M)=L
      I=K 
      M=M+1
      GOTO380
360 	IL(M)=K
      IU(M)=J
      J=L 
      M=M+1
      GOTO380
370 	M=M-1
      IF(M.EQ.0)RETURN
      I=IL(M)
      J=IU(M)
380 	JMI=J-I
      IF(JMI.GE.11)GOTO310
      IF(I.EQ.1)GOTO305
      I=I-1
390 	I=I+1
      IF(I.EQ.J)GOTO370
      AMED=Y(I+1)
      IF(Y(I).LE.AMED)GOTO390 
      K=I 
395 	Y(K+1)=Y(K)
      K=K-1
      IF(AMED.LT.Y(K))GOTO395 
      Y(K+1)=AMED
      GOTO390
      END 
