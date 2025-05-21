	subroutine basis
c	Determine basis functions (e.g. eqn 59 of Vosse96b) and their
c	derivatives at the GLL points.  Assume n=5 of wikipedia's
c	convention (N=4 Vosse's convention).  That is, GLL points are at
c	three interior points plus +-1
	parameter (lmax=3)
	real*8 x0(lmax+1),dx0(lmax+1)
	real*8 legi,y
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
	call gl
	do i=1,lmax+2
	y=x(i)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	do j=1,lmax+2
	if(i.ne.j) phi(i,j)=0.d0
	if(i.eq.j) phi(i,j)=1.d0
	if(i.ne.j) then
	y=x(j)
	call legend(y,x0,dx0)
	dphi(i,j)=x0(lmax+1)/((x(j)-x(i))*legi)
c		write(6,*)'x(j)-x(i)=',x(j)-x(i),'legi=',legi
	endif
	if(i.eq.j) then
	y=x(j)
	call legend(y,x0,dx0)
	dphi(i,j)=0.5d0*dx0(lmax+1)/x0(lmax+1)
	endif
	enddo
	enddo
c	3D basis functions
	do i1=1,lmax+2
	do i2=1,lmax+2
	do i3=1,lmax+2
	k=(lmax+2)*((lmax+2)*(i3-1)+i1-1)+i2
	do j1=1,lmax+2
	do j2=1,lmax+2
	do j3=1,lmax+2
	l=(lmax+2)*((lmax+2)*(j3-1)+j1-1)+j2
c	phi2(k,l)=phi_i1 (x_j1) * phi_i2 (x_j2) * phi_i3 (x_j3)
	phi2(k,l)=phi(i1,j1)*phi(i2,j2)*phi(i3,j3)
c	dphi2x(k,l) = (d/d(x1)) [ phi_i1 (x1) * phi_i2 (x2) * phi_i3 (x3) ] | x1=x_j1,x2=x_j2,x3=x_j3
	dphi2x(k,l)=dphi(i1,j1)*phi(i2,j2)*phi(i3,j3)
c	dphi2y(k,l) = (d/d(x2)) [ phi_i1 (x1) * phi_i2 (x2) * phi_i3 (x3) ] | x1=x_j1,x2=x_j2,x3=x_j3
	dphi2y(k,l)=phi(i1,j1)*dphi(i2,j2)*phi(i3,j3)
c	dphi2z(k,l) = (d/d(x3)) [ phi_i1 (x1) * phi_i2 (x2) * phi_i3 (x3) ] | x1=x_j1,x2=x_j2,x3=x_j3
	dphi2z(k,l)=phi(i1,j1)*phi(i2,j2)*dphi(i3,j3)
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
c	3D weights
	do j1=1,lmax+2
	do j2=1,lmax+2
	do j3=1,lmax+2
	l=(lmax+2)*((lmax+2)*(j3-1)+j1-1)+j2
	wt2(l)=wt(j1)*wt(j2)*wt(j3)
c		write(6,*)'BASIS: wt2(',l,')=',wt2(l)
	enddo
	enddo
	enddo
	return
	end	

	subroutine legend(x,x0,dx0)
c	Returns the Legendre function P sub n (x)  
c	for n=1 through lmax+1. 
c	Array x0 contains P sub n (x)
c	Array dx0 contains (d/dx) P sub n (x).
	parameter (lmax=3)
	real*8 x,tox,bj,bjm,bjp,x0(lmax+1),dx0(lmax+1)
cOLD	if(x.lt.-1.0d0.or.x.gt.1.0d0) pause 
cOLD     &	'x must be nonnegative in legend'
		tox=2.d0*x
		bjm=1.d0
		bj=x
		x0(1)=x
		dx0(1)=1.d0
		n=lmax+1
		do 11 j=1,n-1
		  bjp=((dble(j)+0.5d0)*tox*bj-dble(j)*bjm)/dble(j+1)
		  bjm=bj
		  bj=bjp
		  x0(j+1)=bj
	  dx0(j+1)=dble(j+1)*bjm+x*dx0(j)
11		continue
	return
	end  

c	subroutine gl
cc	Determine Gauss-Lobatto quadrature points for n=6 (Vosse's N=5)
cc	There are four interior points plus +-1
c	parameter (lmax=4)
c	real*8 wt,x
c	common/glpts/wt(lmax+2),x(lmax+2)
cc	See file SpectralElementMethod/Lobatto-quadrature-pt2
c	wt(1)=1.d0/15.d0
c	wt(2)=(14.d0-dsqrt(7.d0))/30.d0
c	wt(3)=(14.d0+dsqrt(7.d0))/30.d0
c	wt(4)=wt(3)
c	wt(5)=wt(2)
c	wt(6)=wt(1)
c	x(1)=-1.d0
c	x(2)=-dsqrt((7.d0+2.d0*sqrt(7.d0))/21.d0)
c	x(3)=-dsqrt((7.d0-2.d0*sqrt(7.d0))/21.d0)
c	x(4)=-x(3)
c	x(5)=-x(2)
c	x(6)=-x(1)
c	return
c	end

	subroutine gl
c	Determine Gauss-Lobatto quadrature points for n=5 (Vosse's N=4)
c	There are three interior points plus +-1
	parameter (lmax=3)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c	See file SpectralElementMethod/Lobatto-quadrature-pt2
	wt(1)=1.d0/10.d0
	wt(2)=49.d0/90.d0
	wt(3)=32.d0/45.d0
	wt(4)=wt(2)
	wt(5)=wt(1)
	x(1)=-1.d0
	x(2)=-dsqrt(3.d0/7.d0)
	x(3)=0.d0
	x(4)=-x(2)
	x(5)=-x(1)
	return
	end

	function phiv(ilx,ily,ilz,xval,yval,zval)
c	Evaluate phi_[ilx](xval) * phi_[llz](zval) at (xval,zval)
c	not necessarily at nodal points.
	real*8 y,legi
	real*8 phiv,phix,phiy,phiz
	parameter (lmax=3)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
	real*8 x0(lmax+1),dx0(lmax+1)
c
	y=x(ilx)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(xval)
	  if((y-x(ilx)).lt.0.001d0.and.(y-x(ilx)).ge.0.d0) y=x(ilx)+0.001
	  if((y-x(ilx)).gt.-0.001d0.and.(y-x(ilx)).le.0.d0) y=x(ilx)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phix=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilx))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ily)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(yval)
	  if((y-x(ily)).lt.0.001d0.and.(y-x(ily)).ge.0.d0) y=x(ily)+0.001
	  if((y-x(ily)).gt.-0.001d0.and.(y-x(ily)).le.0.d0) y=x(ily)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiy=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ily))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ilz)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(zval)
	  if((y-x(ilz)).lt.0.001d0.and.(y-x(ilz)).ge.0.d0) y=x(ilz)+0.001
	  if((y-x(ilz)).gt.-0.001d0.and.(y-x(ilz)).le.0.d0) y=x(ilz)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiz=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilz))) * (1.d0 -y*y)*dx0(lmax+1) 
	phiv=phix*phiy*phiz
	return
	end

	function dphivx(ilx,ily,ilz,xval,yval,zval)
c	Evaluate d[phi_[ilx](x)/dx]|x=xval * phi_[ily](yval) * phi_[ilz](zval) at (xval,yval,zval)
c	not necessarily at nodal points.
c	Do not include multiplication by (2./dx)
	real*8 y,legi
	real*8 dphivx,phix,dphix,phiy,phiz
	parameter (lmax=3)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
	real*8 x0(lmax+1),dx0(lmax+1)
c
	y=x(ilx)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(xval)
	  if((y-x(ilx)).lt.0.001d0.and.(y-x(ilx)).ge.0.d0) y=x(ilx)+0.001
	  if((y-x(ilx)).gt.-0.001d0.and.(y-x(ilx)).le.0.d0) y=x(ilx)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phix=-1.d0/(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilx))) * (1.d0 -y*y)*dx0(lmax+1) 
	dphix=(-phix+x0(lmax+1)/legi)/(y-x(ilx))
	y=x(ily)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(yval)
	  if((y-x(ily)).lt.0.001d0.and.(y-x(ily)).ge.0.d0) y=x(ily)+0.001
	  if((y-x(ily)).gt.-0.001d0.and.(y-x(ily)).le.0.d0) y=x(ily)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiy=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ily))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ilz)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(zval)
	  if((y-x(ilz)).lt.0.001d0.and.(y-x(ilz)).ge.0.d0) y=x(ilz)+0.001
	  if((y-x(ilz)).gt.-0.001d0.and.(y-x(ilz)).le.0.d0) y=x(ilz)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiz=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilz))) * (1.d0 -y*y)*dx0(lmax+1) 
	dphivx=dphix*phiy*phiz
	return
	end

	function dphivy(ilx,ily,ilz,xval,yval,zval)
c	Evaluate phi_[ilx](xval) * d[phi_[ily](y)/dy]|y=yval * phi_[ilz](zval) at (xval,yval,zval)
c	not necessarily at nodal points.
c	Do not include multiplication by (2./dy)
	real*8 y,legi
	real*8 dphivy,phix,phiy,dphiy,phiz
	parameter (lmax=3)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
	real*8 x0(lmax+1),dx0(lmax+1)
c
	y=x(ilx)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(xval)
	  if((y-x(ilx)).lt.0.001d0.and.(y-x(ilx)).ge.0.d0) y=x(ilx)+0.001
	  if((y-x(ilx)).gt.-0.001d0.and.(y-x(ilx)).le.0.d0) y=x(ilx)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phix=-1.d0/(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilx))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ily)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(yval)
	  if((y-x(ily)).lt.0.001d0.and.(y-x(ily)).ge.0.d0) y=x(ily)+0.001
	  if((y-x(ily)).gt.-0.001d0.and.(y-x(ily)).le.0.d0) y=x(ily)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiy=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ily))) * (1.d0 -y*y)*dx0(lmax+1) 
	dphiy=(-phiy+x0(lmax+1)/legi)/(y-x(ily))
	y=x(ilz)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(zval)
	  if((y-x(ilz)).lt.0.001d0.and.(y-x(ilz)).ge.0.d0) y=x(ilz)+0.001
	  if((y-x(ilz)).gt.-0.001d0.and.(y-x(ilz)).le.0.d0) y=x(ilz)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiz=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilz))) * (1.d0 -y*y)*dx0(lmax+1) 
	dphivy=phix*dphiy*phiz
	return
	end

	function dphivz(ilx,ily,ilz,xval,yval,zval)
c	Evaluate phi_[ilx](xval) * phi_[ily](yval) * d[phi_[ilz](z)/dz]|z=zval at (xval,zval)
c	not necessarily at nodal points.
c	Do not include multiplication by (2./dz)
	real*8 y,legi
	real*8 dphivz,phix,phiy,dphiz,phiz
	parameter (lmax=3)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
	real*8 x0(lmax+1),dx0(lmax+1)
c
	y=x(ilx)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(xval)
	  if((y-x(ilx)).lt.0.001d0.and.(y-x(ilx)).ge.0.d0) y=x(ilx)+0.001
	  if((y-x(ilx)).gt.-0.001d0.and.(y-x(ilx)).le.0.d0) y=x(ilx)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phix=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilx))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ily)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(yval)
	  if((y-x(ily)).lt.0.001d0.and.(y-x(ily)).ge.0.d0) y=x(ily)+0.001
	  if((y-x(ily)).gt.-0.001d0.and.(y-x(ily)).le.0.d0) y=x(ily)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiy=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ily))) * (1.d0 -y*y)*dx0(lmax+1) 
	y=x(ilz)
	call legend(y,x0,dx0)
	legi=x0(lmax+1)
	y=dble(zval)
	  if((y-x(ilz)).lt.0.001d0.and.(y-x(ilz)).ge.0.d0) y=x(ilz)+0.001
	  if((y-x(ilz)).gt.-0.001d0.and.(y-x(ilz)).le.0.d0) y=x(ilz)-0.001
	  if(y.gt.1.d0) y=1.d0
	  if(y.lt.-1.d0) y=-1.d0
	call legend(y,x0,dx0)
	phiz=-1./(dble(lmax+1)*dble(lmax+2)*legi*(y-x(ilz))) * (1.d0 -y*y)*dx0(lmax+1) 
	dphiz=(-phiz+x0(lmax+1)/legi)/(y-x(ilz))
	dphivz=phix*phiy*dphiz
	return
	end

