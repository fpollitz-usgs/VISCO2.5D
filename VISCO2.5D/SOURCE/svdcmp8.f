	subroutine svdcmp(a,m,n,mp,np,w,v)
c       Given a matrix a, with logical dimensions m by n and physical
c       dimensions mp by np, this routine computes its singular value
c       decomposition, a=u.w.v(transpose).  The matrix u replaces a on
c       output.  The diagonal matrix of singular values w is output as
c       a vector w.  The matrix v (not the transpose of v) is output
c       as v.  m must be greater than or equal to n; if it is smaller,
c       then a should be filled up to square with zero rows.
	parameter (nmax=2528)
	implicit real*8 (a-h,o-z)
	dimension a(mp,np),w(np),v(np,np),rv1(nmax)
c        if(m.lt.n) pause 'You must augment a with extra zeros'
c       Householder reduction to bidiagonal form.
	g=0.0d0
	scale=0.0d0
	anorm=0.0d0
	do 25 i=1,n
	l=i+1
	rv1(i)=scale*g
	g=0.0d0
	s=0.0d0
	scale=0.0d0
	if(i.le.m) then
		do 11 k=i,m
		scale=scale+abs(a(k,i))
11              continue
		if(scale.ne.0.0d0) then
		  do 12 k=i,m
		  a(k,i)=a(k,i)/scale
		  s=s+a(k,i)*a(k,i)
12                continue
		  f=a(i,i)
		  g=-sign(dsqrt(s),f)
		  h=f*g-s
		  a(i,i)=f-g
		  if(i.ne.n) then
			do 15 j=l,n
			s=0.0d0
			do 13 k=i,m
			s=s+a(k,i)*a(k,j)
13                      continue
			f=s/h
			do 14 k=i,m
			a(k,j)=a(k,j)+f*a(k,i)
14                      continue
15                      continue
		  endif
		  do 16 k=i,m
		  a(k,i)=scale*a(k,i)
16                continue
		endif
	endif
	w(i)=scale*g
	g=0.0d0
	s=0.0d0
	scale=0.0d0
	if(i.le.m.and.i.ne.n) then
		do 17 k=l,n
		scale=scale+abs(a(i,k))
17              continue
		if(scale.ne.0.0d0) then
		  do 18 k=l,n
		  a(i,k)=a(i,k)/scale
		  s=s+a(i,k)*a(i,k)
18                continue
		  f=a(i,l)
		  g=-sign(dsqrt(s),f)
		  h=f*g-s
		  a(i,l)=f-g
		  do 19 k=l,n
		  rv1(k)=a(i,k)/h
19                continue
		  if(i.ne.m) then
			do 23 j=l,m
			s=0.0d0
			do 21 k=l,n
			s=s+a(j,k)*a(i,k)
21                      continue
			do 22 k=l,n
			a(j,k)=a(j,k)+s*rv1(k)
22                      continue
23                      continue
		  endif
		  do 24 k=l,n
		  a(i,k)=scale*a(i,k)
24                continue
		endif
	endif
	anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25      continue
c       Accumulation of right-hand transformations.
	do 32 i=n,1,-1
	if(i.lt.n) then
		if(g.ne.0.0d0) then
		  do 26 j=l,n
		  v(j,i)=(a(i,j)/a(i,l))/g
26                continue
		  do 29 j=l,n
		  s=0.0d0
		  do 27 k=l,n
		  s=s+a(i,k)*v(k,j)
27                continue
		  do 28 k=l,n
		  v(k,j)=v(k,j)+s*v(k,i)
28                continue
29                continue
		endif
		do 31 j=l,n
		v(i,j)=0.0d0
		v(j,i)=0.0d0
31              continue
	endif
	v(i,i)=1.0d0
	g=rv1(i)
	l=i
32      continue
c       Accumulation of left-hand transformations
	do 39 i=n,1,-1
	l=i+1
	g=w(i)
	if(i.lt.n) then
		do 33 j=l,n
		a(i,j)=0.0d0
33              continue
	endif
	if(g.ne.0.0d0) then
		g=1.0d0/g
		if(i.ne.n) then
		  do 36 j=l,n
		  s=0.0d0
		  do 34 k=l,m
		  s=s+a(k,i)*a(k,j)
34                continue
		  f=(s/a(i,i))*g
		  do 35 k=i,m
		  a(k,j)=a(k,j)+f*a(k,i)
35                continue
36                continue
		endif
		do 37 j=i,m
		a(j,i)=a(j,i)*g
37              continue
	else
		do 38 j=i,m
		a(j,i)=0.0d0
38              continue
	endif
	a(i,i)=a(i,i)+1.0d0
39      continue
c       Diagonalization of the bidiagonal form.
	do 49 k=n,1,-1
	do 48 its=1,40
	do 41 l=k,1,-1
	nm=l-1
	if((abs(rv1(l))+anorm).eq.anorm) go to 2
	if((abs(w(nm))+anorm).eq.anorm) go to 1
41      continue
1       c=0.0d0
	s=1.0d0
	do 43 i=l,k
	f=s*rv1(i)
	if((abs(f)+anorm).ne.anorm) then
		g=w(i)
		h=dsqrt(f*f+g*g)
		w(i)=h
		h=1.0d0/h
		c=(g*h)
		s=-(f*h)
		do 42 j=1,m
		y=a(j,nm)
		z=a(j,i)
		a(j,nm)=(y*c)+(z*s)
		a(j,i)=-(y*s)+(z*c)
42              continue
	endif
43      continue
2       z=w(k)
	if(l.eq.k) then
		if(z.lt.0.0d0) then
		  w(k)=-z
		  do 44 j=1,n
		  v(j,k)=-v(j,k)
44                continue
		endif
		go to 3
	endif
	if(its.eq.40) write(6,*)'No convergence in 40 iterations'
	if(its.eq.40) stop
	x=w(l)
	nm=k-1
	y=w(nm)
	g=rv1(nm)
	h=rv1(k)
	f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
	g=dsqrt(f*f+1.0d0)
	f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
c       Next QR transformation.
	c=1.0d0
	s=1.0d0
	do 47 j=l,nm
	i=j+1
	g=rv1(i)
	y=w(i)
	h=s*g
	g=c*g
	z=dsqrt(f*f+h*h)
	rv1(j)=z
	c=f/z
	s=h/z
	f=(x*c)+(g*s)
	g=-(x*s)+(g*c)
	h=y*s
	y=y*c
	do 45 nm=1,n
	x=v(nm,j)
	z=v(nm,i)
	v(nm,j)=(x*c)+(z*s)
	v(nm,i)=-(x*s)+(z*c)
45      continue
	z=dsqrt(f*f+h*h)
	w(j)=z
	if(z.ne.0.0d0) then
		z=1.0d0/z
		c=f*z
		s=h*z
	endif
	f=(c*g)+(s*y)
	x=-(s*g)+(c*y)
	do 46 nm=1,m
	y=a(nm,j)
	z=a(nm,i)
	a(nm,j)=(y*c)+(z*s)
	a(nm,i)=-(y*s)+(z*c)
46      continue
47      continue
	rv1(l)=0.0d0
	rv1(k)=f
	w(k)=x
48      continue
3       continue
49      continue
	return
	end
	subroutine svdksb(u,w,v,m,n,mp,np,b,x)
c       Solves (a).x=b for a vector x, where (a) is specified by the
c       arays u,w,v as returned by svdcmp.  m and n are the logical
c       dimensions of (a), and will be equal for square matrices.
c       mp and np are the physical dimensions of (a).  b is the input
c       right-hand side.  x is the output solution vector.  No input
c       quantities are destroyed, so the routine may be called sequentially
c       with different b's.
	parameter (nmax=2528)
	implicit real*8 (a-h,o-z)
	dimension u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(nmax)
	common/smooth/wm(370)
	do 12 j=1,n
	s=0.d0
        if(w(j).ne.0.) then
		do 11 i=1,m
		s=s+u(i,j)*b(i)
11              continue
c		s=s*w(j)/(w(j)*w(j)+wm(j))
                s=s/w(j)
        endif
	tmp(j)=s
12      continue
	do 14 j=1,n
	s=0.d0
	do 13 jj=1,n
	s=s+v(j,jj)*tmp(jj)
13      continue
	x(j)=s
14      continue
	return
	end
	subroutine svdksc(u,w,v,m,n,mp,np,b,x)
c       Solves (a).x=b for a vector x, where (a) is specified by the
c       arays u,w,v as returned by svdcmp.  m and n are the logical
c       dimensions of (a), and will be equal for square matrices.
c       mp and np are the physical dimensions of (a).  b is the input
c       right-hand side.  x is the output solution vector.  No input
c       quantities are destroyed, so the routine may be called sequentially
c       with different b's.
	parameter (nmax=2528)
	implicit real*8 (a-h,o-z)
	dimension u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(nmax)
	do 12 j=1,n
	s=0.d0
	if(w(j).ne.0.) then
		do 11 i=1,m
		s=s+u(i,j)*b(i)
11              continue
		s=s/w(j)
	endif
	tmp(j)=s
12      continue
	do 14 j=1,n
	s=0.d0
	do 13 jj=1,n
	s=s+v(j,jj)*tmp(jj)
13      continue
	x(j)=s
14      continue
	return
	end
