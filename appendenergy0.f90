program addColumn
implicit none
!real*4 a,b,c,d,e,f,g,h,j,k,l,m,n,o
character(len=100) energyfile
character(len=100) anglefile
character(len=4) collision(6), collisionNew(6)
data collision/'SI','DI','TI','DCAI','SS','DS'/
data collisionNew/'si','di','ti','dcai','ss','ds'/

integer, dimension(6,2) :: q_limit
DATA q_limit/&
1,1,2,3,1,1,&
9,9,9,9,8,7/
real junk,func
real,allocatable,dimension(:) :: a,b,cc,d,e,f,g,h,j,k,l,m,n,o
real,allocatable,dimension(:) :: a1,b1,cc1,d1,e1,f1,g1,h1,j1,k1,l1,m1,nn1,o1
integer i,p,maxlines,nlines,nlines1,nlines2,c,q,n1,n2,err,n3

do c=1,6 !Collision types, see indexing comments above
  do q=q_limit(c,1),q_limit(c,2) !Loop through charge states
!* Open energy file for a given collision type (c) and charge-state (q)
    write(energyfile, '("./ElectronDist1/",A,"_",I0,"_nrg.dat")')& !* energy file name
    trim(collision(c)),q-1 !* Specifies collision type (c) & charge-state (q)
    open(1,file=trim(energyfile),status='old') !* Open energy file

    if(c.eq.4)goto 3000
    write(energyfile, '("./UpdatedElectronDist/sp-",A,"E-O",I0,".dat")')& !* energy file name
    trim(collisionNew(c)),q-1 !* Specifies collision type (c) & charge-state (q)
    open(2,file=trim(energyfile),status='old') !* Open energy file
    3000 continue
    write(energyfile, '("./ElectronDist/",A,"_",I0,"_nrg.dat")')& !* energy file name
    trim(collision(c)),q-1 !* Specifies collision type (c) & charge-state (q)
    open(3,file=trim(energyfile),status='unknown') !* Open energy file


    maxlines=1000
    nlines=0
    do p=1,maxlines
      read(1,*,END=1000) junk
      if(p.eq.maxlines) write(*,*) 'Error: maximum number of lines exceeded. Increase maxlines.'
      nlines=nlines+1
    end do
    1000 continue
    rewind(1)
    write(*,*) c,q
    allocate(a(nlines))
    allocate(b(nlines))
    allocate(cc(nlines))
    allocate(d(nlines))
    allocate(e(nlines))
    allocate(f(nlines))
    allocate(g(nlines))
    allocate(h(nlines))
    allocate(j(nlines))
    allocate(k(nlines))
    do p=1,nlines
      read(1,*) a(p),b(p),cc(p),d(p),e(p),f(p),g(p),h(p),j(p),k(p)
    end do
    if(c.eq.4)then
      nlines1=nlines
      allocate(l(nlines1))
      allocate(m(nlines1))
      allocate(n(nlines1))
      allocate(o(nlines1))
      l=100000.0
      goto 4000
    end if
    nlines1=0
    do p=1,maxlines
      read(2,*,END=2000) junk
      if(p.eq.maxlines) write(*,*) 'Error: maximum number of lines exceeded. Increase maxlines.'
      nlines1=nlines1+1
    end do
    2000 continue
    allocate(l(nlines1))
    allocate(m(nlines1))
    allocate(n(nlines1))
    allocate(o(nlines1))
    rewind(2)
    do p=1,nlines1
      read(2,*) l(p),m(p),n(p),o(p)
    end do
    4000 continue
    n1=1
    n2=1
    n3=1
    err=0
    do i = 1,10000
      if(a(n1).eq.l(n2).and.err.lt.1)then
        n1=n1+1
        n2=n2+1
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 5000
        if(n1.gt.nlines)err=1
        if(n2.gt.nlines1)err=2
      elseif(a(n1).lt.l(n2).and.err.ne.1.or.err.eq.2)then
        n1=n1+1
        n2=n2
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 5000
        if(n1.gt.nlines.and.c.eq.4)goto 5000
        if(n1.gt.nlines)err=1
      elseif(a(n1).gt.l(n2).and.err.ne.2.or.err.eq.1)then
        n1=n1
        n2=n2+1
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 5000
        if(n2.gt.nlines1)err=2
      endif
    end do
    5000 continue
    nlines2=n3
    allocate(a1(nlines2))
    allocate(b1(nlines2),cc1(nlines2),d1(nlines2),e1(nlines2),&
    f1(nlines2),g1(nlines2),h1(nlines2),j1(nlines2),k1(nlines2),l1(nlines2),&
    m1(nlines2),nn1(nlines2),o1(nlines2))
    a1=0.0;b1=0.0;cc1=0.0;d1=0.0;e1=0.0;f1=0.0;g1=0.0;h1=0.0;j1=0.0;k1=0.0
    l1=0.0;m1=0.0;nn1=0.0;o1=0.0
    n1=1
    n2=1
    n3=1
    err=0
    do i = 1,10000
      if(a(n1).eq.l(n2).and.err.lt.1)then
        a1(n3)=a(n1)
        b1(n3)=b(n1)
        cc1(n3)=cc(n1)
        d1(n3)=d(n1)
        e1(n3)=e(n1)
        f1(n3)=f(n1)
        g1(n3)=g(n1)
        h1(n3)=h(n1)
        j1(n3)=j(n1)
        k1(n3)=k(n1)
        l1(n3)=l(n1)
        m1(n3)=m(n1)
        nn1(n3)=n(n1)
        o1(n3)=o(n1)
        n1=n1+1
        n2=n2+1
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 6000
        if(n1.gt.nlines)err=1
        if(n2.gt.nlines1)err=2
      elseif(a(n1).lt.l(n2).and.err.ne.1.or.err.eq.2)then
        a1(n3)=a(n1)
        b1(n3)=b(n1)
        cc1(n3)=cc(n1)
        d1(n3)=d(n1)
        e1(n3)=e(n1)
        f1(n3)=f(n1)
        g1(n3)=g(n1)
        h1(n3)=h(n1)
        j1(n3)=j(n1)
        k1(n3)=k(n1)
        if(n2.le.nlines1)then
          func=(a1(n3)-l(n2-1))/(l(n2)-l(n2-1))
          if(n2.eq.1)then
            m1(n3)=1.0e-30
            nn1(n3)=1.0e-30
            o1(n3)=1.0e-30
          else
            m1(n3)=(1-func)*m(n2-1)+func*m(n2)
            nn1(n3)=(1-func)*n(n2-1)+func*n(n2)
            o1(n3)=(1-func)*o(n2-1)+func*o(n2)
          end if
        else
          m1(n3)=1.0e-30
          nn1(n3)=1.0e-30
          o1(n3)=1.0e-30
        end if
        n1=n1+1
        n2=n2
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 6000
        if(n1.gt.nlines.and.c.eq.4)goto 6000
        if(n1.gt.nlines)err=1
      elseif(a(n1).gt.l(n2).and.err.ne.2.or.err.eq.1)then
        a1(n3)=l(n2)
        m1(n3)=m(n2)
        nn1(n3)=n(n2)
        o1(n3)=o(n2)
        if(n1.le.nlines)then
          func=(a1(n3)-a(n1-1))/(a(n1)-a(n1-1))
          if(n1.eq.1)then
            b1(n3)=1.0e-30
            cc1(n3)=1.0e-30
            d1(n3)=1.0e-30
            e1(n3)=1.0e-30
            f1(n3)=1.0e-30
            g1(n3)=1.0e-30
            h1(n3)=1.0e-30
            j1(n3)=1.0e-30
            k1(n3)=1.0e-30
          else
            b1(n3)=(1-func)*b(n1-1)+func*b(n1)
            cc1(n3)=(1-func)*cc(n1-1)+func*cc(n1)
            d1(n3)=(1-func)*d(n1-1)+func*d(n1)
            e1(n3)=(1-func)*e(n1-1)+func*e(n1)
            f1(n3)=(1-func)*f(n1-1)+func*f(n1)
            g1(n3)=(1-func)*g(n1-1)+func*g(n1)
            h1(n3)=(1-func)*h(n1-1)+func*h(n1)
            j1(n3)=(1-func)*j(n1-1)+func*j(n1)
            k1(n3)=(1-func)*k(n1-1)+func*k(n1)
          endif
        else
          b1(n3)=1.0e-30
          cc1(n3)=1.0e-30
          d1(n3)=1.0e-30
          e1(n3)=1.0e-30
          f1(n3)=1.0e-30
          g1(n3)=1.0e-30
          h1(n3)=1.0e-30
          j1(n3)=1.0e-30
          k1(n3)=1.0e-30
        end if
        n1=n1
        n2=n2+1
        n3=n3+1
        if(n1.gt.nlines.and.n2.gt.nlines1)goto 6000
        if(n2.gt.nlines1)err=2
      endif
      !write(*,*) a(i),b1(i),l(i),func
    end do
    6000 continue
    !do i=1,7!nlines2
    !end do
    do i=1,nlines2-1
      write(3,'(F10.3,12(ES9.2E2))') a1(i),b1(i),cc1(i),d1(i),e1(i),f1(i),&
      g1(i),h1(i),j1(i),k1(i),m1(i),nn1(i),o1(i)
    end do
    deallocate(a,b,cc,d,e,f,g,h,j,k,l,m,n,o)
    deallocate(a1,b1,cc1,d1,e1,f1,g1,h1,j1,k1,l1,m1,nn1,o1)
    !  elseif(a(i).lt.l(i))
  end do
end do
close(1)
close(2)
close(3)
end
