program ProbDist
!*******************************************************************************
!* Created by Stephen J. Houston 3.7.18
!*******************************************************************************
!* Interpolates the ejected electron energy/angle distribution files
!* with a log-cubic-spline interpolator.
!* Number of interpolation points can be adjusted with "no(a)interp".
!* Uses D.R. Schultz's SDXS (2016) to determine the probability
!* distribution functions for energy (2 stream bins) and angle (0-180).
!*
!* ProbDist.f08 ReadElectDist.f08 spline.f08 splineinterp.f08
!*
!* Collision Type:
!*    1) Single Ionization
!*    2) Double Ionization
!*    3) Transfer Ionization
!*    4) Double-Capture Auto Ionization
!*    5) Single Stripping
!*    6) Double stripping
!*
!* Possible Ion Energies:
!*    1, 10, 50, 75, 100, 200, 500, 1000, 2000, 5000, 10000, 25000
!*
!*******************************************************************************

implicit none

!**************************** Variable Declaration *****************************
integer i, j, k, l, m, n, ne, na, nointerp, noainterp, ii, jj, kk
integer noebins, noabins, stepse, stepsa, noainterp2
real*8 stepsize,pi
!parameter(nointerp=19400000, noainterp=18000, stepsize=0.01d0, noebins=149,&
!noabins=54, ne=2600, na=1800,pi=4.0*atan(1.0d0))
parameter(nointerp=260000, noainterp=18000, stepsize=0.01d0,ne=2600, na=1800,&
          pi=4.0*atan(1.0d0))
integer ionE(12),ODC(6) !Conversion for process numbers from OXSDist
real*8, dimension(6,9,13,149) :: edist
real*8, dimension(6,9,13,54) :: adist
real*8 stepsizeb(noainterp)
real*8 OXSDist(38,10,12) !No. of processes, charge states, energies
real*8 probfunce,probfunca
real*8 newenergybins(nointerp), newanglebins(noainterp) !interpolated bin values
real*8 stepsizebins(ne),ssb(nointerp)
real*8 intesdxs(nointerp) !interpolated energy and angle
real*8 intasdxs(noainterp) ! singly differential cross sections
!real*8 twostreambins(ne) !two stream energy bins
real*8 finalabins(na)
real*8 Es, As !Ejected electron energy Es for probability distribution function
real*8 eb, f
real*8 EsXS, totXS, AsXS, totAXS !Ejected electron cross sections up to energy Es/infinity
real*8 integralasdxs,integralesdxs
real*8, dimension(6,9,13,ne) :: eprobfunc
real*8, dimension(6,9,13,na) :: aprobfunc
real*8,allocatable :: varystepsize(:),delAbins(:),delEbins(:),anglebinsrad(:)
real*8,allocatable :: asdxs(:),energybins(:),esdxs(:),anglebins(:),asdxsp(:)
real*8,allocatable :: esdxs2(:),asdxs2(:) !2nd derivative of sdxs
!data twostreambins/20*0.5,70*1.0,10*2.0,20*5.0,10*10.0,20*10.0,10*50.0,10*100.0,40*200.0,&
!        10*400,10*1000,10*2000,10*5000,10*10000.0/
data stepsizebins/200*0.05,700*0.1,100*0.2,200*0.5,100*1.0,200*1.0,100*5.0,100*10.0,400*20.0,&
        100*40,100*100,100*200,100*500,100*1000.0/
data ssb/20000*0.0005,70000*0.001,10000*0.002,20000*0.005,10000*0.01,&
        20000*0.01,10000*0.05,10000*0.1,40000*0.2,10000*0.4,10000*1,10000*2,&
        10000*5,10000*10.0/
data stepsizeb/18000*0.01/
data finalabins/1800*0.1/
data ionE/1, 10, 50, 75, 100, 200, 500, 1000, 2000, 5000, 10000, 25000/
data ODC/1,6,11,16,36,37/ !Oxygen xs Distribution Conversion

! open(10, FILE='./NewElectronDist/ProbDistFunc/eprobfunc.dat', STATUS='unknown')
!open(20, FILE='./NewElectronDist/ProbDistFunc/aprobfuncSinTheta.dat', STATUS='unknown')
!open(30, FILE='./ElectronDist/ProbDistFunc/eprobfuncinterp.dat', STATUS='unknown')
!*****************************************************************************
!         Used when trouble shooting for making changes faster
ii = 6 !Collision
jj = 1 !Charge state
kk = 6 !Ion Energy
!*****************************************************************************
call ReadElectDist(edist,adist) !Read in all electron distribution information
call ReadOXSDist(OXSDist)
n=0
newenergybins=0.0
do i=1,nointerp
  eb=eb+ssb(i)
  newenergybins(i)=eb
end do
!~do i=1,nointerp
!~  newenergybins(i)=i*stepsize
!~end do
!do i=1,nointerp
!  if(i.eq.1)stepsizebins(i)=twostreambins(i)
!  if(i.gt.1)stepsizebins(i)=twostreambins(i)-twostreambins(i-1)
!  write(*,*) stepsizebins(i)
!end do
do i=1,noainterp
  newanglebins(i)=i*stepsize*pi/180.0 !In radians
end do
!Calculating the second derivative of each point first
do i=ii,ii!1,6 !Do it for every collision type
  write(*,*) '------------------NEW COLLISION------------------'
  do j=jj,jj!1,9 !For every charge state
    n=0
    write(*,*)'------------NEW CHARGE STATE------------'
    do k=kk,kk!2,13 !For every initial ion energy (the first column is energy bins)
      intesdxs = 0.0
      intasdxs = 0.0
      totxs = 0.0
      totaxs = 0.0
      if(i.eq.1)then !Schultz had different binning depending on the collision type.
        noebins=133
        noabins=46
      elseif(i.eq.2)then
        noebins=115
        noabins=46
      elseif(i.eq.3)then
        noebins=72
        noabins=54
      elseif(i.eq.4)then
        noebins=64
        noabins=45
      elseif(i.eq.5)then
        noebins=149
        noabins=48
      elseif(i.eq.6)then
        noebins=138
        noabins=43
      end if
      integralesdxs=0.0
      integralasdxs=0.0
      allocate(esdxs(noebins),energybins(noebins),delEbins(noebins))
      allocate(esdxs2(noebins))
      esdxs=0.0;energybins=0.0;delEbins=0.0;esdxs2=0.0
      do l=1,noebins
        esdxs(l) = edist(i,j,k,l)
        energybins(l) = edist(i,j,1,l)
        if(l.eq.1)then
          delEbins(l)=energybins(l)-0.0
        else
          delEbins(l)=energybins(l)-energybins(l-1)
        end if
        write(*,'(F10.3,2x,ES10.2E2,2x,F10.4,2x,ES10.2E2)') energybins(l), esdxs(l),delEbins(l),delEbins(l)*esdxs(l)
      end do
      call simp(noebins,delEbins,esdxs,integralesdxs)
      write(*,*) integralesdxs
      stop
      !noabins=4
      allocate(delAbins(noabins),anglebinsrad(noabins),asdxs(noabins))
      allocate(anglebins(noabins),asdxs2(noabins),asdxsp(noabins))
      delAbins=0.0;anglebinsrad=0.0;asdxs=0.0;anglebins=0.0;asdxs2=0.0;asdxsp=0.0
      anglebinsrad = 0.0
      delAbins = 0.0
      do l=1,noabins
        asdxs(l) = adist(i,j,k,l)
        anglebins(l) = adist(i,j,1,l)
        anglebinsrad(l) = adist(i,j,1,l)*pi/180.0
        asdxsp(l)=asdxs(l)*sin(anglebinsrad(l))!*2*pi
      end do
      do l=1,noabins
        if(l.eq.noabins)then
          delAbins(l)=(pi-anglebinsrad(l))!*pi/180.0
        elseif(l.eq.1)then
          delAbins(l)=(anglebinsrad(l)/2)
        else
          delAbins(l)=(anglebinsrad(l+1)-anglebinsrad(l))!*pi/180.0
        end if
!        write(*,'(F7.2,2x,ES10.2E3,2F8.4)') anglebins(l), asdxsp(l),delAbins(l)!*sin(anglebins(l)*3.14/180.0)
      end do
      call simp(noabins,delAbins,asdxsp,integralasdxs)
!      write(*,*) integralasdxs
!      write(*,*) '1'
!      call spline(log(energybins),log(esdxs),noebins,esdxs2) !Let esdxs=0
!      call spline(log(anglebins),log(asdxsp),noabins,asdxs2)
!      write(*,*) '2'
!goto 1000
      l=1
      do m=1,nointerp
       call splineinterp(log(newenergybins(m)),log(energybins),log(esdxs),&
       noebins,esdxs2,intesdxs(m))
       intesdxs(m)=exp(intesdxs(m))
       !Prevents the function from "blowing up" to infinity
       if(newenergybins(m).gt.energybins(noebins)) intesdxs(m)=0.0
       !Prevents the functions from starting at infinity
       if(newenergybins(m).lt.energybins(1)) intesdxs(m)=0.0
       if(isnan(intesdxs(m))) intesdxs(m)=0.0
       if(intesdxs(m).lt.0) intesdxs(m)=0.0
!       if(m.lt.3700)write(*,*) newenergybins(m), intesdxs(m)
!       if(mod(m,1000).eq.0)write(*,*) newenergybins(m), intesdxs(m), m
!       if(m.ge.nointerp-100)write(*,*) newenergybins(m), intesdxs(m)
!       if(m.gt.100000.and.m.lt.400000.and.mod(m,1000).eq.0)&
!        write(*,*)newenergybins(m), intesdxs(m)
      end do
!      stop
!1000 continue
goto 2000
      l=1
      f=0
      !asdxs=asdxs*1e19
      !noainterp2=150
      do m=1,noainterp
!        call splineinterp(log(newanglebins(m)),log(anglebins),log(asdxsp),&
!        noabins,asdxs2,intasdxs(m))
!        intasdxs(m)=exp(intasdxs(m))
        !If the new angle bin is closer to the next bin, go to it.
        9000 continue
        if(newanglebins(m).lt.anglebinsrad(1))then
          intasdxs(m)=1.00E-30*pi/180.0
          goto 8000
        elseif(l.eq.noabins)then
          intasdxs(m)=asdxsp(l)!*sin(newanglebins(m))*2*pi
          goto 8000
        elseif(newanglebins(m).lt.anglebinsrad(l+1))then
          f=(newanglebins(m)-anglebinsrad(l))/(anglebinsrad(l+1)-anglebinsrad(l))
        else
          l=l+1
          goto 9000
        end if
        intasdxs(m)=((1-f)*asdxsp(l)+f*asdxsp(l+1))!*sin(newanglebins(m))*2*pi
        8000 continue
!        if(m.lt.50)write(*,*) m,l,asdxsp(l),newanglebins(m),anglebins(l), intasdxs(m)
        !Prevents the function from "blowing up" to infinity
      !  if(newanglebins(m).gt.anglebins(noabins)) intasdxs(m)=0.0
        !Prevents the functions from starting at infinity
      !  if(newanglebins(m).lt.anglebins(1)) intasdxs(m)=0.0
      !  if(isnan(intasdxs(m))) intasdxs(m)=0.0
      !  if(intasdxs(m).lt.0) intasdxs(m)=0.0
      !  if(m.gt.0.and.mod(m,1).eq.0)then
      !  if(l.gt.1)  write(*,*) newanglebins(m)*180.0/pi, intasdxs(m)!, asdxsp(l+1),anglebins(l+1)*180.0/pi,&
        !                                              asdxsp(l), anglebins(l)*180.0/pi
      !  end if
  !      if(m.le.100)write(*,*) newanglebins(m), intasdxs(m), asdxs(l-1),anglebins(l-1),&
  !      asdxs(l), anglebins(l)
!        if(m.gt.16000.and.mod(m,100).eq.50)write(*,'(F7.2,2x,ES10.2E3)') newanglebins(m), intasdxs(m)
      end do
2000 continue
      Es=0
      As=0
!      call simp(nointerp,ssb,intesdxs,totxs)
     call simp1(nointerp,stepsize,intesdxs,totxs)
!~      call simp(noainterp,stepsizeb*pi/180,intasdxs,totaxs)
      !call simp1(noainterp2,stepsize,intasdxs,totaxs)
!      write(*,*) totAXS, integralasdxs, OXSDist(ODC(i),j+1,k-1)
!      if(n.eq.0)then
!        write(*,*)'Coll  ChSt  IonEng  iEXSratio  EXSratio  iAXSratio  AXSratio'
!        n=1
!      end if
!      write(*,101) i,j-1,ionE(k-1),1-totxs/OXSDist(ODC(i),j+1,k-1),&
!      1-integralesdxs/OXSDist(ODC(i),j+1,k-1),1-totaxs/OXSDist(ODC(i),j+1,k-1),&
!      1-integralasdxs/OXSDist(ODC(i),j+1,k-1)
      write(*,100) 'Collision: ', i, '  Charge State: ', j-1, ' Ion Energy: ', &
      ionE(k-1), ' Total XS:  ', totxs, integralesdxs, OXSDist(ODC(i),j+1,k-1),&
      1-totxs/OXSDist(ODC(i),j+1,k-1),1-integralesdxs/OXSDist(ODC(i),j+1,k-1)
      ! write(*,100) 'Collision: ', i, '  Charge State: ', j-1, ' Ion Energy: ', &
      ! ionE(k-1), ' Total AXS: ', totaxs, integralasdxs, OXSDist(ODC(i),j+1,k-1),&
      ! 1-totaxs/OXSDist(ODC(i),j+1,k-1), 1-integralasdxs/OXSDist(ODC(i),j+1,k-1)

!      write(*,*)'       bin #      steps   2stream bin     Es Cross Section     TotalXS        Probfunc'

!************ Probability Distribution Function for all interpolated data ************
!***      do l=1,50000!nointerp
!***        probfunc=0.0
!***        EsXS=0.0
!***        call simp(l,stepsize,intesdxs,EsXS)
!***        probfunc=1-(EsXS/totXS)
!***        if(mod(l,10000).eq.0) write(*,*) l, steps, newenergybins(l), EsXS, totXS, probfunc
!***        eprobfunc(i,j,1,l) = newenergybins(l)
!***        eprobfunc(i,j,k,l) = probfunc
!***      end do

!************ Probability Distribution Function for 2 stream binning ************
!goto 3000
      do l=1,ne
        probfunce=0.0
        EsXS=0.0
        Es=0.0
        stepse=0
!~       Es=Es+stepsizebins(l)!twostreambins(l)
        stepse=int(Es/(stepsizebins(l)/100))
!~       stepse=int(Es/stepsize)
        allocate(varystepsize(l*100))
        varystepsize=0.0
        do m=1,l*100
          varystepsize(m)=ssb(m)
          Es=Es+varystepsize(m)
        end do
        if(totxs.eq.0)then
          probfunce=0.0
          goto 5000
        end if
!        call simp(l*100,varystepsize,intesdxs,EsXS)
       if(mod(stepse,2).eq.0)call simp1(stepse,stepsize,intesdxs,EsXS)
       if(mod(stepse,2).eq.1)call simp2(stepse,stepsize,intesdxs,EsXS)
!       write(*,*) '3', l, ne
        if(EsXS.gt.totXS)EsXS=totXS
        if(EsXS.lt.0)EsXS=totXS
        probfunce=1-(EsXS/totXS)
        if(probfunce.gt.1)probfunce=0.0
        if(probfunce.lt.0)probfunce=0.0
        if(isnan(probfunce))probfunce=0.0
!       write(*,*) l, l*100, Es, EsXS, totXS, probfunce
       write(*,*) l, stepse, Es, EsXS, totXS, probfunce
        5000 continue
        eprobfunc(i,j,1,l) = Es
        eprobfunc(i,j,k,l) = probfunce
        deallocate(varystepsize)
      end do
!3000 continue
      goto 4000
      do l=1,na
        probfunca=0.0
        AsXS=0.0
        stepsa=0
        As=As+finalabins(l)!*pi/180
        if(totaxs.eq.0)then
          probfunca=0.0
          goto 6000
        end if
        stepsa=int(As/stepsize)
        call simp1(stepsa,stepsize*pi/180,intasdxs,AsXS)
        probfunca=1-(AsXS/totAXS)
        if(probfunca.lt.0)probfunca=0.0
!        if(l.lt.10)write(*,3001) l, stepsa, As, AsXS, totAXS, probfunca
!        if(mod(l,10).eq.0)write(*,3001) l, stepsa, As, AsXS, totAXS, probfunca
        6000 continue
        aprobfunc(i,j,1,l) = As
        aprobfunc(i,j,k,l) = probfunca
      end do
      deallocate(esdxs,energybins,delEbins,esdxs2)
      deallocate(delAbins,anglebinsrad,asdxs,anglebins,asdxs2,asdxsp)
      4000 continue
    end do
  end do
end do
Es = 0.0
As = 0.0
!do i=1,ne
!  Es=Es+twostreambins(i)
!  write(30,*) Es, (eprobfunc(6,j,1,i),j=1,9) !eprobfunc(Collision,Charge State,Ion Energy,ne)
!end do
! write(10,*) eprobfunc
! write(20,*) aprobfunc
! close(10)
! close(20)
100 format(A11,I1,A16,I1,A13,I5,A12,ES9.3E2,2x,ES9.3E2,2x,ES9.3E2,2x,F8.3,2x,F8.3)
101 format(2x,I1,4x,I1,5x,I5,4x,F8.3,2x,F8.3,3x,F8.3,2x,F8.3)
3001 format (I6,2x,I8,2x,F8.2,3x,ES8.2E2,3x,ES8.2E2,3x,F6.4)
end program

!**********************************************************************
subroutine simp (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx (can't be variable)
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N,i
real*8 F(N), F0, F2(N/2), H(N), INTEG, SF1,SF2
real*8,allocatable,dimension(:) :: F1

if(mod(N,2).eq.0)then
  allocate(F1(N/2-1))
else
  allocate(F1(N/2))
end if
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)
SF1 = 0.0
SF2 = 0.0
!!!write(*,*) F
F1=F(3:N:2)
F2=F(2:N:2)
do i=1,size(F1)
  SF1=SF1+F1(i)*H((i*2)+1)
end do
do i=1,size(F2)
  SF2=SF2+F2(i)*H(i*2)
end do
!SF1=sum(F1) !Summation of all odd elements
!SF2=sum(F2) !Summation of all even elements
F0 = H(1)*F(1) + H(N)*F(N)

INTEG = (F0 + 2*SF2 + 4*SF1)/3
deallocate(F1)

return
end

!**********************************************************************
subroutine simp1 (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx (can't be variable)
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N
real*8 F(N), F0, F1(N/2-1), F2(N/2), H, INTEG, SF1,SF2
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)
!!!write(*,*) F
F1=F(3:N:2)
F2=F(2:N:2)
SF1=sum(F1) !Summation of all odd elements
SF2=sum(F2) !Summation of all even elements
F0 = F(1) + F(N)

INTEG = H * (F0 + 2*SF2 + 4*SF1)/3

return
end

!**********************************************************************
subroutine simp2 (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx (can't be variable)
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N
real*8 F(N), F0, F1(N/2), F2(N/2), H, INTEG, SF1,SF2
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)
!!!write(*,*) F
F1=F(3:N:2)
F2=F(2:N:2)
SF1=sum(F1) !Summation of all odd elements
SF2=sum(F2) !Summation of all even elements
F0 = F(1) + F(N)

INTEG = H * (F0 + 2*SF2 + 4*SF1)/3

return
end
