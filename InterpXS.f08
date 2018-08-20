subroutine InterpXS(process,processI,Ioffset,ChS)
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18
!*******************************************************************************
implicit real*8(a-h,o-z)

integer lin1, lin2 !Linear interpolation boundaries.
!lin1 is the energy at which the FIRST linear interpolator STOPS.
!lin2 is the energy at which the SECOND linear interpolator STARTS.
integer, intent(in) :: ChS
integer nEng !Number of charge states and energies.

parameter (nEng=12)

real*8, dimension(ChS,nEng), intent(in) :: process
real*8, dimension(ChS,25000), intent(out) :: ProcessI
real*8, dimension(nEng) :: processVector,pV2
real*8, dimension(25000) :: eVector !Charge states and energy vect.
integer, dimension(ChS,2), intent(in) :: Ioffset !Interpolator offset
real*8 EngData(nEng)

data EngData /1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,&
              10000.0,25000.0/
!*******************************************************************************
!open(unit=200,file='dcdsI.txt',status='unknown')
processI=0.0 !Final interpolated matrix
do i=1,ChS
  E=1.0
  processVector=0.0 !Vector of XS data (One charge state)
  pV2=0.0 !Second derivative of process vector
  eVector=0.0 !Energy vector
  lin1 = Ioffset(i,1)
  lin2 = Ioffset(i,2)
  do j=1,nEng
    processVector(j)=process(i,j)
  end do
  do j=1,lin1
    if (j.le.10) then
      pVI = log(processVector(1))+(log(E)-log(EngData(1)))*&
            (log(processVector(2))-log(processVector(1)))/&
            (log(EngData(2))-log(EngData(1)))
    elseif(j.le.50) then
      pVI = log(processVector(2))+(log(E)-log(EngData(2)))*&
            (log(processVector(3))-log(processVector(2)))/&
            (log(EngData(3))-log(EngData(2)))
    elseif(j.le.75) then
      pVI = log(processVector(3))+(log(E)-log(EngData(3)))*&
            (log(processVector(4))-log(processVector(3)))/&
            (log(EngData(4))-log(EngData(3)))
    elseif(j.le.100) then
      pVI = log(processVector(4))+(log(E)-log(EngData(4)))*&
            (log(processVector(5))-log(processVector(4)))/&
            (log(EngData(5))-log(EngData(4)))
    elseif(j.le.200) then
      pVI = log(processVector(5))+(log(E)-log(EngData(5)))*&
            (log(processVector(6))-log(processVector(5)))/&
            (log(EngData(6))-log(EngData(5)))
    elseif(j.le.500) then
      pVI = log(processVector(6))+(log(E)-log(EngData(6)))*&
            (log(processVector(7))-log(processVector(6)))/&
            (log(EngData(7))-log(EngData(6)))
    elseif(j.le.1000) then
      pVI = log(processVector(7))+(log(E)-log(EngData(7)))*&
            (log(processVector(8))-log(processVector(7)))/&
            (log(EngData(8))-log(EngData(7)))
    elseif(j.le.2000) then
      pVI = log(processVector(8))+(log(E)-log(EngData(8)))*&
            (log(processVector(9))-log(processVector(8)))/&
            (log(EngData(9))-log(EngData(8)))
    elseif(j.le.5000) then
      pVI = log(processVector(9))+(log(E)-log(EngData(9)))*&
            (log(processVector(10))-log(processVector(9)))/&
            (log(EngData(10))-log(EngData(9)))
    elseif(j.le.10000) then
      pVI = log(processVector(10))+(log(E)-log(EngData(10)))*&
            (log(processVector(11))-log(processVector(10)))/&
            (log(EngData(11))-log(EngData(10)))
    elseif(j.le.25000) then
      pVI = log(processVector(11))+(log(E)-log(EngData(11)))*&
            (log(processVector(12))-log(processVector(11)))/&
            (log(EngData(12))-log(EngData(11)))
    end if
    if(i.eq.1) ProcessI(i,j)=exp(pVI) !Interpolated values of process vector
    if(i.eq.2) ProcessI(i,j)=exp(pVI)
    if(i.eq.3) ProcessI(i,j)=exp(pVI)
    if(i.eq.4) ProcessI(i,j)=exp(pVI)
    if(i.eq.5) ProcessI(i,j)=exp(pVI)
    if(i.eq.6) ProcessI(i,j)=exp(pVI)
    if(i.eq.7) ProcessI(i,j)=exp(pVI)
    if(i.eq.8) ProcessI(i,j)=exp(pVI)
    if(i.eq.9) ProcessI(i,j)=exp(pVI)
    if(i.eq.10) ProcessI(i,j)=exp(pVI)
    E = E+1.0
!    write(*,*) pVI, exp(pVI), EngData(1), EngData(2), E, processVector(1), processVector(2),&
!    processVector(1)-processVector(2)
  end do
!  stop
  do j=lin1+1,lin2
    call spline(log(EngData),log(processVector),nEng,pV2)
    call splineinterp(log(E),log(EngData),log(processVector),nEng,pV2,pVI)
    if(exp(pVI).lt.0) write(*,*) "(XS9)WARNING: CROSS-SECTION LESS THAN ZERO"
    if(i.eq.1) ProcessI(i,j)=exp(pVI) !Interpolated values of process vector
    if(i.eq.2) ProcessI(i,j)=exp(pVI)
    if(i.eq.3) ProcessI(i,j)=exp(pVI)
    if(i.eq.4) ProcessI(i,j)=exp(pVI)
    if(i.eq.5) ProcessI(i,j)=exp(pVI)
    if(i.eq.6) ProcessI(i,j)=exp(pVI)
    if(i.eq.7) ProcessI(i,j)=exp(pVI)
    if(i.eq.8) ProcessI(i,j)=exp(pVI)
    if(i.eq.9) ProcessI(i,j)=exp(pVI)
    if(i.eq.10) ProcessI(i,j)=exp(pVI)
    eVector(j)=E
    E=E+1.0
  end do
  do j=lin2+1,25000
    if (j.le.10) then
      pVI = log(processVector(1))+(log(E)-log(EngData(1)))*&
            (log(processVector(2))-log(processVector(1)))/&
            (log(EngData(2))-log(EngData(1)))
    elseif(j.le.50) then
      pVI = log(processVector(2))+(log(E)-log(EngData(2)))*&
            (log(processVector(3))-log(processVector(2)))/&
            (log(EngData(3))-log(EngData(2)))
    elseif(j.le.75) then
      pVI = log(processVector(3))+(log(E)-log(EngData(3)))*&
            (log(processVector(4))-log(processVector(3)))/&
            (log(EngData(4))-log(EngData(3)))
    elseif(j.le.100) then
      pVI = log(processVector(4))+(log(E)-log(EngData(4)))*&
            (log(processVector(5))-log(processVector(4)))/&
            (log(EngData(5))-log(EngData(4)))
    elseif(j.le.200) then
      pVI = log(processVector(5))+(log(E)-log(EngData(5)))*&
            (log(processVector(6))-log(processVector(5)))/&
            (log(EngData(6))-log(EngData(5)))
    elseif(j.le.500) then
      pVI = log(processVector(6))+(log(E)-log(EngData(6)))*&
            (log(processVector(7))-log(processVector(6)))/&
            (log(EngData(7))-log(EngData(6)))
    elseif(j.le.1000) then
      pVI = log(processVector(7))+(log(E)-log(EngData(7)))*&
            (log(processVector(8))-log(processVector(7)))/&
            (log(EngData(8))-log(EngData(7)))
    elseif(j.le.2000) then
      pVI = log(processVector(8))+(log(E)-log(EngData(8)))*&
            (log(processVector(9))-log(processVector(8)))/&
            (log(EngData(9))-log(EngData(8)))
    elseif(j.le.5000) then
      pVI = log(processVector(9))+(log(E)-log(EngData(9)))*&
            (log(processVector(10))-log(processVector(9)))/&
            (log(EngData(10))-log(EngData(9)))
    elseif(j.le.10000) then
      pVI = log(processVector(10))+(log(E)-log(EngData(10)))*&
            (log(processVector(11))-log(processVector(10)))/&
            (log(EngData(11))-log(EngData(10)))
    elseif(j.le.25000) then
      pVI = log(processVector(11))+(log(E)-log(EngData(11)))*&
            (log(processVector(12))-log(processVector(11)))/&
            (log(EngData(12))-log(EngData(11)))
    end if
    if(i.eq.1) ProcessI(i,j)=exp(pVI) !Interpolated values of process vector
    if(i.eq.2) ProcessI(i,j)=exp(pVI)
    if(i.eq.3) ProcessI(i,j)=exp(pVI)
    if(i.eq.4) ProcessI(i,j)=exp(pVI)
    if(i.eq.5) ProcessI(i,j)=exp(pVI)
    if(i.eq.6) ProcessI(i,j)=exp(pVI)
    if(i.eq.7) ProcessI(i,j)=exp(pVI)
    if(i.eq.8) ProcessI(i,j)=exp(pVI)
    if(i.eq.9) ProcessI(i,j)=exp(pVI)
    if(i.eq.10) ProcessI(i,j)=exp(pVI)
    E = E+1.0
!    write(*,*) pVI, exp(pVI), EngData(1), EngData(2), E, processVector(1), processVector(2),&
!    processVector(1)-processVector(2)
  end do
end do
return
end subroutine
