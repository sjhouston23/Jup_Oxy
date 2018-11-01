subroutine EjectedElectron(E,c,qtmp,eprobfunc,aprobfunc,electron_energy,&
  electron_angle,ebin)
!program EjectedElectron
!*******************************************************************************
!* Created by Stephen J. Houston 1.23.17
!*******************************************************************************
!* This subroutine will calculate the energy and angle of an ejected
!* electron for a given ion energy, ion charge-state, & collision-type.
!* Taking the probability distribution calculated (see ProbDist.f08)
!* using the singly differential cross sections as a function of the
!* ejected electron energy or angle, a monte carlo method similar to
!* that used in CollisionSim.f08 is used to calculate the ejected
!* electron's energy and angle.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of the beam
!*		Type: Real
!*		Units: keV/u
!*
!*		q --> Charge state
!*		Type: Integer
!*		Units: None
!*
!*		c --> Process or collision-type (see edist & adist indexing comments below)
!*		Type: Integer
!*		Units: None
!*
!*    eprobfunc --> Energy probability distribution functions,
!*    returned from ProbDist.f08
!*    Type: Real rank 4
!*    Units: None
!*
!*    aprobfunc --> Angle probability distribution functions,
!*    returned from ProbDist.f08
!*    Type: Real rank 4
!*    Units: None
!*
!*  Returns:
!*    electron_energy --> energy of ejected electron (2 stream binning)
!*    Type: Real
!*    Units: eV
!*
!*    electron_angle --> angle of ejected electron
!*    Type: Real
!*    Units: degrees (0-180)
!
!*******************************************************************************

implicit none

!**************************** Variable Declaration *****************************

integer, intent(in) :: c, qtmp !* c=process, q=IQQ
integer, intent(out) :: ebin
integer i,j,k,t,ne,na,numb,q
!* i=ion energy bin
!* j=corresponding electron energy/angle
!* ne=number of energy bins
parameter(ne=2600, na=1800, numb=4)
real*8, intent(in) :: E
real*8 electron_energytmp1, electron_energytmp2
real*8 electron_angletmp1, electron_angletmp2
real*8, intent(out) :: electron_energy, electron_angle
real, dimension(numb) :: rvec
real randomNumber !Rando number for angle
real*8 f !fraction of ion energy bins
real*8, intent(in), dimension(6,9,13,ne) :: eprobfunc
real*8, intent(in), dimension(6,9,13,na) :: aprobfunc
real*8, dimension(12) :: ion_energies
DATA ion_energies/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,&
                  10000.0,25000.0/

!******************************** Main Program *********************************

ebin=0
i = 0
k = 0
f = 0
q = qtmp
if(q.gt.1)then !Treat negative ion electron ejection same as neutral
  q = q-1 !Need to adjust the arrays since they don't include negative ions
end if !But if q equals 1, then treat it as a neutral oxy electron distribution

do t=12,1,-1 !* Loop through energy bins of differential cross sections
  if(E.ge.real(ion_energies(t-1)+(ion_energies(t)-ion_energies(t-1))/2))then
    i = t
    k = 1
    goto 10
  elseif(E.ge.real(ion_energies(t-1)))then
    i = t-1
    k = 2
    goto 10
  endif
enddo
10 continue
if (i.eq.0) write(206,*) 'EjectedElectron.f08: Error! i=0'
!* Want to use f to somewhat interpolate the cross-section for ion energies that
!* lie between energy bins.
if (k.eq.1) f=(E-ion_energies(i-1))/(ion_energies(i)-ion_energies(i-1))
if (k.eq.2) f=(E-ion_energies(i))/(ion_energies(i+1)-ion_energies(i))

!write(*,*) E,ion_energies(i),i,k,f
!write(*,*) E, i, Q, f

call ranlux(rvec,numb)
!write(*,*) rvec(3)

electron_energy = 0.0 !* Initialize ejected electron energy
electron_energytmp1 = 0.0 !* Will calculate two different electron energies for
electron_energytmp2 = 0.0 !* a simple linear interpolation
electron_angle = 0.0 !* Initialize ejected electron angle
electron_angletmp1 = 0.0 !* Will calculate two different electron angles for
electron_angletmp2 = 0.0 !* a simple linear interpolation
if(c.eq.6.and.eprobfunc(c,q,i+1,2300).gt.0.1)then
  electron_energytmp1=10.0
  goto 20
end if
do j=1,ne
!  if(c.eq.5)write(*,*) eprobfunc(c,q,i+1,j),eprobfunc(c,q,1,j)
  if(rvec(2).ge.eprobfunc(c,q,i+1,j)) then !First i variable is the energy bin
    electron_energytmp1=eprobfunc(c,q,1,j)
    goto 20
  end if
end do
20  continue
!* Loop through all of the energies probabilities to find the electron ejection
!* energy.
do j=1,ne
  if(f.lt.0.5)then !Need to know which energy to multiply by (1-f) and f
    if(c.eq.6.and.eprobfunc(c,q,i+2,2300).gt.0.1)then
      electron_energytmp2=10.0
      goto 30
    end if
    if(rvec(2).ge.eprobfunc(c,q,i+2,j))then
      electron_energytmp2=eprobfunc(c,q,1,j)
      if(c.eq.6.and.electron_energytmp2.gt.500)then
        if(eprobfunc(c,q,i+2,2282).ge.0.2)then
          electron_energytmp2=10.0
        end if
      end if
      electron_energy=(1-f)*electron_energytmp1+f*electron_energytmp2
      goto 30 !Once the electron energy is found, get out of the loop
    end if
  elseif(f.ge.0.5)then !Need to know which energy to multiply by (1-f) and f
    if(c.eq.6.and.eprobfunc(c,q,i,2300).gt.0.1)then
      electron_energytmp2=10.0
      goto 30
    end if
    if(rvec(2).ge.eprobfunc(c,q,i,j))then
      electron_energytmp2=eprobfunc(c,q,1,j)
      if(c.eq.6.and.electron_energytmp2.gt.500)then
        if(eprobfunc(c,q,i+2,2282).ge.0.2)then
          electron_energytmp2=10.0
        end if
      end if
      electron_energy=(1-f)*electron_energytmp2+f*electron_energytmp1
      goto 30 !Once the electron energy is found, get out of the loop
    end if
  end if
end do
30  continue
!* Loop through all of the angle probabilities to find the electron ejection
!* angle.
randomNumber=rvec(3)
!if(c.eq.5.or.c.eq.6)randomNumber=rvec(2) !Correlation for SS and DS
do j=1,na
  if(rvec(3).ge.aprobfunc(c,q,i+1,j)) then
    electron_angletmp1=aprobfunc(c,q,1,j)
    goto 40 !Once the electron angle is found, get out of the loop
  end if
end do
40  continue
do j=1,na
  if(f.lt.0.5)then !Need to know which energy to multiply by (1-f) and f
    if(rvec(3).ge.aprobfunc(c,q,i+2,j))then
      electron_angletmp2=aprobfunc(c,q,1,j)
      electron_angle=(1-f)*electron_angletmp1+f*electron_angletmp2
      goto 50 !Once the electron energy is found, get out of the loop
    end if
  elseif(f.ge.0.5)then !Need to know which energy to multiply by (1-f) and f
    if(rvec(3).ge.aprobfunc(c,q,i,j))then
      electron_angletmp2=aprobfunc(c,q,1,j)
      electron_angle=(1-f)*electron_angletmp2+f*electron_angletmp1
      goto 50 !Once the electron energy is found, get out of the loop
    end if
  end if
end do
50  continue
do j=1,ne
  if(electron_energy.le.eprobfunc(c,q,1,j)) then
    !Find to which 2-stream energy bin the electron energy corresponds
    ebin=floor(real(j/10.0))+1
    goto 60
  end if
end do
60  continue
! if(E.lt.6.05.and.E.gt.6.03.and.c.eq.6)write(*,*) E, c,i, q, f, ebin, electron_energy, k, electron_energytmp1, electron_energytmp2,&! eprobfunc(c,q,1,ebin)!, eprobfunc(c,q,1,j-1),&
!  eprobfunc(c,q,1,j), eprobfunc(c,q,1,j+1),rvec(2)
!write(*,*) electron_energy
!end do
return
end subroutine
