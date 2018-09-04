program JEDIOxyPrecip
!*******************************************************************************
!* Created by Stephen J. Houston 9.4.18
!*******************************************************************************
!* This program reads in a JEDI spectrum (*.d2s) and normalizes it based on the
!* energy and the energy bin width. It then multiplies the new input ion flux
!* by normalized (1 input ion/cm^2/s), previously calculated results that
!* correspond to the JEDI energy bins ~(11,15,20,30,47,60,78,121,218,456 keV).
!* Some of the higher energy JEDI bins overlap - I treat them as if they don't.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV !Used for kind=int64
use formatting !Used for formatting.f08
implicit real*8(a-h,o-z) !i,j,k,l,m,n are assumed to be integers

!*******************************************************************************
integer energy,atmosLen,run

parameter(atmosLen=1544,nProc=36,nChS=10) !Atmosphere, processes, charge states
parameter(nE2strBins=260) !Number of 2-stream energy bins
parameter(nOutputFiles=15) !Number of data files from ion precip code
parameter(number_of_energies=10) !Number of JEDI energy bins

real*8 Eion(number_of_energies) !Ion energies
real*8,dimension(atmosLen) :: altitude !Altitude array
real*8,dimension(number_of_energies,atmosLen) :: totalHp,totalH2p,H2Ex
real*8,dimension(number_of_energies,atmosLen,nProc) :: Hp,H2p
real*8,dimension(number_of_energies,atmosLen,nProc,nChS) :: oxygen
real*8,dimension(number_of_energies,atmosLen,nE2strBins) :: prode2stF,prode2stB
!* JEDI variables
real*8,dimension(number_of_energies) :: Jenergy,Jintensity,Jebins,Jflux
!*   Jenergy - JEDI energy bin [keV]
!*   Jintensity - JEDI ion flux [c/s/ster/cm^2/keV]
!*   Jebins - Size of JEDI energy bins [keV]
!*   Jflux - Jintensity values converted to [counts/cm^2/s]

!Same as previous variables except has a leading "J"
real*8 IntegratedXRay(2,nChS) !1 is DE, 2 is CX
real*8,dimension(atmosLen) :: JtotalHp,JtotalH2p,JH2Ex
real*8,dimension(atmosLen,nProc) :: JHp,JH2p
real*8,dimension(atmosLen,nProc,nChS) :: Joxygen
real*8,dimension(atmosLen,nE2strBins) :: Jprode2stF,Jprode2stB

character(len=10) date,version
character(len=12) time
character(len=100) filename,filenames(nOutputFiles)
character(len=1000) HpHeader,Hp2Header

!****************************** Data Declaration *******************************
!* Initial ion enegy input:
data Eion/10.625,15.017,20.225,29.783,46.653,59.770,77.522,120.647,218.125,&
          456.250/ !Juno energy bins from JEDI.
data filenames/'H+_Prod','H2+_Prod','H2_Excite_Prod','Oxy_Neg','Oxy0_','Oxy1_',&
'Oxy2_','Oxy3_','Oxy4_','Oxy5_','Oxy6_','Oxy7_','Oxy8_','2Str_Elect_Fwd',&
'2Str_Elect_Bwd'/
!* Width of JEDI energy bins: (May eventually need to be adjusted)
data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
pi=4.0d0*atan(1.0d0);Jenergy=0.0;Jintensity=0.0;Jbins=0.0;Jflux=0.0
!*************************** Open JEDI Ion Spectrum ****************************
write(version,'("PJ7-1")') !Filename of a JEDI spectrum (.d2s file)
write(filename,'("./JunoData/Spectra/",A,".d2s")') trim(version)
open(unit=100,file=trim(filename),status='old')
write(*,*)
write(*,1000) 'File:','Date:','Time:' !Write out general information
do i=1,25 !Reading in the data measured by JEDI
  if(i.le.2.or.i.ge.4.and.i.le.15)read(100,*)
  if(i.eq.3)read(100,1001) date, time !Read the date and time of the flyby
  if(i.eq.3)write(*,*) trim(version),'  ',date,'  ',time !Write to screen
  if(i.ge.16)read(100,1002) Jenergy(i-15),Jintensity(i-15)
end do
write(*,*)
write(*,1003)'Energy Bin:','JEDI Intensity:','Energy Bin Width:',&
             'Normalized Flux:' !Write out general information
do run=1,number_of_energies !Convert to [counts/cm^2/s]
!* The first 3 energy bins include both sulfur and oxygen. I'm assuming a 1:2
!* ratio of sulfur:oxygen (from SO_2)
  if(run.le.3)Jflux(run)=Jintensity(run)*2*pi*Jebins(run)*2/3
  if(run.ge.4)Jflux(run)=Jintensity(run)*2*pi*Jebins(run)
  write(*,1004)Jenergy(run),Jintensity(run),Jebins(run),Jflux(run)
end do
do i=1,number_of_energies
  write(*,*) Jflux(i)/Jflux(number_of_energies)
end do
close(100) !Close JEDI measurement file
!********************************* Initialize **********************************
energy=0;altitude=0.0;Hp=0.0;totalHp=0.0;H2p=0.0;totalH2p=0.0;H2Ex=0.0
oxygen=0.0;prode2stF=0.0;prode2stB=0.0
!********** Open output data files for each set of initial energies ************
do run=1,number_of_energies !Loop through each initial ion energy
  energy=nint(Eion(run))
  do i=1,nOutputFiles !Open all of the files
    write(filename,'("./Output/Juno/",I0,"keV/",A,"_Comb.dat")') &
          energy,trim(filenames(i))
    filename=trim(filename)
    open(unit=100+i,file=filename,status='old')
  end do
  do i=1,2
    read(101,*) !Read the headers of the first four files
    read(102,*)
    read(103,*)
    read(103,*)
  end do
  read(101,'(1x,A)') HpHeader !Save these headers to output later
  read(102,'(1x,A)') Hp2Header
  do i=1,atmosLen !Hydrogen Ionization/Excitation vs. altitude
    read(101,F01) altitude(i),(Hp(run,i,j),j=1,31),totalHp(run,i)
    read(102,F01) dum,(H2p(run,i,j),j=1,11),totalH2p(run,i)
    read(103,F02) dum,H2Ex(run,i)
  end do
  do i=1,nChS !Loop through every charge state
    read(103+i,*) !Oxygen header
    do j=1,atmosLen
      read(103+i,F01) dum,(oxygen(run,j,k,i),k=1,nProc)
    end do
  end do
  do j=1,nE2strBins !Electron production rate for 2-stream
    read(114,F2Str) (prode2stF(run,i,j),i=atmosLen,1,-1)
    read(115,F2str) (prode2stB(run,i,j),i=atmosLen,1,-1)
  end do
end do !End of do-loop for each energy
do i=1,nOutputFiles
  close(100+i) !close files
end do
!********************************* Initialize **********************************
JHp=0.0;JtotalHp=0.0;JH2p=0.0;JtotalH2p=0.0;JH2Ex=0.0;Joxygen=0.0
Jprode2stF=0.0;Jprode2stB=0.0
!************************* Calculate JEDI Productions **************************
write(*,*) 'Calculating JEDI production rates...'
do run=1,number_of_energies !Loop through every energy bin
  do i=1,atmosLen !Loop through entire atmosphere
    JtotalHp(i)=JtotalHp(i)+totalHp(run,i)*Jflux(run) !Total H+
    JtotalH2p(i)=JtotalH2p(i)+totalH2p(run,i)*Jflux(run) !Total H2+
    JH2Ex(i)=JH2Ex(i)+H2Ex(run,i)*Jflux(run) !Total H2*
    do j=1,nProc !Loop through every processes
      JHp(i,j)=JHp(i,j)+Hp(run,i,j)*Jflux(run) !H+ by process
      JH2p(i,j)=JH2p(i,j)+H2p(run,i,j)*Jflux(run) !H2+ by process
      do k=1,nChS !Loop through every charge state
        Joxygen(i,j,k)=Joxygen(i,j,k)+Oxygen(run,i,j,k)*Jflux(run) !Oxygen
      end do !End charge state loop - k
    end do !End processes loop - j
    do j=1,nE2strBins !Loop through 2-stream energy bins
      Jprode2stF(i,j)=Jprode2stF(i,j)+prode2stF(run,i,j)*Jflux(run) !e- forward
      Jprode2stB(i,j)=Jprode2stB(i,j)+prode2stB(run,i,j)*Jflux(run) !e- backward
    end do !End 2-stream energy bins loop
  end do !End atmsophere loop
end do !End energy bin loop
!************************* Write Out JEDI Productions **************************
write(*,*) 'Writing output files...'
do i=1,nOutputFiles
  write(filename,'("./Output/Juno/",A,"/",A,".dat")') &
        trim(version),trim(filenames(i))
  filename=trim(filename)
  open(unit=200+i,file=filename,status='unknown')
end do
write(201,H01) !H+ header
write(201,*) trim(HpHeader)
write(202,H02) !H2+ header
write(202,*) trim(Hp2Header)
write(203,H03) !H2* header
do i=1,atmosLen !Ionization/Excitation vs. altitude
  write(201,F01) altitude(i),(JHp(i,j),j=1,31),JtotalHp(i)
  write(202,F01) altitude(i),(JH2p(i,j),j=1,11),JtotalH2p(i)
  write(203,F02) altitude(i),JH2Ex(i)
end do
do i=1,nChS !Oxygen production
  write(203+i,*) "Alt [km] ", (HProc(k),k=1,nProc)
  do j=1,atmosLen
    write(203+i,F01) altitude(j),(Joxygen(j,k,i),k=1,nProc)
  end do
end do
!****************************** X-Ray Production *******************************
write(filename,'("./Output/Juno/",A,"/XRay_DE.dat")') trim(version)
open(unit=301,file=trim(filename),status='unknown') !Open X-Ray DE
write(filename,'("./Output/Juno/",A,"/XRay_CX.dat")') trim(version)
open(unit=302,file=trim(filename),status='unknown') !Open X-Ray CX
altDelta=2.0e5 !2 km = 200,000 cm
do i=1,2
  write(300+i,N01) !X-Ray Note
  write(300+i,*) !Extra space
  write(300+i,H11) !Extra space
  write(300+i,*) !Extra space
  write(300+i,H08) !Altitude integrated X-ray production header
  write(300+i,H09) !Charge state header
end do
do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
  do j=1,nChS
    IntegratedXRay(1,j)=IntegratedXRay(1,j)+&
      real(Joxygen(i,27,j)+Joxygen(i,29,j)+Joxygen(i,32,j))*altDelta !DE
    IntegratedXRay(2,j)=IntegratedXRay(2,j)+&
      real(Joxygen(i,19,j)+Joxygen(i,25,j)+Joxygen(i,30,j))*altDelta !CX
  end do
end do
do i=1,2
  write(300+i,F05) altDelta/1e5,(IntegratedXRay(i,j),j=1,nChS)
  write(300+i,*) !Extra space
  write(300+i,H10) !X-Ray production vs. altitude header
  write(300+i,H06) !Charge state header
end do
do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
  write(301,F05) altitude(i),& !X-Ray production from direct excitation
    ((Joxygen(i,27,j)+Joxygen(i,29,j)+Joxygen(i,32,j)),j=1,nChs)
  write(302,F05) altitude(i),& !X-Ray production from charge exchange
    ((Joxygen(i,19,j)+Joxygen(i,25,j)+Joxygen(i,30,j)),j=1,nChs)
end do
close(301)
close(302)
!***************************** Secondary Electrons *****************************
do j=1,nE2strBins !2-Stream electrons, forward and backward
  write(214,F2Str) (Jprode2stF(i,j),i=atmosLen,1,-1)
  write(215,F2Str) (Jprode2stB(i,j),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles !Close the combine output files
  close(200+i)
end do

1000 format(1x,A5,10x,A5,7x,A5)
1001 format(39X,A10,1X,A12)
1002 format(5X,ES12.9,1X,ES13.10)
1003 format(3x,A11,2x,A15,2x,A17,2x,A16)
1004 format(F11.3,2x,F15.3,2x,F17.3,2x,F16.3)
end program