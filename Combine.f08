program comb_files
!*******************************************************************************
!* Created by Stephen J. Houston 6.20.18
!*******************************************************************************
!* This program combines multiple files of the same information into one.
!* Used after making a large run on multiple cores (CRC) to simulate as if a
!* single run had taken place.
!*******************************************************************************

use formatting
implicit real*8(a-h,o-z)

!*******************************************************************************
integer energy,atmosLen

parameter(nOutputFiles=19,MaxnTrials=1000,atmosLen=1544,nProc=36,nChS=10)
parameter(nOxEngBins=5000,nStopPowerEBins=295,nE2strBins=260)

integer trial(MaxnTrials),collisions(8,5),Ccollisions(8,5)
integer,dimension(nStopPowerEBins) :: nSPions,CnSPions

real*8 oxEngBins(nOxEngBins)
real*8,dimension(nChS,nOxEngBins) :: OxyVsEng,COxyVsEng
real*8,dimension(atmosLen) :: altitude,totalHp,totalH2p,H2Ex,CtotalHp,&
                              CtotalH2p,CH2Ex
real*8,dimension(nProc,atmosLen) :: totHp,totH2p,CtotHp,CtotH2p
real*8,dimension(nStopPowerEBins) :: stopPowerEBins,SPvsEng,SigTotvsEng,&
  dEvsEng,dNvsEng,SigdEvsEng,CSPvsEng,CSigTotvsEng,CdEvsEng,CdNvsEng,CSigdEvsEng
real*8,dimension(atmosLen,nE2strBins) :: prode2stF,prode2stB,Cprode2stF,&
  Cprode2stB
real*8,dimension(nProc,atmosLen,nChS) :: oxygen,Coxygen
real*8,dimension(atmosLen,nChS) :: oxygenCX,CoxygenCX

character(len=100) filename,filenames(nOutputFiles)
character(len=1000) HpHeader,Hp2Header

!****************************** Data Declaration *******************************
data filenames/'H+_Prod','H2+_Prod','H2_Excite_Prod','Oxy_Vs_Energy',&
'Stopping_Power','Processes','Oxy_Neg','Oxy0_','Oxy1_','Oxy2_','Oxy3_','Oxy4_',&
'Oxy5_','Oxy6_','Oxy7_','Oxy8_','Oxy_CX','2Str_Elect_Fwd','2Str_Elect_Bwd'/
!********************************* Initialize **********************************
nTrials=0;CtotHp=0.0;CtotH2p=0.0;CtotalHp=0.0;CtotalH2p=0.0;CH2Ex=0.0
trial=0;COxyVsEng=0.0;CSPvsEng=0.0;CSigTotvsEng=0.0;CdEvsEng=0.0;CSigdEvsEng=0.0
CnSPions=0;Ccollisions=0;Coxygen=0.0;CoxygenCX=0.0;CdNvsEng=0.0
!********** Open output data files for each set of initial energies ************
energy=2000
write(filename,'("./Output/",I0,"keV/Elapsed_Times.dat")') energy
open(unit=100,file=filename,status='old')
do i=1,MaxnTrials
  read(100,*,end=1000) trial(i)
  nTrials=nTrials+1
end do
1000 continue
close(100)
write(*,*) 'Reading in and combining files...'
write(*,*) 'Number of files: ',nTrials,'At an energy of: ',energy,'keV/u.'
do n=1,nTrials
!********************************* Initialize **********************************
  totHp=0.0;totH2p=0.0;totalHp=0.0;totalH2p=0.0;H2Ex=0.0
  OxyVsEng=0.0;SPvsEng=0.0;SigTotvsEng=0.0;dEvsEng=0.0;SigdEvsEng=0.0
  nSPions=0;collisions=0;oxygen=0.0;oxygenCX=0.0;dNvsEng=0.0
!*******************************************************************************
  m=trial(n)
  write(*,*) 'File:',n,'Trial:',m
  do i=1,nOutputFiles
    write(filename,'("./Output/",I0,"keV/",A,I0,".dat")') &
          energy,trim(filenames(i)),trial(n)
    filename=trim(filename)
    open(unit=100+i,file=filename,status='old')
  end do
  do i=1,2
    read(101,*)
    read(102,*)
    read(103,*)
    read(104,*)
    read(105,*)
    read(105,*)
  end do
  read(101,'(1x,A)') HpHeader
  read(102,'(1x,A)') Hp2Header
  read(103,*)
  read(103,*)
  do i=1,atmosLen !Ionization/Excitation vs. altitude
    read(101,F01) altitude(i),(totHp(j,i),j=1,31),totalHp(i)
    read(102,F01) altitude(i),(totH2p(j,i),j=1,11),totalH2p(i)
    read(103,F02) altitude(i),H2Ex(i)
  end do
  CtotHp=CtotHp+totHp
  CtotH2p=CtotH2p+totH2p
  CtotalHp=CtotalHp+totalHp
  CtotalH2p=CtotalH2p+totalH2p
  CH2Ex=CH2Ex+H2Ex
  do i=1,energy/5
    read(104,F03) oxEngBins(i),(OxyVsEng(j,i),j=1,nChS)
  end do
  COxyVsEng=COxyVsEng+OxyVsEng
  do i=1,energy/10
    read(105,F04) stopPowerEbins(i),SPvsEng(i),SigTotvsEng(i),dEvsEng(i),&
                  dNvsEng(i),SigdEvsEng(i),nSPions(i)
  end do
  CSPvsEng=CSPvsEng+SPvsEng
  CSigTotvsEng=CSigTotvsEng+SigTotvsEng
  CdEvsEng=CdEvsEng+dEvsEng
  CdNvsEng=CdNvsEng+dNvsEng
  CSigdEvsEng=CSigdEvsEng+SigdEvsEng
  CnSPions=CnSPions+nSPions
  do i=1,8
    read(106,*) (collisions(i,j),j=1,5)
  end do
  Ccollisions=Ccollisions+collisions
  do i=1,nChS
    read(106+i,*) !Oxygen header
    do j=1,atmosLen
      read(106+i,F01) altitude(j),(oxygen(k,j,i),k=1,nProc)
    end do
  end do
  Coxygen=Coxygen+oxygen
  read(117,*) !Oxygen charge exchange header
  do i=1,atmosLen
    read(117,F05) altitude(i),(oxygenCX(i,j),j=1,nChS)
  end do
  CoxygenCX=CoxygenCX+oxygenCX
  do j=1,nE2strBins
    read(118,F2Str) (prode2stF(i,j),i=atmosLen,1,-1)
    read(119,F2str) (prode2stB(i,j),i=atmosLen,1,-1)
  end do
  Cprode2stF=Cprode2stF+prode2stF
  Cprode2stB=Cprode2stB+prode2stF
  do i=1,nOutputFiles
    close(100+i)
  end do
end do
!*******************************************************************************
write(*,*) 'Writing output files...'
do i=1,nOutputFiles
  write(filename,'("./Output/",I0,"keV/",A,"_Comb.dat")') &
        energy,trim(filenames(i))
  filename=trim(filename)
  open(unit=200+i,file=filename,status='unknown')
end do
write(201,H01) !H+ header
write(201,*) trim(HpHeader)
write(202,H02) !H2+ header
write(202,*) trim(Hp2Header)
write(203,H03) !H2* header
do i=1,atmosLen !Ionization/Excitation vs. altitude
  write(201,F01) altitude(i),(CtotHp(j,i)/real(nTrials),j=1,31),CtotalHp(i)
  write(202,F01) altitude(i),(CtotH2p(j,i)/real(nTrials),j=1,11),CtotalH2p(i)
  write(203,F02) altitude(i),CH2Ex(i)/real(nTrials)
end do
write(204,H04) !Oxy vs energy header
do i=1,nOxEngBins !Oxygen charge state distribution
  write(204,F03) oxEngBins(i),(COxyVsEng(j,i)/real(nTrials),j=1,nChS)
end do
write(205,H05) !Stopping power header
do i=1,nStopPowerEBins !Stopping power vs. ion energy
  write(205,F04) stopPowerEBins(i), &
                 CSPvsEng(i)/real(nTrials), &
                 CSigTotvsEng(i)/real(nTrials), &
                 CdEvsEng(i)/real(nTrials), &
                 CdNvsEng(i)/real(nTrials), &
                 (CSigTotvsEng(i)*CdEvsEng(i))/(real(nTrials)**2), &
                 CnSPions(i)
end do
do i=1,8 !Total number of each type of collision
  write(206,*) (Ccollisions(i,j),j=1,5)
end do
do i=1,nChS !Oxygen production
  write(206+i,*) "Alt [km] ", (HProc(k),k=1,nProc)
  do j=1,atmosLen
    write(206+i,F01) altitude(j),(Coxygen(k,j,i)/real(nTrials),k=1,nProc)
  end do
end do
write(217,H06)
do i=1,atmosLen
  write(217,F05) altitude(i),(CoxygenCX(i,j)/real(nTrials),j=1,nChS)
end do
!***************************** Secondary Electrons *****************************
do i=1,atmosLen
  do j=1,nE2strBins
    Cprode2stF(i,j)=real(Cprode2stF(i,j))/real(nTrials)
    Cprode2stB(i,j)=real(Cprode2stB(i,j))/real(nTrials)
  end do
end do
do j=1,nE2strBins !2-Stream electrons, forward and backward
  write(218,F2Str) (prode2stF(i,j),i=atmosLen,1,-1)
  write(219,F2Str) (prode2stB(i,j),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles
  close(200+i)
end do
end program
