program comb_files
!*******************************************************************************
!* Created by Stephen J. Houston 6.20.18
!*******************************************************************************
!* This program combines multiple files of the same information into one.
!* Used after making a large run on multiple cores (CRC) to simulate as if a
!* single run had taken place.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV
use formatting
implicit real*8(a-h,o-z)

!*******************************************************************************
integer energy,atmosLen,run

parameter(nOutputFiles=19,MaxnTrials=1000,atmosLen=1544,nProc=36,nChS=10)
parameter(nOxEngBins=5000,nStopPowerEBins=295,nE2strBins=260,MaxnLines=100000)
parameter(number_of_energies=12)

integer err,start,nerr !Error upon opening files
integer trial(MaxnTrials),nLines(nOutputFiles) !Number of trials/lines in a file
integer nL(nOutputFiles) !Number of lines there should be in a file
integer(kind=int64) :: collisions(8,5),Ccollisions(8,5),collSUM(8),CcollSUM(8)
integer(kind=int64) :: SIM(5)
integer(kind=int64),dimension(nStopPowerEBins) :: nSPions,CnSPions

real*8 collPerc(8,5),CcollPerc(8,5),collPSUM(8),CcollPSUM(8)
real*8 oxEngBins(nOxEngBins),Eion(number_of_energies)
real*8 IntegratedXRay(2,nChS) !1 is DE, 2 is CX
real*8 altDelta !Width of atmosphere bins (2.0 km)
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
!data Eion/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,10000.0,&
!          25000.0/
data Eion/10.625,15.017,20.225,29.783,46.653,59.770,77.522,120.647,218.125,&
          456.250/ !Juno energy bins from JEDI.
data filenames/'H+_Prod','H2+_Prod','H2_Excite_Prod','Oxy_Vs_Energy',&
'Stopping_Power','Processes','Oxy_Neg','Oxy0_','Oxy1_','Oxy2_','Oxy3_','Oxy4_',&
'Oxy5_','Oxy6_','Oxy7_','Oxy8_','Oxy_CX','2Str_Elect_Fwd','2Str_Elect_Bwd'/
data nL/1547,1547,1548,5002,300,23,1545,1545,1545,1545,1545,1545,1545,1545,&
        1545,1545,1545,40300,40300/
!******************************** Main Program *********************************
do run=1,number_of_energies
!********************************* Initialize **********************************
  nTrials=0;CtotHp=0.0;CtotH2p=0.0;CtotalHp=0.0;CtotalH2p=0.0;CH2Ex=0.0
  trial=0;COxyVsEng=0.0;CSPvsEng=0.0;CSigTotvsEng=0.0;CdEvsEng=0.0
  CnSPions=0;Ccollisions=0;Coxygen=0.0;CoxygenCX=0.0;CdNvsEng=0.0
  CSigdEvsEng=0.0;start=1;nerr=0;CcollSUM=0;CcollPerc=0.0;CcollPSUM=0
!********** Open output data files for each set of initial energies ************
  energy=Eion(run)
  write(filename,'("./Output/Juno/",I0,"keV/Elapsed_Times.dat")') energy
  open(unit=100,file=filename,status='old')
  do i=1,MaxnTrials
    read(100,*,end=1000) trial(i)
    nTrials=nTrials+1
  end do
  1000 continue
  close(100)
  write(*,*) 'Reading in and combining files...'
  write(*,*) 'Number of files: ',nTrials,'At an energy of: ',energy,'keV/u.'
  1002 continue
  do n=1,nTrials
!********************************* Initialize **********************************
    totHp=0.0;totH2p=0.0;totalHp=0.0;totalH2p=0.0;H2Ex=0.0
    OxyVsEng=0.0;SPvsEng=0.0;SigTotvsEng=0.0;dEvsEng=0.0;SigdEvsEng=0.0
    nSPions=0;collisions=0;oxygen=0.0;oxygenCX=0.0;dNvsEng=0.0
    nLines=0;IntegratedXRay=0.0
!*******************************************************************************
    m=trial(n)
    do i=1,nOutputFiles !Open all of the files
      write(filename,'("./Output/Juno/",I0,"keV/",A,I0,".dat")') &
            energy,trim(filenames(i)),trial(n)
      filename=trim(filename)
      open(unit=100+i,file=filename,status='old',iostat=err)
      if(err.gt.0)then
        write(*,*) 'File:',n,'Trial:',m,'ERROR! File:',filename
        start=n+1
        nerr=nerr+1
        goto 1002
      end if
    end do
    write(*,*) 'File:',n,'Trial:',m
    do i=1,nOutputFiles
      do j=1,MaxnLines
        read(100+i,*,end=1001)
        nLines(i)=nLines(i)+1
      end do
    1001 continue
    if(nLines(i).eq.0.0)then
      start=n+1
      goto 1002
    end if
    if(nLines(i).ne.nL(i))write(*,*)'nLines different for unit: ',100+i
    rewind(100+i)
    end do
    do i=1,2
      read(101,*) !Read the headers of the first four files
      read(102,*)
      read(103,*)
      read(103,*)
      read(104,*)
    end do
    read(101,'(1x,A)') HpHeader !Save these headers to output later
    read(102,'(1x,A)') Hp2Header
    do i=1,atmosLen !Ionization/Excitation vs. altitude
      if(nLines(1).eq.nL(1))read(101,F01) altitude(i),&
        (totHp(j,i),j=1,31),totalHp(i)
      if(nLines(2).eq.nL(2))read(102,F01) altitude(i),&
        (totH2p(j,i),j=1,11),totalH2p(i)
      if(nlines(3).eq.nL(3))read(103,F02) altitude(i),H2Ex(i)
    end do
    CtotHp=CtotHp+totHp !Add up all the hydrogen ionizations/excitation
    CtotH2p=CtotH2p+totH2p
    CtotalHp=CtotalHp+totalHp
    CtotalH2p=CtotalH2p+totalH2p
    CH2Ex=CH2Ex+H2Ex
    do i=1,nL(4)-2 !Oxygen charge state distribution
      read(104,F03) oxEngBins(i),(OxyVsEng(j,i),j=1,nChS)
    end do
    COxyVsEng=COxyVsEng+OxyVsEng !Add up all the oxygen vs energy
    do i=1,5 !Stopping power header
      read(105,*)
    end do
    do i=1,nLines(5)-5
      read(105,F04) stopPowerEbins(i),SPvsEng(i),SigTotvsEng(i),dEvsEng(i),&
                    dNvsEng(i),SigdEvsEng(i),nSPions(i)
    end do
    CSPvsEng=CSPvsEng+SPvsEng !Add up all the stopping power variables
    CSigTotvsEng=CSigTotvsEng+SigTotvsEng
    CdEvsEng=CdEvsEng+dEvsEng
    CdNvsEng=CdNvsEng+dNvsEng
    CSigdEvsEng=CSigdEvsEng+SigdEvsEng
    CnSPions=CnSPions+nSPions
    read(106,*) !Processes header
    do i=1,8 !Process counts
      if(nlines(6).eq.nL(6))read(106,F06) dum,(collisions(i,j),j=1,5),collSUM(i)
    end do
    do i=1,4
      read(106,*) !Additional lines to skip
    end do
    do i=1,8 !Processes by percentage
      if(nlines(6).eq.nL(6))read(106,F07) dum,(collPerc(i,j),j=1,5),collPSUM(i)
    end do
    Ccollisions=Ccollisions+collisions !Number of collisions of each type
    CcollSUM=CcollSUM+collSUM !Total number of collisions
    CcollPerc=CcollPerc+collPerc !Percentage of each type of collision
    CcollPSUM=CcollPSUM+collPSUM !Total percentage of each type of collision
    do i=1,nChS
      read(106+i,*) !Oxygen header
      do j=1,atmosLen !oxygen production rates
        if(nlines(6+i).eq.nL(6+i))then
          read(106+i,F01) altitude(j),(oxygen(k,j,i),k=1,nProc)
        end if
      end do
    end do
    Coxygen=Coxygen+oxygen !Adding up oxygen production rates
    read(117,*) !Oxygen charge exchange header
    do i=1,atmosLen !Oxygen production from charge exchange
      if(nlines(17).eq.nL(17))read(117,F05) altitude(i),(oxygenCX(i,j),j=1,nChS)
    end do
    CoxygenCX=CoxygenCX+oxygenCX !Adding up oxygen charge exchange
    do j=1,nE2strBins !Electron production rate for 2-stream
      if(nlines(18).eq.nL(18))read(118,F2Str) (prode2stF(i,j),i=atmosLen,1,-1)
      if(nlines(19).eq.nL(19))read(119,F2str) (prode2stB(i,j),i=atmosLen,1,-1)
    end do
    Cprode2stF=Cprode2stF+prode2stF !Adding up electrons for 2-stream
    Cprode2stB=Cprode2stB+prode2stF
    do i=1,nOutputFiles !Close all of the files
      close(100+i)
    end do
  end do
  nTrials=nTrials-nerr !If there was an error opening a file, subtract it out
  !*******************************************************************************
  write(*,*) 'Writing output files...'
  do i=1,nOutputFiles
    write(filename,'("./Output/Juno/",I0,"keV/",A,"_Comb.dat")') &
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
  SIM=sum(Ccollisions,dim=1)
  write(206,H07) !Collisions header
  do i=1,8 !Total number of each type of collision
    write(206,F06) Coll(i), (Ccollisions(i,j),j=1,5),CcollSUM(i)
  end do
  write(206,*) '------------------------------------------------------------&
                --------------------------'
  write(206,F06) 'Sum     ',SIM(1),SIM(2),SIM(3),SIM(4),SIM(5)
  write(206,*) ''
  write(206,H07) !Collisions percentage header
  do i=1,8 !Total percentage of each type of collision
    write(206,F07) Coll(i), (CcollPerc(i,j)/real(nTrials),j=1,5),&
    CcollPSUM(i)/real(nTrials)
  end do
  write(206,*) '------------------------------------------------------------&
                --------------------------'
  write(206,F07) 'Sum     ',real(SIM(1))/real(sum(SIM))*100,&
   real(SIM(2))/real(sum(SIM))*100,real(SIM(3))/real(sum(SIM))*100,&
   real(SIM(4))/real(sum(SIM))*100,real(SIM(5))/real(sum(SIM))*100
  do i=1,nChS !Oxygen production
    write(206+i,*) "Alt [km] ", (HProc(k),k=1,nProc)
    do j=1,atmosLen
      write(206+i,F01) altitude(j),(Coxygen(k,j,i)/real(nTrials),k=1,nProc)
    end do
  end do
  write(217,H06) !Oxygen charge exchange header
  do i=1,atmosLen !Oxygen production from charge exchange
    write(217,F05) altitude(i),(CoxygenCX(i,j)/real(nTrials),j=1,nChS)
  end do
  !****************************** X-Ray Production *******************************
  write(filename,'("./Output/Juno/",I0,"keV/XRay_DE_Comb.dat")') energy
  open(unit=300,file=trim(filename),status='unknown') !Open X-Ray DE
  write(filename,'("./Output/Juno/",I0,"keV/XRay_CX_Comb.dat")') energy
  open(unit=301,file=trim(filename),status='unknown') !Open X-Ray CX
  altDelta=2.0e5 !2 km = 200,000 cm
  write(300,N01) !X-Ray DE Note
  write(301,N02) !X-Ray CX Note
  write(300,*) !Extra space
  write(301,*) !Extra space
  write(300,H11) !Extra space
  write(301,H11) !Extra space
  write(300,*) !Extra space
  write(301,*) !Extra space
  write(300,H08) !Altitude integrated X-ray production header
  write(301,H08) !Altitude integrated X-ray production header
  write(300,H09) !Charge state header
  write(301,H09) !Charge state header
  do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
    do j=1,nChS
      IntegratedXRay(1,j)=IntegratedXRay(1,j)+&
        real(Coxygen(27,i,j)+Coxygen(29,i,j)+Coxygen(32,i,j))*altDelta !DE
      IntegratedXRay(2,j)=IntegratedXRay(2,j)+&
        real(Coxygen(19,i,j)+Coxygen(25,i,j)+Coxygen(30,i,j))*altDelta !CX
    end do
  end do
  write(300,F05) altDelta/1e5,(IntegratedXRay(1,j)/real(nTrials),j=1,nChS) !DE
  write(301,F05) altDelta/1e5,(IntegratedXRay(2,j)/real(nTrials),j=1,nChS) !CX
  write(300,*) !Extra space
  write(301,*) !Extra space
  write(300,H10) !X-Ray production vs. altitude header
  write(301,H10) !X-Ray production vs. altitude header
  write(300,H06) !Charge state header
  write(301,H06) !Charge state header
  do i=1,atmosLen !DE - TEX+SPEX,SI+SPEX,DI+SPEX, CX - SC+SS,TI,SC
    write(300,F05) altitude(i),& !X-Ray production from direct excitation
      (real(Coxygen(27,i,j)+Coxygen(29,i,j)+Coxygen(32,i,j))/real(nTrials),&
      j=1,nChs)
    write(301,F05) altitude(i),& !X-Ray production from charge exchange
      (real(Coxygen(19,i,j)+Coxygen(25,i,j)+Coxygen(30,i,j))/real(nTrials),&
      j=1,nChs)
  end do
  !***************************** Secondary Electrons *****************************
  do i=1,atmosLen !Normalize electrons for 2-stream code
    do j=1,nE2strBins
      Cprode2stF(i,j)=real(Cprode2stF(i,j))/real(nTrials)
      Cprode2stB(i,j)=real(Cprode2stB(i,j))/real(nTrials)
    end do
  end do
  do j=1,nE2strBins !2-Stream electrons, forward and backward
    write(218,F2Str) (prode2stF(i,j),i=atmosLen,1,-1)
    write(219,F2Str) (prode2stB(i,j),i=atmosLen,1,-1)
  end do
  do i=1,nOutputFiles !Close the combine output files
    close(200+i)
  end do
  close(300)
  close(301)
end do
end program
