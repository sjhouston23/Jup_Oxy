program ProductionRates
!*******************************************************************************
!* Created by Stephen J. Houston 2.4.19
!*******************************************************************************
!* This program loops through all of the output files and combines the
!* production rates of various energies into output files by charge state for
!* easier plotting.
!*******************************************************************************

use formatting !Formatting module to avoid cluttering the end of the program
implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer Eng,ChS,Alt,AtmosLen

parameter(nFiles=10,nChS=10,AtmosLen=1544,nEnergies=27)

integer,dimension(nEnergies) :: Eion

real*8 Altitude(AtmosLen)
real*8,dimension(AtmosLen,nEnergies,nChS,2) :: Production

character(len=100) filename,files(2)
!****************************** Data Declaration *******************************
! data Eion/10,50,75,100,200,300,500,1000,2000,5000,10000,25000/
data Eion/10,25,50,75,100,125,150,175,200,250,300,350,400,450,500,600,700,800,&
         900,1000,1250,1500,1750,2000,5000,10000,25000/
!******************************** Main Program *********************************
do Eng=1,nEnergies
  write(filename,'("./Output/",I0,"keV/XRay_CX_Comb.dat")') Eion(Eng)
  open(unit=100,file=trim(filename),status='unknown') !Open X-Ray DE
  write(filename,'("./Output/",I0,"keV/XRay_DE_Comb.dat")') Eion(Eng)
  open(unit=101,file=trim(filename),status='unknown') !Open X-Ray CX
  do i=1,13
    read(100,*) !Skip the header
    read(101,*) !Skip the header
  end do
  read(100,*) !Has an additional header line
  do i=1,AtmosLen
    read(100,*) Altitude(i),(Production(i,Eng,ChS,1),ChS=1,nChS) !CX production
    read(101,*) Altitude(i),(Production(i,Eng,ChS,2),ChS=1,nChS) !DE production
  end do
  open(unit=102,file="./Output/TotalProductions/TotLatex.dat",access='append')
  write(102,1003) Eion(eng),&
    (sum(Production(:,Eng,8,1))+sum(Production(:,Eng,8,2)))*2e5,&
    (sum(Production(:,Eng,9,1))+sum(Production(:,Eng,9,2)))*2e5
1003 format(1x,I5,2(' & ',ES8.2),' \\ \hline')
  close(100)
  close(101)
  close(102)
end do
stop
do ChS=1,nChS
  write(filename,'("./Output/TotalProductions/XRay_CX_",I0,".dat")') ChS-2
  open(unit=100,file=trim(filename),status='unknown') !Open X-Ray DE
  write(filename,'("./Output/TotalProductions/XRay_DE_",I0,".dat")') ChS-2
  open(unit=101,file=trim(filename),status='unknown') !Open X-Ray CX
  write(100,1000) ' ΔAlt [km]', (Eion(Eng),Eng=1,nEnergies) !CX
  write(101,1000) ' ΔAlt [km]', (Eion(Eng),Eng=1,nEnergies) !DE
  write(100,1001) 2.0,(sum(Production(:,Eng,ChS,1))*2e5,Eng=1,nEnergies) !CX
  write(101,1001) 2.0,(sum(Production(:,Eng,ChS,2))*2e5,Eng=1,nEnergies) !DE
  write(100,*)
  write(101,*)
  write(100,1002) ' Alt [km] ', (Eion(Eng),Eng=1,nEnergies) !CX
  write(101,1002) ' Alt [km] ', (Eion(Eng),Eng=1,nEnergies) !DE
  do Alt=1,AtmosLen
    write(100,1001) Altitude(Alt),(Production(Alt,Eng,ChS,1),Eng=1,nEnergies)
    write(101,1001) Altitude(Alt),(Production(Alt,Eng,ChS,2),Eng=1,nEnergies)
  end do
  close(100)
  close(101)
end do


1000 format(A11,12(3x,I8))
1001 format(F10.2,12(1x,'&',1x,ES8.2E2))
1002 format(A10,12(2x,I8))



end program
