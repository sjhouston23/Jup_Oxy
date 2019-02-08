program AddXRayProcesses
!*******************************************************************************
!* Created by Stephen J Houston 10.17.18
!*******************************************************************************
!* This program is designed to read in oxygen integral cross-section data
!* that produces photons and add them together to look at the total cross-
!* section from those processes.
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS,Eng

parameter(nChS=9,nEng=12,nProc=3,nIter=25000)

real*8,dimension(nChS,nIter) :: TotalXS
real*8,dimension(nProc,nChS,nIter) :: ProcessXS

character(len=100) :: files(nProc),filename

!****************************** Data Declaration *******************************

! data files/'siss','sids','tids','diss','dids','dcaids','scds','texss','texds'/
data files/'sispex','dispex','texspex'/

!******************************** Main Program *********************************
ProcessXS=0.0;TotalXS=0.0
do i=1,nProc
  write(filename,'("./SIMXSInterp/",A,".dat")') trim(files(i))
  open(unit=100+i,file=trim(filename),status='old')
end do
do i=1,nIter
  read(101,1008) dum,(ProcessXS(1,j,i),j=2,9) !SI+SS
  read(102,1008) dum,(ProcessXS(2,j,i),j=2,9) !SI+DS
  read(103,1008) dum,(ProcessXS(3,j,i),j=2,9) !TI+DS
  ! read(104,1009) dum,(ProcessXS(4,j,i),j=1,9) !DI+SS
  ! read(105,1007) dum,(ProcessXS(5,j,i),j=2,8) !DI+DS
  ! read(106,1005) dum,(ProcessXS(6,j,i),j=4,8) !DCAI+DS
  ! read(107,1006) dum,(ProcessXS(7,j,i),j=3,8) !SC+DS
  ! read(108,1009) dum,(ProcessXS(8,j,i),j=1,9) !TEX+SS
  ! read(109,1007) dum,(ProcessXS(9,j,i),j=2,8) !TEX+DS
end do
do i=1,nProc
  close(100+i)
end do
TotalXS=sum(ProcessXS,dim=1)
open(unit=200,file='./SIMXSInterp/DEProcesses.dat')
do i=1,nIter
  write(200,1008) real(i),(TotalXS(j,i),j=2,9)
end do
close(200)

1005 format (F8.1,3x,5(ES11.5E2,2x))
1006 format (F8.1,3x,6(ES11.5E2,2x))
1007 format (F8.1,3x,7(ES11.5E2,2x))
1008 format (F8.1,3x,8(ES11.5E2,2x))
1009 format (F8.1,3x,9(ES11.5E2,2x))


end program
