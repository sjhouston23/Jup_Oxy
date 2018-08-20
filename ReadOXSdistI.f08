program ReadOXSDistI
!*******************************************************************************
!* Created by Stephen J. Houston 7.12.18
!*******************************************************************************
!* This program reads in all of the interpolated integral cross-section data
!* for oxygen and adds the cross-sections for charge exchange and electron
!* removal processes.
!* CX1 - Charge exchange  resulting in O^(q-1) (Charge states  0 - 8) - ChS9
!* CX2 - Charge exchange  resulting in O^(q-2) (Charge states  2 - 8) - ChS7
!* ER1 - Electron removal resulting in O^(q+1) (Charge states -1 - 7) - ChS9
!* ER2 - Electron removal resulting in O^(q+2) (Charge states -1 - 6) - ChS8
!* Note:
!*   The negative ion channel can only be accesed from O neutral to O-;
!*   therefore, O- cannot be reached for CX2, giving CX2 only 7 charge state
!*   possibilities
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS1,ChS5,ChS6,ChS7,ChS8,ChS9,ChS10 !Number of charge states
integer nEng !Number of energies.
integer nFiles !Number of cross-section files.
!ChS10 includes negative ion channel

parameter (ChS1=1,ChS5=5,ChS6=6,ChS7=7,ChS8=8,ChS9=9,ChS10=10,nEng=25000,nFiles=36)

real*8,dimension(ChS1,nEng) :: neg !Negative ion channel
real*8,dimension(ChS5,nEng) :: dcds,dcdpex,dcaids,dcaidpex
real*8,dimension(ChS6,nEng) :: tids,tidpex,scds,scdpex,dcss,dcspex,dcaiss,&
                               dcaispex
real*8,dimension(ChS7,nEng) :: tiss,tispex,texds,texdpex,sids,sidpex,scss,&
                               scspex,dids,didpex,dcai,dc,ds,CX2
real*8,dimension(ChS8,nEng) :: ti,texspex,sispex,sc,dispex,ss,ER2
real*8,dimension(ChS9,nEng) :: texss,siss,diss,CX1,ER1
real*8,dimension(ChS10,nEng):: tex,si,di

character(len=100) filename !Will be open() filename.
character(len=20) XSFile(nFiles) !Cross-section filenames

!****************************** Data Declaration *******************************

data XSFile /'dcds','dcdpex','dcaids','dcaidpex','tids','tidpex','scds',&
             'scdpex','dcss','dcspex','dcaiss','dcaispex','tiss','tispex',&
             'texds','texdpex','sids','sidpex','scss','scspex','dids','didpex',&
             'dcai','dc','ti','texss','texspex','siss','sispex','sc','diss',&
             'dispex','tex','si','di','neg'/

!******************************** Main Program *********************************

do i=1,nFiles
  write(filename,'("./SIMXSInterp/",A,".dat")'), trim(XSfile(i))
  open(unit=99+i,file=filename,status='old') !Open all files
!  read(99+i,*) !Skip the first line of the files.
end do

do j=1,nEng !Negative ion channel
  read(135,1001) e, (neg(i,j),i=1,ChS1)
!  write(*,*) e,j
!end do
!do i=1,ChS5 !Process with 5 charge states
  read(100,1005) e, (dcds(i,j),i=1,ChS5) !DC+DS, O2+ - O6+
  read(101,1005) e, (dcdpex(i,j),i=1,ChS5)
  read(102,1005) e, (dcaids(i,j),i=1,ChS5)
  read(103,1005) e, (dcaidpex(i,j),i=1,ChS5)
!end do
!do i=1,ChS6 !Process with 6 charge states
  read(104,1006) e, (tids(i,j),i=1,ChS6)
  read(105,1006) e, (tidpex(i,j),i=1,ChS6)
  read(106,1006) e, (scds(i,j),i=1,ChS6)
  read(107,1006) e, (scdpex(i,j),i=1,ChS6)
  read(108,1006) e, (dcss(i,j),i=1,ChS6)
  read(109,1006) e, (dcspex(i,j),i=1,ChS6)
  read(110,1006) e, (dcaiss(i,j),i=1,ChS6)
  read(111,1006) e, (dcaispex(i,j),i=1,ChS6)
!end do
!do i=1,ChS7 !Process with 7 charge states
  read(112,1007) e, (tiss(i,j),i=1,ChS7)
  read(113,1007) e, (tispex(i,j),i=1,ChS7)
  read(114,1007) e, (texds(i,j),i=1,ChS7)
  read(115,1007) e, (texdpex(i,j),i=1,ChS7)
  read(116,1007) e, (sids(i,j),i=1,ChS7)
  read(117,1007) e, (sidpex(i,j),i=1,ChS7)
  read(118,1007) e, (scss(i,j),i=1,ChS7)
  read(119,1007) e, (scspex(i,j),i=1,ChS7)
  read(120,1007) e, (dids(i,j),i=1,ChS7)
  read(121,1007) e, (didpex(i,j),i=1,ChS7)
  read(122,1007) e, (dcai(i,j),i=1,ChS7)
  read(123,1007) e, (dc(i,j),i=1,ChS7)
!  read(136,1007) (ds(i,j),i=1,ChS7)
!end do
!do i=1,ChS8 !Process with 8 charge states
!  write(*,*) j
  read(124,1008) e, (ti(i,j),i=1,ChS8)
  read(126,1008) e, (texspex(i,j),i=1,ChS8)
  read(128,1008) e, (sispex(i,j),i=1,ChS8)
  read(129,1008) e, (sc(i,j),i=1,ChS8)
  read(131,1008) e, (dispex(i,j),i=1,ChS8)
!  read(135,1008) (ss(i,j),i=1,ChS8)
!end do
!do i=1,ChS9 !Process with 9 charge states
  read(125,1009) e, (texss(i,j),i=1,ChS9)
  read(127,1009) e, (siss(i,j),i=1,ChS9)
  read(130,1009) e, (diss(i,j),i=1,ChS9)
!end do
!do i=1,ChS10 !Process with 9 charge states
  read(132,1010) e, (tex(i,j),i=1,ChS10)
  read(133,1010) e, (si(i,j),i=1,ChS10)
  read(134,1010) e, (di(i,j),i=1,ChS10)
end do
do i=1,nFiles
  close(99+i)
end do
!******************************** CX processes *********************************
CX1=0.0;CX2=0.0;ER1=0.0;ER2=0.0
do i=1,ChS9 !Loop through every charge state, set XS=0 if CS doesn't exist
  do j=1,nEng !Loop through every energy
    if(i.le.5)then !Charge states  2 - 6
      CX2(i,j)=dc(i,j)+dcspex(i,j)+dcdpex(i,j) !CX2 O2+ thru O6+
    elseif(i.eq.6)then
      CX2(i,j)=dc(i,j)+dcspex(i,j) !CX2 O7+
    elseif(i.eq.7)then
      CX2(i,j)=dc(i,j) !CX2 O8+
    end if
    if(i.le.7)then !Charge states 0 - 6
      ER2(i,j)=sids(i,j)+dids(i,j)+texds(i,j) !ER2 O thru O6+
    end if
    select case(i)
    case(1) !O -> O-,O2+ -> O,O- -> O,O -> O2+
      CX1(i,j)=neg(i,j) !CX1 O neutral
      ER1(i,j)=siss(i,j)+diss(i,j)+texss(i,j) !ER1 O-
    case(2) !O+ -> O
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j) !CX1 O+
      ER1(i,j)=siss(i,j)+diss(i,j)+texss(i,j) !ER1 O
    case(3)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j)+dcai(i-2,j)+dcaispex(i-2,j)+dcaidpex(i-2,j)&
              +dcss(i-2,j) !CX1 O2+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+texss(i,j)+scds(i-2,j) !ER1 O+
    case(4)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j)+dcai(i-2,j)+dcaispex(i-2,j)+dcaidpex(i-2,j)&
              +dcss(i-2,j) !CX1 O3+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+dcaids(i-3,j)+texss(i,j)&
              +scds(i-2,j) !ER1 O2+
    case(5)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j)+dcai(i-2,j)+dcaispex(i-2,j)+dcaidpex(i-2,j)&
              +dcss(i-2,j) !CX1 O4+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+dcaids(i-3,j)+texss(i,j)&
              +scds(i-2,j) !ER1 O3+
    case(6)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j)+dcai(i-2,j)+dcaispex(i-2,j)+dcaidpex(i-2,j)&
              +dcss(i-2,j) !CX1 O5+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+dcaids(i-3,j)+texss(i,j)&
              +scds(i-2,j) !ER1 O4+
    case(7)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+tidpex(i-1,j)+sc(i-1,j)+scspex(i-1,j)&
              +scdpex(i-1,j)+dcai(i-2,j)+dcaispex(i-2,j)+dcaidpex(i-2,j)&
              +dcss(i-2,j) !CX1 O6+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+dcaids(i-3,j)+texss(i,j)&
              +scds(i-2,j) !ER1 O5+
    case(8)
      CX1(i,j)=ti(i-1,j)+tispex(i-1,j)+sc(i-1,j)+scspex(i-1,j)+dcai(i-2,j)&
              +dcaispex(i-2,j)+dcss(i-2,j) !CX1 O7+
      ER1(i,j)=siss(i,j)+tids(i-2,j)+diss(i,j)+dcaids(i-3,j)+texss(i,j)&
              +scds(i-2,j) !ER1 O6+
    case(9)
      CX1(i,j)=ti(i-1,j)+sc(i-1,j)+dcai(i-2,j) !CX1 O8+
      ER1(i,j)=siss(i,j)+diss(i,j)+texss(i,j) !ER1 O7+
    end select
  end do
end do
open(unit=200,file='./SIMXSInterp/CX1.dat')
open(unit=201,file='./SIMXSInterp/CX2.dat')
open(unit=202,file='./SIMXSInterp/ER1.dat')
open(unit=203,file='./SIMXSInterp/ER2.dat')
do j=1,nEng
  write(200,1009) real(j), (CX1(i,j),i=1,ChS9)
  write(201,1009) real(j), (CX2(i,j),i=1,ChS7)
  write(202,1009) real(j), (ER1(i,j),i=1,ChS9)
  write(203,1009) real(j), (ER2(i,j),i=1,ChS8)
end do
close(200)
close(201)
close(202)
close(203)
1001 format (F8.1,3x,ES11.5E2)
1005 format (F8.1,3x,5(ES11.5E2,2x))
1006 format (F8.1,3x,6(ES11.5E2,2x))
1007 format (F8.1,3x,7(ES11.5E2,2x))
1008 format (F8.1,3x,8(ES11.5E2,2x))
1009 format (F8.1,3x,9(ES11.5E2,2x))
1010 format (F8.1,3x,10(ES11.5E2,2x))
end program
