program ReadCrossSections
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS1,ChS5,ChS6,ChS7,ChS8,ChS9,ChS10,nEng !Number of charge states and energies.
integer nFiles !Number of cross-section files.
!ChS10 includes negative ion channel

parameter (ChS1=1,ChS5=5,ChS6=6,ChS7=7,ChS8=8,ChS9=9,ChS10=10,nEng=12,nFiles=36)

integer, dimension(ChS10,2):: Gen10 !Generic interpolator offset for ChS10
integer, dimension(ChS10,2):: OffDI,OffSI,OffTEX !Offset for Double Ionization..
integer, dimension(ChS9,2) :: Gen9
integer, dimension(ChS9,2) :: OffDISS,OffSISS,OffTEXSS
integer, dimension(ChS8,2) :: Gen8,OffDISPEX,OffSISPEX,OffTEXSPEX,OffSC,OffTI
integer, dimension(ChS7,2) :: Gen7,OffSIDS,OffSIDPEX,OffDIDS,OffDIDPEX,OffTEXDS
integer, dimension(ChS7,2) :: OffTEXDPEX,OffDC,OffSCSS,OffSCSPEX
integer, dimension(ChS6,2) :: Gen6,OffDCSS,OffDCSPEX,OffSCDS,OffSCDPEX
integer, dimension(ChS5,2) :: Gen5,OffDCDS,OffDCDPEX
integer, dimension(ChS1,2) :: Gen1

real*8, dimension(ChS1,nEng) :: neg !Negative ion channel
real*8, dimension(ChS5,nEng) :: dcds,dcdpex,dcaids,dcaidpex
real*8, dimension(ChS6,nEng) :: tids,tidpex,scds,scdpex,dcss,dcspex,dcaiss,&
                                dcaispex
real*8, dimension(ChS7,nEng) :: tiss,tispex,texds,texdpex,sids,sidpex,scss,&
                                scspex,dids,didpex,dcai,dc
real*8, dimension(ChS8,nEng) :: ti,texspex,sispex,sc,dispex
real*8, dimension(ChS9,nEng) :: texss,siss,diss
real*8, dimension(ChS10,nEng):: tex,si,di

real*8, dimension(ChS1,25000) :: negI !Interpolated
real*8, dimension(ChS5,25000) :: dcdsI,dcdpexI,dcaidsI,dcaidpexI
real*8, dimension(ChS6,25000) :: tidsI,tidpexI,scdsI,scdpexI,dcssI,dcspexI,&
                                 dcaissI,dcaispexI
real*8, dimension(ChS7,25000) :: tissI,tispexI,texdsI,texdpexI,sidsI,sidpexI,&
                                 scssI,scspexI,didsI,didpexI,dcaiI,dcI
real*8, dimension(ChS8,25000) :: tiI,texspexI,sispexI,scI,dispexI
real*8, dimension(ChS9,25000) :: texssI,sissI,dissI
real*8, dimension(ChS10,25000):: texI,siI,diI,sigTot
real*8 EngData(nEng)

character(len=100) filename !Will be open() filename.
character(len=20) XSFile(nFiles) !Cross-section filenames

!****************************** Data Declaration *******************************

data XSFile /'dcds','dcdpex','dcaids','dcaidpex','tids','tidpex','scds',&
             'scdpex','dcss','dcspex','dcaiss','dcaispex','tiss','tispex',&
             'texds','texdpex','sids','sidpex','scss','scspex','dids','didpex',&
             'dcai','dc','ti','texss','texspex','siss','sispex','sc','diss',&
             'dispex','tex','si','di','neg'/
data EngData /1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,5000.0,&
              10000.0,25000/
! CS:       -1    0     1     2     3     4     5     6     7     8
data Gen10 /1,    1,    1,    1,    1,    1,    1,    1,    1,    1,&
            10000,10000,10000,10000,10000,10000,10000,10000,10000,10000/
data OffSI /10, 10, 10, 10, 10, 10, 10, 10, 10, 50,&
            500,500,500,500,500,500,500,500,500,500/
data OffDI /500,500,1,   1,    10, 10,   10,   50,   50,   1,&
            500,500,5000,10000,500,10000,10000,10000,10000,10000/
data OffTEX /50, 50, 75,   100,  1,    1,    1,    10,   10,   10,&
             500,500,10000,10000,10000,10000,10000,10000,10000,10000/
! CS:         -1    0     1   2    3    4     5    6     7
data OffSISS /1,    1,    1,    1,  1,   1,   50,   50,  100,&
              10000,10000,10000,500,1000,5000,10000,1000,500/
data OffDISS /1,    1,    1,    1,    1,   50,  50,   50,  50,&
              10000,10000,10000,10000,1000,1000,10000,1000,500/
data OffTEXSS /1,    1,    1,   500,50,  50,  50,   500,  500,&
               10000,10000,1000,500,1000,1000,10000,10000,10000/
! CS:      0   1   2   3   4   5   6   7
data Gen8 /1,    1,    1,    1,    1,    1,    1,    1,  &
           10000,10000,10000,10000,10000,10000,10000,10000/
data OffSISPEX /1,  1,  1,   100,10, 50,   50, 100,&
                500,500,2000,200,200,10000,200,200/
data OffDISPEX /10,   1,    1,    1,    50, 50, 50, 200,&
                10000,10000,10000,10000,500,500,500,500/
data OffTEXSPEX /1,    1,  500,100,  500,500, 500, 500,&
                 10000,200,500,10000,500,2000,5000,10000/
data OffSC /1,    1,  10,  10,  10, 10, 10, 10,&
            10000,200,1000,2000,500,500,500,1000/
data OffTI /1,  10, 1,  1,  1,  1, 50, 10,&
            100,100,200,200,200,500,10000,1000/
! CS:      0   1   2   3   4   5   6
data Gen7 /1,    1,    1,    1,    1,    1,    1,&
           10000,10000,10000,10000,10000,10000,10000/
data OffSIDS /1,    1,    100,75,   50,  50, 50,&
              10000,10000,500,10000,5000,500,10000/
data OffSIDPEX /100, 1,    1,  100,100, 1,    200,&
                2000,10000,100,500,1000,10000,10000/
data OffDIDS /1,   1,   1,   1,    50,  50,  50,&
              5000,5000,5000,2000,2000,1000,2000/
data OffDIDPEX /10,  1,    1,  10,  100, 100,  200,&
                2000,10000,500,2000,1000,10000,2000/
data OffTEXDS /1,  100, 100,  100,  100,  1000, 1000,&
               500,2000,10000,10000,10000,10000,10000/
data OffTEXDPEX /500,500,500, 500,500,200,  500,&
                 500,500,1000,500,500,10000,10000/
data OffDC /10, 1,  1,  1,  1,   50, 200,&
            200,200,500,500,1000,500,1000/
data OffSCSS /1,    1,    1,    1,    1,    100,  100,&
              10000,10000,10000,10000,10000,10000,10000/
data OffSCSPEX /1,    1,    1,    1,    1,    50,   1,&
                10000,10000,10000,10000,10000,10000,10000/
! CS:      0   1   2   3   4   5
data Gen6 /1,    1,    1,    1,    1,    1,&
           10000,10000,10000,10000,10000,10000/
data OffDCSS /10,  1,   1,   1, 50, 50,&
              1000,2000,1000,75,500,1000/
data OffDCSPEX /1, 1,  1,   1, 50, 75,&
                200,200,1000,75,500,500/
data OffSCDS /1,    10,   1,   1,     50,   100,&
              10000,10000,10000,10000,10000,10000/
data OffSCDPEX /1,    1,    50,   50,   100,  1,&
                10000,10000,10000,10000,10000,10000/
! CS:      0   1   2   3   4
data Gen5 /1,    1,    1,    1,    1,&
           10000,10000,10000,10000,10000/
data OffDCDS /10,1,  50,  200,50,&
              50,500,1000,500,500/
data OffDCDPEX /1,  10,  10, 50, 75,&
                200,2000,500,500,5000/
! CS:      -1
data Gen1 /1,&
           10/

!******************************** Main Program *********************************

do i=1,nFiles
  write(filename,'("./SIMXS/",A,".txt")'), trim(XSfile(i))
  open(unit=99+i,file=filename,status='old') !Open all files
  read(99+i,*) !Skip the first line of the files.
end do

do i=1,ChS1 !Negative ion channel
  read(135,1000) (neg(i,j),j=1,nEng)
end do
!write(*,*) neg
do i=1,ChS5 !Process with 5 charge states
  read(100,1000) (dcds(i,j),j=1,nEng) !DC+DS, O2+ - O6+
  read(101,1000) (dcdpex(i,j),j=1,nEng)
  read(102,1000) (dcaids(i,j),j=1,nEng)
  read(103,1000) (dcaidpex(i,j),j=1,nEng)
end do
do i=1,ChS6 !Process with 6 charge states
  read(104,1000) (tids(i,j),j=1,nEng)
  read(105,1000) (tidpex(i,j),j=1,nEng)
  read(106,1000) (scds(i,j),j=1,nEng)
  read(107,1000) (scdpex(i,j),j=1,nEng)
  read(108,1000) (dcss(i,j),j=1,nEng)
  read(109,1000) (dcspex(i,j),j=1,nEng)
  read(110,1000) (dcaiss(i,j),j=1,nEng)
  read(111,1000) (dcaispex(i,j),j=1,nEng)
end do
do i=1,ChS7 !Process with 7 charge states
  read(112,1000) (tiss(i,j),j=1,nEng)
  read(113,1000) (tispex(i,j),j=1,nEng)
  read(114,1000) (texds(i,j),j=1,nEng)
  read(115,1000) (texdpex(i,j),j=1,nEng)
  read(116,1000) (sids(i,j),j=1,nEng)
  read(117,1000) (sidpex(i,j),j=1,nEng)
  read(118,1000) (scss(i,j),j=1,nEng)
  read(119,1000) (scspex(i,j),j=1,nEng)
  read(120,1000) (dids(i,j),j=1,nEng)
  read(121,1000) (didpex(i,j),j=1,nEng)
  read(122,1000) (dcai(i,j),j=1,nEng)
  read(123,1000) (dc(i,j),j=1,nEng)
end do
do i=1,ChS8 !Process with 8 charge states
  read(124,1000) (ti(i,j),j=1,nEng)
  read(126,1000) (texspex(i,j),j=1,nEng)
  read(128,1000) (sispex(i,j),j=1,nEng)
  read(129,1000) (sc(i,j),j=1,nEng)
  read(131,1000) (dispex(i,j),j=1,nEng)
end do
do i=1,ChS9 !Process with 9 charge states
  read(125,1000) (texss(i,j),j=1,nEng)
  read(127,1000) (siss(i,j),j=1,nEng)
  read(130,1000) (diss(i,j),j=1,nEng)
end do
do i=1,ChS10 !Process with 9 charge states
  read(132,1000) (tex(i,j),j=1,nEng)
  read(133,1000) (si(i,j),j=1,nEng)
  read(134,1000) (di(i,j),j=1,nEng)
end do
do i=1,nFiles
  close(99+i)
end do
call InterpXS(neg,negI,Gen1,ChS1)
call InterpXS(dcds,dcdsI,OffDCDS,ChS5)
call InterpXS(dcdpex,dcdpexI,OffDCDPEX,ChS5)
call InterpXS(dcaids,dcaidsI,Gen5,ChS5)
call InterpXS(dcaidpex,dcaidpexI,Gen5,ChS5)
call InterpXS(tids,tidsI,Gen6,ChS6)
call InterpXS(tidpex,tidpexI,Gen6,ChS6)
call InterpXS(scds,scdsI,OffSCDS,ChS6)
call InterpXS(scdpex,scdpexI,OffSCDPEX,ChS6)
call InterpXS(dcss,dcssI,OffDCSS,ChS6)
call InterpXS(dcspex,dcspexI,OffDCSPEX,ChS6)
call InterpXS(dcaiss,dcaissI,Gen6,ChS6)
call InterpXS(dcaispex,dcaispexI,Gen6,ChS6)
call InterpXS(tiss,tissI,Gen7,ChS7)
call InterpXS(tispex,tispexI,Gen7,ChS7)
call InterpXS(texds,texdsI,OffTEXDS,ChS7)
call InterpXS(texdpex,texdpexI,OffTEXDPEX,ChS7)
call InterpXS(sids,sidsI,OffSIDS,ChS7)
call InterpXS(sidpex,sidpexI,OffSIDPEX,ChS7)
call InterpXS(scss,scssI,OffSCSS,ChS7)
call InterpXS(scspex,scspexI,OffSCSPEX,ChS7)
call InterpXS(dids,didsI,OffDIDS,ChS7)
call InterpXS(didpex,didpexI,OffDIDPEX,ChS7)
call InterpXS(dcai,dcaiI,Gen7,ChS7)
call InterpXS(dc,dcI,OffDC,ChS7)
call InterpXS(ti,tiI,OffTI,ChS8)
call InterpXS(texspex,texspexI,OffTEXSPEX,ChS8)
call InterpXS(sispex,sispexI,OffSISPEX,ChS8)
call InterpXS(sc,scI,OffSC,ChS8)
call InterpXS(dispex,dispexI,OffDISPEX,ChS8)
call InterpXS(texss,texssI,OffTEXSS,ChS9)
call InterpXS(siss,sissI,OffSISS,ChS9)
call InterpXS(diss,dissI,OffDISS,ChS9)
call InterpXS(tex,texI,OffTEX,ChS10)
call InterpXS(si,siI,OffSI,ChS10)
call InterpXS(di,diI,OffDI,ChS10)
do i=1,nFiles!24!nFiles
  write(filename,'("./SIMXSInterp/",A,".dat")'), trim(XSfile(i))
  open(unit=199+i,file=filename,status='unknown') !Open all files
  write(filename,'("./SIMXSPoints/",A,"p.dat")'), trim(XSfile(i)) !Write data points out
  open(unit=299+i,file=filename,status='unknown') !Open all files
end do
do i=1,25000 !CAN WRITE OUT I AS ENERGY VECTOR BECAUSE EVECTOR IS IN STEPS OF 1!
  write(200,3005) real(i), (dcdsI(j,i),j=1,ChS5)
  write(201,3005) real(i), (dcdpexI(j,i),j=1,ChS5)
  write(202,3005) real(i), (dcaidsI(j,i),j=1,ChS5)
  write(203,3005) real(i), (dcaidpexI(j,i),j=1,ChS5)
  write(204,3006) real(i), (tidsI(j,i),j=1,ChS6)
  write(205,3006) real(i), (tidpexI(j,i),j=1,ChS6)
  write(206,3006) real(i), (scdsI(j,i),j=1,ChS6)
  write(207,3006) real(i), (scdpexI(j,i),j=1,ChS6)
  write(208,3006) real(i), (dcssI(j,i),j=1,ChS6)
  write(209,3006) real(i), (dcspexI(j,i),j=1,ChS6)
  write(210,3006) real(i), (dcaissI(j,i),j=1,ChS6)
  write(211,3006) real(i), (dcaispexI(j,i),j=1,ChS6)
  write(212,3007) real(i), (tissI(j,i),j=1,ChS7)
  write(213,3007) real(i), (tispexI(j,i),j=1,ChS7)
  write(214,3007) real(i), (texdsI(j,i),j=1,ChS7)
  write(215,3007) real(i), (texdpexI(j,i),j=1,ChS7)
  write(216,3007) real(i), (sidsI(j,i),j=1,ChS7)
  write(217,3007) real(i), (sidpexI(j,i),j=1,ChS7)
  write(218,3007) real(i), (scssI(j,i),j=1,ChS7)
  write(219,3007) real(i), (scspexI(j,i),j=1,ChS7)
  write(220,3007) real(i), (didsI(j,i),j=1,ChS7)
  write(221,3007) real(i), (didpexI(j,i),j=1,ChS7)
  write(222,3007) real(i), (dcaiI(j,i),j=1,ChS7)
  write(223,3007) real(i), (dcI(j,i),j=1,ChS7)
  write(224,3008) real(i), (tiI(j,i),j=1,ChS8)
  write(226,3008) real(i), (texspexI(j,i),j=1,ChS8)
  write(228,3008) real(i), (sispexI(j,i),j=1,ChS8)
  write(229,3008) real(i), (scI(j,i),j=1,ChS8)
  write(231,3008) real(i), (dispexI(j,i),j=1,ChS8)
  write(225,3009) real(i), (texssI(j,i),j=1,ChS9)
  write(227,3009) real(i), (sissI(j,i),j=1,ChS9)
  write(230,3009) real(i), (dissI(j,i),j=1,ChS9)
  write(232,3010) real(i), (texI(j,i),j=1,ChS10)
  write(233,3010) real(i), (siI(j,i),j=1,ChS10)
  write(234,3010) real(i), (diI(j,i),j=1,ChS10)
  write(235,3010) real(i), (negI(j,i),j=1,ChS1)
end do
do i=1,nFiles
  close(199+i)
end do
!goto 5000
sigTot=1.0e-44
do i=1,ChS10
  do j=1,25000
    if(i.eq.1) then !O-
      sigTot(i,j)=siI(i,j)+sissI(i,j)+diI(i,j)+dissI(i,j)+texI(i,j)+texssI(i,j)
    else if(i.eq.2)then !O
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)+&
        negI(i-1,j)
    else if(i.eq.3)then !O+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.4)then !O++
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaidsI(i-3,j)+dcaispexI(i-3,j)+dcaidpexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcdsI(i-3,j)+dcspexI(i-3,j)+dcdpexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.5)then !O3+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaidsI(i-3,j)+dcaispexI(i-3,j)+dcaidpexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcdsI(i-3,j)+dcspexI(i-3,j)+dcdpexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.6)then !O4+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaidsI(i-3,j)+dcaispexI(i-3,j)+dcaidpexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcdsI(i-3,j)+dcspexI(i-3,j)+dcdpexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.7)then !O5+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaidsI(i-3,j)+dcaispexI(i-3,j)+dcaidpexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcdsI(i-3,j)+dcspexI(i-3,j)+dcdpexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.8)then !O6+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sidsI(i-1,j)+sispexI(i-1,j)+sidpexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+didsI(i-1,j)+dispexI(i-1,j)+didpexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tidsI(i-2,j)+tispexI(i-2,j)+tidpexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaidsI(i-3,j)+dcaispexI(i-3,j)+dcaidpexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scdsI(i-2,j)+scspexI(i-2,j)+scdpexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcdsI(i-3,j)+dcspexI(i-3,j)+dcdpexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texdsI(i-1,j)+texspexI(i-1,j)+texdpexI(i-1,j)
    else if(i.eq.9)then !O7+
      sigTot(i,j)=siI(i,j)+sissI(i,j)+sispexI(i-1,j)+&
        diI(i,j)+dissI(i,j)+dispexI(i-1,j)+&
        tiI(i-2,j)+tissI(i-2,j)+tispexI(i-2,j)+&
        dcaiI(i-3,j)+dcaissI(i-3,j)+dcaispexI(i-3,j)+&
        scI(i-2,j)+scssI(i-2,j)+scspexI(i-2,j)+&
        dcI(i-3,j)+dcssI(i-3,j)+dcspexI(i-3,j)+&
        texI(i,j)+texssI(i,j)+texspexI(i-1,j)
    else if(i.eq.10)then !O8+
      sigTot(i,j)=siI(i,j)+diI(i,j)+tiI(i-2,j)+dcaiI(i-3,j)+scI(i-2,j)+dcI(i-3,j)+&
                  texI(i,j)
    end if
  end do
end do
!5000 continue
open(unit=300,file='./SIMXSInterp/TotalXS.dat',status='unknown')
!write(*,*) siI(9,25000),diI(9,25000),tiI(8,25000),dcaiI(7,25000),scI(8,25000),dcI(7,25000),&
!            texI(9,25000)
!write(*,*) siI(9,25000)+diI(9,25000)+tiI(8,25000)+dcaiI(7,25000)+scI(8,25000)+dcI(7,25000)+&
!            texI(9,25000)
do i=1,25000
  write(300,3011) i, (sigTot(j,i),j=1,ChS10)
end do
close(300)
do i=1,12
  write(300,3005) EngData(i), (dcds(j,i),j=1,Chs5)
  write(301,3005) EngData(i), (dcdpex(j,i),j=1,ChS5)
  write(302,3005) EngData(i), (dcaids(j,i),j=1,ChS5)
  write(303,3005) EngData(i), (dcaidpex(j,i),j=1,ChS5)
  write(304,3006) EngData(i), (tids(j,i),j=1,ChS6)
  write(305,3006) EngData(i), (tidpex(j,i),j=1,ChS6)
  write(306,3006) EngData(i), (scds(j,i),j=1,ChS6)
  write(307,3006) EngData(i), (scdpex(j,i),j=1,ChS6)
  write(308,3006) EngData(i), (dcss(j,i),j=1,ChS6)
  write(309,3006) EngData(i), (dcspex(j,i),j=1,ChS6)
  write(310,3006) EngData(i), (dcaiss(j,i),j=1,ChS6)
  write(311,3006) EngData(i), (dcaispex(j,i),j=1,ChS6)
  write(312,3007) EngData(i), (tiss(j,i),j=1,ChS7)
  write(313,3007) EngData(i), (tispex(j,i),j=1,ChS7)
  write(314,3007) EngData(i), (texds(j,i),j=1,ChS7)
  write(315,3007) EngData(i), (texdpex(j,i),j=1,ChS7)
  write(316,3007) EngData(i), (sids(j,i),j=1,ChS7)
  write(317,3007) EngData(i), (sidpex(j,i),j=1,ChS7)
  write(318,3007) EngData(i), (scss(j,i),j=1,ChS7)
  write(319,3007) EngData(i), (scspex(j,i),j=1,ChS7)
  write(320,3007) EngData(i), (dids(j,i),j=1,ChS7)
  write(321,3007) EngData(i), (didpex(j,i),j=1,ChS7)
  write(322,3007) EngData(i), (dcai(j,i),j=1,ChS7)
  write(323,3007) EngData(i), (dc(j,i),j=1,ChS7)
  write(324,3008) EngData(i), (ti(j,i),j=1,ChS8)
  write(325,3009) EngData(i), (texss(j,i),j=1,ChS9)
  write(326,3008) EngData(i), (texspex(j,i),j=1,ChS8)
  write(327,3009) EngData(i), (siss(j,i),j=1,ChS9)
  write(328,3008) EngData(i), (sispex(j,i),j=1,ChS8)
  write(329,3008) EngData(i), (sc(j,i),j=1,ChS8)
  write(330,3009) EngData(i), (diss(j,i),j=1,ChS9)
  write(331,3008) EngData(i), (dispex(j,i),j=1,ChS8)
  write(332,3010) EngData(i), (tex(j,i),j=1,ChS10)
  write(333,3010) EngData(i), (si(j,i),j=1,ChS10)
  write(334,3010) EngData(i), (di(j,i),j=1,ChS10)
  write(335,3001) EngData(i), (neg(j,i),j=1,ChS1)
end do
do i=1,nFiles
  close(299+i)
end do
1000 format (9x,12ES9.7)
2000 format (12(ES13.7,2x))
3001 format (F8.1,3x,ES11.5E2)
3005 format (F8.1,3x,5(ES11.5E2,2x))
3006 format (F8.1,3x,6(ES11.5E2,2x))
3007 format (F8.1,3x,7(ES11.5E2,2x))
3008 format (F8.1,3x,8(ES11.5E2,2x))
3009 format (F8.1,3x,9(ES11.5E2,2x))
3010 format (F8.1,3x,10(ES11.5E2,2x))
3011 format (I5,3x,10(ES11.5E2,2x))
end program
