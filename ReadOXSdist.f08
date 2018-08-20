subroutine ReadOXSDist(OXSDist)
!*******************************************************************************
!* Created by Stephen J. Houston 4.25.18
!*******************************************************************************
!* This subroutine reads in all of the integral cross-section data for oxygen
!* and puts it into a single variable.
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS1,ChS5,ChS6,ChS7,ChS8,ChS9,ChS10 !Number of charge states
integer nEng !Number of energies.
integer nFiles !Number of cross-section files.
!ChS10 includes negative ion channel

parameter (ChS1=1,ChS5=5,ChS6=6,ChS7=7,ChS8=8,ChS9=9,ChS10=10,nEng=12,nFiles=38)

real*8,dimension(ChS1,nEng) :: neg !Negative ion channel
real*8,dimension(ChS5,nEng) :: dcds,dcdpex,dcaids,dcaidpex
real*8,dimension(ChS6,nEng) :: tids,tidpex,scds,scdpex,dcss,dcspex,dcaiss,&
                               dcaispex
real*8,dimension(ChS7,nEng) :: tiss,tispex,texds,texdpex,sids,sidpex,scss,&
                               scspex,dids,didpex,dcai,dc,ds
real*8,dimension(ChS8,nEng) :: ti,texspex,sispex,sc,dispex,ss
real*8,dimension(ChS9,nEng) :: texss,siss,diss
real*8,dimension(ChS10,nEng):: tex,si,di

real*8 OXSDist(nFiles,ChS10,nEng) !No. of processes, charge states, energies

character(len=100) filename !Will be open() filename.
character(len=20) XSFile(nFiles) !Cross-section filenames

!****************************** Data Declaration *******************************

data XSFile /'dcds','dcdpex','dcaids','dcaidpex','tids','tidpex','scds',&
             'scdpex','dcss','dcspex','dcaiss','dcaispex','tiss','tispex',&
             'texds','texdpex','sids','sidpex','scss','scspex','dids','didpex',&
             'dcai','dc','ti','texss','texspex','siss','sispex','sc','diss',&
             'dispex','tex','si','di','ss','ds','neg'/

!******************************** Main Program *********************************

do i=1,nFiles
  write(filename,'("./SIMXS/",A,".txt")'), trim(XSfile(i))
  open(unit=99+i,file=filename,status='old') !Open all files
  read(99+i,*) !Skip the first line of the files.
end do

do i=1,ChS1 !Negative ion channel
  read(137,1000) (neg(i,j),j=1,nEng)
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
  read(136,1000) (ds(i,j),j=1,nEng)
end do
do i=1,ChS8 !Process with 8 charge states
  read(124,1000) (ti(i,j),j=1,nEng)
  read(126,1000) (texspex(i,j),j=1,nEng)
  read(128,1000) (sispex(i,j),j=1,nEng)
  read(129,1000) (sc(i,j),j=1,nEng)
  read(131,1000) (dispex(i,j),j=1,nEng)
  read(135,1000) (ss(i,j),j=1,nEng)
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
OXSDist=0.0
do i=1,nFiles !Loop through every process
  do j=1,ChS10 !Loop through every charge state, set XS=0 if CS doesn't exist
    do k=1,nEng !Loop through every energy
      select case(i)
      case(1)
        OXSDist(i,j,k)=si(j,k)
      case(2)
        if(j.le.9)OXSDist(i,j,k)=siss(j,k)
      case(3)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=sids(j-1,k)
      case(4)
        if(j.gt.1.and.j.le.9)OXSDist(i,j,k)=sispex(j-1,k)
      case(5)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=sidpex(j-1,k)
      case(6)
        OXSDist(i,j,k)=di(j,k)
      case(7)
        if(j.le.9)OXSDist(i,j,k)=diss(j,k)
      case(8)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=dids(j-1,k)
      case(9)
        if(j.gt.1.and.j.le.9)OXSDist(i,j,k)=dispex(j-1,k)
      case(10)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=didpex(j-1,k)
      case(11)
        if(j.gt.2)OXSDist(i,j,k)=ti(j-2,k)
      case(12)
        if(j.gt.2.and.j.le.9)OXSDist(i,j,k)=tiss(j-2,k)
      case(13)
        if(j.gt.2.and.j.le.8)OXSDist(i,j,k)=tids(j-2,k)
      case(14)
        if(j.gt.2.and.j.le.9)OXSDist(i,j,k)=tispex(j-2,k)
      case(15)
        if(j.gt.2.and.j.le.8)OXSDist(i,j,k)=tidpex(j-2,k)
      case(16)
        if(j.gt.3)OXSDist(i,j,k)=dcai(j-3,k)
      case(17)
        if(j.gt.3.and.j.le.9)OXSDist(i,j,k)=dcaiss(j-3,k)
      case(18)
        if(j.gt.3.and.j.le.8)OXSDist(i,j,k)=dcaids(j-3,k)
      case(19)
        if(j.gt.3.and.j.le.9)OXSDist(i,j,k)=dcaispex(j-3,k)
      case(20)
        if(j.gt.3.and.j.le.8)OXSDist(i,j,k)=dcaidpex(j-3,k)
      case(21)
        if(j.gt.2)OXSDist(i,j,k)=sc(j-2,k)
      case(22)
        if(j.gt.2.and.j.le.9)OXSDist(i,j,k)=scss(j-2,k)
      case(23)
        if(j.gt.2.and.j.le.8)OXSDist(i,j,k)=scds(j-2,k)
      case(24)
        if(j.gt.2.and.j.le.9)OXSDist(i,j,k)=scspex(j-2,k)
      case(25)
        if(j.gt.2.and.j.le.8)OXSDist(i,j,k)=scdpex(j-2,k)
      case(26)
        if(j.gt.3)OXSDist(i,j,k)=dc(j-3,k)
      case(27)
        if(j.gt.3.and.j.le.9)OXSDist(i,j,k)=dcss(j-3,k)
      case(28)
        if(j.gt.3.and.j.le.8)OXSDist(i,j,k)=dcds(j-3,k)
      case(29)
        if(j.gt.3.and.j.le.9)OXSDist(i,j,k)=dcspex(j-3,k)
      case(30)
        if(j.gt.3.and.j.le.8)OXSDist(i,j,k)=dcdpex(j-3,k)
      case(31)
        OXSDist(i,j,k)=tex(j,k)
      case(32)
        if(j.le.9)OXSDist(i,j,k)=texss(j,k)
      case(33)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=texds(j-1,k)
      case(34)
        if(j.gt.1.and.j.le.9)OXSDist(i,j,k)=texspex(j-1,k)
      case(35)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=texdpex(j-1,k)
      case(36)
        if(j.gt.1.and.j.le.9)OXSDist(i,j,k)=ss(j-1,k)
      case(37)
        if(j.gt.1.and.j.le.8)OXSDist(i,j,k)=ds(j-1,k)
      case(38)
        if(j.eq.2)OXSDist(i,j,k)=neg(j-1,k) !Only applies to O neutral
      end select
    end do
  end do
end do
1000 format (9x,12ES9.7)
end subroutine
