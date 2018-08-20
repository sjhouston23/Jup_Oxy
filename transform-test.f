      implicit real*8 (a-h,o-z)

      pi=4.d0*atan(1.d0)
      degrad=pi/180.

      write(6,*) 'Test of target-to-projectile frame transformation'
      write(6,*) 

      write(6,*) 'Test 1: Example from the appendix of the new paper'
      write(6,*) 'SS, O3+ 50 keV/u'
      write(6,*) 'target frame E_e = 92.71 eV'
      write(6,*) 'target frame theta_e = 9.116 deg'
      vrel=1.41471
      en=92.71
      an=9.116

      vz=(sqrt(2.*en/27.2116)*cos(an*degrad))-vrel
      v2=(2.*en/27.2116) - (2*vz*vrel) - (vrel*vrel) 
      ex=v2/2.*27.2116
      tx=acos(vz/sqrt(v2))/degrad
 
      write(6,*) 'projectile frame E_e = ',ex,' eV'
      write(6,*) 'projectile frame theta_e = ',tx,' deg'

      write(6,*) 

      write(6,*) 'Test 2: Dave result for'
      write(6,*) 'SS, O6+ 10 MeV/u'
      write(6,*) 'target frame E_e = 8590.354 eV'
      write(6,*) 'target frame theta_e = 5.86625 deg'
      vrel=20.00071
      en=8590.254
      an=5.86625

      vz=(sqrt(2.*en/27.2116)*cos(an*degrad))-vrel
      v2=(2.*en/27.2116) - (2*vz*vrel) - (vrel*vrel) 
      ex=v2/2.*27.2116
      tx=acos(vz/sqrt(v2))/degrad
 
      write(6,*) 'projectile frame E_e = ',ex,' eV'
      write(6,*) 'projectile frame theta_e = ',tx,' deg'

      write(6,*) 

      write(6,*) 'Test 3: Stephen result for'
      write(6,*) 'SS, O6+ 10 MeV/u'
      write(6,*) 'target frame E_e = 8726 eV'
      write(6,*) 'target frame theta_e = 7.28 deg'
      vrel=20.00071
      en=8726.
      an=7.28

      vz=(sqrt(2.*en/27.2116)*cos(an*degrad))-vrel
      v2=(2.*en/27.2116) - (2*vz*vrel) - (vrel*vrel) 
      ex=v2/2.*27.2116
      tx=acos(vz/sqrt(v2))/degrad
 
      write(6,*) 'projectile frame E_e = ',ex,' eV'
      write(6,*) 'projectile frame theta_e = ',tx,' deg'

      stop
      end
