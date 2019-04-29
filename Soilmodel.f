
      Program Soilmodel

      Implicit none

      integer*2, parameter :: nlayer = 40             ! Number of subsurface layers
      real*8, parameter :: zmax = 6.                ! Maximum depth (m)
      real*8, parameter :: fac = 1.15                  ! factor for the definition of subsurface layers (1.15)
      real*8, parameter :: pi = 3.14159
      real*8, parameter :: S0 = 600.                  ! Solar constant at 1.52 AU (W/m²)
      real*8, parameter :: eccen = 0.0933             ! Eccentricity
      real*8, parameter :: albedo0 = 0.15
      real*8, parameter :: albedoCO2 = 0.65
      real*8, parameter :: sigma = 5.67E-08           ! Stefan-Boltzmann constant (W.m⁻².K⁻⁴)
      real*8, parameter :: emi = 1.                   ! Emissivity
      real*8, parameter :: marsday = 88775.244        ! Period of a Martian day (s)
      real*8, parameter :: hour = marsday/24.5        ! Period of a Martian hour (s)
      real*8, parameter :: marsyear = 59355128.1384   ! Period of a Martian year (s)
      real*8, parameter :: dt = 1849.4845           ! Time step (s), 1849.4845 = 1/2 Martian hour
      integer, parameter :: nsteps = 2*24*668*100       ! Number of time steps, multiple of 48 = 1 Mars day (= Sol)
      real*8, parameter :: qflux = 0.03               ! Heat flux (W/m²)
      real*8, parameter :: LHCO2 = 5.9E05             ! Latent heat of vaporization of CO2 (J/kg)
      real*8, parameter :: hs = 5.E-11                ! Radiogenic heat production (W/kg)
      real*8, parameter :: psurf = 610.               ! Surface pressure (Pa)
      real*8, parameter :: gmars = 3.72               ! Gravity on Mars (m s⁻²)
      real*8, parameter :: f_IR = 0.04
      real*8, parameter :: f_scat = 0.02

      real*8 :: obliq = 25.20                         ! Obliquity (degrees)
      real*8 :: lat = -40.0                            ! Latitude (degrees)
      real*8 :: omega = 286.5                         ! Argument of perihelion measured from the vernal equinox of Mars (degrees)
      real*8 :: TCO2frost = 148.

      real*8 :: slope                                 ! inclination of the slope (rad)
      real*8, parameter :: az = 0.*pi/180.	        	! slope orientation westward from south (NH), north(SH)

      real*8, dimension(nlayer) :: z                  ! Depth of soil layers (m)
      real*8, dimension(nlayer) :: tn                 ! Thickness of layers (m)
      real*8, dimension(nlayer) :: tsoil              ! Soil temperature (K)
      real*8, dimension(nlayer) :: cpsoil             ! Specific heat of soil layers (J kg⁻¹ K⁻¹)
      real*8, dimension(nlayer) :: rhosoil            ! Density of soil layers (kg/m³)
      real*8, dimension(nlayer) :: ksoil              ! Thermal conductivity of soil layers (W m⁻¹ K⁻¹)
      real*8, dimension(nlayer) :: diffu              ! Thermal diffusivity of soil layers (m²/s)
      real*8, dimension(nlayer+1) :: temp, tempnew    ! Vectors of surface temperature + soil temperatures
      real*8, dimension(nlayer+1) :: pres             ! Vectors of surface pressure + soil pressures (Pa)
      real*8, dimension(nlayer) :: alpha
      real*8, dimension(nlayer) :: gamma
      real*8, dimension(nlayer+1) :: u, v, w
      real*8, dimension(1000) :: gam

      integer :: i, j, k, m, ik,compteur
      real*8 :: dz                                    ! Space step (m)
      real*8 :: provis = 0
      real*8 :: delta                                 ! Skin depth (m)
      real*8 :: Q_IR, Q_land, Q_sun, Q_scat, Q_tot
      real*8 :: Smean, Smean_flat, Smean_year, Tmean_year         ! Average daily/yearly insolation (W/m²) and temperature (K)
      real*8 :: Tsurf, Tsurf_flat                                 ! Surface temperature (K)
      real*8 :: psy
      real*8 :: Tmax
      real*8 :: bet, test
      real*8 :: time
      real*8 :: day, hours
      real*8 :: MCO2, DMCO2, MCO2_fl, DMCO2_fl                      ! Mass of CO2 (kg m⁻²)
      real*8 :: mean_anom                             ! Mean anomaly
      real*8 :: eccen_anom                            ! Eccentric anomaly
      real*8 :: true_anom                             ! True anomaly
      real*8 :: r_sunmars
      logical :: flag1, flag2, flag3
      real*8 :: albedo
      real*8 :: ea_new, fct
      real*8 :: ea_new_prime, fct_prime
      real*8 :: guess1, guess2, guess3
      real*8 :: dec                                   ! Declination of sun in the sky
      integer :: rev
      real*8 :: ls                                    ! Solar longitude
      real*8 :: heta
      real*8 :: Fatm                                  ! Atmospheric thermal radiation in the downward direction
      real*8 :: H                                     ! Hour angle of the Sun
      real*8 :: cos_zh, cos_zh_noon                   ! cos of solar zenithal angle (flat surface)
      real*8 :: cos_zh_slope			                    ! cos of solar zenithal angle (over the slope)
      real*8 :: cos_psy                               ! cos of solar azimuth angle


      !==============================================================================
      !     Initialization of cpsoil, rhosoil, ksoil and diffu (values from Sizemore and Mellon, 2006)
      !==============================================================================

      do i=1, nlayer
        cpsoil(i) = 837.0                ! soil properties for a dry soil
        rhosoil(i) = 1255.0              ! Chevrier : rho*c = 1,05 E+6
        ksoil(i) = 0.085
        diffu(i) = ksoil(i)/(rhosoil(i)*cpsoil(i))
      enddo

!      do i=10, nlayer
!         cpsoil(i) = .4*(1540.)+.6*(837.)              ! soil properties for an ice-cemented soil (40%)
!         rhosoil(i) = .4*(927.)+.6*(1255.)
!         ksoil(i) = .4*(3.2)+.6*(.085)
!         diffu(i) = ksoil(i)/(rhosoil(i)*cpsoil(i))
!      enddo

      delta = sqrt(diffu(1)*marsday/pi)   ! diurnal skin depth

      !==============================================================================
      !     Definition of subsurface layers similar to M-SIM (Schorghofer)
      !==============================================================================

      if(nlayer .lt. 6) then
        write(*,*) 'Inappropriate number of layers (nlayer)'
      endif

      do i=1,nlayer
        dz = zmax/nlayer
        z(i) = (i-0.5)*dz
      enddo

      if(fac .gt. 1.) then
        dz=zmax/(3.+2.*fac*(fac**(nlayer-2)-1.)/(fac-1.))
        z(1)=dz
        z(2) = 3*z(1)
        do i=3,nlayer
          z(i)=(1.+fac)*z(i-1)-fac*z(i-2)
        enddo
      else
        write(*,*) 'Inappropriate input to fac'
      endif
      if(z(1) .lt. 1.e-5) then
        write(*,*) 'WARNING: first grid point is too shallow'
      endif

      if(z(6) .gt. delta) then
        write(*,*) 'WARNING: less than 6 points within diurnal',
     &                                            ' skin depth'
      endif

      tn(1) = z(1)
      do i=2,nlayer
        tn(i) = z(i)-z(i-1)
      enddo

      !      do i=1,nlayer
      !          write(*,*)i , z(i), tn(i), delta
      !      enddo


      !==============================================================================
      ! Layer temperatures are initialized by calculating the average annual equilibrium
      ! surface temperature.
      !==============================================================================


      lat = lat*pi/180.
      obliq = obliq*pi/180.
      omega = omega*pi/180.

      Smean_year=(S0/4.*(1.-eccen**2)**(0.5))*((1.5-(2.*SIN(obliq))/pi)
     &                -(1.5-(6.*SIN(obliq))/pi)*(SIN(lat))**2.)

      Tmean_year=(Smean_year*(1.-albedo0)/(sigma*emi))**0.25

      !      write(*,*) Smean_year, Tmean_year

      Tsurf = Tmean_year
      Tsurf_flat = Tmean_year

      do i=1,nlayer
        tsoil(i) = Tmean_year     ! tsoil holds the initial temperature profile
      enddo

      !==============================================================================

      time = 0.0
      MCO2 = 0.0       ! Initial mass of CO2
      MCO2_fl = 0.0

      open(4,file='Temperature.dat',status ='old')
      do k = 35,35
        slope=k*pi/180.
        Tmax=148.

        do j=1,nsteps
          time = j*dt
!          write(*,*) j, int(100.*(j-1)/nsteps), '% complete...'

          day = int(time/marsday)
          hours = int((time-day*marsday)/hour)

          mean_anom = 2.*pi*time/marsyear
          eccen_anom = mean_anom
          flag1 = .true.
          flag2 = .true.
          flag3 = .true.

          if(Tsurf .le. 148.0 .or. MCO2 .gt. 0.0) then
            albedo = albedoCO2
            flag2 = .false.

          else
            albedo = albedo0

          endif

          do while(flag1)

            ea_new = -mean_anom + eccen_anom - eccen*sin(eccen_anom)
            ea_new_prime = 1-eccen*cos(eccen_anom)
            guess1 = eccen_anom - ea_new/ea_new_prime

            if(abs(guess1-eccen_anom) .le. 1.0E-03) then
              flag1 = .false.
            endif

            eccen_anom = guess1

          enddo

        true_anom=2.*atan(sqrt((1+eccen)/(1-eccen))*tan(eccen_anom/2.))

          dec = asin(sin(obliq)*sin(true_anom + omega))

          rev = int((true_anom + omega)/(2.*pi))

          ls = (true_anom + omega)*180./pi - rev*360.

          H = pi*(abs(1-(time-day*marsday)/(12.*hour)))

          cos_zh = sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(H)
          cos_zh_noon = sin(lat)*sin(dec)+cos(lat)*cos(dec)

          cos_psy = (-sin(dec)*cos(lat)+cos(dec)*sin(lat)*cos(H))/
     &                sqrt(1-cos_zh**2)
          psy=acos(cos_psy)

          if (cos_zh .le. 0.0) then
            cos_zh_slope = 0.0
            cos_zh = 0
          else
            cos_zh_slope=cos(slope)*cos_zh-sqrt(1-cos_zh**2)*
     &            sin(slope)*cos(psy-az)
            provis=1
          endif

          if (cos_zh_slope .le. 0.0) then
            cos_zh_slope = 0.0
          endif

          !====================================================================
          !        For average daily insolation
          !====================================================================

          !         Smean = ((S0*(1+eccen*cos(true_anom))**2)/(1-eccen**2)**2)/pi*
          !     &            (heta*sin(dec)*sin(lat)+sin(heta)*cos(dec)*cos(lat))*

          !====================================================================
          !        For insolation at each time step
          !====================================================================

          Smean=(S0*(1+eccen*cos(true_anom))**2)/(1-eccen**2)**2

          !====================================================================

          Q_sun = Smean*cos_zh*(1.-albedo)
          Q_IR = f_IR*Smean *cos_zh_noon
          Q_scat = 1/2*f_scat*Smean*provis

          do while(flag3)
            do ik = 1, 50

            fct = -(emi*sigma*Tsurf_flat**4.)+Q_sun+Q_scat+Q_IR+
     &      ksoil(1)*(tsoil(1)-Tsurf_flat)/dz
            fct_prime = -(4.*emi*sigma*Tsurf_flat**3.)-ksoil(1)/dz

            guess3 = Tsurf_flat - fct/fct_prime

!            if(abs(guess3 - Tsurf_flat) .le. 1.0E-4) then
!               flag3 = .false.
!            endif

            Tsurf_flat = guess3
            end do
            flag3 = .false.
          enddo

          if(Tsurf_flat .le. TCO2frost .or. MCO2_fl .gt. 0.0) then
            Tsurf_flat = TCO2frost

            DMCO2_fl = (emi*sigma*Tsurf_flat**4.-(Q_sun+Q_scat+Q_IR)
     &              - ksoil(1)*(tsoil(1)-Tsurf_flat)/dz)*(dt/LHCO2)
            MCO2_fl = MCO2_fl + DMCO2_fl

            if(MCO2_fl .le. 0.0) then
              MCO2_fl = 0.0
              Tsurf_flat = 148.1
            endif

          endif

          Q_sun = Smean*cos_zh_slope*(1.-albedo)
          Q_IR = f_IR*Smean *cos_zh_noon*cos(slope/2)**2
          Q_scat = 1/2*f_scat*Smean*provis*cos(slope/2)**2
          Q_land = sigma*emi*Tsurf_flat**4*sin(slope/2)**2
          Q_tot=Q_sun+Q_scat+Q_IR+Q_land

          do while(flag2)
            do ik=0,50
            fct = -(emi*sigma*Tsurf**4.)+Q_tot+
     &      ksoil(1)*(tsoil(1)-Tsurf)/dz
            fct_prime = -(4.*emi*sigma*Tsurf**3.)-ksoil(1)/dz

            guess2 = Tsurf - fct/fct_prime

!            if(abs(guess2 - Tsurf) .le. 1.0E-4) then
!              flag2 = .false.
!            endif

            Tsurf = guess2
            end do
            flag2 = .false.
          enddo

          if(Tsurf .le. TCO2frost .or. MCO2 .gt. 0.0) then
            Tsurf = TCO2frost

            DMCO2 = (emi*sigma*Tsurf**4.-(Q_sun+Q_scat+Q_IR+Q_land)
     &              - ksoil(1)*(tsoil(1)-Tsurf)/dz)*(dt/LHCO2)
            MCO2 = MCO2 + DMCO2

            if(MCO2 .le. 0.0) then
              MCO2 = 0.0
              Tsurf = 148.1
            endif

          endif

          !==============================================================================
          ! Crank-Nicolson algorihtm to find soil temperatures.
          !==============================================================================

          do i=1, nlayer-1
            alpha(i) = diffu(i+1)*dt/(tn(i+1)*(tn(i+1)+tn(i)))
            gamma(i) = diffu(i)*dt/(tn(i)*(tn(i+1)+tn(i)))
          enddo
          gamma(nlayer) = diffu(i)*dt/(2.*tn(nlayer)**2)

          temp(1) = Tsurf
          temp(2) = gamma(1)*temp(1)+(1.-alpha(1)-gamma(1))*tsoil(1)
     &             +alpha(1)*tsoil(2)+dt*hs/cpsoil(1)
          do i=3, nlayer
            temp(i) = gamma(i-1)*tsoil(i-2)+(1.-alpha(i-1)-gamma(i-1))
     &               *tsoil(i-1)+alpha(i-1)*tsoil(i)+dt*hs/cpsoil(i-1)
          enddo
        temp(nlayer+1)=gamma(nlayer)*tsoil(nlayer-1)+(1.-gamma(nlayer))
     &          *tsoil(nlayer)+dt*hs/cpsoil(nlayer)+dt*(qflux-hs*
     &           rhosoil(nlayer)*z(nlayer))/(rhosoil(nlayer)*
     &           cpsoil(nlayer)*tn(nlayer))


          !------Resolution of the matrix inversion problem A * tempnew = temp-----------
          !           A is a tridiagonal matrix with the vectors u,v,w

          u(1) = 0.
          do i=2, nlayer+1
            u(i) = -gamma(i-1)
          enddo

          v(1) = 1.
          do i=2, nlayer
            v(i) = 1.+alpha(i-1)+gamma(i-1)
          enddo
          v(nlayer+1) = 1.+gamma(nlayer)

          w(1) = 0.
          do i =2, nlayer
            w(i) = -alpha(i-1)
          enddo
          w(nlayer+1) = 0.

          bet = v(1)
          tempnew(1) = temp(1)/bet
          do i=2, nlayer+1
            gam(i) = w(i-1)/bet
            bet = v(i)-u(i)*gam(i)
            tempnew(i) = (temp(i)-u(i)*tempnew(i-1))/bet
          enddo

          do i=nlayer+1,1,-1
            tempnew(i)=tempnew(i)-gam(i+1)*tempnew(i+1)
          enddo

          temp = tempnew
          pres(1) = psurf
          do i=1, nlayer
            tsoil(i) = temp(i+1)
            pres(i+1)=psurf+rhosoil(i)*gmars*z(i)
          enddo

          if (Tsurf .gt. Tmax) then
            Tmax=Tsurf
          end if

          if(time .gt. 1849.4845*48*668*99) then                ! 99 ans
            compteur = compteur+1
            if(mod(compteur,48) .eq. 0.0) then
              do i=1, nlayer
               if(
!     &           (min(abs(ls-0),abs(ls-360)) .lt. 0.5) .or.
     &           (abs(ls-30) .lt. 0.5) .or.
!     &           (abs(ls-60) .lt. 0.5) .or.
     &           (abs(ls-90) .lt. 0.5) .or.
!     &           (abs(ls-120) .lt. 0.5) .or.
     &           (abs(ls-150) .lt. 0.5) .or.
!     &           (abs(ls-180) .lt. 0.5) .or.
     &           (abs(ls-210) .lt. 0.5) .or.
!     &           (abs(ls-240) .lt. 0.5) .or.
     &           (abs(ls-270) .lt. 0.5) .or.
!     &           (abs(ls-300) .lt. 0.5) .or.
     &           (abs(ls-330) .lt. 0.5)
     &           ) then

                  write(4,*) ls, z(i), temp(i)

               endif
              end do
            endif
          endif

        enddo
!        write(4,*) k, Tmax
      enddo

      close(4)

      End program Soilmodel
