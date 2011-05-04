!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_aerosol
!
      use mod_runparams
      use mod_message
      use mod_trachem
      use mod_slice
      use mod_main
      use mod_mainchem
!
      private
!
      public :: nspi
      public :: aermm , aerlwtr
      public :: dgmix , dssamix , dextmix
      public :: tauxar_mix , tauasc_mix
      public :: gtota_mix , ftota_mix
      public :: tauxar_mix_cs , tauasc_mix_cs
      public :: gtota_mix_cs , ftota_mix_cs
      public :: allocate_mod_aerosol , aermix , aeroppt , aerout
!
!     MODIF 16/09/2005 IBRAH Internal mixing
!
!      Num of spctrl intervals across solar spectrum
!
      integer , parameter :: nspi = 19   ! Spectral index
!
      integer , parameter :: ncoefs = 5  ! Number of coefficients
!
      real(8) , parameter :: d10e5  = 1.0D+05
      real(8) , parameter :: d10e4  = 1.0D+04
      real(8) , parameter :: nearone  = 0.99D+00
      real(8) , parameter :: minimum_aerosol = 1.0D-14
      real(8) , parameter :: minimum_utaer   = 1.0D-10
      real(8) , parameter :: minimum_waer   = 1.0D-30
      real(8) , parameter :: minimum_gaer   = 1.0D-20
      real(8) , parameter :: fiveothree  = d_five/d_three
!
! kscoef  - specific extinction (m2/g)
! wscoef  - single partical albedo
! gscoef  - asymmetry parameter
! ksbase  - specific extinction (m2/g) base
! wsbase  - single partical albedo base
! gsbase  - asymmetry parameter base
! ksdust  - specific extinction (m2/g) dust
! wsdust  - single partical albedo dust
! gsdust  - asymmetry parameter dust
!
      integer , private :: ii , jj ! coefficient index
!
      real(8) , dimension(nspi) :: gsbase , gsbc_hb , gsbc_hl ,   &
                                   gsoc_hb , gsoc_hl , ksbase ,   &
                                   ksbc_hb , ksbc_hl , ksoc_hb ,  &
                                   ksoc_hl , wsbase , wsbc_hb ,   &
                                   wsbc_hl , wsoc_hb , wsoc_hl
      real(8) , dimension(nspi,ncoefs) :: gscoef , kscoef , wscoef
      real(8) , dimension(nspi,4) :: gsdust , ksdust , wsdust
      real(4) , dimension(4,19,11,11,11,11) :: dextmix , dgmix , dssamix
!
!     Aerosol mass mixing ratio
!
      real(8) , allocatable , dimension(:,:,:) :: aermm
!
!     Background aerosol mass mixing ratio
!
      real(8) , allocatable , dimension(:,:) :: aermmb
!
!     Radiation level aerosol mass mixing ratio
!
      real(8) , allocatable , dimension(:,:,:) :: aermmr
!
!     Aerosol optical properties (for the mixing) 
!
      real(8) , allocatable , dimension(:,:,:) :: ftota_mix ,   &
                  gtota_mix , tauasc_mix , tauxar_mix
!
      real(8) , allocatable , dimension(:,:) :: ftota_mix_cs ,  &
                  gtota_mix_cs , tauasc_mix_cs , tauxar_mix_cs
!
!     Work arrays for aeroppt (aerosol individual optical properties SW)
!
      real(8) , allocatable , dimension(:,:) :: aermtot , aervtot
      real(8) , allocatable , dimension(:,:,:) :: fa , ga , tauxar , &
                                                 uaer , wa
      real(8) , allocatable , dimension(:,:) :: faer , gaer , tauaer , &
                                 utaer , waer
      real(8) , dimension(4) :: prop
!
!   Aersol LW optical properties
      real(8) ,allocatable, dimension(:,:,:) ::  aerlwtr 
!
!------------------------------------------------------------------------------
!                  DATA SECTION
!------------------------------------------------------------------------------
!
      data ksbase/5.206D0 , 5.206D0 , 5.206D0 , 5.206D0 , 5.206D0 ,    &
          5.206D0 , 5.206D0 , 3.203D0 , 3.203D0 , 1.302D0 , 5.992D-01 ,&
          2.948D-01 , 1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 ,  &
          5.770D-02 , 2.290D-01 , 2.270D-01/
!
      data ((kscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/1.126D+01 ,      &
           -2.502D-01 , -1.087D+00 , -1.794D+02 , 1.556D+01 ,       &
           1.126D+01 , -2.502D-01 , -1.087D+00 , -1.794D+02 ,       &
           1.556D+01 , 1.126D+01 , -2.502D-01 , -1.087D+00 ,        &
           -1.794D+02 , 1.556D+01 , 1.126D+01 , -2.502D-01 ,        &
           -1.087D+00 , -1.794D+02 , 1.556D+01 , 1.126D+01 ,        &
           -2.502D-01 , -1.087D+00 , -1.794D+02 , 1.556D+01 ,       &
           1.126D+01 , -2.502D-01 , -1.087D+00 , -1.794D+02 ,       &
           1.556D+01 , 1.126D+01 , -2.502D-01 , -1.087D+00 ,        &
           -1.794D+02 , 1.556D+01 , 1.124D+01 , -3.040D-01 ,        &
           -1.088D+00 , -1.776D+02 , 1.537D+01 , 1.124D+01 ,        &
           -3.040D-01 , -1.088D+00 , -1.776D+02 , 1.537D+01 ,       &
           1.222D+01 , -3.770D-01 , -1.089D+00 , -1.898D+02 ,       &
           1.504D+01 , 1.357D+01 , -4.190D-01 , -1.087D+00 ,        &
           -2.070D+02 , 1.478D+01 , 1.557D+01 , -4.353D-01 ,        &
           -1.083D+00 , -2.382D+02 , 1.486D+01 , 1.758D+01 ,        &
           -4.389D-01 , -1.078D+00 , -2.716D+02 , 1.505D+01 ,       &
           1.597D+01 , -4.337D-01 , -1.073D+00 , -2.510D+02 ,       &
           1.527D+01 , 2.107D+01 , -3.041D-01 , -1.067D+00 ,        &
           -2.494D+02 , 1.166D+01 , -2.424D-01 , -1.770D-01 ,       &
           -1.032D+00 , 1.469D-01 , 1.947D+00 , 2.535D+01 ,         &
           -2.270D-01 , -1.052D+00 , -2.528D+02 , 9.888D+00 ,       &
           -1.545D-01 , -1.661D-01 , -1.030D+00 , -4.698D-04 ,      &
           7.275D-02 , 8.835D-01 , -1.590D-01 , -1.029D+00 ,        &
           -2.838D+01 , 2.734D+01/
!
      data wsbase/7.371D-08 , 7.371D-08 , 7.371D-08 , 7.371D-08 ,     &
          7.371D-08 , 7.371D-08 , 7.371D-08 , 6.583D-08 , 6.583D-08 , &
          3.656D-06 , 4.919D-05 , 3.539D-03 , 2.855D-02 , 2.126D-01 , &
          8.433D-01 , 9.653D-01 , 6.198D-01 , 9.642D-01 , 9.699D-01/
!
      data ((wscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/2.492D+00 ,        &
           -5.210D-02 , -1.036D+00 , -4.398D+01 , 1.724D+01 ,         &
           2.492D+00 , -5.210D-02 , -1.036D+00 , -4.398D+01 ,         &
           1.724D+01 , 2.492D+00 , -5.210D-02 , -1.036D+00 ,          &
           -4.398D+01 , 1.724D+01 , 2.492D+00 , -5.210D-02 ,          &
           -1.036D+00 , -4.398D+01 , 1.724D+01 , 2.492D+00 ,          &
           -5.210D-02 , -1.036D+00 , -4.398D+01 , 1.724D+01 ,         &
           2.492D+00 , -5.210D-02 , -1.036D+00 , -4.398D+01 ,         &
           1.724D+01 , 2.492D+00 , -5.210D-02 , -1.036D+00 ,          &
           -4.398D+01 , 1.724D+01 , 1.139D+00 , -1.110D-02 ,          &
           -1.011D+00 , -7.754D+00 , 6.737D+00 , 1.139D+00 ,          &
           -1.110D-02 , -1.011D+00 , -7.754D+00 , 6.737D+00 ,         &
           1.848D+00 , -3.920D-04 , -9.924D-01 , -1.607D+00 ,         &
           8.587D-01 , 5.459D+00 , 9.357D-01 , -1.626D+00 ,           &
           -5.282D+00 , 1.066D+00 , 1.187D+00 , 2.241D-01 ,           &
           -1.226D+00 , 1.442D+01 , -1.402D+01 , -3.640D+00 ,         &
           2.552D-01 , -1.168D+00 , 4.458D+01 , 1.152D+01 ,           &
           -5.634D+00 , 2.068D-01 , -1.122D+00 , 7.528D+01 ,          &
           1.290D+01 , 1.826D-01 , 6.588D-02 , -1.098D+00 ,           &
           -1.996D-02 , 1.618D-01 , 2.164D+00 , 1.194D-01 ,           &
           -1.044D+00 , -3.221D+01 , 1.564D+01 , 2.268D-01 ,          &
           3.266D-02 , -1.064D+00 , -2.677D-02 , 1.309D-01 ,          &
           2.178D+00 , 1.151D-01 , -1.042D+00 , -3.325D+01 ,          &
           1.600D+01 , 1.713D+00 , 9.166D-02 , -1.039D+00 ,           &
           -2.660D+01 , 1.629D+01/
!
      data gsbase/6.899D-01 , 6.899D-01 , 6.899D-01 , 6.899D-01 ,     &
          6.899D-01 , 6.899D-01 , 6.899D-01 , 6.632D-01 , 6.632D-01 , &
          5.912D-01 , 5.111D-01 , 4.269D-01 , 3.321D-01 , 2.197D-01 , &
          1.305D-01 , 7.356D-02 , 1.602D-01 , 6.883D-02 , 6.304D-02/
!
      data ((gscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/ -9.874D-01 ,   &
           -3.033D+01 , -2.138D+01 , -2.265D+00 , 5.238D+00 ,      &
           -9.874D-01 , -3.033D+01 , -2.138D+01 , -2.265D+00 ,     &
           5.238D+00 , -9.874D-01 , -3.033D+01 , -2.138D+01 ,      &
           -2.265D+00 , 5.238D+00 , -9.874D-01 , -3.033D+01 ,      &
           -2.138D+01 , -2.265D+00 , 5.238D+00 , -9.874D-01 ,      &
           -3.033D+01 , -2.138D+01 , -2.265D+00 , 5.238D+00 ,      &
           -9.874D-01 , -3.033D+01 , -2.138D+01 , -2.265D+00 ,     &
           5.238D+00 , -9.874D-01 , -3.033D+01 , -2.138D+01 ,      &
           -2.265D+00 , 5.238D+00 , -3.666D-01 , -1.319D+00 ,      &
           -3.311D+00 , -2.821D-02 , 8.844D-01 , -3.666D-01 ,      &
           -1.319D+00 , -3.311D+00 , -2.821D-02 , 8.844D-01 ,      &
           5.824D-01 , -1.875D-01 , -1.567D+00 , -4.402D+00 ,      &
           6.268D+00 , 1.238D+00 , -1.550D-01 , -1.368D+00 ,       &
           -1.127D+01 , 8.334D+00 , 2.299D+00 , -1.686D-01 ,       &
           -1.304D+00 , -2.677D+01 , 1.101D+01 , 3.037D+00 ,       &
           -1.447D-01 , -1.223D+00 , -2.609D+01 , 8.267D+00 ,      &
           4.683D+00 , -2.307D-01 , -1.241D+00 , -4.312D+01 ,      &
           8.838D+00 , 3.842D+00 , -6.301D-01 , -1.367D+00 ,       &
           -4.144D+01 , 9.620D+00 , 3.237D+00 , -4.530D-01 ,       &
           -1.204D+00 , -3.234D+01 , 8.946D+00 , 4.181D+00 ,       &
           -4.140D-01 , -1.284D+00 , -4.489D+01 , 9.950D+00 ,      &
           3.378D+00 , -4.334D-01 , -1.188D+00 , -3.664D+01 ,      &
           9.786D+00 , 3.943D+00 , -3.952D-01 , -1.170D+00 ,       &
           -4.415D+01 , 1.031D+01/
!
!-----------------------------------------------------------------------
!
      data ksbc_hb /20.7830D0 , 17.2120D0 , 15.8640D0 , 15.0530D0 , &
                    14.3040D0 , 13.6130D0 , 11.9660D0 ,  6.5782D0 , &
                     4.3961D0 ,  4.3800D0 ,  2.1100D0 ,  2.1100D0 , &
                     2.1100D0 ,  2.1100D0 ,  1.4000D0 ,  1.4000D0 , &
                     1.4000D0 ,  1.4000D0 ,  1.4000D0/
 
      data wsbc_hb /0.245240D0 , 0.209620D0 , 0.195000D0 , 0.185860D0 , &
                    0.177190D0 , 0.168970D0 , 0.148490D0 , 0.071748D0 , &
                    0.037536D0 , 0.089000D0 , 0.025000D0 , 0.025000D0 , &
                    0.025000D0 , 0.025000D0 , 0.009000D0 , 0.009000D0 , &
                    0.009000D0 , 0.009000D0 , 0.009000D0/
 
      data gsbc_hb /0.213870D0 , 0.171540D0 , 0.155640D0 , 0.146100D0 , &
                    0.137320D0 , 0.129230D0 , 0.110060D0 , 0.049175D0 , &
                    0.026638D0 , 0.220000D0 , 0.123000D0 , 0.123000D0 , &
                    0.123000D0 , 0.123000D0 , 0.073000D0 , 0.073000D0 , &
                    0.073000D0 , 0.073000D0 , 0.073000D0/
! 
!----------------------------------------------------------------------
!
      data ksbc_hl /14.8510D0 , 14.2580D0 , 13.9430D0 , 13.7240D0 , &
                    13.5070D0 , 13.2950D0 , 12.7220D0 ,  9.4434D0 , &
                     6.9653D0 ,  4.3800D0 ,  2.1100D0 ,  2.1100D0 , &
                     2.1100D0 ,  2.1100D0 ,  1.4000D0 ,  1.4000D0 , &
                     1.4000D0 ,  1.4000D0 ,  1.4000D0/
 
      data wsbc_hl /0.46081D0 , 0.44933D0 , 0.44397D0 , 0.44065D0 , &
                    0.43737D0 , 0.43394D0 , 0.42346D0 , 0.35913D0 , &
                    0.29579D0 , 0.08900D0 , 0.02500D0 , 0.02500D0 , &
                    0.02500D0 , 0.02500D0 , 0.00900D0 , 0.00900D0 , &
                    0.00900D0 , 0.00900D0 , 0.00900D0/
 
      data gsbc_hl /0.69038D0 , 0.65449D0 , 0.63711D0 , 0.62542D0 , &
                    0.61407D0 , 0.60319D0 , 0.57467D0 , 0.42050D0 , &
                    0.30660D0 , 0.22000D0 , 0.12300D0 , 0.12300D0 , &
                    0.12300D0 , 0.12300D0 , 0.07300D0 , 0.07300D0 , &
                    0.07300D0 , 0.07300D0 , 0.07300D0/
! 
!---------------------------------------------------------------------
!
      data ksoc_hb /6.0584D0 ,  6.0654D0 ,  6.1179D0 ,  6.0102D0 ,  &
                    5.8000D0 ,  5.6957D0 ,  5.6494D0 ,  4.3283D0 ,  &
                    3.1485D0 ,  1.3020D0 ,  5.992D-01 , 2.948D-01 , &
                    1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 , &
                    5.770D-02 , 2.290D-01 , 2.270D-01/
 
      data wsoc_hb /0.91735D0 , 0.92365D0 , 0.92941D0 , 0.93067D0 , &
                    0.93311D0 , 0.93766D0 , 0.94042D0 , 0.95343D0 , &
                    0.95480D0 , 0.70566D0 , 0.70566D0 , 0.70566D0 , &
                    0.70566D0 , 0.70566D0 , 0.70566D0 , 0.70566D0 , &
                    0.70566D0 , 0.70566D0 , 0.70566D0/
 
      data gsoc_hb /0.67489D0 ,  0.67003D0 ,  0.67725D0 ,  0.65487D0 ,  &
                    0.65117D0 ,  0.66116D0 ,  0.64547D0 ,  0.60033D0 ,  &
                    0.55389D0 ,  5.9120D-01 , 5.1110D-01 , 4.2690D-01 , &
                    3.3210D-01 , 2.1970D-01 , 1.3050D-01 , 7.3560D-02 , &
                    1.6020D-01 , 6.8830D-02 , 6.3040D-02/
! 
!---------------------------------------------------------------------
!
      data ksoc_hl / 3.5430D0 ,  3.6230D0 ,  3.7155D0  , 3.7120D0 ,  &
                     3.6451D0 ,  3.6376D0 ,  3.6844D0  , 3.4588D0 ,  &
                     2.9846D0 ,  1.3020D0 ,  5.992D-01 , 2.948D-01 , &
                     1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 , &
                     5.770D-02 , 2.290D-01 , 2.270D-01/
! 
      data wsoc_hl /0.87931D0 , 0.88292D0 , 0.89214D0 , 0.89631D0 , &
                    0.89996D0 , 0.90540D0 , 0.90805D0 , 0.93423D0 , &
                    0.95012D0 , 0.85546D0 , 0.85546D0 , 0.85546D0 , &
                    0.85546D0 , 0.85546D0 , 0.85546D0 , 0.85546D0 , &
                    0.85546D0 , 0.85546D0 , 0.85546D0/
! 
      data gsoc_hl /0.73126D0 , 0.71089D0 , 0.72042D0 , 0.69924D0 , &
                    0.69908D0 , 0.70696D0 , 0.68479D0 , 0.64879D0 , &
                    0.63433D0 , 5.912D-01 , 5.111D-01 , 4.269D-01 , &
                    3.321D-01 , 2.197D-01 , 1.305D-01 , 7.356D-02 , &
                    1.602D-01 , 6.883D-02 , 6.304D-02/
! 
!     DUST OP data base for external mixing : maximum of 4 bin for the
!     momeent , determined from Zender et al.
! 
      data ((ksdust(ii,jj),jj=1,4),ii=1,nspi) / 1.88010D0 , 0.76017D0 , &
            0.36681D0 , 0.16933D0 , 2.02540D0 , 0.78378D0 , 0.36845D0 , &
            0.17002D0 , 1.95470D0 , 0.76389D0 , 0.37119D0 , 0.17032D0 , &
            1.89960D0 , 0.74916D0 , 0.37220D0 , 0.17052D0 , 1.79460D0 , &
            0.74044D0 , 0.36984D0 , 0.17072D0 , 1.71490D0 , 0.75401D0 , &
            0.36892D0 , 0.17090D0 , 1.54310D0 , 0.82322D0 , 0.37505D0 , &
            0.17143D0 , 2.44820D0 , 0.85680D0 , 0.38078D0 , 0.17396D0 , &
            3.10670D0 , 0.74488D0 , 0.43688D0 , 0.18104D0 , 0.50391D0 , &
            0.67245D0 , 0.53605D0 , 0.20599D0 , 0.50391D0 , 0.67245D0 , &
            0.53605D0 , 0.20599D0 , 0.50391D0 , 0.67245D0 , 0.53605D0 , &
            0.20599D0 , 0.50391D0 , 0.67245D0 , 0.53605D0 , 0.20599D0 , &
            0.50391D0 , 0.67245D0 , 0.53605D0 , 0.20599D0 , 0.50391D0 , &
            0.67245D0 , 0.53605D0 , 0.20599D0 , 0.50391D0 , 0.67245D0 , &
            0.53605D0 , 0.20599D0 , 0.50391D0 , 0.67245D0 , 0.53605D0 , &
            0.20599D0 , 0.50391D0 , 0.67245D0 , 0.53605D0 , 0.20599D0 , &
            0.50391D0 , 0.67245D0 , 0.53605D0 , 0.20599D0/
! 
      data ((wsdust(ii,jj),jj=1,4),ii=1,nspi) / 0.64328D0 , 0.55196D0 , &
            0.53748D0 , 0.54342D0 , 0.67757D0 , 0.56909D0 , 0.53639D0 , &
            0.54232D0 , 0.67316D0 , 0.56027D0 , 0.53875D0 , 0.54181D0 , &
            0.66245D0 , 0.55338D0 , 0.53947D0 , 0.54149D0 , 0.68132D0 , &
            0.57440D0 , 0.54328D0 , 0.54143D0 , 0.67960D0 , 0.58467D0 , &
            0.54242D0 , 0.54113D0 , 0.72679D0 , 0.68744D0 , 0.58564D0 , &
            0.54576D0 , 0.94730D0 , 0.88181D0 , 0.80761D0 , 0.70455D0 , &
            0.97536D0 , 0.89161D0 , 0.86378D0 , 0.75800D0 , 0.89568D0 , &
            0.96322D0 , 0.95008D0 , 0.89293D0 , 0.89568D0 , 0.96322D0 , &
            0.95008D0 , 0.89293D0 , 0.89568D0 , 0.96322D0 , 0.95008D0 , &
            0.89293D0 , 0.89568D0 , 0.96322D0 , 0.95008D0 , 0.89293D0 , &
            0.89568D0 , 0.96322D0 , 0.95008D0 , 0.89293D0 , 0.89568D0 , &
            0.96322D0 , 0.95008D0 , 0.89293D0 , 0.89568D0 , 0.96322D0 , &
            0.95008D0 , 0.89293D0 , 0.89568D0 , 0.96322D0 , 0.95008D0 , &
            0.89293D0 , 0.89568D0 , 0.96322D0 , 0.95008D0 , 0.89293D0 , &
            0.89568D0 , 0.96322D0 , 0.95008D0 , 0.89293D0/
! 
      data ((gsdust(ii,jj),jj=1,4),ii=1,nspi) / 0.87114D0 , 0.92556D0 , &
            0.94542D0 , 0.94831D0 , 0.86127D0 , 0.92100D0 , 0.94355D0 , &
            0.94813D0 , 0.83800D0 , 0.91194D0 , 0.94304D0 , 0.94803D0 , &
            0.81760D0 , 0.90442D0 , 0.94239D0 , 0.94796D0 , 0.77088D0 , &
            0.88517D0 , 0.93710D0 , 0.94775D0 , 0.73925D0 , 0.88364D0 , &
            0.93548D0 , 0.94763D0 , 0.60695D0 , 0.86086D0 , 0.91824D0 , &
            0.94473D0 , 0.64393D0 , 0.76457D0 , 0.81330D0 , 0.87784D0 , &
            0.74760D0 , 0.62041D0 , 0.80665D0 , 0.85974D0 , 0.26761D0 , &
            0.56045D0 , 0.68897D0 , 0.68174D0 , 0.26761D0 , 0.56045D0 , &
            0.68897D0 , 0.68174D0 , 0.26761D0 , 0.56045D0 , 0.68897D0 , &
            0.68174D0 , 0.26761D0 , 0.56045D0 , 0.68897D0 , 0.68174D0 , &
            0.26761D0 , 0.56045D0 , 0.68897D0 , 0.68174D0 , 0.26761D0 , &
            0.56045D0 , 0.68897D0 , 0.68174D0 , 0.26761D0 , 0.56045D0 , &
            0.68897D0 , 0.68174D0 , 0.26761D0 , 0.56045D0 , 0.68897D0 , &
            0.68174D0 , 0.26761D0 , 0.56045D0 , 0.68897D0 , 0.68174D0 , &
            0.26761D0 , 0.56045D0 , 0.68897D0 , 0.68174D0/
!
      contains
! 
      subroutine allocate_mod_aerosol
      implicit none   
#ifdef MPP1
      allocate(aermm(iym1,kz,jxp))
#else
#ifdef BAND
      allocate(aermm(iym1,kz,jx))
#else
      allocate(aermm(iym1,kz,jxm1))
#endif 
#endif 
      allocate(aermmb(iym1,kz))
      allocate(ftota_mix(iym1,0:kz,nspi))
      allocate(gtota_mix(iym1,0:kz,nspi))
      allocate(tauasc_mix(iym1,0:kz,nspi))
      allocate(tauxar_mix(iym1,0:kz,nspi))
      allocate(ftota_mix_cs(iym1,nspi))
      allocate(gtota_mix_cs(iym1,nspi))
      allocate(tauasc_mix_cs(iym1,nspi))
      allocate(tauxar_mix_cs(iym1,nspi))
      allocate(aermtot(iym1,kz))
      allocate(aervtot(iym1,kz))
      allocate(aerlwtr(iym1,kzp1,kzp1))
      if ( ichem == 1 ) then
        allocate(aermmr(iym1,kz,ntr))
        allocate(fa(iym1,0:kz,ntr))
        allocate(ga(iym1,0:kz,ntr))
        allocate(tauxar(iym1,0:kz,ntr))
        allocate(uaer(iym1,0:kz,ntr))
        allocate(wa(iym1,0:kz,ntr))
        allocate(faer(iym1,ntr))
        allocate(gaer(iym1,ntr))
        allocate(tauaer(iym1,ntr))
        allocate(utaer(iym1,ntr))
        allocate(waer(iym1,ntr))
      end if
      aermm = d_zero
      aermmb = d_zero
      ftota_mix = d_zero
      gtota_mix = d_zero
      tauasc_mix = d_zero
      tauxar_mix = d_zero
      ftota_mix_cs = d_zero
      gtota_mix_cs = d_zero
      tauasc_mix_cs = d_zero
      tauxar_mix_cs = d_zero
      aermtot = d_zero
      aervtot = d_zero
      aerlwtr = d_zero
      if ( ichem == 1 ) then
        aermmr = d_zero
        fa = d_zero
        ga = d_zero
        tauxar = d_zero
        uaer = d_zero
        wa = d_zero
        faer = d_zero
        gaer = d_zero
        tauaer = d_zero
        utaer = d_zero
        waer = d_zero
      end if
      end subroutine allocate_mod_aerosol
!
!-----------------------------------------------------------------------
!
! SUBROUTINE AERMIX
!
! Set global mean tropospheric aerosol
!
! Specify aerosol mixing ratio and compute relative humidity for later
! adjustment of aerosol optical properties. Aerosol mass mixing ratio
! is specified so that the column visible aerosol optical depth is a
! specified global number (tauvis). This means that the actual mixing
! ratio depends on pressure thickness of the lowest three atmospheric
! layers near the surface.
!
! Optical properties and relative humidity parameterization are from:
!
! J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
! and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
! 16 April 1993
!
! Visible (vis) here means 0.5-0.7 micro-meters
! Forward scattering fraction is taken as asymmetry parameter squared
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb  March 1995
! Standarized:       L. Buja,     Feb 1996
! Reviewed:          B. Briegleb, Mar 1996
!
!-----------------------------------------------------------------------
!
      subroutine aermix(pint , rh , j , istart , iend , nx , nk , ntrac)
! 
      implicit none
!
      integer :: j , istart , iend , nx , nk , ntrac
!     Radiation level interface pressures (dynes/cm2)
      real(8) , dimension(nx,nk+1) :: pint
!     Radiation level relative humidity (fraction)
      real(8) , dimension(nx,nk) :: rh
!
      intent (in) j , pint , istart , iend , nk , ntrac
      intent (out) rh
!
! i      - longitude index
! k      - level index
! mxaerl - max nmbr aerosol levels counting up from surface
! tauvis - visible optical depth
! kaervs - visible extinction coefficiant of aerosol (m2/g)
! omgvis - visible omega0
! gvis   - visible forward scattering asymmetry parameter
!
!-----------------------------------------------------------------------
! 
      real(8) :: gvis , kaervs , omgvis , rhfac , tauvis
      integer :: i , itr , k , mxaerl
!
      data kaervs /5.3012D0/        ! multiplication factor for kaer
      data omgvis /0.999999D0/
      data gvis   /0.694889D0/
      data rhfac  /1.6718D0/        ! EES added for efficiency
!
!--------------------------------------------------------------------------
!
      mxaerl = 4
!
!fil  tauvis = 0.01D0
      tauvis = 0.04D0
!
!     Set relative humidity and factor; then aerosol amount:
!
      do k = 1 , nk
        do i = istart , iend
 
!added    July 13, 2000: needed for aerosols in radiation
          rh(i,k) = dmin1(rhb3d(i,k,j),nearone)
!EES:     do not change to 1.00:  wscoef(3,10) in radcsw = .9924 and is
!         divided by RH.  rh is limited to .99 to avoid dividing by zero
!added
!
!         Define background aerosol
!         Find constant aerosol mass mixing ratio for specified levels
!         in the column, converting units where appropriate
!         for the moment no more used
!
          if ( k >= nk + 1 - mxaerl ) then
            aermmb(i,k) = egravgts*tauvis/ &
                         (d10e4*kaervs*rhfac*(d_one-omgvis*gvis*gvis) &
                         *(pint(i,kzp1)-pint(i,kz + 1 - mxaerl)))
          else
            aermmb(i,k) = d_zero
          end if
! 
          if ( ichem == 1 ) then
            do itr = 1 , ntrac
              aermmr(i,k,itr) = chia(i,k,j,itr)/sps1%ps(i,j)
            end do
          end if

        end do
      end do
!
      end subroutine aermix
!
! SUBROUTINE AEROPPT
!
      subroutine aeroppt(rh,pint)
!
      implicit none
!
!     Interface pressure, relative humidity
!
      real(8) , dimension(iym1,kzp1) :: pint
      real(8) , dimension(iym1,kz) :: rh
      intent (in) pint , rh
!
      integer :: i , i1 , i2 , i3 , i4 , ibin , itr , k , k1, k2 , ns
      real(8) :: path , uaerdust , qabslw
!
! uaer, tauxar  - Aerosol radiative properties (local arrays)
! wa            - Aerosol single scattering albedo
! ga            - Aerosol asimmetry parameter
! fa            - Aerosol forward scattered fraction
! utaer, tauaer - Total column aerosol extinction
! waer          - Aerosol single scattering albedo
! gaer          - Aerosol asymmetry parameter
! faer          - Aerosol forward scattered fraction
!
      if ( ichem /= 1 ) then
        tauxar_mix_cs(:,:) = d_zero
        tauasc_mix_cs(:,:) = d_zero
        gtota_mix_cs(:,:) = d_zero
        ftota_mix_cs(:,:) = d_zero
        tauxar_mix(:,:,:) = d_zero
        tauasc_mix(:,:,:) = d_zero
        gtota_mix(:,:,:) = d_zero
        ftota_mix(:,:,:) = d_zero
        aerlwtr (:,:,:) = d_one
        return
      end if
!
      if ( idirect >= 1 ) then
        tauxar = d_zero
        wa = d_zero
        ga = d_zero
        fa = d_zero
      else
!
!       Nothing to do for this
!
        tauxar_mix_cs(:,:) = d_zero
        tauasc_mix_cs(:,:) = d_zero
        gtota_mix_cs(:,:) = d_zero
        ftota_mix_cs(:,:) = d_zero
        tauxar_mix(:,:,:) = d_zero
        tauasc_mix(:,:,:) = d_zero
        gtota_mix(:,:,:) = d_zero
        ftota_mix(:,:,:) = d_zero
        aerlwtr (:,:,:) = d_one
        return
      end if
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!       Option 1  melange externe
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if ( mixtype == 1 ) then
!
!       Spectral loop
!
        do ns = 1 , nspi

          tauxar_mix_cs(:,ns) = d_zero
          tauasc_mix_cs(:,ns) = d_zero
          gtota_mix_cs(:,ns) = d_zero
          ftota_mix_cs(:,ns) = d_zero
!
          tauxar_mix(:,:,ns) = d_zero
          tauasc_mix(:,:,ns) = d_zero
          gtota_mix(:,:,ns) = d_zero
          ftota_mix(:,:,ns) = d_zero
!
          uaer(:,0,:) = d_zero
          tauxar(:,0,:) = d_zero
          wa(:,0,:) = d_zero
          ga(:,0,:) = d_zero
          fa(:,0,:) = d_zero
          utaer(:,:) = d_zero
          tauaer(:,:) = d_zero
          waer(:,:) = d_zero
          gaer(:,:) = d_zero
          faer(:,:) = d_zero
! 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!         calculate optical properties of each aerosol component
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
          do k = 1 , kz
            do i = 1 , iym1
              path = (pint(i,k+1)-pint(i,k))*regravgts
              ibin = 0
              do itr = 1 , ntr
                uaer(i,k,itr) = d_zero
                if ( rh(i,k) < d_zero .or. rh(i,k) > d_one ) &
                  print * , i , k , rh(i,k) , '  RH WARNING !!!!!'
!
                if ( chtrname(itr) == 'XXXXX') then
                  continue
                end if
!
                if ( chtrname(itr) == 'DUST' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 gaffe au facteur !!
!                 gaffe au ntr/bins
                  ibin = ibin + 1
                  if ( ibin > 4 ) print * , 'DUST OP PBLEME !!!!'
!
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksdust(ns,ibin)
                  wa(i,k,itr) = wsdust(ns,ibin)
                  ga(i,k,itr) = gsdust(ns,ibin)
                  fa(i,k,itr) = gsdust(ns,ibin)*gsdust(ns,ibin)
! 
                else if ( chtrname(itr) == 'SO4' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksbase(ns)       &
                                   *dexp(kscoef(ns,1)+kscoef(ns,2)     &
                                   /(rh(i,k)+kscoef(ns,3))             &
                                   +kscoef(ns,4)                       &
                                   /(rh(i,k)+kscoef(ns,5)))
!
                  wa(i,k,itr) = d_one - wsbase(ns)                      &
                               *dexp(wscoef(ns,1)+wscoef(ns,2)         &
                               /(rh(i,k)+wscoef(ns,3))+wscoef(ns,4)    &
                               /(rh(i,k)+wscoef(ns,5)))
!
                  ga(i,k,itr) = gsbase(ns)                              &
                               *dexp(gscoef(ns,1)+gscoef(ns,2)         &
                               /(rh(i,k)+gscoef(ns,3))+gscoef(ns,4)    &
                               /(rh(i,k)+gscoef(ns,5)))
!
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
! 
                else if ( chtrname(itr) == 'OC_HL' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Humidity effect !
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksoc_hl(ns)      &
                                   *(d_one-rh(i,k))**(-0.2D0)
                  wa(i,k,itr) = wsoc_hl(ns)
                  ga(i,k,itr) = gsoc_hl(ns)
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                else if ( chtrname(itr) == 'BC_HL' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Humidity effect !
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksbc_hl(ns)      &
                                   *(d_one-rh(i,k))**(-0.25D0)
                  wa(i,k,itr) = wsbc_hl(ns)
                  ga(i,k,itr) = gsbc_hl(ns)
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                else if ( chtrname(itr) == 'OC_HB' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksoc_hb(ns)
                  wa(i,k,itr) = wsoc_hb(ns)
                  ga(i,k,itr) = gsoc_hb(ns)
                  fa(i,k,itr) = gsoc_hb(ns)*gsoc_hb(ns)
                else if ( chtrname(itr) == 'BC_HB' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Absorbing aerosols (soot type)
                  tauxar(i,k,itr) = d10e5*uaer(i,k,itr)*ksbc_hb(ns)
                  wa(i,k,itr) = wsbc_hb(ns)
                  ga(i,k,itr) = gsbc_hb(ns)
                  fa(i,k,itr) = gsbc_hb(ns)*gsbc_hb(ns)
                end if
              end do  ! end tracer loop
            end do
          end do
!
!         optical properties for the clear sky diagnostic
!
          do i = 1 , iym1
            do itr = 1 , ntr
              do k = 1 , kz
                utaer(i,itr) = utaer(i,itr) + uaer(i,k,itr)
                tauaer(i,itr) = tauaer(i,itr) + tauxar(i,k,itr)
                waer(i,itr) = waer(i,itr) + wa(i,k,itr)*uaer(i,k,itr)
                gaer(i,itr) = gaer(i,itr) + ga(i,k,itr)*uaer(i,k,itr)
                faer(i,itr) = faer(i,itr) + fa(i,k,itr)*uaer(i,k,itr)
              end do
              if ( utaer(i,itr) <= minimum_utaer ) &
                utaer(i,itr) = minimum_utaer
              waer(i,itr) = waer(i,itr)/utaer(i,itr)
              gaer(i,itr) = gaer(i,itr)/utaer(i,itr)
              faer(i,itr) = faer(i,itr)/utaer(i,itr)
            end do
          end do
!
!         Calculate the EXTERNAL Mixing of aerosols
!         melange externe
!
          do i = 1 , iym1
            do itr = 1 , ntr
!             only for climatic feedback allowed
              do k = 0 , kz
                tauxar_mix(i,k,ns) = tauxar_mix(i,k,ns)               &
                                    + tauxar(i,k,itr)
                tauasc_mix(i,k,ns) = tauasc_mix(i,k,ns)               &
                                    + tauxar(i,k,itr)*wa(i,k,itr)
                gtota_mix(i,k,ns) = gtota_mix(i,k,ns) + ga(i,k,itr)   &
                                   *tauxar(i,k,itr)*wa(i,k,itr)
                ftota_mix(i,k,ns) = ftota_mix(i,k,ns) + fa(i,k,itr)   &
                                   *tauxar(i,k,itr)*wa(i,k,itr)
              end do
!
!             Clear sky (always calcuated if idirect >=1 for
!             diagnostic radiative forcing)
!
              tauxar_mix_cs(i,ns) = tauxar_mix_cs(i,ns)               &
                                   + tauaer(i,itr)
              if (waer(i,itr) > minimum_waer) then
                tauasc_mix_cs(i,ns) = tauasc_mix_cs(i,ns) +           &
                                   tauaer(i,itr)*waer(i,itr)
              end if
              if (gaer(i,itr) > minimum_gaer .and.  &
                  waer(i,itr) > minimum_gaer) then
                gtota_mix_cs(i,ns) = gtota_mix_cs(i,ns) + gaer(i,itr) * &
                                  tauaer(i,itr)*waer(i,itr)
                ftota_mix_cs(i,ns) = ftota_mix_cs(i,ns) + faer(i,itr) * &
                                  tauaer(i,itr)*waer(i,itr)
              end if
            end do
          end do
!
        end do ! end spectral loop
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!         Option 2  melange interne
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      else if ( mixtype == 2 ) then
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Spectral loop
!
        do ns = 1 , nspi

          tauxar_mix_cs(:,ns) = d_zero
          tauasc_mix_cs(:,ns) = d_zero
          gtota_mix_cs(:,ns)  = d_zero
          ftota_mix_cs(:,ns)  = d_zero
          tauxar_mix(:,:,ns)  = d_zero
          tauasc_mix(:,:,ns)  = d_zero
          gtota_mix(:,:,ns)   = d_zero
          ftota_mix(:,:,ns)   = d_zero
! 
          uaer(:,0,:)   = d_zero
          tauxar(:,0,:) = d_zero
          wa(:,0,:)     = d_zero
          ga(:,0,:)     = d_zero
          fa(:,0,:)     = d_zero
          utaer(:,:)    = d_zero
          tauaer(:,:)   = d_zero
          waer(:,:)     = d_zero
          gaer(:,:)     = d_zero
          faer(:,:)     = d_zero
!
!         calculate optical properties of each aerosol component
!
          do i = 1 , iym1
            do k = 1 , kz
              path = (pint(i,k+1)-pint(i,k))*regravgts
              if ( rh(i,k) < d_zero .or. rh(i,k) > d_one ) &
                write ( 6,* ) 'WARNING RH : ' , i , k , rh(i,k)
!             sum of hydrophilic aerosols
              aervtot(i,k) = d_zero
              aermtot(i,k) = d_zero
!
              if ( iso4 /= 0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iso4)/rhoso4
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iso4)
              end if
! 
              if ( ibchl /= 0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,ibchl)/rhobc
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,ibchl)
              end if
!
              if ( iochl /= 0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iochl)/rhooc
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iochl)
              end if
! 
              if ( idust(1) /= 0 ) then
                aervtot(i,k) = aervtot(i,k) + &
                               aermmr(i,k,idust(1))/rhodust
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,idust(1))
              end if
!             minimum quantity of total aerosol
              if ( aermtot(i,k) > minimum_aerosol ) then
!               indexes in the internal mixing table
                prop(1) = (aermmr(i,k,iso4)/rhoso4)/aervtot(i,k)
                prop(2) = (aermmr(i,k,ibchl)/rhobc)/aervtot(i,k)
                prop(3) = (aermmr(i,k,iochl)/rhooc)/aervtot(i,k)
                prop(4) = (aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
!               FIND THE GREATEST FRACTIONAL PART
                if ( iso4 /= 0 ) then
                  i1 = idnint(d_10*prop(1)) + 1
                else
                  i1 = 0 + 1
                end if
                if ( ibchl /= 0 ) then
                  i2 = idnint(d_10*prop(2)) + 1
                else
                  i2 = 0 + 1
                end if
                if ( iochl /= 0 ) then
                  i3 = idnint(d_10*prop(3)) + 1
                else
                  i3 = 0 + 1
                end if
                if ( idust(1) /= 0 ) then
                  i4 = idnint(d_10*prop(4)) + 1
                else
                  i4 = 0 + 1
                end if
!
!               final optical parameters
!
                if ( i1+i2+i3+i4 == 13 ) i4 = i4 + 1
                if ( i1+i2+i3+i4 == 15 ) then
                  if ( i4 /= 1 ) i4 = i4 - 1
                end if
! 
                if ( i1+i2+i3+i4 == 15 ) then
                  if ( i1 /= 1 ) i1 = i1 - 1
                end if
! 
                if ( i1+i2+i3+i4 == 15 ) then
                  if ( i3 /= 1 ) i3 = i3 - 1
                end if
! 
                if ( i1+i2+i3+i4 == 15 ) then
                  call fatal(__FILE__,__LINE__,&
                    'WRONG COMBINATION. SHOULD NEVER HAPPEN')
                end if
! 
                if ( i1+i2+i3+i4 /= 14 ) then
                  print * , i1 , i2 , i3 , i4 , i1 + i2 + i3 + i4
                  print * , idust(1) , iochl , ibchl , iso4
                  print * , 'OC HL' , aermmr(i,k,iochl)/rhooc
                  print * , 'BC HL' , aermmr(i,k,ibchl)/rhobc
                  print * , 'SO4' , aermmr(i,k,iso4)/rhoso4
                  print * , 'DUST' , aermmr(i,k,idust(1))/rhodust
                  print * , 'VOL TOT' , aervtot(i,k)
                  print * , 'OC HL%' , d_10*(aermmr(i,k,iochl)/rhooc) &
                       /aervtot(i,k)
                  print * , 'BC HL%' , d_10*(aermmr(i,k,ibchl)/rhobc) &
                       /aervtot(i,k)
                  print * , 'SO4 %' , d_10*(aermmr(i,k,iso4)/rhoso4)  &
                       /aervtot(i,k)
                  print * , 'SO4 %' ,                              &
                       idnint(d_10*(aermmr(i,k,iso4)/rhoso4)/      &
                       aervtot(i,k))
                  print * , 'DUST %' ,                                &
                       d_10*(aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
                  print * , 'DUST %' ,                             &
                       idnint(d_10*(aermmr(i,k,idust(1))/rhodust)  &
                       /aervtot(i,k))
                  call fatal(__FILE__,__LINE__,                       &
                            'SOMETHING WRONG ON SPECIES ABUNDANCE')
                end if
!
                tauxar_mix(i,k,ns) = dextmix(1,ns,i4,i2,i3,i1)     &
                                    *aermtot(i,k)*path*d10e5
                tauasc_mix(i,k,ns) = dssamix(1,ns,i4,i2,i3,i1)     &
                                    *tauxar_mix(i,k,ns)
                gtota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)        &
                                   *tauasc_mix(i,k,ns)             &
                                   *tauxar_mix(i,k,ns)
                ftota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)        &
                                   *dgmix(1,ns,i4,i2,i3,i1)        &
                                   *tauasc_mix(i,k,ns)             &
                                   *tauxar_mix(i,k,ns)
! 
!               clear sky dignostic
                utaer(i,1) = utaer(i,1) + aermtot(i,k)*path
 
                tauaer(i,1) = tauaer(i,1) + dextmix(1,ns,i4,i2,i3,i1) &
                             *aermtot(i,k)*path*d10e5
                waer(i,1) = waer(i,1) + dssamix(1,ns,i4,i2,i3,i1)     &
                           *aermtot(i,k)*path
                gaer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                           *aermtot(i,k)*path
                faer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                           *dgmix(1,ns,i4,i2,i3,i1)*aermtot(i,k)*path
! 
              end if ! end minimum concentration conditions
            end do ! end k loop
! 
            if ( utaer(i,1) > minimum_utaer ) then
              waer(i,1) = waer(i,1)/utaer(i,1)
              gaer(i,1) = gaer(i,1)/utaer(i,1)
              faer(i,1) = faer(i,1)/utaer(i,1)
            end if
!           clear sky final effective optical properties
            tauxar_mix_cs(i,ns) = tauaer(i,1)
            tauasc_mix_cs(i,ns) = waer(i,1)*tauaer(i,1)
            gtota_mix_cs(i,ns) = gaer(i,1)*waer(i,1)*tauaer(i,1)
            ftota_mix_cs(i,ns) = faer(i,1)*waer(i,1)*tauaer(i,1)
          end do ! end i loop
        end do ! end spectral loop
!
      else
        write ( 6,* ) 'MIXTYPE = ', mixtype
        call fatal(__FILE__,__LINE__,'UNSUPPORTED MIXTYPE IN AEROPPT')
      end if
!
! FAB 
! DUST LW emissivity 
! qabslw = absorption coefficient between k1 and  k2 (m2.g-1) in the LW : 
      qabslw = d_r10
!     initialisation à 1 = perfect transmittivity
      aerlwtr (:,:,:) = d_one
!
      if ( idirect >= 1 ) then
!
        do k1 = 1 , kzp1
          do k2 = 1 , kzp1
            do i = 1 , iym1
              if ( k1==k2 ) aerlwtr(i,k1,k2) = d_one
!             aerosol path btw k1 and k2 flux level
              ibin = 0
              uaerdust = d_zero
              do itr = 1 , ntr     
                if ( chtrname(itr) == 'DUST' ) then
                  ibin = ibin+1
                  if ( k1<k2 ) then
                    uaerdust =  uaerdust + d10e5 *         &
                               (sum(uaer(i,k1:k2-1,itr)))
                    aerlwtr(i,k1,k2) = &
                             dexp(-fiveothree * qabslw * uaerdust)
                  else if ( k1>k2 ) then
                    uaerdust =  uaerdust + d10e5 *         &
                               (sum(uaer(i,k2:k1-1,itr)))
                    aerlwtr(i,k1,k2) = &
                             dexp(-fiveothree * qabslw * uaerdust)
                  end if
                end if
              end do
            end do
          end do
        end do
      end if
!     
      end subroutine aeroppt
!
! SUBROUTINE AEROUT
!
      subroutine aerout(jslc,aeradfo,aeradfos,aerlwfo,aerlwfos)
!
      implicit none
!
      integer :: jslc
      real(8) , dimension(iym1) :: aeradfo , aeradfos, aerlwfo ,        &
                                  aerlwfos
      intent (in) aeradfo , aeradfos, aerlwfo , aerlwfos
!
      real(8) :: rntim
      integer :: i , k
! 
      do k = 1 , kz
        do i = 2 , iym1
#ifdef MPP1
          aerext(i-1,k,jslc) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc) = gtota_mix(i,k,8)
#else
#ifdef BAND
          aerext(i-1,k,jslc) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc) = gtota_mix(i,k,8)
#else
          aerext(i-1,k,jslc-1) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc-1) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc-1) = gtota_mix(i,k,8)
#endif
#endif
        end do
      end do
!
!     CARE :Average the radiative forcing between chem output time
!     steps (in hour) according to radfrq (in min), aertarf is reset to
!     0 at each chem output (cf output.f)
!
      rntim = minph*(chemfrq/radfrq)
!
!     aersol radative forcing (care cgs to mks after radiation scheme !)
!
      do i = 2 , iym1
#ifdef MPP1
        aertarf(i-1,jslc) = aertarf(i-1,jslc) +            &
                            aeradfo(i)*d_r1000/rntim
        aersrrf(i-1,jslc) = aersrrf(i-1,jslc) +            &
                            aeradfos(i)*d_r1000/rntim
        aertalwrf(i-1,jslc) = aertalwrf(i-1,jslc) +        &
                             aerlwfo(i)*d_r1000/rntim
        aersrlwrf(i-1,jslc) = aersrlwrf(i-1,jslc) +        &
                             aerlwfos(i)*d_r1000/rntim
#else
#ifdef BAND
        aertarf(i-1,jslc) = aertarf(i-1,jslc) + aeradfo(i)     &
                             *d_r1000/rntim
        aersrrf(i-1,jslc) = aersrrf(i-1,jslc) + aeradfos(i)    &
                             *d_r1000/rntim
        aertalwrf(i-1,jslc) = aertalwrf(i-1,jslc) +            &
                               aerlwfo(i) * d_r1000/rntim
        aersrlwrf(i-1,jslc) = aersrlwrf(i-1,jslc) +            &
                               aerlwfos(i) * d_r1000/rntim
#else
        aertarf(i-1,jslc-1) = aertarf(i-1,jslc-1) + aeradfo(i)   &
                             *d_r1000/rntim
        aersrrf(i-1,jslc-1) = aersrrf(i-1,jslc-1) + aeradfos(i)  &
                             *d_r1000/rntim
        aertalwrf(i-1,jslc-1) = aertalwrf(i-1,jslc-1) +          &
                               aerlwfo(i) * d_r1000/rntim
        aersrlwrf(i-1,jslc-1) = aersrlwrf(i-1,jslc-1) +          &
                               aerlwfos(i) * d_r1000/rntim
#endif
#endif
      end do
! 
      end subroutine aerout
!
      end module mod_aerosol
