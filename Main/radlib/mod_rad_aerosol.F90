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

module mod_rad_aerosol
!
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_memutil
  use mod_mpmessage
  use mod_rad_common
!
  private
!
  public :: tauxar3d , tauasc3d , gtota3d , ftota3d
  public :: tauxar , tauasc , gtota , ftota
  public :: aertrlw
  public :: allocate_mod_rad_aerosol , aermix , aeroppt
!
! MODIF 16/09/2005 IBRAH Internal mixing
!
  integer , parameter :: ncoefs = 5  ! Number of coefficients
  integer , parameter :: nwav = 19
  integer , parameter :: nih = 8
!
  real(dp) , parameter :: d10e5  = 1.0D+05
  real(dp) , parameter :: d10e4  = 1.0D+04
  real(dp) , parameter :: minimum_aerosol = 1.0D-14
  real(dp) , parameter :: minimum_utaer   = 1.0D-10
  real(dp) , parameter :: minimum_waer   = 1.0D-30
  real(dp) , parameter :: minimum_gaer   = 1.0D-20
  real(dp) , parameter :: fiveothree  = d_five/d_three
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
  real(dp) , dimension(nspi) :: gsbase , gsbc_hb , gsbc_hl , gsoc_hb , &
            gsoc_hl , ksbase , ksbc_hb , ksbc_hl , ksoc_hb , ksoc_hl , &
            wsbase , wsbc_hb , wsbc_hl , wsoc_hb , wsoc_hl
  real(dp) , dimension(nspi,ncoefs) :: gscoef , kscoef , wscoef
  real(dp) , dimension(nwav,2,nih) :: ksslt , wsslt , gsslt
  real(dp) , dimension(nspi,4) :: gsdust , ksdust , wsdust
  real(dp) , dimension(nspi,2) :: gssslt , kssslt , wssslt
  !
  real(dp) , dimension(8) :: rhp
!
! Aerosol mass mixing ratio
!
  real(dp) , pointer , dimension(:,:,:) :: aermm
!
! Background aerosol mass mixing ratio
!
  real(dp) , pointer , dimension(:,:) :: aermmb
!
! Aerosol optical properties (for the mixing) 
!
  real(dp) , pointer , dimension(:,:,:) :: ftota3d ,   &
                 gtota3d , tauasc3d , tauxar3d
!
  real(dp) , pointer , dimension(:,:) :: ftota ,  &
         gtota , tauasc , tauxar
!
! Work arrays for aeroppt (aerosol individual optical properties SW)
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
!
  real(dp) , pointer , dimension(:,:) :: aermtot , aervtot
  real(dp) , pointer , dimension(:,:,:) :: fa , ga , tx , uaer , wa
  real(dp) , pointer , dimension(:,:) :: faer , gaer , tauaer , utaer , waer
  real(dp) , dimension(4) :: prop
  integer :: ll , mm , nn
  integer :: npoints
!
! Aersol LW optical properties
!
  real(dp) , pointer , dimension(:,:,:) ::  aertrlw 
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
  data (((ksslt(nn,ll,mm),ll=1,2),nn=1,nwav),mm=1,nih) /                    &
    1.10670D0 , 0.24350D0 , 1.13470D0 , 0.24490D0 , 1.14490D0 , 0.24530D0 , &
    1.15360D0 , 0.24570D0 , 1.16140D0 , 0.24610D0 , 1.16810D0 , 0.24640D0 , &
    1.18630D0 , 0.24730D0 , 1.26130D0 , 0.25120D0 , 1.28830D0 , 0.25610D0 , &
    1.11740D0 , 0.26480D0 , 1.07340D0 , 0.26760D0 , 1.01450D0 , 0.27070D0 , &
    0.93780D0 , 0.27400D0 , 0.88700D0 , 0.27610D0 , 0.85710D0 , 0.27750D0 , &
    0.77260D0 , 0.28170D0 , 0.55190D0 , 0.29620D0 , 0.23690D0 , 0.31880D0 , &
    0.23690D0 , 0.31880D0 , 2.76050D0 , 0.62500D0 , 2.79270D0 , 0.62740D0 , &
    2.80620D0 , 0.62820D0 , 2.81920D0 , 0.62900D0 , 2.83250D0 , 0.62970D0 , &
    2.84560D0 , 0.63040D0 , 2.88990D0 , 0.63220D0 , 3.08680D0 , 0.64020D0 , &
    3.20500D0 , 0.64830D0 , 2.99350D0 , 0.66540D0 , 2.91440D0 , 0.67000D0 , &
    2.83310D0 , 0.67460D0 , 2.69720D0 , 0.68040D0 , 2.60410D0 , 0.68440D0 , &
    2.53200D0 , 0.68740D0 , 2.31690D0 , 0.69580D0 , 1.74510D0 , 0.71490D0 , &
    0.96070D0 , 0.77950D0 , 0.96070D0 , 0.77950D0 , 3.43580D0 , 0.78970D0 , &
    3.48660D0 , 0.79210D0 , 3.50620D0 , 0.79300D0 , 3.52400D0 , 0.79370D0 , &
    3.54090D0 , 0.79440D0 , 3.55630D0 , 0.79510D0 , 3.60340D0 , 0.79700D0 , &
    3.84490D0 , 0.80600D0 , 4.01160D0 , 0.81490D0 , 3.83520D0 , 0.83560D0 , &
    3.75150D0 , 0.84090D0 , 3.67460D0 , 0.84630D0 , 3.53150D0 , 0.85240D0 , &
    3.43520D0 , 0.85700D0 , 3.35000D0 , 0.86030D0 , 3.09100D0 , 0.86980D0 , &
    2.36570D0 , 0.89080D0 , 1.39570D0 , 0.96980D0 , 1.39570D0 , 0.96980D0 , &
    4.13430D0 , 0.95220D0 , 4.18210D0 , 0.95690D0 , 4.20170D0 , 0.95850D0 , &
    4.21980D0 , 0.95980D0 , 4.23730D0 , 0.96080D0 , 4.25330D0 , 0.96160D0 , &
    4.30420D0 , 0.96320D0 , 4.59390D0 , 0.97370D0 , 4.80120D0 , 0.98300D0 , &
    4.68220D0 , 1.00670D0 , 4.59660D0 , 1.01270D0 , 4.52830D0 , 1.01870D0 , &
    4.38560D0 , 1.02600D0 , 4.29240D0 , 1.03140D0 , 4.19710D0 , 1.03510D0 , &
    3.90190D0 , 1.04550D0 , 3.03000D0 , 1.06800D0 , 1.89030D0 , 1.15860D0 , &
    1.89030D0 , 1.15860D0 , 5.84240D0 , 1.36500D0 , 5.88300D0 , 1.36820D0 , &
    5.90140D0 , 1.36940D0 , 5.91970D0 , 1.37050D0 , 5.93840D0 , 1.37140D0 , &
    5.95670D0 , 1.37230D0 , 6.02000D0 , 1.37490D0 , 6.38780D0 , 1.38870D0 , &
    6.69990D0 , 1.40120D0 , 6.73060D0 , 1.43000D0 , 6.65350D0 , 1.43720D0 , &
    6.61700D0 , 1.44470D0 , 6.49700D0 , 1.45260D0 , 6.42870D0 , 1.45840D0 , &
    6.31860D0 , 1.46280D0 , 5.96400D0 , 1.47540D0 , 4.77770D0 , 1.50540D0 , &
    3.29640D0 , 1.61800D0 , 3.29640D0 , 1.61800D0 , 8.47370D0 , 2.03340D0 , &
    8.57490D0 , 2.03540D0 , 8.61230D0 , 2.03640D0 , 8.64390D0 , 2.03730D0 , &
    8.67140D0 , 2.03810D0 , 8.69340D0 , 2.03880D0 , 8.74920D0 , 2.04110D0 , &
    9.21250D0 , 2.06220D0 , 9.61140D0 , 2.07770D0 , 9.94990D0 , 2.11370D0 , &
    9.90040D0 , 2.12310D0 , 9.92670D0 , 2.13290D0 , 9.87340D0 , 2.14430D0 , &
    9.86950D0 , 2.15220D0 , 9.75570D0 , 2.15760D0 , 9.36830D0 , 2.17230D0 , &
    7.80820D0 , 2.20780D0 , 5.96950D0 , 2.34700D0 , 5.96950D0 , 2.34700D0 , &
   14.68840D0 , 3.60330D0 ,14.82000D0 , 3.60900D0 ,14.86800D0 , 3.61130D0 , &
   14.90890D0 , 3.61330D0 ,14.94520D0 , 3.61520D0 ,14.97530D0 , 3.61700D0 , &
   15.05470D0 , 3.62240D0 ,15.63030D0 , 3.64850D0 ,16.16270D0 , 3.66570D0 , &
   17.07870D0 , 3.72400D0 ,17.13010D0 , 3.73710D0 ,17.29570D0 , 3.75280D0 , &
   17.43470D0 , 3.76940D0 ,17.60120D0 , 3.78100D0 ,17.52590D0 , 3.78800D0 , &
   17.21790D0 , 3.80740D0 ,15.24990D0 , 3.85640D0 ,13.22780D0 , 4.03860D0 , &
   13.22780D0 , 4.03860D0 ,22.48880D0 , 5.62270D0 ,22.62770D0 , 5.61380D0 , &
   22.67970D0 , 5.61090D0 ,22.72570D0 , 5.60950D0 ,22.76830D0 , 5.60980D0 , &
   22.80590D0 , 5.61170D0 ,22.91520D0 , 5.62400D0 ,23.61060D0 , 5.66380D0 , &
   24.24360D0 , 5.69990D0 ,25.74500D0 , 5.76950D0 ,25.92050D0 , 5.78650D0 , &
   26.22020D0 , 5.80620D0 ,26.57360D0 , 5.82770D0 ,26.92400D0 , 5.84240D0 , &
   26.92910D0 , 5.85120D0 ,26.83940D0 , 5.87590D0 ,24.85840D0 , 5.94260D0 , &
   23.35360D0 , 6.16670D0 ,23.35360D0 , 6.16670D0/

  data (((wsslt(nn,ll,mm),ll=1,2),nn=1,nwav),mm=1,nih) /                    & 
    0.99970D0 , 0.99780D0 , 0.99980D0 , 0.99870D0 , 0.99980D0 , 0.99900D0 , &
    0.99990D0 , 0.99920D0 , 0.99990D0 , 0.99940D0 , 0.99990D0 , 0.99950D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    0.99780D0 , 0.98760D0 , 0.99700D0 , 0.98450D0 , 0.99480D0 , 0.97590D0 , &
    0.99260D0 , 0.96710D0 , 0.99250D0 , 0.96520D0 , 0.99110D0 , 0.96080D0 , &
    0.98690D0 , 0.94790D0 , 0.95600D0 , 0.85510D0 , 0.98970D0 , 0.97820D0 , &
    0.98970D0 , 0.97820D0 , 0.99980D0 , 0.99930D0 , 0.99990D0 , 0.99960D0 , &
    1.00000D0 , 0.99960D0 , 1.00000D0 , 0.99970D0 , 1.00000D0 , 0.99980D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 0.99890D0 , 0.99100D0 , 0.99610D0 , 0.98380D0 , &
    0.97580D0 , 0.96040D0 , 0.96060D0 , 0.94340D0 , 0.96860D0 , 0.94900D0 , &
    0.96730D0 , 0.94590D0 , 0.95970D0 , 0.93470D0 , 0.70640D0 , 0.71750D0 , &
    0.93370D0 , 0.85470D0 , 0.93370D0 , 0.85470D0 , 0.99980D0 , 0.99940D0 , &
    0.99990D0 , 0.99970D0 , 1.00000D0 , 0.99970D0 , 1.00000D0 , 0.99980D0 , &
    1.00000D0 , 0.99990D0 , 1.00000D0 , 0.99990D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 0.99900D0 , 0.99100D0 , &
    0.99600D0 , 0.98320D0 , 0.97510D0 , 0.95980D0 , 0.95960D0 , 0.94280D0 , &
    0.96790D0 , 0.94820D0 , 0.96670D0 , 0.94510D0 , 0.95920D0 , 0.93370D0 , &
    0.70160D0 , 0.71770D0 , 0.92940D0 , 0.83360D0 , 0.92940D0 , 0.83360D0 , &
    1.00000D0 , 0.99940D0 , 1.00000D0 , 0.99970D0 , 1.00000D0 , 0.99970D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , 1.00000D0 , 0.99990D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    0.99900D0 , 0.99070D0 , 0.99590D0 , 0.98250D0 , 0.97480D0 , 0.95920D0 , &
    0.95920D0 , 0.94210D0 , 0.96760D0 , 0.94730D0 , 0.96640D0 , 0.94410D0 , &
    0.95910D0 , 0.93230D0 , 0.70100D0 , 0.71720D0 , 0.92720D0 , 0.81720D0 , &
    0.92720D0 , 0.81720D0 , 1.00000D0 , 0.99960D0 , 1.00000D0 , 0.99970D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , &
    1.00000D0 , 0.99990D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 0.99900D0 , 0.99010D0 , 0.99580D0 , 0.98100D0 , &
    0.97440D0 , 0.95790D0 , 0.95880D0 , 0.94050D0 , 0.96730D0 , 0.94520D0 , &
    0.96620D0 , 0.94170D0 , 0.95920D0 , 0.92880D0 , 0.70350D0 , 0.71500D0 , &
    0.92460D0 , 0.78760D0 , 0.92460D0 , 0.78760D0 , 1.00000D0 , 0.99970D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , &
    1.00000D0 , 0.99990D0 , 1.00000D0 , 0.99990D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 0.99890D0 , 0.98890D0 , &
    0.99550D0 , 0.97880D0 , 0.97400D0 , 0.95610D0 , 0.95860D0 , 0.93790D0 , &
    0.96700D0 , 0.94180D0 , 0.96600D0 , 0.93770D0 , 0.95930D0 , 0.92320D0 , &
    0.70870D0 , 0.71080D0 , 0.92240D0 , 0.75570D0 , 0.92240D0 , 0.75570D0 , &
    1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    0.99860D0 , 0.98660D0 , 0.99500D0 , 0.97530D0 , 0.97330D0 , 0.95280D0 , &
    0.95800D0 , 0.93350D0 , 0.96630D0 , 0.93590D0 , 0.96530D0 , 0.93080D0 , &
    0.95900D0 , 0.91320D0 , 0.71670D0 , 0.70160D0 , 0.91790D0 , 0.71290D0 , &
    0.91790D0 , 0.71290D0 , 1.00000D0 , 0.99980D0 , 1.00000D0 , 0.99990D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , 1.00000D0 , &
    1.00000D0 , 1.00000D0 , 0.99830D0 , 0.98460D0 , 0.99430D0 , 0.97210D0 , &
    0.97240D0 , 0.94980D0 , 0.95720D0 , 0.92930D0 , 0.96540D0 , 0.93050D0 , &
    0.96440D0 , 0.92440D0 , 0.95820D0 , 0.90420D0 , 0.72210D0 , 0.69280D0 , &
    0.91180D0 , 0.68270D0 , 0.91180D0 , 0.68270D0/

  data (((gsslt(nn,ll,mm),ll=1,2), nn=1,nwav),mm=1,nih) /                   &
    0.73080D0 , 0.80840D0 , 0.71820D0 , 0.81040D0 , 0.71420D0 , 0.81080D0 , &
    0.71090D0 , 0.81090D0 , 0.70820D0 , 0.81060D0 , 0.70600D0 , 0.81000D0 , &
    0.70120D0 , 0.80630D0 , 0.69510D0 , 0.79950D0 , 0.69590D0 , 0.79020D0 , &
    0.70140D0 , 0.77760D0 , 0.69940D0 , 0.77390D0 , 0.69710D0 , 0.76920D0 , &
    0.69470D0 , 0.76430D0 , 0.69600D0 , 0.76180D0 , 0.69530D0 , 0.76130D0 , &
    0.69220D0 , 0.75980D0 , 0.64050D0 , 0.74180D0 , 0.59490D0 , 0.70590D0 , &
    0.59490D0 , 0.70590D0 , 0.78530D0 , 0.84490D0 , 0.78490D0 , 0.84630D0 , &
    0.78460D0 , 0.84680D0 , 0.78400D0 , 0.84720D0 , 0.78310D0 , 0.84760D0 , &
    0.78200D0 , 0.84800D0 , 0.77700D0 , 0.84930D0 , 0.77180D0 , 0.84890D0 , &
    0.77310D0 , 0.84390D0 , 0.77860D0 , 0.83450D0 , 0.77770D0 , 0.83330D0 , &
    0.77760D0 , 0.83590D0 , 0.77890D0 , 0.83640D0 , 0.78110D0 , 0.83260D0 , &
    0.78220D0 , 0.83280D0 , 0.78470D0 , 0.83430D0 , 0.77590D0 , 0.88690D0 , &
    0.70530D0 , 0.80650D0 , 0.70530D0 , 0.80650D0 , 0.80300D0 , 0.84880D0 , &
    0.79650D0 , 0.84790D0 , 0.79440D0 , 0.84790D0 , 0.79260D0 , 0.84810D0 , &
    0.79120D0 , 0.84850D0 , 0.79000D0 , 0.84900D0 , 0.78750D0 , 0.85180D0 , &
    0.78010D0 , 0.85380D0 , 0.77990D0 , 0.84840D0 , 0.78650D0 , 0.84180D0 , &
    0.78580D0 , 0.84050D0 , 0.78610D0 , 0.84320D0 , 0.78810D0 , 0.84410D0 , &
    0.79060D0 , 0.84040D0 , 0.79200D0 , 0.84070D0 , 0.79570D0 , 0.84230D0 , &
    0.79390D0 , 0.89370D0 , 0.72450D0 , 0.81840D0 , 0.72450D0 , 0.81840D0 , &
    0.80440D0 , 0.83590D0 , 0.80240D0 , 0.85000D0 , 0.80150D0 , 0.85400D0 , &
    0.80070D0 , 0.85670D0 , 0.79980D0 , 0.85840D0 , 0.79900D0 , 0.85900D0 , &
    0.79620D0 , 0.85590D0 , 0.78600D0 , 0.85570D0 , 0.78510D0 , 0.85480D0 , &
    0.79040D0 , 0.84740D0 , 0.79010D0 , 0.84670D0 , 0.79080D0 , 0.84990D0 , &
    0.79340D0 , 0.85060D0 , 0.79610D0 , 0.84700D0 , 0.79780D0 , 0.84720D0 , &
    0.80260D0 , 0.84880D0 , 0.80620D0 , 0.89840D0 , 0.73900D0 , 0.82710D0 , &
    0.73900D0 , 0.82710D0 , 0.81710D0 , 0.84450D0 , 0.81550D0 , 0.85040D0 , &
    0.81470D0 , 0.85230D0 , 0.81380D0 , 0.85370D0 , 0.81290D0 , 0.85460D0 , &
    0.81200D0 , 0.85500D0 , 0.80870D0 , 0.85420D0 , 0.79550D0 , 0.86100D0 , &
    0.79060D0 , 0.85680D0 , 0.79590D0 , 0.85560D0 , 0.79590D0 , 0.85500D0 , &
    0.79740D0 , 0.85830D0 , 0.80080D0 , 0.86020D0 , 0.80350D0 , 0.85730D0 , &
    0.80560D0 , 0.85770D0 , 0.81190D0 , 0.85960D0 , 0.82680D0 , 0.90580D0 , &
    0.76040D0 , 0.84160D0 , 0.76040D0 , 0.84160D0 , 0.82960D0 , 0.83970D0 , &
    0.82610D0 , 0.84290D0 , 0.82470D0 , 0.84430D0 , 0.82360D0 , 0.84560D0 , &
    0.82270D0 , 0.84680D0 , 0.82200D0 , 0.84800D0 , 0.82050D0 , 0.85190D0 , &
    0.80440D0 , 0.86460D0 , 0.79840D0 , 0.86330D0 , 0.79950D0 , 0.86230D0 , &
    0.79980D0 , 0.86220D0 , 0.80170D0 , 0.86620D0 , 0.80580D0 , 0.86830D0 , &
    0.80820D0 , 0.86600D0 , 0.81070D0 , 0.86650D0 , 0.81830D0 , 0.86880D0 , &
    0.84470D0 , 0.91370D0 , 0.78000D0 , 0.85830D0 , 0.78000D0 , 0.85830D0 , &
    0.83790D0 , 0.83920D0 , 0.83800D0 , 0.83790D0 , 0.83780D0 , 0.83810D0 , &
    0.83740D0 , 0.83880D0 , 0.83680D0 , 0.84020D0 , 0.83600D0 , 0.84210D0 , &
    0.83240D0 , 0.85090D0 , 0.81740D0 , 0.86400D0 , 0.80800D0 , 0.87030D0 , &
    0.80450D0 , 0.86700D0 , 0.80400D0 , 0.86850D0 , 0.80560D0 , 0.87300D0 , &
    0.80970D0 , 0.87640D0 , 0.81120D0 , 0.87540D0 , 0.81390D0 , 0.87650D0 , &
    0.82260D0 , 0.88010D0 , 0.86370D0 , 0.92410D0 , 0.79820D0 , 0.88210D0 , &
    0.79820D0 , 0.88210D0 , 0.84460D0 , 0.80750D0 , 0.84500D0 , 0.82000D0 , &
    0.84500D0 , 0.82460D0 , 0.84480D0 , 0.82850D0 , 0.84450D0 , 0.83200D0 , &
    0.84400D0 , 0.83500D0 , 0.84160D0 , 0.84330D0 , 0.82800D0 , 0.86270D0 , &
    0.81760D0 , 0.87020D0 , 0.80900D0 , 0.87030D0 , 0.80800D0 , 0.87240D0 , &
    0.80920D0 , 0.87730D0 , 0.81270D0 , 0.88150D0 , 0.81290D0 , 0.88140D0 , &
    0.81550D0 , 0.88280D0 , 0.82420D0 , 0.88730D0 , 0.87450D0 , 0.93080D0 , &
    0.80830D0 , 0.90060D0 , 0.80830D0 , 0.90060D0/
!
  data rhp /0.0D0,0.5D0,0.7D0,0.8D0,0.9D0,0.95D0,0.98D0,0.99D0/
!
  contains
! 
  subroutine allocate_mod_rad_aerosol(ichem)
    implicit none
    integer , intent(in) :: ichem
    npoints = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem2d(aermmb,1,npoints,1,kz,'aerosol:aermmb')
    call getmem3d(ftota3d,1,npoints,0,kz,1,nspi,'aerosol:ftota3d')
    call getmem3d(gtota3d,1,npoints,0,kz,1,nspi,'aerosol:gtota3d')
    call getmem3d(tauasc3d,1,npoints,0,kz,1,nspi,'aerosol:tauasc3d')
    call getmem3d(tauxar3d,1,npoints,0,kz,1,nspi,'aerosol:tauxar3d')
    call getmem2d(ftota,1,npoints,1,nspi,'aerosol:ftota')
    call getmem2d(gtota,1,npoints,1,nspi,'aerosol:gtota')
    call getmem2d(tauasc,1,npoints,1,nspi,'aerosol:tauasc')
    call getmem2d(tauxar,1,npoints,1,nspi,'aerosol:tauxar')
    call getmem2d(aermtot,1,npoints,1,kz,'aerosol:aermtot')
    call getmem2d(aervtot,1,npoints,1,kz,'aerosol:aervtot')
    call getmem3d(aertrlw,1,npoints,1,kzp1,1,kzp1,'aerosol:aertrlw')
    if ( ichem == 1 ) then
      call getmem3d(fa,1,npoints,0,kz,1,ntr,'aerosol:fa')
      call getmem3d(ga,1,npoints,0,kz,1,ntr,'aerosol:ga')
      call getmem3d(tx,1,npoints,0,kz,1,ntr,'aerosol:tx')
      call getmem3d(uaer,1,npoints,0,kz,1,ntr,'aerosol:uaer')
      call getmem3d(wa,1,npoints,0,kz,1,ntr,'aerosol:wa')
      call getmem2d(faer,1,npoints,1,ntr,'aerosol:faer')
      call getmem2d(gaer,1,npoints,1,ntr,'aerosol:gaer')
      call getmem2d(tauaer,1,npoints,1,ntr,'aerosol:tauaer')
      call getmem2d(utaer,1,npoints,1,ntr,'aerosol:utaer')
      call getmem2d(waer,1,npoints,1,ntr,'aerosol:waer')
    end if
  end subroutine allocate_mod_rad_aerosol
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
  subroutine aermix(pint,n1,n2)
! 
    implicit none
!
    integer , intent(in) :: n1 , n2
!   Radiation level interface pressures (dynes/cm2)
    real(dp) , intent(in) , pointer , dimension(:,:) :: pint
!
!-----------------------------------------------------------------------
!
! mxaerl - max nmbr aerosol levels counting up from surface
! tauvis - visible optical depth
! kaervs - visible extinction coefficiant of aerosol (m2/g)
! omgvis - visible omega0
! gvis   - visible forward scattering asymmetry parameter
!
!-----------------------------------------------------------------------
! 
    real(dp) :: gvis , kaervs , omgvis , rhfac , tauvis
    integer :: n , k , mxaerl
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
!   Set relative humidity and factor; then aerosol amount:
!
    do k = 1 , kz
      do n = n1 , n2
!
!       Define background aerosol
!       Find constant aerosol mass mixing ratio for specified levels
!       in the column, converting units where appropriate
!       for the moment no more used
!
        if ( k >= kz + 1 - mxaerl ) then
          aermmb(n,k) = egravgts * tauvis / &
                       (d10e4*kaervs*rhfac*(d_one-omgvis*gvis*gvis) * &
                       (pint(n,kzp1)-pint(n,kzp1-mxaerl)))
        else
          aermmb(n,k) = d_zero
        end if
      end do
    end do
!
  end subroutine aermix
!
! SUBROUTINE AEROPPT
!
  subroutine aeroppt(rh,aermmr,pint,n1,n2)
!
    implicit none
!
!   Interface pressure, relative humidity
!
    integer , intent(in) :: n1 , n2
    real(dp) , intent(in) , pointer , dimension(:,:) :: pint
    real(dp) , intent(in) , pointer , dimension(:,:,:) :: aermmr
    real(dp) , intent(in) , pointer , dimension(:,:) :: rh
!
    integer :: n , l , ibin , jbin , itr , k , k1, k2 , ns
    real(dp) :: path , uaerdust , qabslw , rh0
!
    if ( .not. lchem ) then
      tauxar(:,:) = d_zero
      tauasc(:,:) = d_zero
      gtota(:,:) = d_zero
      ftota(:,:) = d_zero
      tauxar3d(:,:,:) = d_zero
      tauasc3d(:,:,:) = d_zero
      gtota3d(:,:,:) = d_zero
      ftota3d(:,:,:) = d_zero
      aertrlw (:,:,:) = d_one
      return
    end if
!
    tx = d_zero
    wa = d_zero
    ga = d_zero
    fa = d_zero
!   
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!   Melange externe
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!   Spectral loop
!
    do ns = 1 , nspi

      tauxar(:,ns) = d_zero
      tauasc(:,ns) = d_zero
      gtota(:,ns) = d_zero
      ftota(:,ns) = d_zero
!
      tauxar3d(:,:,ns) = d_zero
      tauasc3d(:,:,ns) = d_zero
      gtota3d(:,:,ns) = d_zero
      ftota3d(:,:,ns) = d_zero
!
      uaer(:,0,:) = d_zero
      tx(:,0,:) = d_zero
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
!     calculate optical properties of each aerosol component
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      do k = 1 , kz
        do n = n1 , n2
          path = (pint(n,k+1)-pint(n,k))*regravgts
          ibin = 0
          jbin = 0
          do itr = 1 , ntr
            uaer(n,k,itr) = d_zero
            if ( rh(n,k) < d_zero .or. rh(n,k) > d_one ) then
              print * , n , k , rh(n,k) , '  RH WARNING !!!!!'
            end if
            if ( tracname(itr) == 'XXXXX') then
              continue
            end if
            if ( tracname(itr)(1:4) == 'DUST' ) then
              uaer(n,k,itr) = aermmr(n,k,itr)*path
              ibin = ibin + 1
              if ( ibin > 4 ) then
                call fatal(__FILE__,__LINE__,'DUST BINS MAX is 4')
              end if
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksdust(ns,ibin)
              wa(n,k,itr) = wsdust(ns,ibin)
              ga(n,k,itr) = gsdust(ns,ibin)
              fa(n,k,itr) = gsdust(ns,ibin)*gsdust(ns,ibin)
            else if ( tracname(itr) == 'SO4' ) then
              ! maximum limit for effect on sulfate extinction 
              rh0 = dmin1(0.97D0,dmax1(d_zero,rh(n,k)))
              uaer(n,k,itr) = aermmr(n,k,itr)*path
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksbase(ns) *               &
                   dexp(kscoef(ns,1)+kscoef(ns,2)/(rh0+kscoef(ns,3)) + &
                   kscoef(ns,4)/(rh0+kscoef(ns,5)))
              wa(n,k,itr) = d_one - wsbase(ns) * &
                   dexp(wscoef(ns,1)+wscoef(ns,2) / &
                   (rh(n,k)+wscoef(ns,3))+wscoef(ns,4)/(rh(n,k)+wscoef(ns,5)))
              ga(n,k,itr) = gsbase(ns) * dexp(gscoef(ns,1)+gscoef(ns,2) / &
                   (rh(n,k)+gscoef(ns,3))+gscoef(ns,4)/(rh(n,k)+gscoef(ns,5)))
              fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
            else if ( tracname(itr)(1:4) == 'SSLT' ) then
              rh0 = dmin1(0.99D0,dmax1(d_zero,rh(n,k)))
              jbin = jbin+1
              if ( jbin > 2 ) then
                call fatal(__FILE__,__LINE__,'SEA SALT BINS MAX is 2')
              end if
              do l = 1 , 7
                if ( rh0 > rhp(1) .and. rh0 <= rhp(l+1) ) then
! FAB : test according to li et al., ksslt cannot exceed 1.3
! quick fix for now, update parameterisation to LI et al, ACP 2008 in a near future
                  kssslt(ns,jbin) = dmin1(ksslt(ns,jbin,l),1.2D0)
                  gssslt(ns,jbin) = gsslt(ns,jbin,l)
                  wssslt(ns,jbin) = wsslt(ns,jbin,l)
                end if
              end do
              uaer(n,k,itr) = aermmr(n,k,itr)*path
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*kssslt(ns,jbin)
              wa(n,k,itr) = wssslt(ns,jbin)
              ga(n,k,itr) = gssslt(ns,jbin)
              fa(n,k,itr) = gssslt(ns,jbin)*gssslt(ns,jbin)
            else if ( tracname(itr) == 'OC_HL' ) then
              uaer(n,k,itr) = aermmr(n,k,itr)*path
!             Humidity effect !
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksoc_hl(ns) * &
                            (d_one-rh(n,k))**(-0.2D0)
              wa(n,k,itr) = wsoc_hl(ns)
              ga(n,k,itr) = gsoc_hl(ns)
              fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
            else if ( tracname(itr) == 'BC_HL' ) then
              uaer(n,k,itr) = aermmr(n,k,itr)*path
!             Humidity effect !
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksbc_hl(ns) * &
                            (d_one-rh(n,k))**(-0.25D0)
              wa(n,k,itr) = wsbc_hl(ns)
              ga(n,k,itr) = gsbc_hl(ns)
              fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
            else if ( tracname(itr) == 'OC_HB' ) then
              uaer(n,k,itr) = aermmr(n,k,itr)*path
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksoc_hb(ns)
              wa(n,k,itr) = wsoc_hb(ns)
              ga(n,k,itr) = gsoc_hb(ns)
              fa(n,k,itr) = gsoc_hb(ns)*gsoc_hb(ns)
            else if ( tracname(itr) == 'BC_HB' ) then
              uaer(n,k,itr) = aermmr(n,k,itr)*path
!             Absorbing aerosols (soot type)
              tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksbc_hb(ns)
              wa(n,k,itr) = wsbc_hb(ns)
              ga(n,k,itr) = gsbc_hb(ns)
              fa(n,k,itr) = gsbc_hb(ns)*gsbc_hb(ns)
            end if
          end do  ! end tracer loop
        end do
      end do
!
!     optical properties for the clear sky diagnostic
!
      do itr = 1 , ntr
        do k = 1 , kz
          do n = n1 , n2
            utaer(n,itr) = utaer(n,itr) + uaer(n,k,itr)
            tauaer(n,itr) = tauaer(n,itr) + tx(n,k,itr)
            waer(n,itr) = waer(n,itr) + wa(n,k,itr)*uaer(n,k,itr)
            gaer(n,itr) = gaer(n,itr) + ga(n,k,itr)*uaer(n,k,itr)
            faer(n,itr) = faer(n,itr) + fa(n,k,itr)*uaer(n,k,itr)
          end do
        end do
      end do
      do itr = 1 , ntr
        do n = n1 , n2
          if ( utaer(n,itr) <= minimum_utaer ) utaer(n,itr) = minimum_utaer
          waer(n,itr) = waer(n,itr)/utaer(n,itr)
          gaer(n,itr) = gaer(n,itr)/utaer(n,itr)
          faer(n,itr) = faer(n,itr)/utaer(n,itr)
        end do
      end do
!
!     Calculate the EXTERNAL Mixing of aerosols
!     melange externe
!
!     only for climatic feedback allowed
      do itr = 1 , ntr
        do k = 0 , kz
          do n = n1 , n2
            tauxar3d(n,k,ns) = tauxar3d(n,k,ns) + tx(n,k,itr)
            tauasc3d(n,k,ns) = tauasc3d(n,k,ns) + tx(n,k,itr)*wa(n,k,itr)
            gtota3d(n,k,ns) = gtota3d(n,k,ns) + ga(n,k,itr) * &
                                tx(n,k,itr)*wa(n,k,itr)
            ftota3d(n,k,ns) = ftota3d(n,k,ns) + fa(n,k,itr) * &
                                tx(n,k,itr)*wa(n,k,itr)
          end do
        end do
      end do
!
!     Clear sky (always calcuated if ichdir >=1 for
!     diagnostic radiative forcing)
!
      do itr = 1 , ntr
        do n = n1 , n2
          tauxar(n,ns) = tauxar(n,ns) + tauaer(n,itr)
          if (waer(n,itr) > minimum_waer) then
            tauasc(n,ns) = tauasc(n,ns) + tauaer(n,itr)*waer(n,itr)
          end if
          if (gaer(n,itr) > minimum_gaer .and.  &
              waer(n,itr) > minimum_gaer) then
            gtota(n,ns) = gtota(n,ns) + gaer(n,itr) * &
                              tauaer(n,itr)*waer(n,itr)
            ftota(n,ns) = ftota(n,ns) + faer(n,itr) * &
                              tauaer(n,itr)*waer(n,itr)
          end if
        end do
      end do
!
    end do ! end spectral loop
!
!   FAB 
!   DUST LW emissivity 
!   qabslw = absorption coefficient between k1 and  k2 (m2.g-1) in the LW : 
    qabslw = d_r10
!   initialisation à 1 = perfect transmittivity
    aertrlw (:,:,:) = d_one
!
    do k1 = 1 , kzp1
      do k2 = 1 , kzp1
        do n = n1 , n2
          if ( k1==k2 ) aertrlw(n,k1,k2) = d_one
!         aerosol path btw k1 and k2 flux level
          ibin = 0
          uaerdust = d_zero
          do itr = 1 , ntr     
            if ( tracname(itr) == 'DUST' ) then
              ibin = ibin+1
              if ( k1<k2 ) then
                uaerdust =  uaerdust + d10e5 * (sum(uaer(n,k1:k2-1,itr)))
                aertrlw(n,k1,k2) = dexp(-fiveothree * qabslw * uaerdust)
              else if ( k1>k2 ) then
                uaerdust =  uaerdust + d10e5 * (sum(uaer(n,k2:k1-1,itr)))
                aertrlw(n,k1,k2) = dexp(-fiveothree * qabslw * uaerdust)
              end if
            end if
          end do
        end do
      end do
    end do
!     
  end subroutine aeroppt
!
end module mod_rad_aerosol
