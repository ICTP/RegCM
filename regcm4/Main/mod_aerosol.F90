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

      use mod_dynparam

      implicit none

!     MODIF 16/09/2005 IBRAH Internal mixing

!
! PARAMETER definitions
!
!      Num of spctrl intervals across solar spectrum
!
      integer , parameter :: nspi = 19   ! Spectral index
!
      integer , parameter :: ncoefs = 5  ! Number of coefficients
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
      real(8) , dimension(nspi) :: gsbase , gsbc_hb , gsbc_hl ,         &
                                 & gsoc_hb , gsoc_hl , ksbase ,         &
                                 & ksbc_hb , ksbc_hl , ksoc_hb ,        &
                                 & ksoc_hl , wsbase , wsbc_hb ,         &
                                 & wsbc_hl , wsoc_hb , wsoc_hl
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
!     Aerosol extinction optical depth
!
      real(8) , allocatable , dimension(:,:,:) :: ftota_mix ,           &
                 & gtota_mix , tauasc_mix , tauxar_mix
!
      real(8) , allocatable , dimension(:,:) :: ftota_mix_cs ,          &
                 & gtota_mix_cs , tauasc_mix_cs , tauxar_mix_cs
!
!     Work arrays for aeroppt
!
      real(8) , allocatable , dimension(:,:) :: aermtot , aervtot
      real(8) , allocatable , dimension(:,:,:) :: fa , ga , tauxar ,    &
                              &                   uaer , wa
      real(8) , allocatable , dimension(:,:) :: faer , gaer , tauaer ,  &
                              &   utaer , waer
      real(8) , dimension(4) :: frac , prop


!------------------------------------------------------------------------------
!                  DATA SECTION
!------------------------------------------------------------------------------

      data ksbase/5.206D0 , 5.206D0 , 5.206D0 , 5.206D0 , 5.206D0 ,     &
         & 5.206D0 , 5.206D0 , 3.203D0 , 3.203D0 , 1.302D0 , 5.992D-01 ,&
         & 2.948D-01 , 1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 ,  &
         & 5.770D-02 , 2.290D-01 , 2.270D-01/
!
      data ((kscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/1.126D+01 ,          &
          & -2.502D-01 , -1.087D+00 , -1.794D+02 , 1.556D+01 ,          &
          & 1.126D+01 , -2.502D-01 , -1.087D+00 , -1.794D+02 ,          &
          & 1.556D+01 , 1.126D+01 , -2.502D-01 , -1.087D+00 ,           &
          & -1.794D+02 , 1.556D+01 , 1.126D+01 , -2.502D-01 ,           &
          & -1.087D+00 , -1.794D+02 , 1.556D+01 , 1.126D+01 ,           &
          & -2.502D-01 , -1.087D+00 , -1.794D+02 , 1.556D+01 ,          &
          & 1.126D+01 , -2.502D-01 , -1.087D+00 , -1.794D+02 ,          &
          & 1.556D+01 , 1.126D+01 , -2.502D-01 , -1.087D+00 ,           &
          & -1.794D+02 , 1.556D+01 , 1.124D+01 , -3.040D-01 ,           &
          & -1.088D+00 , -1.776D+02 , 1.537D+01 , 1.124D+01 ,           &
          & -3.040D-01 , -1.088D+00 , -1.776D+02 , 1.537D+01 ,          &
          & 1.222D+01 , -3.770D-01 , -1.089D+00 , -1.898D+02 ,          &
          & 1.504D+01 , 1.357D+01 , -4.190D-01 , -1.087D+00 ,           &
          & -2.070D+02 , 1.478D+01 , 1.557D+01 , -4.353D-01 ,           &
          & -1.083D+00 , -2.382D+02 , 1.486D+01 , 1.758D+01 ,           &
          & -4.389D-01 , -1.078D+00 , -2.716D+02 , 1.505D+01 ,          &
          & 1.597D+01 , -4.337D-01 , -1.073D+00 , -2.510D+02 ,          &
          & 1.527D+01 , 2.107D+01 , -3.041D-01 , -1.067D+00 ,           &
          & -2.494D+02 , 1.166D+01 , -2.424D-01 , -1.770D-01 ,          &
          & -1.032D+00 , 1.469D-01 , 1.947D+00 , 2.535D+01 ,            &
          & -2.270D-01 , -1.052D+00 , -2.528D+02 , 9.888D+00 ,          &
          & -1.545D-01 , -1.661D-01 , -1.030D+00 , -4.698D-04 ,         &
          & 7.275D-02 , 8.835D-01 , -1.590D-01 , -1.029D+00 ,           &
          & -2.838D+01 , 2.734D+01/
!
      data wsbase/7.371D-08 , 7.371D-08 , 7.371D-08 , 7.371D-08 ,       &
         & 7.371D-08 , 7.371D-08 , 7.371D-08 , 6.583D-08 , 6.583D-08 ,  &
         & 3.656D-06 , 4.919D-05 , 3.539D-03 , 2.855D-02 , 2.126D-01 ,  &
         & 8.433D-01 , 9.653D-01 , 6.198D-01 , 9.642D-01 , 9.699D-01/
!
      data ((wscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/2.492D+00 ,          &
          & -5.210D-02 , -1.036D+00 , -4.398D+01 , 1.724D+01 ,          &
          & 2.492D+00 , -5.210D-02 , -1.036D+00 , -4.398D+01 ,          &
          & 1.724D+01 , 2.492D+00 , -5.210D-02 , -1.036D+00 ,           &
          & -4.398D+01 , 1.724D+01 , 2.492D+00 , -5.210D-02 ,           &
          & -1.036D+00 , -4.398D+01 , 1.724D+01 , 2.492D+00 ,           &
          & -5.210D-02 , -1.036D+00 , -4.398D+01 , 1.724D+01 ,          &
          & 2.492D+00 , -5.210D-02 , -1.036D+00 , -4.398D+01 ,          &
          & 1.724D+01 , 2.492D+00 , -5.210D-02 , -1.036D+00 ,           &
          & -4.398D+01 , 1.724D+01 , 1.139D+00 , -1.110D-02 ,           &
          & -1.011D+00 , -7.754D+00 , 6.737D+00 , 1.139D+00 ,           &
          & -1.110D-02 , -1.011D+00 , -7.754D+00 , 6.737D+00 ,          &
          & 1.848D+00 , -3.920D-04 , -9.924D-01 , -1.607D+00 ,          &
          & 8.587D-01 , 5.459D+00 , 9.357D-01 , -1.626D+00 ,            &
          & -5.282D+00 , 1.066D+00 , 1.187D+00 , 2.241D-01 ,            &
          & -1.226D+00 , 1.442D+01 , -1.402D+01 , -3.640D+00 ,          &
          & 2.552D-01 , -1.168D+00 , 4.458D+01 , 1.152D+01 ,            &
          & -5.634D+00 , 2.068D-01 , -1.122D+00 , 7.528D+01 ,           &
          & 1.290D+01 , 1.826D-01 , 6.588D-02 , -1.098D+00 ,            &
          & -1.996D-02 , 1.618D-01 , 2.164D+00 , 1.194D-01 ,            &
          & -1.044D+00 , -3.221D+01 , 1.564D+01 , 2.268D-01 ,           &
          & 3.266D-02 , -1.064D+00 , -2.677D-02 , 1.309D-01 ,           &
          & 2.178D+00 , 1.151D-01 , -1.042D+00 , -3.325D+01 ,           &
          & 1.600D+01 , 1.713D+00 , 9.166D-02 , -1.039D+00 ,            &
          & -2.660D+01 , 1.629D+01/
!
      data gsbase/6.899D-01 , 6.899D-01 , 6.899D-01 , 6.899D-01 ,       &
         & 6.899D-01 , 6.899D-01 , 6.899D-01 , 6.632D-01 , 6.632D-01 ,  &
         & 5.912D-01 , 5.111D-01 , 4.269D-01 , 3.321D-01 , 2.197D-01 ,  &
         & 1.305D-01 , 7.356D-02 , 1.602D-01 , 6.883D-02 , 6.304D-02/
!
      data ((gscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/ - 9.874D-01 ,       &
          & -3.033D+01 , -2.138D+01 , -2.265D+00 , 5.238D+00 ,          &
          & -9.874D-01 , -3.033D+01 , -2.138D+01 , -2.265D+00 ,         &
          & 5.238D+00 , -9.874D-01 , -3.033D+01 , -2.138D+01 ,          &
          & -2.265D+00 , 5.238D+00 , -9.874D-01 , -3.033D+01 ,          &
          & -2.138D+01 , -2.265D+00 , 5.238D+00 , -9.874D-01 ,          &
          & -3.033D+01 , -2.138D+01 , -2.265D+00 , 5.238D+00 ,          &
          & -9.874D-01 , -3.033D+01 , -2.138D+01 , -2.265D+00 ,         &
          & 5.238D+00 , -9.874D-01 , -3.033D+01 , -2.138D+01 ,          &
          & -2.265D+00 , 5.238D+00 , -3.666D-01 , -1.319D+00 ,          &
          & -3.311D+00 , -2.821D-02 , 8.844D-01 , -3.666D-01 ,          &
          & -1.319D+00 , -3.311D+00 , -2.821D-02 , 8.844D-01 ,          &
          & 5.824D-01 , -1.875D-01 , -1.567D+00 , -4.402D+00 ,          &
          & 6.268D+00 , 1.238D+00 , -1.550D-01 , -1.368D+00 ,           &
          & -1.127D+01 , 8.334D+00 , 2.299D+00 , -1.686D-01 ,           &
          & -1.304D+00 , -2.677D+01 , 1.101D+01 , 3.037D+00 ,           &
          & -1.447D-01 , -1.223D+00 , -2.609D+01 , 8.267D+00 ,          &
          & 4.683D+00 , -2.307D-01 , -1.241D+00 , -4.312D+01 ,          &
          & 8.838D+00 , 3.842D+00 , -6.301D-01 , -1.367D+00 ,           &
          & -4.144D+01 , 9.620D+00 , 3.237D+00 , -4.530D-01 ,           &
          & -1.204D+00 , -3.234D+01 , 8.946D+00 , 4.181D+00 ,           &
          & -4.140D-01 , -1.284D+00 , -4.489D+01 , 9.950D+00 ,          &
          & 3.378D+00 , -4.334D-01 , -1.188D+00 , -3.664D+01 ,          &
          & 9.786D+00 , 3.943D+00 , -3.952D-01 , -1.170D+00 ,           &
          & -4.415D+01 , 1.031D+01/
!
!-----------------------------------------------------------------------
!
      data ksbc_hb/20.783 , 17.212 , 15.864 , 15.053 , 14.304 , 13.613 ,&
         & 11.966 , 6.5782 , 4.3961 , 4.38 , 2.11 , 2.11 , 2.11 , 2.11 ,&
         & 1.40 , 1.40 , 1.40 , 1.40 , 1.40/
 
      data wsbc_hb/0.24524 , 0.20962 , 0.195 , 0.18586 , 0.17719 ,      &
         & 0.16897 , 0.14849 , 0.071748 , 0.037536 , .089 , .025 ,      &
         & .025 , .025 , .025 , .009 , .009 , .009 , .009 , .009/
 
      data gsbc_hb/0.21387 , 0.17154 , 0.15564 , 0.1461 , 0.13732 ,     &
         & 0.12923 , 0.11006 , 0.049175 , 0.026638 , .220 , .123 ,      &
         & .123 , .123 , .123 , .073 , .073 , .073 , .073 , .073/
 
!----------------------------------------------------------------------

      data ksbc_hl/14.851 , 14.258 , 13.943 , 13.724 , 13.507 , 13.295 ,&
         & 12.722 , 9.4434 , 6.9653 , 4.38 , 2.11 , 2.11 , 2.11 , 2.11 ,&
         & 1.40 , 1.40 , 1.40 , 1.40 , 1.40/
 
      data wsbc_hl/0.46081 , 0.44933 , 0.44397 , 0.44065 , 0.43737 ,    &
         & 0.43394 , 0.42346 , 0.35913 , 0.29579 , .089 , .025 , .025 , &
         & .025 , .025 , .009 , .009 , .009 , .009 , .009/
 
      data gsbc_hl/0.69038 , 0.65449 , 0.63711 , 0.62542 , 0.61407 ,    &
         & 0.60319 , 0.57467 , 0.4205 , 0.3066 , .220 , .123 , .123 ,   &
         & .123 , .123 , .073 , .073 , .073 , .073 , .073/
 
!---------------------------------------------------------------------

      data ksoc_hb/6.0584 , 6.0654 , 6.1179 , 6.0102 , 5.8 , 5.6957 ,   &
         & 5.6494 , 4.3283 , 3.1485 , 1.302 , 5.992D-01 , 2.948D-01 ,   &
         & 1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 , 5.770D-02 ,  &
         & 2.290D-01 , 2.270D-01/
 
      data wsoc_hb/0.91735 , 0.92365 , 0.92941 , 0.93067 , 0.93311 ,    &
         & 0.93766 , 0.94042 , 0.95343 , 0.9548 , 0.70566 , 0.70566 ,   &
         & 0.70566 , 0.70566 , 0.70566 , 0.70566 , 0.70566 , 0.70566 ,  &
         & 0.70566 , 0.70566/
 
      data gsoc_hb/0.67489 , 0.67003 , 0.67725 , 0.65487 , 0.65117 ,    &
         & 0.66116 , 0.64547 , 0.60033 , 0.55389 , 5.912D-01 ,          &
         & 5.111D-01 , 4.269D-01 , 3.321D-01 , 2.197D-01 , 1.305D-01 ,  &
         & 7.356D-02 , 1.602D-01 , 6.883D-02 , 6.304D-02/
 
!---------------------------------------------------------------------

      data ksoc_hl/3.543 , 3.623 , 3.7155 , 3.712 , 3.6451 , 3.6376 ,   &
         & 3.6844 , 3.4588 , 2.9846 , 1.302 , 5.992D-01 , 2.948D-01 ,   &
         & 1.475D-01 , 7.387D-02 , 1.683D-01 , 2.655D-01 , 5.770D-02 ,  &
         & 2.290D-01 , 2.270D-01/
 
      data wsoc_hl/0.87931 , 0.88292 , 0.89214 , 0.89631 , 0.89996 ,    &
         & 0.9054 , 0.90805 , 0.93423 , 0.95012 , 0.85546 , 0.85546 ,   &
         & 0.85546 , 0.85546 , 0.85546 , 0.85546 , 0.85546 , 0.85546 ,  &
         & 0.85546 , 0.85546/
 
      data gsoc_hl/0.73126 , 0.71089 , 0.72042 , 0.69924 , 0.69908 ,    &
         & 0.70696 , 0.68479 , 0.64879 , 0.63433 , 5.912D-01 ,          &
         & 5.111D-01 , 4.269D-01 , 3.321D-01 , 2.197D-01 , 1.305D-01 ,  &
         & 7.356D-02 , 1.602D-01 , 6.883D-02 , 6.304D-02/
 
!     DUST OP data base for external mixing : maximum of 4 bin for the
!     momeent , determined from Zender et al.
 
      data ((ksdust(ii,jj),jj=1,4),ii=1,nspi)/1.8801 , 0.76017 ,        &
          & 0.36681 , 0.16933 , 2.0254 , 0.78378 , 0.36845 , 0.17002 ,  &
          & 1.9547 , 0.76389 , 0.37119 , 0.17032 , 1.8996 , 0.74916 ,   &
          & 0.3722 , 0.17052 , 1.7946 , 0.74044 , 0.36984 , 0.17072 ,   &
          & 1.7149 , 0.75401 , 0.36892 , 0.1709 , 1.5431 , 0.82322 ,    &
          & 0.37505 , 0.17143 , 2.4482 , 0.8568 , 0.38078 , 0.17396 ,   &
          & 3.1067 , 0.74488 , 0.43688 , 0.18104 , 0.50391 , 0.67245 ,  &
          & 0.53605 , 0.20599 , 0.50391 , 0.67245 , 0.53605 , 0.20599 , &
          & 0.50391 , 0.67245 , 0.53605 , 0.20599 , 0.50391 , 0.67245 , &
          & 0.53605 , 0.20599 , 0.50391 , 0.67245 , 0.53605 , 0.20599 , &
          & 0.50391 , 0.67245 , 0.53605 , 0.20599 , 0.50391 , 0.67245 , &
          & 0.53605 , 0.20599 , 0.50391 , 0.67245 , 0.53605 , 0.20599 , &
          & 0.50391 , 0.67245 , 0.53605 , 0.20599 , 0.50391 , 0.67245 , &
          & 0.53605 , 0.20599/
 
      data ((wsdust(ii,jj),jj=1,4),ii=1,nspi)/0.64328 , 0.55196 ,       &
          & 0.53748 , 0.54342 , 0.67757 , 0.56909 , 0.53639 , 0.54232 , &
          & 0.67316 , 0.56027 , 0.53875 , 0.54181 , 0.66245 , 0.55338 , &
          & 0.53947 , 0.54149 , 0.68132 , 0.5744 , 0.54328 , 0.54143 ,  &
          & 0.6796 , 0.58467 , 0.54242 , 0.54113 , 0.72679 , 0.68744 ,  &
          & 0.58564 , 0.54576 , 0.9473 , 0.88181 , 0.80761 , 0.70455 ,  &
          & 0.97536 , 0.89161 , 0.86378 , 0.758 , 0.89568 , 0.96322 ,   &
          & 0.95008 , 0.89293 , 0.89568 , 0.96322 , 0.95008 , 0.89293 , &
          & 0.89568 , 0.96322 , 0.95008 , 0.89293 , 0.89568 , 0.96322 , &
          & 0.95008 , 0.89293 , 0.89568 , 0.96322 , 0.95008 , 0.89293 , &
          & 0.89568 , 0.96322 , 0.95008 , 0.89293 , 0.89568 , 0.96322 , &
          & 0.95008 , 0.89293 , 0.89568 , 0.96322 , 0.95008 , 0.89293 , &
          & 0.89568 , 0.96322 , 0.95008 , 0.89293 , 0.89568 , 0.96322 , &
          & 0.95008 , 0.89293/
 
      data ((gsdust(ii,jj),jj=1,4),ii=1,nspi)/0.87114 , 0.92556 ,       &
          & 0.94542 , 0.94831 , 0.86127 , 0.921 , 0.94355 , 0.94813 ,   &
          & 0.838 , 0.91194 , 0.94304 , 0.94803 , 0.8176 , 0.90442 ,    &
          & 0.94239 , 0.94796 , 0.77088 , 0.88517 , 0.9371 , 0.94775 ,  &
          & 0.73925 , 0.88364 , 0.93548 , 0.94763 , 0.60695 , 0.86086 , &
          & 0.91824 , 0.94473 , 0.64393 , 0.76457 , 0.8133 , 0.87784 ,  &
          & 0.7476 , 0.62041 , 0.80665 , 0.85974 , 0.26761 , 0.56045 ,  &
          & 0.68897 , 0.68174 , 0.26761 , 0.56045 , 0.68897 , 0.68174 , &
          & 0.26761 , 0.56045 , 0.68897 , 0.68174 , 0.26761 , 0.56045 , &
          & 0.68897 , 0.68174 , 0.26761 , 0.56045 , 0.68897 , 0.68174 , &
          & 0.26761 , 0.56045 , 0.68897 , 0.68174 , 0.26761 , 0.56045 , &
          & 0.68897 , 0.68174 , 0.26761 , 0.56045 , 0.68897 , 0.68174 , &
          & 0.26761 , 0.56045 , 0.68897 , 0.68174 , 0.26761 , 0.56045 , &
          & 0.68897 , 0.68174/
!
      contains

! 
      subroutine allocate_mod_aerosol
      implicit none   
#ifdef MPP1
      allocate(aermm(iym1,kz,jxp))
#else
      allocate(aermm(iym1,kz,jxm1))
#endif 
      allocate(aermmb(iym1,kz))
      allocate(aermmr(iym1,kz,ntr))
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
      end subroutine allocate_mod_aerosol

!
! SUBROUTINE AERMIX
!
      subroutine aermix(pint , rh , j , istart , iend , nx , nk , ntrac)
 
!-----------------------------------------------------------------------
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
 
      use mod_param2 , only : ichem
      use mod_slice , only : rhb3d
      use mod_main , only : psa
      use mod_mainchem , only : chia
      use mod_constants , only : gtigts
      implicit none
!
! Dummy arguments
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
!---------------------------Local variables-----------------------------
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
 
      real(8) :: gvis , kaervs , omgvis , rhfac , tauvis
      integer :: i , itr , k , mxaerl
!
      data kaervs/5.3012D0/         ! multiplication factor for kaer
      data omgvis/0.999999D0/
      data gvis/0.694889D0/
      data rhfac/1.6718D0/          ! EES added for efficiency
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
          rh(i,k) = dmin1(rhb3d(i,k,j),0.99D0)
!EES:     do not change to 1.00:  wscoef(3,10) in radcsw = .9924 and is
!         divided by RH.  rh is limited to .99 to avoid dividing by zero
!added
 
!
!         Define background aerosol
!         Find constant aerosol mass mixing ratio for specified levels
!         in the column, converting units where appropriate
!         for the moment no more used
!
          if ( k.ge.nk + 1 - mxaerl ) then
            aermmb(i,k) = gtigts*tauvis/(1.0D4*kaervs*rhfac*(1.-omgvis* &
                        & gvis*gvis)                                    &
                        & *(pint(i,kzp1)-pint(i,kz + 1 - mxaerl)))
          else
            aermmb(i,k) = 0.0D0
          end if
 
!         if(ichem .eq. 1 .and. idirect .eq. 1) then
          if ( ichem.eq.1 ) then
            do itr = 1 , ntrac
!             aermmr(i,k,itr)= dmax1( chia(i,k,j,itr)/psa(i,j)
!             $                               ,aermmb(i,k) )
              aermmr(i,k,itr) = chia(i,k,j,itr)/psa(i,j)
            end do
          else if ( ehso4 ) then
            do itr = 1 , ntrac
!             aermmr(i,k,itr) = aermmb(i,k) + aermm(i,k,j)
!Dec.11       aermmr(i,k,itr) = aermm(i,k,j)
              aermmr(i,k,itr) = 0.0D0
!Dec.11_
            end do
          else
            do itr = 1 , ntrac
!             aermmr(i,k,itr)= 0.0D0
              aermmr(i,k,itr) = aermmb(i,k)
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
      use mod_param2 , only : ichem , idirect
      use mod_trachem , only : mixtype , chtrname , iso4 , ibchl ,      &
                    &          iochl , idust
      use mod_message , only : fatal
      use mod_constants , only : rhoso4 , rhobc, rhooc, rhodust , gtigts
      implicit none
!
! Dummy arguments
!
!     Interface pressure, relative humidity
!
      real(8) , dimension(iym1,kzp1) :: pint
      real(8) , dimension(iym1,kz) :: rh
      intent (in) pint , rh
!
! Local variables
!
      integer :: i , i1 , i2 , i3 , i4 , ibin , itr , k , ns
      real(8) :: path
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
      if ( ichem.ne.1 ) then
        tauxar_mix_cs(:,:) = 0.0
        tauasc_mix_cs(:,:) = 0.0
        gtota_mix_cs(:,:) = 0.0
        ftota_mix_cs(:,:) = 0.0
        tauxar_mix(:,:,:) = 0.0
        tauasc_mix(:,:,:) = 0.0
        gtota_mix(:,:,:) = 0.0
        ftota_mix(:,:,:) = 0.0
        return
      end if

      if ( idirect.ge.1 ) then
        tauxar = 0.0
        wa = 0.0
        ga = 0.0
        fa = 0.0
      else
!       Nothing to do for this
        tauxar_mix_cs(:,:) = 0.0
        tauasc_mix_cs(:,:) = 0.0
        gtota_mix_cs(:,:) = 0.0
        ftota_mix_cs(:,:) = 0.0
        tauxar_mix(:,:,:) = 0.0
        tauasc_mix(:,:,:) = 0.0
        gtota_mix(:,:,:) = 0.0
        ftota_mix(:,:,:) = 0.0
        return
      end if
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!       Option 1  melange externe
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if ( mixtype.eq.1 ) then
!
!       Spectral loop
!
        do ns = 1 , nspi

          tauxar_mix_cs(:,ns) = 0.0
          tauasc_mix_cs(:,ns) = 0.0
          gtota_mix_cs(:,ns) = 0.0
          ftota_mix_cs(:,ns) = 0.0

          tauxar_mix(:,:,ns) = 0.0
          tauasc_mix(:,:,ns) = 0.0
          gtota_mix(:,:,ns) = 0.0
          ftota_mix(:,:,ns) = 0.0

          uaer(:,0,:) = 0.0
          tauxar(:,0,:) = 0.0
          wa(:,0,:) = 0.0
          ga(:,0,:) = 0.0
          fa(:,0,:) = 0.0
          utaer(:,:) = 0.0
          tauaer(:,:) = 0.0
          waer(:,:) = 0.0
          gaer(:,:) = 0.0
          faer(:,:) = 0.0
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
!
!         calculate optical properties of each aerosol component
!

          do k = 1 , kz
            do i = 1 , iym1
              path = (pint(i,k+1)-pint(i,k))/gtigts
              ibin = 0
              do itr = 1 , ntr
                uaer(i,k,itr) = 0.
                if ( rh(i,k).lt.0.0 .or. rh(i,k).gt.1.0 ) print * ,     &
                   & i , k , rh(i,k) , '  RH WARNING !!!!!'

                if ( chtrname(itr).eq.'XXXXX') then
                  continue
                end if

                if ( chtrname(itr).eq.'DUST' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 gaffe au facteur !!
!                 gaffe au ntr/bins
                  ibin = ibin + 1
                  if ( ibin.gt.4 ) print * , 'DUST OP PBLEME !!!!'

                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksdust(ns,ibin)
                  wa(i,k,itr) = wsdust(ns,ibin)
                  ga(i,k,itr) = gsdust(ns,ibin)
                  fa(i,k,itr) = gsdust(ns,ibin)*gsdust(ns,ibin)
 
                else if ( chtrname(itr).eq.'SO4' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbase(ns)       &
                                  & *exp(kscoef(ns,1)+kscoef(ns,2)      &
                                  & /(rh(i,k)+kscoef(ns,3))             &
                                  & +kscoef(ns,4)                       &
                                  & /(rh(i,k)+kscoef(ns,5)))
!
                  wa(i,k,itr) = 1.0 - wsbase(ns)                        &
                              & *exp(wscoef(ns,1)+wscoef(ns,2)          &
                              & /(rh(i,k)+wscoef(ns,3))+wscoef(ns,4)    &
                              & /(rh(i,k)+wscoef(ns,5)))
!
                  ga(i,k,itr) = gsbase(ns)                              &
                              & *exp(gscoef(ns,1)+gscoef(ns,2)          &
                              & /(rh(i,k)+gscoef(ns,3))+gscoef(ns,4)    &
                              & /(rh(i,k)+gscoef(ns,5)))
!
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
 
                else if ( chtrname(itr).eq.'OC_HL' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Humidity effect !
                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksoc_hl(ns)      &
                                  & *(1-rh(i,k))**(-0.2)
                  wa(i,k,itr) = wsoc_hl(ns)
                  ga(i,k,itr) = gsoc_hl(ns)
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                else if ( chtrname(itr).eq.'BC_HL' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Humidity effect !
                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbc_hl(ns)      &
                                  & *(1-rh(i,k))**(-0.25)
                  wa(i,k,itr) = wsbc_hl(ns)
                  ga(i,k,itr) = gsbc_hl(ns)
                  fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                else if ( chtrname(itr).eq.'OC_HB' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksoc_hb(ns)
                  wa(i,k,itr) = wsoc_hb(ns)
                  ga(i,k,itr) = gsoc_hb(ns)
                  fa(i,k,itr) = gsoc_hb(ns)*gsoc_hb(ns)
                else if ( chtrname(itr).eq.'BC_HB' ) then
                  uaer(i,k,itr) = aermmr(i,k,itr)*path
!                 Absorbing aerosols (soot type)
                  tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbc_hb(ns)
                  wa(i,k,itr) = wsbc_hb(ns)
                  ga(i,k,itr) = gsbc_hb(ns)
                  fa(i,k,itr) = gsbc_hb(ns)*gsbc_hb(ns)
                else
                end if
              end do  ! end tracer loop
            end do
          end do

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
              if ( utaer(i,itr).le.1.D-10 ) utaer(i,itr) = 1.D-10
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
                                   & + tauxar(i,k,itr)
                tauasc_mix(i,k,ns) = tauasc_mix(i,k,ns)               &
                                   & + tauxar(i,k,itr)*wa(i,k,itr)
                gtota_mix(i,k,ns) = gtota_mix(i,k,ns) + ga(i,k,itr)   &
                                  & *tauxar(i,k,itr)*wa(i,k,itr)
                ftota_mix(i,k,ns) = ftota_mix(i,k,ns) + fa(i,k,itr)   &
                                  & *tauxar(i,k,itr)*wa(i,k,itr)
              end do
!
!             Clear sky (always calcuated if idirect >=1 for
!             diagnostic radiative forcing)
!
              tauxar_mix_cs(i,ns) = tauxar_mix_cs(i,ns)               &
                                  & + tauaer(i,itr)
              tauasc_mix_cs(i,ns) = tauasc_mix_cs(i,ns)               &
                                  & + tauaer(i,itr)*waer(i,itr)
              gtota_mix_cs(i,ns) = gtota_mix_cs(i,ns) + gaer(i,itr)   &
                                 & *tauaer(i,itr)*waer(i,itr)
              ftota_mix_cs(i,ns) = ftota_mix_cs(i,ns) + faer(i,itr)   &
                                 & *tauaer(i,itr)*waer(i,itr)
            end do
          end do
!
        end do ! end spectral loop

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!         Option 2  melange interne
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      else if ( mixtype.eq.2 ) then

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Spectral loop
!
        do ns = 1 , nspi

          tauxar_mix_cs(:,ns) = 0.0
          tauasc_mix_cs(:,ns) = 0.0
          gtota_mix_cs(:,ns) = 0.0
          ftota_mix_cs(:,ns) = 0.0
          tauxar_mix(:,:,ns) = 0.0
          tauasc_mix(:,:,ns) = 0.0
          gtota_mix(:,:,ns) = 0.0
          ftota_mix(:,:,ns) = 0.0
 
          uaer(:,0,:) = 0.0
          tauxar(:,0,:) = 0.0
          wa(:,0,:) = 0.0
          ga(:,0,:) = 0.0
          fa(:,0,:) = 0.0
          utaer(:,:) = 0.0
          tauaer(:,:) = 0.0
          waer(:,:) = 0.0
          gaer(:,:) = 0.0
          faer(:,:) = 0.0
 
!
!         calculate optical properties of each aerosol component
!
          do i = 1 , iym1
            do k = 1 , kz
              path = (pint(i,k+1)-pint(i,k))/gtigts
              if ( rh(i,k).lt.0.0 .or. rh(i,k).gt.1.0 )                 &
                write ( 6,* ) 'WARNING RH : ' , i , k , rh(i,k)
 
!             sum of hydrophilic aerosols
              aervtot(i,k) = 0.
              aermtot(i,k) = 0.
 
              if ( iso4.ne.0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iso4)/rhoso4
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iso4)
              end if
 
              if ( ibchl.ne.0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,ibchl)/rhobc
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,ibchl)
              end if
 
              if ( iochl.ne.0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iochl)/rhooc
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iochl)
              end if
 
              if ( idust(1).ne.0 ) then
                aervtot(i,k) = aervtot(i,k) + aermmr(i,k,idust(1))    &
                             & /rhodust
                aermtot(i,k) = aermtot(i,k) + aermmr(i,k,idust(1))
              end if
 
!             minimum quantity of total aerosol
 
              if ( aermtot(i,k).gt.1.D-14 ) then
!               indexes in the internal mixing table
                prop(1) = (aermmr(i,k,iso4)/rhoso4)/aervtot(i,k)
                prop(2) = (aermmr(i,k,ibchl)/rhobc)/aervtot(i,k)
                prop(3) = (aermmr(i,k,iochl)/rhooc)/aervtot(i,k)
                prop(4) = (aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
                frac(1) = fraction(prop(1))
                frac(2) = fraction(prop(2))
                frac(3) = fraction(prop(3))
                frac(4) = fraction(prop(4))
!               FIND THE GREATEST FRACTIONAL PART
 
                if ( iso4.ne.0 ) then
                  i1 = nint(10*prop(1)) + 1
                else
                  i1 = 0 + 1
                end if
 
                if ( ibchl.ne.0 ) then
                  i2 = nint(10*prop(2)) + 1
                else
                  i2 = 0 + 1
                end if
 
                if ( iochl.ne.0 ) then
                  i3 = nint(10*prop(3)) + 1
                else
                  i3 = 0 + 1
                end if
 
                if ( idust(1).ne.0 ) then
                  i4 = nint(10*prop(4)) + 1
                else
                  i4 = 0 + 1
                end if
!
!               final optical parameters
!
                if ( i1+i2+i3+i4.eq.13 ) i4 = i4 + 1
                if ( i1+i2+i3+i4.eq.15 ) then
                  if ( i4.ne.1 ) i4 = i4 - 1
                end if
 
                if ( i1+i2+i3+i4.eq.15 ) then
                  if ( i1.ne.1 ) i1 = i1 - 1
                end if
 
                if ( i1+i2+i3+i4.eq.15 ) then
                  if ( i3.ne.1 ) i3 = i3 - 1
                end if
 
                if ( i1+i2+i3+i4.eq.15 ) call fatal(__FILE__,__LINE__,&
                    &'WRONG COMBINATION. SHOULD NEVER HAPPEN')
 
                if ( i1+i2+i3+i4.ne.14 ) then
                  print * , i1 , i2 , i3 , i4 , i1 + i2 + i3 + i4
                  print * , idust(1) , iochl , ibchl , iso4
                  print * , 'OC HL' , aermmr(i,k,iochl)/rhooc
                  print * , 'BC HL' , aermmr(i,k,ibchl)/rhobc
                  print * , 'SO4' , aermmr(i,k,iso4)/rhoso4
                  print * , 'DUST' , aermmr(i,k,idust(1))/rhodust
                  print * , 'VOL TOT' , aervtot(i,k)
                  print * , 'OC HL%' , 10*(aermmr(i,k,iochl)/rhooc)   &
                      & /aervtot(i,k)
                  print * , 'BC HL%' , 10*(aermmr(i,k,ibchl)/rhobc)   &
                      & /aervtot(i,k)
                  print * , 'SO4 %' , 10*(aermmr(i,k,iso4)/rhoso4)    &
                      & /aervtot(i,k)
                  print * , 'SO4 %' ,                                 &
                      & nint(10*(aermmr(i,k,iso4)/rhoso4)/aervtot(i,k)&
                      & )
                  print * , 'DUST %' ,                                &
                      & 10*(aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
                  print * , 'DUST %' ,                                &
                      & nint(10*(aermmr(i,k,idust(1))/rhodust)        &
                      & /aervtot(i,k))
                  call fatal(__FILE__,__LINE__,                       &
                            &'SOMETHING WRONG ON SPECIES ABUNDANCE')
                end if

                tauxar_mix(i,k,ns) = dextmix(1,ns,i4,i2,i3,i1)        &
                                   & *aermtot(i,k)*path*1D5
                tauasc_mix(i,k,ns) = dssamix(1,ns,i4,i2,i3,i1)        &
                                   & *tauxar_mix(i,k,ns)
                gtota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)           &
                                  & *tauasc_mix(i,k,ns)               &
                                  & *tauxar_mix(i,k,ns)
                ftota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)           &
                                  & *dgmix(1,ns,i4,i2,i3,i1)          &
                                  & *tauasc_mix(i,k,ns)               &
                                  & *tauxar_mix(i,k,ns)
 
!               clear sky dignostic
                utaer(i,1) = utaer(i,1) + aermtot(i,k)*path
 
                tauaer(i,1) = tauaer(i,1) + dextmix(1,ns,i4,i2,i3,i1) &
                            & *aermtot(i,k)*path*1D5
                waer(i,1) = waer(i,1) + dssamix(1,ns,i4,i2,i3,i1)     &
                          & *aermtot(i,k)*path
                gaer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                          & *aermtot(i,k)*path
                faer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                          & *dgmix(1,ns,i4,i2,i3,i1)*aermtot(i,k)*path
 
              end if ! end minimum concentration conditions
            end do ! end k loop
 
            if ( utaer(i,1).gt.1.D-12 ) then
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      else
        write ( 6,* ) 'MIXTYPE = ', mixtype
        call fatal(__FILE__,__LINE__,'UNSUPPORTED MIXTYPE IN AEROPPT')
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      end subroutine aeroppt
!
! SUBROUTINE AEROUT
!
      subroutine aerout(jslc,aeradfo,aeradfos)
!
      use mod_dynparam , only : chemfrq
      use mod_param2 , only : radfrq
      use mod_trachem , only : aerext , aerssa , aerasp , aertarf ,     &
                      &        aersrrf
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(iym1) :: aeradfo , aeradfos
      intent (in) aeradfo , aeradfos
!
! Local variables
!
      integer :: i , k , ntim
! 
      do k = 1 , kz
        do i = 2 , iym1
#ifdef MPP1
          aerext(i-1,k,jslc) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc) = gtota_mix(i,k,8)
#else
          aerext(i-1,k,jslc-1) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc-1) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc-1) = gtota_mix(i,k,8)
#endif
        end do
      end do
!
!     CARE :Average the radiative forcing between chem output time
!     steps (in hour) according to radfrq (in min), aertarf is reset to
!     0 at each chem output (cf output.f)
!
      ntim = 60*chemfrq/radfrq
!
!     aersol radative forcing (care cgs to mks after radiation scheme !)
!
      do i = 2 , iym1
#ifdef MPP1
        aertarf(i-1,jslc) = aertarf(i-1,jslc) + aeradfo(i)*1.E-3/ntim
        aersrrf(i-1,jslc) = aersrrf(i-1,jslc) + aeradfos(i)*1.E-3/ntim
#else
        aertarf(i-1,jslc-1) = aertarf(i-1,jslc-1) + aeradfo(i)          &
                            & *1.E-3/ntim
        aersrrf(i-1,jslc-1) = aersrrf(i-1,jslc-1) + aeradfos(i)         &
                            & *1.E-3/ntim
#endif
      end do
 
      end subroutine aerout

      end module mod_aerosol
