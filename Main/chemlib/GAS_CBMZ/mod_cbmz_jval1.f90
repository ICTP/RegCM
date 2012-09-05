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

module mod_cbmz_jval1

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam

  private

  public :: jvalpro , readhv

  contains

    subroutine jvalpro(nhv,hvmat,jarray,jparam,jval)
      implicit none
!
      real(rk8) , dimension(22,40) :: hvmat
      real(rk8) , dimension(22) :: jparam
      real(rk8) , dimension(80,510,56) :: jarray
      real(rk8) , dimension(56) :: jval
      integer(ik4) , dimension(22) :: nhv
      intent (in) hvmat , jarray , jparam , nhv
      intent (inout) jval
!
      real(rk8) , dimension(56) :: cfac , jfaerz , jfsur
      real(rk8) :: falt , fzen , x
      real(rk8) :: fkn
      integer(ik4) :: i , ialt , id , ig01 , ig02 , ig11 , ig12 , ig21 ,     &
                 ig22 , ij , im , iwri , iy , izen , j , jaer , jalb ,  &
                 jc , jcld , jct , jtem , k , kn
      real(rk8) , dimension(20,56) :: jfrac
      real(rk8) , dimension(20) :: jfx
!
!     This Subroutine takes jparams (zenith, altitude, etc.)
!     and generates jvals from table, for up to 56 species.
!
!     Critical INDEX OPTION:  jaer index (6) jcld index (9)
!     identifies aerosol optical depth index (always followed by SSA)
!     and cloud-above optical depth (followed by cloud-below and alt.
 
!     adju NOTE:  for future VECTORIZATION:  variables listed as (   n)
!     may be switched to ( kk,n)
!     Look out for indices:  ig11, etc.
 
!     INPUTS (parameters for specified case)
!     INPUT ARRAY DATA
!     (note: current scale jarray (k,ig,jc):  k=1,56)
!     OUTPUTS
 
!     INTERNAL
 
!     cfac = 2CLOUD ABOVE-BELOW FACTOR  (link to OPTION below)
 
!     CORRECT FACTORS BASED ON JTAB AND (jclb*jclbm-1)*(1-jcla*jclam)
      data cfac/0.000D+00 , 6.175D-01 , 2.079D+00 , 1.774D+00 , 2.407D+00 , &
                2.479D+00 , 2.365D+00 , 1.495D+00 , 0.000D+00 , 0.000D+00 , &
                1.424D+00 , 1.732D+00 , 1.180D+00 , 1.202D+00 , 1.373D+00 , &
                1.538D+00 , 9.732D-01 , 0.000D+00 , 0.000D+00 , 1.228D+00 , &
                1.911D+00 , 1.831D+00 , 8.667D-01 , 1.481D+00 , 1.170D+00 , &
                1.336D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , 0.000D+00 , &
                0.000D+00/
 
!     ORIGINAL (OLD) FACTORS CALCULATED BASED  ON JTAB (old correct)
!     (4=1.319, 1.699,  1.839 2.147)
!     data  cfac/ 0.,      5.973E-01,2.025E+00,1.839E+00,2.322E+00,
!     *           2.397E+00,       0.,1.557E+00,       0.,       0.,
!     *           1.470E+00,1.809E+00,1.187E+00,1.204E+00,1.400E+00,
!     *           1.602E+00,9.530E-01,       0.,0.000E+00,1.234E+00,
!     *           1.926E+00,1.877E+00,8.461E-01,1.531E+00,1.177E+00,
!     *           1.373E+00, 0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0./
 
!     FACTORS CALCULATED BASED  ON JBASE (OLD)
!     data  cfac/ 0., 4.896E-01, 1.056E+00, 1.319E+00, 1.134E+00,
!     *           1.215E+00,   0., 1.255E+00,  0.     ,   0.      ,
!     *           1.189E+00,1.343E+00,9.741E-01,9.885E-01,1.145E+00,
!     *           1.281E+00, 7.782E-01,   0.   ,0.000E+00,1.011E+00,
!     *           1.218E+00,1.281E+00,6.906E-01,1.227E+00, 9.677E-01,
!     *           1.122E+00, 0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!     *           0.,0.,0.,0.,0.,0./
 
!     -----------------------------------------------------------------
 
!     VECTOR PARAMETER - MAYBE ADD LATER
!     kk=1
 
!     TEST WRITE INDEX (=1 to write)
!     (=2 for SPECIAL WRITE TO 58-for 2CLOUD ABOVE-BELOW
!     PARAMETERIZATION)

      iwri = 0
 
!     INDEX OPTION:  jaer (=6), jcld (originally =9) jalb (=8)
      jaer = 6
      jcld = 9
      jalb = 8
      jtem = 13
 
!     OPTION:  TOTAL NUMBER OF J-VALUES TO BE READ. (jct=4 test 26 trop
!     56 f jct=4
      jct = 26
 
!     PRELIMINARY:  ENTER JVALS AS ZERO.
!     RETURN IF ZENITH>94.  (CUT IF VECTORIZED)
      jval(:) = d_zero
      if ( jparam(1) >= hvmat(1,nhv(1)) ) return
 
!     LOOP:  Establish index and fractions for each j-parameter
!     Note:  allow for hvmat intervals monotonically increasing or
!     decreasi Note:  SKIP (but do not exit) for nhv=1.
!     NOTE:  Loop was 1,20; before TEMP.  Does DATE use this?  Probably
!     not do i=1,20
      jparamloop: &
      do i = 1 , 19
        jfx(i) = d_zero
        if ( nhv(i) <= 0 ) exit
        if ( nhv(i) /= 1 ) then
 
!         SPECIAL TEMPERATURE:
!         (not here!  AFTER jfx() set, reset.)
!         jfx(jtem)=2.  Only  one value.
 
!         Test to exit if matrix intervals are zero
!         (note:  matrix must be monotonically increasing or decreasing)

          if ( dabs(hvmat(i,1)-hvmat(i,nhv(i))) > dlowval ) then
            do j = 1 , (nhv(i)-1)
              if ( dabs(hvmat(i,j)-hvmat(i,j+1)) < dlowval ) then
                exit jparamloop
              end if
            end do
 
!           Enter fraction for parameter outside matrix range
            if ( hvmat(i,1) < hvmat(i,nhv(i)) ) then
              if ( jparam(i) <= hvmat(i,1) ) jfx(i) = d_one
              if ( jparam(i) >= hvmat(i,nhv(i)) ) jfx(i) = dble(nhv(i))
            end if
            if ( hvmat(i,1) > hvmat(i,nhv(i)) ) then
              if ( jparam(i) >= hvmat(i,1) ) jfx(i) = d_one
              if ( jparam(i) <= hvmat(i,nhv(i)) ) jfx(i) = dble(nhv(i))
            end if
 
!           Enter fraction for parameter inside matrix range
            do j = 1 , (nhv(i)-1)
              if ( (jparam(i) >= hvmat(i,j) .and. &
                    jparam(i) <= hvmat(i,j+1)) .or. &
                   (jparam(i) <= hvmat(i,j) .and. &
                    jparam(i) >= hvmat(i,j+1)) ) then
                jfx(i) = dble(j) + &
                     (jparam(i)-hvmat(i,j))/(hvmat(i,j+1)-hvmat(i,j))
              end if
            end do
 
!           SPECIAL TEMPERATURE:
!           jfx(jtem)=2.  Only  one value is entered in hvarray matrix
!           for jtem This represents  delta  parameter for t(ialt)+10
!           degrees. jfx(jtem)=2 ensures that  jtem interpolation will 
!           select jfrac(jte as zenith/altitude interpolation from this
!           one data set. (jfx(jtem)=1 would correspond to base case
 
!           parameter value = hvmatz
            if ( i == jtem ) jfx(i) = d_two
 
!           TEST WRITE
            if ( iwri == 1 ) then
              write (57,99001) i , jparam(i) , jfx(i)
              write (57,99002) (hvmat(i,j),j=1,nhv(i))
            end if
          end if
        end if
 
      end do jparamloop ! establish index and fraction for j-parameters
 
!     Establish ig parameters for base interpolation:  zenith and
!     altitude zenith and altitude are controlled by ig index in
!     jarray(k,ig,jc) ig11=lower zenith index, lower altitude index
!     ig12=upper zenith index, lower altitude index
!     ig21=lower zenith index, upper altitude index
!     ig22=upper zenith index, upper altitude index
 
!     ig01=lower zenith index, surface altitude index
!     ig02=upper zenith index, surface altitude index
 
!     NOTE:  TO VECTORIZE:  fzen, falt, izen, ialt, ig11, etc
!     must all be vectors.  Skip.
 
      izen = idint(jfx(1))
      if ( izen == nhv(1) ) izen = izen - 1
      fzen = d_one + dble(izen) - jfx(1)
      ialt = idint(jfx(2))
      if ( ialt == nhv(2) ) ialt = ialt - 1
      falt = d_one + dble(ialt) - jfx(2)
 
      ig11 = izen + nhv(1)*(ialt-1)
      ig12 = ig11
      if ( izen < nhv(1) ) ig12 = ig11 + 1
      ig21 = ig11
      if ( ialt < nhv(2) ) ig21 = ig11 + nhv(1)
      ig22 = ig21
      if ( izen < nhv(1) ) ig22 = ig21 + 1
 
      ig01 = izen
      ig02 = ig01
      if ( izen < nhv(1) ) ig02 = ig01 + 1
 
!     Establish base j-values.  Enter into jval output.
      do jc = 1 , jct
        jval(jc) = jarray(1,ig11,jc)*fzen*falt + jarray(1,ig12,jc) *        &
                  (d_one-fzen)*falt + jarray(1,ig21,jc)*fzen*(d_one-falt) + &
                   jarray(1,ig22,jc)*(d_one-fzen)*(d_one-falt)
!
!       NONLINEAR CORRECTION for variation with altitude - NO EFFECT
!       (effects HNO3 at 24 km, otherwise zero effect.)
        if ( jval(jc) > 0 .and. ialt <= (nhv(2)-2) ) then
          x = d_two*jarray(1,ig21,jc)/(jarray(1,ig11,jc) + &
                                       jarray(1,(ig21+nhv(1)),jc)) - d_one
          jval(jc) = jval(jc)*(d_one+x*falt*(d_one-falt))
        end if
      end do
 
!     TEST WRITE
      if ( iwri == 1 ) then
        write (57,99003) izen , ialt , ig11 , ig12 , ig21 , ig22 ,      &
                         fzen , falt
        write (57,99004) jval(1) , jarray(1,ig11,1) , jarray(1,ig12,1) ,&
                         jarray(1,ig21,1) , jarray(1,ig22,1)
        write (57,99004) jval(2) , jarray(1,ig11,2) , jarray(1,ig12,2) ,&
                         jarray(1,ig21,2) , jarray(1,ig22,2)
        write (57,99004) jval(4) , jarray(1,ig11,4) , jarray(1,ig12,4) ,&
                         jarray(1,ig21,4) , jarray(1,ig22,4)
      end if
 
!     LOOP: Establish fractional adjustments
!     LOOP for each fraction.  Set k-index
      k = 1
      do i = 3 , 19
 
!       Set initial fraction equal to base value
!       (watch special 11, 12! Make sure base value is OK if nhv=1 and
!       skipp
        do jc = 1 , jct
          jfrac(i,jc) = d_one
        end do
 
!       if(i.eq.6)then
        if ( i == jaer ) then
          do jc = 1 , jct
            jfaerz(jc) = d_one
          end do
        end if
 
!       if(i.eq.9)then
        if ( i == jcld ) then
          do jc = 1 , jct
            jfsur(jc) = d_one
          end do
        end if
 
!       If nhv=1, skip this loop.  Automatically enter jfract=1.
!       Note if nhv=1, values are omitted from jarray(k,ig,jc)
        if ( nhv(i) /= 1 ) then
!         If nhv=0, exit.
          if ( nhv(i) <= 0 ) exit
 
!         SPECIAL TEMPERATURE:
!         Just one array  value is entered for temperature array.
!         Set  jfx(i) so that jfx(i)=nhv(i)=2.
!         This insures that jfrac(jtem)=direct interpol.  from single
 
!         value. Establish adjustment factor and k-index for jarray
!         (VECTORIZE)
          kn = k + idint(jfx(i))
          if ( idint(jfx(i)) == nhv(i) ) kn = kn - 1
          fkn = d_one + dble(kn-k) - jfx(i)
!
!         TEST WRITE:
          if ( fkn < d_zero .or. fkn > d_one ) then
            if ( iwri == 1 ) then
              write (56,99005) i , k , kn , jfx(i) , fkn
            end if
          end if
 
!         Establish fraction from jarray.
!         Use zenith interpolation, itself interpolated between two
!         fractions.
          do jc = 1 , jct
            jfrac(i,jc) = fkn*(jarray(kn,ig11,jc)*fzen*falt + &
                               jarray(kn,ig12,jc)*(d_one-fzen)*falt + &
                               jarray(kn,ig21,jc)*fzen*(d_one-falt) + &
                               jarray(kn,ig22,jc)*(d_one-fzen)*(d_one-falt)) + &
                  (d_one-fkn)*(jarray(kn+1,ig11,jc)*fzen*falt + &
                               jarray(kn+1,ig12,jc)*(d_one-fzen)*falt + &
                               jarray(kn+1,ig21,jc)*fzen*(d_one-falt) + &
                               jarray(kn+1,ig22,jc)*(d_one-fzen)*(d_one-falt))
          end do
 
!         Establish surface fraction-FOR SURFACE IMPACT OF CLOUD BELOW.
!         This uses the CLOUD-BELOW value of cloud optical depth
!         (from the 10th index)
!         but it uses it with the CLOUD-ABOVE value (9th index)
!         Representing the CLOUD-BELOW impact (as CLOUD ABOVE) at the
!         surface. THIS IS USED TO ADJUST THE ALBEDO FRACTION (8)
 
!         ALSO - EXPERIMENTAL - USED TO ADJUST CLOUD-ABOVE
!         TO CORRECT UNDERESTIMATE IN CASE WITH CLOUD-ABOVE AND
!         CLOUD-BELOW. (totally empirical, no reason.)
 
!         if(i.eq.10) then
          if ( i == jcld+1 ) then
!           TEST WRITE
            if ( iwri == 1 ) then
              write (57,99006) jfrac(jalb,1) , jfrac(jalb,2) , jfrac(jalb,4)
            end if
            do jc = 1 , jct
              jfsur(jc) = fkn*(jarray(kn-nhv(jcld),ig01,jc)*fzen + &
                               jarray(kn-nhv(jcld),ig02,jc)*(d_one-fzen)) + &
                  (d_one-fkn)*(jarray(kn+1-nhv(jcld),ig01,jc)*fzen + &
                               jarray(kn+1-nhv(jcld),ig02,jc)*(d_one-fzen))
!             ALBEDO
!             jfrac(8,jc)=1-jfsur(jc)*(1.-jfrac(8,jc))
              jfrac(jalb,jc) = d_one - jfsur(jc)*(d_one-jfrac(jalb,jc))
!             OLD 2CLOUD ABOVE-BELOW CORRECTION FACTOR -   (0, 0.2,
!             0.4*) c    
!             jfrac(9,jc)=1.-(1-jfrac(9,jc))*(0.8+0.2*jfsur(jc)) c    
!             jfrac(9,jc)=1.-(1-jfrac(9,jc))*(0.6+0.4*jfsur(jc))
!             ORIGINAL 2CLOUD ABOVE-BELOW CORRECTION FACTOR -         
!             (see data abo WAS HERE, MOVED BELOW c    
!             jfrac(9,jc)=jfrac(9,jc) c    *       *(1.+
!             cfac(jc)*(1.-jfsur(jc))*(1.-jfrac(9,jc)) )
            end do
!           TEST WRITE
            if ( iwri==1 ) then
              write (57,99007) jfsur(1) , jfsur(2) , jfsur(4) , &
                 jfrac(jalb,1) , jfrac(jalb,2) , jfrac(jalb,4)
            end if
!           SPECIAL WRITE FOR CREATING 2CLOUD ABOVE-BELOW CORRECTION 
!           FACTOR. was here, moved below.
!           END IF - Establish surface fraction and adjust albedo, 
!           cloud above-be
          end if
 
!         Establish AEROSOL ZERO fraction:  case i=6, 1st fraction
!         value. fraction k index is k+1 instead of kn
          if ( i == jaer ) then
            do jc = 1 , jct
              jfaerz(jc) = (jarray(k+1,ig11,jc)*fzen*falt + &
                            jarray(k+1,ig12,jc)*(d_one-fzen)*falt + &
                            jarray(k+1,ig21,jc)*fzen*(d_one-falt) + &
                            jarray(k+1,ig22,jc)*(d_one-fzen)*(d_one-falt))
            end do
          end if
 
!         TEST WRITE
          if ( iwri == 1 ) then
            write (57,99008) i , k , kn , jfx(i) , fkn
            write (57,99009) i , kn , fkn , fzen , falt , jfrac(i,1)
            write (57,99009) i , kn , fkn , fzen , falt , jfrac(i,2)
            write (57,99009) i , kn , fkn , fzen , falt , jfrac(i,4)
            write (57,99010) jarray(kn,ig11,1) , jarray(kn,ig12,1) ,    &
                             jarray(kn,ig21,1) , jarray(kn,ig22,1) ,    &
                             jarray(kn+1,ig11,1) , jarray(kn+1,ig12,1) ,&
                             jarray(kn+1,ig21,1) , jarray(kn+1,ig22,1)
            write (57,99010) jarray(kn,ig11,2) , jarray(kn,ig12,2) ,    &
                             jarray(kn,ig21,2) , jarray(kn,ig22,2) ,    &
                             jarray(kn+1,ig11,2) , jarray(kn+1,ig12,2) ,&
                             jarray(kn+1,ig21,2) , jarray(kn+1,ig22,2)
            write (57,99010) jarray(kn,ig11,4) , jarray(kn,ig12,4) ,    &
                             jarray(kn,ig21,4) , jarray(kn,ig22,4) ,    &
                             jarray(kn+1,ig11,4) , jarray(kn+1,ig12,4) ,&
                             jarray(kn+1,ig21,4) , jarray(kn+1,ig22,4)
          end if
 
!         Special adjustment:  cloud height adjustment factors (j=11,
!         12) represent an adjustment to the cloud-above and
!         cloud-below (j=9,10). These are converted to straight
 
!         fractions here. The adjustment parameter in the array is:
!         Fadj=(1-Ftot/Fcloud)/(1-Fcloud)
!         where Ftot is combined fraction, Fcloud is cloud-alone
!         fraction. Here, the adjustment factor Fadj is replaced with
!         F': Where Ftot=Fc*F'; F'=1-Fadj(1-Fcloud);  F' limited,
 
!         between 0.1 and 1 The same adjustment is applied for AEROSOL
!         SSA (j=7) as an adjustment to the AEROSOL OPTICAL DEPTH
 
!         fraction (j=6) if(i.eq.11.or.i.eq.12          ) then
          if ( i == jcld+2 .or. i == jcld+3 ) then
            ij = i - 2
            do jc = 1 , jct
              jfrac(i,jc) = d_one - jfrac(i,jc)*(d_one-jfrac((ij),jc))
              if ( jfrac(i,jc) < 0.1D0 ) jfrac(i,jc) = 0.1D0
              if ( jfrac(i,jc) > 10.0D0 ) jfrac(i,jc) = 10.0D0
            end do
!           TEST WRITE
            if ( iwri == 1 ) then
              write (57,99011) i , jfrac(i,1)
              write (57,99011) i , jfrac(i,2)
              write (57,99011) i , jfrac(i,4)
            end if
          end if
 
!         Special adjustment:  aerosol SSA (i=7)
!         represent an adjustment to the zero-aerosol fraction (i=6,
!         jfaerz) These are converted to straight fractions here.
 
!         The adjustment parameter in the array is:
!         Fadj=(1-Ftot/Faerbase)/(1-Faerbase/Faerzero)
!         where Ftot is combined fraction, Faerbase is base aerosol and
!         SSA f and Faerzero is fraction for zero aerosol case.
!         (This is necessary because initial value, F=1, does not have
!         zero a Here, the adjustment factor Fadj is replaced with F':
!         Where Ftot=Fc*F'; F'=1-Fadj(1-Faerbase/Faerzero);
!         F' limited, between 0.1 and 10.
 
!         if(i.eq.7) then
          if ( i == jaer+1 ) then
!           ij=6
            ij = jaer
            do jc = 1 , jct
              if ( jfaerz(jc) > 0 ) then
                jfrac(i,jc) = d_one - jfrac(i,jc) * &
                            (d_one-jfrac((ij),jc)/jfaerz(jc))
              else
                jfrac(i,jc) = d_one - jfrac(i,jc)*(d_one-jfrac((ij),jc))
              end if
              if ( jfrac(i,jc) < 0.1D0 ) jfrac(i,jc) = 0.1D0
              if ( jfrac(i,jc) > 10.0D0 ) jfrac(i,jc) = 10.0D0
            end do
!           TEST WRITE
            if ( iwri == 1 ) then
              write (57,99011) i , jfrac(i,1)
              write (57,99011) i , jfrac(i,2) , jfaerz(2)
              write (57,99011) i , jfrac(i,4)
            end if
          end if
 
!         Special adjustment:  albedo fraction in case of cloud-below.
!         Adjust fractional change (relative to 1) based on
!         cloud-below impact on surface:
 
!         (To do this:  establish special jfrac using CLOUD-BELOW jfx,
!         but for CLOUD-ABOVE and for SURFACE.
!         Then adjust albedo fraction:  falb'=1-(1-falb)*fclsurf
!         DONE WITH jfsur() ABOVE.
!         --------------------
 
!         2CLOUD ABOVE-BELOW NONLINEAR ADJUSTMENT:  OPTION.
!         For case with BOTH cloud-above and cloud-below,
!         table fractions underestimate j-values.
!         THIS IS TOTALLY EMPIRICAL-no reason for it.
!
!         Cloud-above fraction is corrected:
!         CLA' = F*(1-CLA*CLAm)*(CLB*CLBm-1)
!         where CLA, CLAm=f9,  f11 (after modification);  CLB,
!         CLBm=f10,f12. MOVED HERE - PREVIOUSLY w/ jfsur() ABOVE.
 
!         if(i.eq.12) then
          if ( i == jcld+3 ) then
!           TEST WRITE
            if ( iwri == 1 ) write (57,99012) jfrac(jcld,1) ,         &
                                  jfrac(jcld,2) , jfrac(jcld,4) ,     &
                                  jfrac(jcld+1,4) , jfrac(jcld+2,4) , &
                                  jfrac(jcld+3,4)
            do jc = 1 , jct
!             2CLOUD AB0VE-BELOW OPTION:  MODIFY HERE..  (use with data
!             OPTION, abov FINAL CORRECT VERSION
!             with INDEX CONTROL:
!             WITHOUT INDEX
              if ( jfrac(jcld,jc)*jfrac(jcld+2,jc) < d_one ) then
                jfrac(jcld,jc) = jfrac(jcld,jc) * &
                  (d_one+cfac(jc)*(jfrac(jcld+1,jc) * &
                   jfrac(jcld+3,jc)-d_one)*(d_one-jfrac(jcld,jc) * &
                   jfrac(jcld+2,jc)))
              end if
!             if(jfrac(9,jc)*jfrac(11,jc).lt.1)
!             *   jfrac(9,jc)=jfrac(9,jc) *(1.+
!             *       cfac(jc)*(jfrac(10,jc)*jfrac(12,jc)-1.)*
!             *                  (1.-jfrac(9,jc)*jfrac(11,jc)) )
!             OLDER-ORIGINAL MODIFICATION.  (use with data OPTION,
!             above) jfrac(9,jc)=jfrac(9,jc) *(1.+
!             *       cfac(jc)*(1.-jfsur(jc))*(1.-jfrac(9,jc)) )
            end do
!           TEST WRITE
            if ( iwri == 1 ) write (57,99013) jfrac(jcld,1) ,  &
                                 jfrac(jcld,2) , jfrac(jcld,4)
          end if
!
!         SPECIAL WRITE FOR CREATING 2CLOUD ABOVE-BELOW CORRECTION 
!         FACTOR. (Save table:  jfsur, jf9,jtab-base, jtuv, jtab for 26
!         species. Then cloud  factor f=
!         (jtuv-jtab)/[jtab*jsurf*(1-jf9)] where j's are summed over
!         table.
          if ( i == jcld+3 ) then
            if ( iwri >= 1 ) write (58,99014) (jfsur(jc),jc=1,26)
            if ( iwri >= 1 ) write (58,99014) (jfrac(jcld,jc),jc=1,26)
            if ( iwri >= 1 ) write (58,99014) (jfrac(jcld+1,jc),jc=1,26)
            if ( iwri >= 1 ) write (58,99014) (jfrac(jcld+2,jc),jc=1,26)
            if ( iwri >= 1 ) write (58,99014) (jfrac(jcld+3,jc),jc=1,26)
            if ( iwri >= 1 ) write (58,99014) (jval(jc),jc=1,26)
            if ( iwri == 1 ) write (57,99014) (jfsur(jc),jc=1,26)
            if ( iwri == 1 ) write (57,99014) (jfrac(jcld,jc),jc=1,26)
            if ( iwri == 1 ) write (57,99014) (jfrac(jcld+1,jc),jc=1,26)
            if ( iwri == 1 ) write (57,99014) (jfrac(jcld+2,jc),jc=1,26)
            if ( iwri == 1 ) write (57,99014) (jfrac(jcld+3,jc),jc=1,26)
            if ( iwri >= 1 ) write (57,99014) (jval(jc),jc=1,26)
          end if
 
!         SPECIAL TEMPERATURE ADJUSTMENT:
!         Just one array  value is entered for temperature array.
!         Above, set  jfx(i) so that jfx(i)=nhv(i)=2.
!         This insures that jfrac(jtem)=direct interpol.  from single
 
!         value. Here, establish temperature parameter from jparam(jtem)
!         jparam(jtem)= dtem if less than  50;
!         else jparam(jtem)=temp, and dtem found as jparam(jtem)-t(alt)
!         where t(alt), std temp, found  from interpolating
!         hvmat(jtem,j),
!         Then  jfrac(jtem) = 1 + 0.1*dtem*jf(jtem); jf(jtem) is factor
!         from t (previously  entered as jfrac(jtem), interpolated for 
!         zenith and (0.1 factor because jf(jtem) is % change for 10
 
!         degree increase.) OPTION - CHANGE JPARAM=dtem FOR TEST
 
          if ( i == jtem ) then
!           TEST WRITE
            if ( iwri == 1 ) write (57,99015) i , jparam(i) , &
                                  jfrac(i,2) , jfrac(i,4)
 
            x = jparam(i)
            if ( jparam(i) > 50 ) then
              x = jparam(i) - (hvmat(jtem,ialt) * &
                      falt+hvmat(jtem,ialt+1)*(d_one-falt))
!             OPTIONAL LINE FOR TESTS
!             jparam(i)=x
!             TEST WRITE
              if ( iwri == 1 ) write (57,99016) ialt , falt , &
                                   hvmat(jtem,ialt) , hvmat(jtem,ialt+1)
            end if
            do jc = 1 , jct
              jfrac(jtem,jc) = d_one + 0.1D0*x*jfrac(jtem,jc)
            end do
 
!           TEST WRITE
            if ( iwri == 1 ) write (57,99017) x , &
                                  jfrac(jtem,2) , jfrac(jtem,4)
 
!           END SPECIAL TEMPERATURE ADJUSTMENT
          end if
 
!         -------
!         ADJUST BASE -JVALUE BY FRACTION - WAS HERE, MOVED TO SEPARATE
!         LOOP BEL ------
 
!         ADVANCE K-COUNTER FOR NEXT LOOP
!         (due to screw-up in TUVGRID1, advance even for NHV=1? No.)
          if ( nhv(i) > 1 ) k = k + nhv(i)
        end if
!       END LOOP:  fractional adjustments
      end do
 
!     TEST J-VALUE  FINAL FRACTIONAL ADJUSTMENT
      if ( iwri == 1 ) write (57,99018) jval(2) , jval(4)
 
!     LOOP TO ADJUST BASE J-VALUE BY FRACTION
 
      do i = 3 , 19
        if ( nhv(i) > 1 ) then
!         ADJUST BASE J-VALUE BY FRACTION
          do jc = 1 , jct
            jval(jc) = jval(jc)*jfrac(i,jc)
          end do
!         TEST J-VALUE  FINAL FRACTIONAL ADJUSTMENT
          if ( iwri == 1 .and. &
               dabs(jfrac(i,1)-d_one) < dlowval .and. &
               dabs(jfrac(i,2)-d_one) < dlowval ) then
            write (57,99019) i , jfrac(i,2) , jval(2) , jfrac(i,4) , jval(4)
          end if
        end if
      end do
 
!     DATE ADJUSTMENT:
!     ASSUME THAT INPUT jparam(20)= DATE  (decimal)
!     either as date factor (-1 to +1)
!     or DAY NUMBER (1-365)
!     or YYMMDD (not 00)
!     DAY FACTOR (X) IS cos(nd*2.*pi/365)
 
      x = jparam(20)
      if ( x > d_one .and. x < 10000.0D0 ) then
        x = dcos(d_two*mathpi*(jparam(20)/dayspy))
      end if
      if ( x >= 10000.0D0 ) then
        iy = idint((x+0.001D0)/10000.0D0)
        im = idint((x+0.001D0-10000.0D0*dble(iy))/100.0D0)
        id = idint((x+0.001D0-10000.0D0*dble(iy)-100.0D0*dble(im)))
        x = dcos(d_two*mathpi*(dble(id+30*(im-1))/dayspy))
      end if
 
!     ADJUSTMENT:  x IS DAY FACTOR, -1. to +1.  NOW MAKE DAY ADJUSTMENT.
      x = d_one + 0.0344D0*x
!     ADJUST BASE J-VALUE FOR DATE
      do jc = 1 , jct
        jval(jc) = jval(jc)*x
      end do

99001 format (' TEST I JPARAM JFX (then hvmat(i,j)):',i5,2F10.3)
99002 format (8F10.3)
99003 format (' TEST izen ialt ig11 ig12 ig21 ig22 fzen falt=',/,6I4,2F10.3)
99004 format (' TEST INITIAL JVAL:',8(1pe10.3))
99005 format ('ERROR fkn(<0 or>1): i k kn jfx fkn(0-1)=',3I3,2F10.4)
99006 format ('ALBEDO ADJUSTMENT FOR CLOUD-BELOW:',/,   &
              'PRIOR    ALBEDO FR (O2 O1D NO2)=',3F12.6)
99007 format ('SURFACE   CLOUD FR (O2 O1D NO2)=',3F12.6,/,              &
             &'MODIFIED ALBEDO FR (O2 O1D NO2)=',3F12.6)
99008 format ('TEST INDEX I K KN JFX(I): kn=k+int(jfx):',3I3,2F10.4)
99009 format (' TEST FRACTION i,kn, fkn fzen falt jfrac =',/,2I5,5F10.4)
99010 format (8(1pe10.3))
99011 format (' TEST ADJ. FRACTION i, jfrac =',i5,2F10.4)
99012 format ('2CLOUD ABOVE-BELOW ADJUSTMENT:    ',/,                   &
              'PRIOR CLOUD-ABOVE  (O2 O1D NO2)=',3F12.6,/,              &
              'CLD-BELOW, CLaltA, CLaltB (NO2)=',3F12.6)
99013 format ('ADJUSTED CLOUD-ABOVE  (O2 O1D NO2)=',3F12.6,/,           &
              'jfrac10, 11, 12 (CLB,m)   (NO2)=',3F12.6)
99014 format (8(1pe10.3))
99015 format (' TEST BEFORE TEMPERATURE ADJ: i, jparam,jfrac O1D,NO2=', &
              i5,4F10.4)
99016 format (' TEST ialt falt hvmat(jt,ialt), hvmat(jt,ia+1)=',i5,     &
              5F10.4)
99017 format (' TEST AFTER  TEMPERATURE ADJ: dtem, fO1D, fNO2=',f10.4,  &
              2(1pe10.3))
99018 format (' BEFORE-FRACTION jO1D,jNO2=',2(1pe10.3))
99019 format (' FINAL J-CALC: F2,J2,F4,J4=',i3,2(f10.3,1pe10.3))
 
    end subroutine jvalpro

! SUCCESSFUL VERSION, 11-22-02 w/MODIFICATIONS.  Use with jvalmain1.f.
!  (use jval2.f, altitude in KM, to test -w/ tuvtest2.f)
!
!  PROGRAM FOR PROCESSING TUV TABLE (TUVGRID)
!  Uses table generated by tuvtab2.f
!   with subroutines for use in testing program and photochemistry.
!
!   MODIFICATIONS:
! WITH SUNSET ZENITH:  change 442 to 510 (FEB 2003)
! WITH ADD FOR TEMPERATURE ADJUSTMENT (jtem=13)
!   DATE, adjustment for AEROSOL SSA, 2CLOUD ABOVE-BELOW. TEMPERATURE.
!   ALSO, ALLOWS NONZERO AEROSOL AS BASE.
!   Aerosol SSA adjustment is relative to aerosol=0, 1st aerosol-not nec
!   Modified ALBEDO, buffered by surface-cloud-below fraction, included.
!    OPTION -  CLOUD ABOVE+BELOW CORRECTION FACTOR - 0.4 INCLUDED (0,0.2
!
! SUCCESSFUL VERSION, 9-26-02, used in TUVtest2.
!  (minor future fix:  cloud-above at alt=8, make  equal to alt=7
!    to allow for high  cloud).
!
! OPTION - altitude as kPa (2) or as km (21) (in readhv-ALTITUDE OPTION)
! OPTION - number of j-values to be read:  4 (test) 26 (trop) 56 (full)
!    (see jct in both subroutines)
! TEST WRITES - see iwri
!
! NOTE:  CHANGES FOR TEST GRID (TUVGRID1) vs FULL GRID (TUVTEST2):
!  TUVGRID1 (original test) HAS JUST 4 SPECIES (jc=1,4)
!  AND ALSO SKIP CLOUD-ABOVE AFTER 7 TROPOSPHERIC LAYERS. (i.le.7)
! TUVGRID2 (full) HAS 26 SPECIES (jc=1,26)
!   AND SKIPS CLOUD-ABOVE AFTER 8 TROPOSPHERIC LAYERS (i.le.8)
!   THESE MUST BE CHANGED IN jtab1.f - see CHANGE GRID OPTION
! -----------------------------------------------------------------
 
    subroutine readhv(lsin,nhv,hvmat,hvmatb,jarray)
      implicit none
!
      integer(ik4) :: lsin
      real(rk8) , dimension(22,40) :: hvmat
      real(rk8) , dimension(22) :: hvmatb
      real(rk8) , dimension(80,510,56) :: jarray
      integer(ik4) , dimension(22) :: nhv
      intent (in) lsin
      intent (inout) hvmat , hvmatb , jarray , nhv
!
      character(4) :: aaa
      real(rk8) , dimension(22) :: hvmatz
      integer(ik4) :: i , ig , iwri , iz , j , jaer , jalb , jc , jcld ,     &
                 jct , jtem , k , m , n , nmax
      real(rk8) :: x , y
!
!     Output matrix:  jarray(k,ig,jc) for
!     k=cases:  k=1 base case, k>1 adjustment factors
!     for cases identified in hvmatrix.
!     ig=matrix of altitudes and zenith angles
!     jc=species j-values
 
!     Includes AUTOMATIC FILL-IN OF LAST ZENITH VALUE (nighttime)
!     AND BASE CASE IN ADJUSTMENT FACTORS.
!     (read from input table.  Base j-values are
!     0 (nighttime) for base case j-values;
!     1 for most adjustment factors, but
!     0 for cloud adjustment factors 11 and 12
!     (which are relative to original cloud factor)
 
!     INPUTS (none)
!     OUTPUTS
!     (note: current scale jarray (k,ig,jc):  k=1,56)
 
!     INTERNAL
!     --------------------------------------------------
!     TEST WRITE:  WRITE IF iwri=1
      write (*,*) 'HVREAD     TUVGRID'
      iwri = 0
 
!     OPTION:  TOTAL NUMBER OF J-VALUES TO BE READ. (jct=4 test 26 trop
!     56 f CHANGE GRID:  TUVGRID1 (jct=4) vs TUVGRID2 (jct=26)
!     see also CHANGE GRID below.
!     jct=4
      jct = 26
 
!     INDEX OPTION:  jaer (=6), jcld (originally =9) jalb (=8)
      jaer = 6
      jcld = 9
      jalb = 8
      jtem = 13
!
!     ZERO INPUT MATRIX
      do i = 1 , 22
        nhv(i) = 0
        hvmatb(i) = d_zero
        do j = 1 , 40
          hvmat(i,j) = d_zero
        end do
      end do
 
!     READ INPUT CONTROLS  (note read index, k, as well)
      rewind lsin
      do i = 1 , 22
        read (lsin,99001) aaa
        if ( iwri == 1 ) write (57,99001) aaa
        read (lsin,99002) k , j , x , y
        if ( iwri == 1 ) write (57,99002) k , j , x , y
        if ( j == 0 ) exit
        if ( k == 0 ) exit
        nhv(k) = j
        hvmatb(k) = x
        hvmatz(k) = y
        read (lsin,99003) (hvmat(k,n),n=1,nhv(k))
        if ( iwri == 1 ) write (57,99003) (hvmat(k,n),n=1,nhv(k))
      end do
 
!     ALT  ALTITUDE OPTION:  kPa or km
      do j = 1 , 40
        hvmat(22,j) = hvmat(2,j)
!       ALTITUDE OPTION;  FOR KM, UNCOMMENT THIS LINE.  FOR KPA,
!       COMMENT OUT hvmat(2,j)=hvmat(21,j)
!       END OPTION
      end do
 
!     SPECIAL TEMPERATURE TREATMENT:
!     INITIAL  READIN HAS nhv(jtem)=nhv(ialt);  and reads hvmat(jtem,j)
!     as standard temperature vs. altitude.
!     DATA READIN (jarray, below) reads  just one data  array,
!     representing parameter for std temp+10 degrees.
!     FOR THIS, RESET nhv(jtem)=2.
!     (hvmatb=hvmat(jtem,1); hvmatz=0.)
 
      if ( jtem > 0 ) then
        if ( nhv(jtem) > 2 ) nhv(jtem) = 2
      end if
 
!     SET INITIAL J-VALUES TO ZERO/ONE
!     J=0 FOR BASE CASE (k=1) AND CLOUD-ADJUSTMENT VALUES (k=11,12)
      do k = 1 , 80
        do ig = 1 , 510
          do jc = 1 , 56
            jarray(k,ig,jc) = d_one
            if ( k == 1 ) jarray(k,ig,jc) = d_zero
          end do
        end do
      end do
 
!     jarray(80,510,56)
 
!     LOOP TO READ ARRAY BASE CASE (k=1) AND ADJUSTMENT CASES (k>1)
!     N REPRESENTS PARAMETER VALUE FROM J-MATRIX.
!     READ FROM ARRAY AND ENTER IMPLICIT VALUES OMITTED FROM ARRAY
      k = 0
      do m = 2 , 20
        nmax = 1
        if ( m >= 3 ) nmax = nhv(m)
        if ( iwri == 1 ) write (57,*) m , nmax , k
!       (Exit if nhv=0)
        if ( nhv(m) <= 0 ) exit
!       (Note: for nhv=1, omit from jarray)
        if ( nhv(m) > 1 ) then
          do n = 1 , nmax
            k = k + 1
            if ( iwri == 1 ) write (57,*) m , n , k
!           ENTER IMPLICIT VALUES FOR ARRAY (nighttime or base case
!           adjustment fac
            do ig = 1 , 510
              do jc = 1 , 56
                jarray(k,ig,jc) = hvmatz(m)
              end do
            end do
 
!           READ ARRAY VALUES (skip for nighttime and base-case
!           adjustment factors (skip for cloud-above at alt>7)
 
!           altitude and zenith loop (skips nighttime zenith)
!           with control to skip base case adjustment factor
!           and to skip for ialt>7 for cloud-above (9, 11:  jcld,
!           jcld+2) NOTE CHANGE with added 1km layer:  ialt>8 skip
 
!           SKIPS BASE CASE ADJUSTMENT FACTOR IN ALL SUBSEQUENT LOOPS
!           BUT INCLUDES IT IN FIRST LOOP (m=2, k=2 was early error)
!           if(k.eq.2.or.(hvmat(m,n).ne.hvmatb(m))) then
            if ( m == 2 .or. dabs((hvmat(m,n)-hvmatb(m))) > dlowval ) then
!             CASE HEADING
              read (lsin,99001) aaa
              if ( iwri == 1 ) write (57,*) k
              if ( iwri == 1 ) write (57,99001) aaa
              do i = 1 , nhv(2)
!               CHANGE GRID OPTION:   SKIP CONTROL - i.le.7 (TUVGRID1).
!               i.le.8 FOR TUV This  skips read for alt>8 for
!               cloud-above or cl-alt-above.
!               if(i.le.7.or.(m.ne.9.and.m.ne.11)) then
                if ( i <= 8 .or. (m /= jcld .and. m /= jcld+2) ) then
                  read (lsin,99001) aaa
                  if ( iwri == 1 ) write (57,*) k , i
                  if ( iwri == 1 ) write (57,99001) aaa
                  do iz = 1 , (nhv(1)-1)
                    ig = iz + nhv(1)*(i-1)
                    read (lsin,99004) (jarray(k,ig,jc),jc=1,jct)
                    if ( iwri == 1 ) write (57,99004) (jarray(k,ig,jc),jc=1,4)
                    if ( ig == 1 .and. iwri == 1 ) write (57,*) k
                    if ( k == 1 .and. iz <= 2 .and. &
                         i <= 2 .and. iwri == 1 ) write (57,*) k , i , iz , ig
                    if ( k == 1 .and. iz <= 2 .and. &
                         i <= 2 .and. iwri == 1 ) then
                      write (57,99004) (jarray(k,ig,jc),jc=1,4)
                    end if
                  end do
                end if
              end do
            end if
!           END LOOP TO READ
          end do
        end if
      end do
99001 format (a4)
99002 format (2I5,2F10.4)
99003 format (8F10.4)
99004 format (8(1pe10.3))
 
    end subroutine readhv

end module mod_cbmz_jval1
