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

module mod_che_sox

  use mod_intkinds
  use mod_realkinds
  use mod_che_bdyco
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  implicit none

  private

  public :: chemsox

  contains

    subroutine chemsox(j,wl,fracloud,fracum,rho,ttb)
     implicit none
!
     integer(ik4) , intent(in) :: j
     real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: ttb , wl , rho
     real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: fracloud , fracum
     real(rkx) :: rxs1 , rxs11 , rxs2 , rxs21 , chimol , cldno , &
                 clmin , remcum , oh1int , so2_rate , so2_avail , krembc
     real(rkx) , dimension(ntr) ::  wetrem , wetrem_cvc
     real(rkx) , dimension(ici1:ici2,kz) :: caircell , so2_snk , concmin
     real(rkx) :: rk_com(ici1:ici2,kz,100)
     real(rkx) :: h2o2mol

     integer(ik4) :: i , k

!FAB
     caircell(:,:) = 1.e-6_rkx * rho(:,:)/amdk * navgdr

     call chemrate(caircell,ttb,rk_com)

     ! clmin = non-precipitating cloud
     ! conversion threshold, clmin=0.01g/m3
     clmin = 0.01_rkx

     ! remcum= removal rate for cumulus
     ! cloud scavenging (s-1) (in-cloud and not grid level)
     remcum = 1.0e-3_rkx

     !---------------------------------------------
     !     SO2 G A Z E O U S    C O N V E R S I O N
     !       =====================================
     !
     ! The sink of SO2 is expressed in term of :
     ! SO2 + OH -----> Products
     ! dC[SO2]/dt = -k C[OH] C[SO2]
     ! C[SO2,t] = C[SO2,t0]*exp(-k C[OH] t)
     ! The variable name of the SO2 sink is snk_so2
     !---------------------------------------------

     do k = 1 , kz
       do i = ici1 , ici2

         cldno = d_one ! no cloud fraction

         if ( ioxclim == 1 ) then
           ! from the oxidant climatology
           oh1int = oxcl(j,i,k,iox_oh)
           if ( czen(j,i) < 0.001_rkx ) then
             oh1int = oh1int * 0.01_rkx
           else
             oh1int = oh1int * 1.99_rkx
           end if
         else
           oh1int = 30.0e5_rkx
           if ( czen(j,i) < 0.001_rkx ) then
             oh1int = oh1int * 0.01_rkx
           else
             oh1int = oh1int * 1.99_rkx
          end if
         end if

         ! Sink & Tendencies
         ! here p1 unit: s^-1  and the ratio of molar mass
         ! of SO4 to SO2 is 96/64 = 1.5
         !---------------------------------------------
         ! rk_com(i,k,12) : rate coef. in cm3/molec-s
         ! rewrite the rate of reaction as first order
         ! rate coefficient  :
         ! SO2_rate = rk_com(i,k,12) * [OH] (S-1)
         !
         ! Note: 1.5 factor accounts for different molecular
         !       weights of  so4 and so2 mw(so4)=96;
         !       mw(so2)=64;  mw(so4)/mw(so2)=3/2=1.5
         !---------------------------------------------
         ! so2_rate = rk_com(i,k,12) * oh1int * d_10
         so2_rate = rk_com(i,k,12) * oh1int
         so2_avail = max(chib(j,i,k,iso2)-mintr,d_zero)
         so2_snk(i,k) = so2_avail*(d_one-exp(-so2_rate*dt))/dt

         chiten(j,i,k,iso2) = chiten(j,i,k,iso2) - so2_snk(i,k) * cldno
         chiten(j,i,k,iso4) = chiten(j,i,k,iso4) + 1.5_rkx*so2_snk(i,k)*cldno

         !  gazeous conversion diagnostic
         if ( ichdiag > 0 ) then
          chemdiag(j,i,k,iso2) = chemdiag(j,i,k,iso2) &
               - so2_snk(i,k) * cldno * cfdout
          chemdiag(j,i,k,iso4) = chemdiag(j,i,k,iso4) &
               +  1.5_rkx*so2_snk(i,k)*cldno  * cfdout
         end if
       end do
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! AQUEOUS CONVERSION IN CLOUDS 
     ! works also when full chestry igaschem == 1
     ! Aqueous conversion from so2 to so4 : control by h2o2
     ! either from climatology / or from gas phase chem
     if (igaschem==0) then
     isulf = iso4
     do k = 1 , kz
       do i = ici1 , ici2
         chimol = 28.9_rkx/64.0_rkx*chib(j,i,k,iso2)/cpsb(j,i) ! kg/kg to mole
         if ( ioxclim == 1 ) then
           h2o2mol =  oxcl(j,i,k,iox_h2o2)
           concmin(i,k) = min(h2o2mol,chimol)*64.0_rkx/28.9_rkx*cpsb(j,i)
         else
           ! cb*kg/kg do tests, suppose h2o2 always enough
           concmin(i,k) = chimol*64._rkx/28.9_rkx*cpsb(j,i)     ! cb*kg/kg
         end if
       end do
     end do
     elseif (igaschem==1 .and. ih2o2 > 0) then 
     isulf = ih2so4
      do k = 1 , kz
       do i = ici1 , ici2
         chimol = 28.9_rkx/64.0_rkx*chib(j,i,k,iso2)/cpsb(j,i) ! kg/kg to mole
         h2o2mol= 28.9_rkx/34.0_rkx*chib(j,i,k,ih2o2)/cpsb(j,i)
         concmin(i,k) = min(h2o2mol,chimol)*64.0_rkx/28.9_rkx*cpsb(j,i)
       end do
      end do
     end if

     ! conversion in   Large scale clouds

     do k = 1 , kz
       do i = ici1 , ici2
         rxs1 = d_zero
         rxs11 = d_zero      ! fraction of conversion, not removed, as SO4 src
         wetrem(iso2) = d_zero
         ! scavenging for SO2, below lsc
         wetrem(isulf) = d_zero
         if ( wl(i,k) > clmin ) then
           ! conversion from so2 to so4
           rxs1 = fracloud(i,k)*chtrsol(iso2)*concmin(i,k) * &
                  (exp(-wl(i,k)/360.0_rkx*dt)-d_one)
           rxs11 = rxs1*1.5_rkx
           ! SO4 src term and the ratio of molar
           ! mass of SO4 to SO2 is 96/64 = 1.5

! FAB TEST : REMOVE WET DEP AS IT IS CALCULATED IN WETDEPA
!!$           if ( cremrat(j,i,k) > d_zero ) then
!!$             wetrem(iso4) = (fracloud(i,k)*chtrsol(iso4)*chib(j,i,k,iso4) - &
!!$                      rxs11)*(exp(-cremrat(j,i,k)/fracloud(i,k)*dt)-d_one)
!!$           end if


         end if
           ! Below cloud scavenging only for SO2 only stratiform precip !
           ! and isthet isnot called
           ! rembc is in calculated in prec, [mm/hr] and converted to
           ! below cloud scavenging rate for SO2 rate, s^-1)
           !     - Levin & Schwatz
           ! s^-1, it is already a grid scale removal rate!
           krembc = 6.5_rkx*1.0e-5_rkx*crembc(j,i,k)**0.68_rkx

         if ( crembc(j,i,k) > d_zero .and. igaschem == 0) then
             wetrem(iso2) =  chtrsol(iso2)*concmin(i,k) * &
                             (exp(-krembc*dt)-d_one)
         end if


         ! Tendancies large scale cloud
         chiten(j,i,k,iso2) = chiten(j,i,k,iso2) + rxs1/dt + &
                              wetrem(iso2)/dt ! zero here
         chiten(j,i,k,isulf) = chiten(j,i,k,isulf) - rxs11/dt + &
                              wetrem(isulf)/dt

         ! and wetdep diagnostics
         ! just for iso2 washout (washout) here.
         ! only the contribution of large scale cloud is accounted for
         washout(j,i,k,iso2) = washout(j,i,k,iso2) - wetrem(iso2)/dt *cfdout
         ! rainout(j,i,k,iso4) = rainout(j,i,k,iso4) - wetrem(iso4)/d_two

         ! chemical aqueous conversion diagnostic

         if ( ichdiag > 0 ) then
           chemdiag(j,i,k,iso2) = chemdiag(j,i,k,iso2) + rxs1/dt * cfdout
           chemdiag(j,i,k,isulf) = chemdiag(j,i,k,isulf) - rxs11/dt * cfdout
         end if
       end do
     end do
     ! cumulus clouds
     ! wet removal by cumulus clouds (over the fraction of grid box
     ! fracum) assume the cloud water content = 2 g/m3  (ref Kasibhatla )
     do i = ici1 , ici2
       if ( kcumtop(j,i) > 0 ) then
         do k = kcumtop(j,i) , kz
           rxs2 = d_zero
           rxs21 = d_zero    ! fraction of conversion, not removed, as SO4 src
           wetrem_cvc(iso2) = d_zero ! scavenging for SO2, below lsc
           wetrem_cvc(isulf) = d_zero

           ! conversion from so2 to so4
           rxs2 = fracum(i,k)*chtrsol(iso2)*concmin(i,k) * &
                  (exp(-d_two/360.0_rkx*dt)-d_one)
           rxs21 = rxs2*1.5_rkx

           ! removal (including theremoval on the rxs21 term)
           ! contratily to LS clouds, remcum is already an in cloud removal rate
!!$ FAB TEST DON'T COSIDER REMOVAL HERE
!!$          wetrem_cvc(iso4) = (fracum(i,k)*chtrsol(iso4)*chib(j,i,k,iso4) - &
!!$                              rxs21)*(exp(-remcum*dt)-d_one)

           ! tendancies due to convective cloud processes
           chiten(j,i,k,iso2) = chiten(j,i,k,iso2) + rxs2/dt

           chiten(j,i,k,isulf) = chiten(j,i,k,isulf) + &
                                wetrem_cvc(isulf)/dt - rxs21/dt

           ! diagnostic of wet deposition: NOT relevant here
           ! only SO2 below large scale cloud washout is considered above
           ! washout(j,i,k,1) = washout(j,i,k,1) - wetrem_cvc(iso2)/2.
           ! washout(j,i,k,iso2) = washout(j,i,k,iso4) - wetrem_cvc(iso4)/d_two

           ! chemical aquesous conversion diagnostic
           ! ( add the contribution of gas + wet lsc + wet cum conversions)
           if ( ichdiag > 0 ) then
             chemdiag(j,i,k,iso2) = chemdiag(j,i,k,iso2) + rxs2/dt * cfdout
             chemdiag(j,i,k,isulf) = chemdiag(j,i,k,isulf) - rxs21/dt  * cfdout
           end if
         end do
       end if
     end do

     ! diagnostic for SO2 durface fluxes

     do i = ici1 , ici2
       wdrout(j,i,iso2) = d_zero
       wdwout(j,i,iso2) = d_zero
       do k = 1 , kz
         ! sum on the vertical to get total surface flux diag fo rain out
         ! and washout (already weighted for time average cfdout !),
         ! also change sign convention normalise by psb to get the right
         ! flux unit
         wdwout(j,i,iso2) = wdwout(j,i,iso2) - &
            washout(j,i,k,iso2)*cdzq(j,i,k) *crhob3d(j,i,k)/cpsb(j,i)
       end do
     end do

!!$!
!!$  if ( (chtrname(itr) == 'DMS' ) .and. iso2 > 0 ) then
!!$
!!$    !---------------------------------------------
!!$    ! rate of reaction for DMS oxidation at daytime using
!!$    ! OH and at nighttime using NO3
!!$    !
!!$    !  - DMS + OH  ---> SO2 + ......
!!$    !  - DMS + OH  ---> 0.6 SO2 + 0.4 DMSO + ....
!!$    !  - DMSO + OH ---> 0.6 SO2 + 0.4 MSA + ......
!!$    !  - DMS + NO3 ---> SO2 +  .....
!!$    !---------------------------------------------
!!$
!!$    do  k = 1 , kx
!!$      do  i = ici1 , ici2
!!$
!!$        no3int = no3(i,k,j)
!!$        oh1int = oh(i,k,j)
!!$
!!$        if ( coszrs(i) < 0.001_rkx ) then
!!$          oh1int = oh1int*0.01_rkx
!!$        end if
!!$
!!$        ratdms_no3 = d_zero
!!$        ratdms_oh  = d_zero
!!$        rattot_dms = d_zero
!!$
!!$        ratdms_no3 = rk_com(i,k,17)  * no3int
!!$        ratdms_oh  = rk_com(i,k,10)  * oh1int
!!$        rattot_dms = ratdms_no3 + ratdms_oh
!!$        ratmsa     = 0.6_rkx * rk_com(i,k,9) * oh1int
!!$
!!$        dmsoh_snk(i,k) = 0.4_rkx*chib(i,k,j,idms) * &
!!$                        (d_one-exp(-ratdms_oh*dt))/dt
!!$        dmsno3_snk(i,k) = chib(i,k,j,idms) * &
!!$                         (d_one-exp(-ratdms_no3*dt))/dt
!!$
!!$        dms_snk(i,k) = chib(i,k,j,idms)*(d_one-exp(-rattot_dms*dt))/dt
!!$        chiten(i,k,j,idms) = chiten(i,k,j,idms) - dms_snk(i,k)
!!$        dms_gas(i,k) = 1.032258_rkx * (dmsoh_snk(i,k) + dmsno3_snk(i,k))
!!$        chiten(i,k,j,iso2) = chiten(i,k,j,iso2) + dms_gas(i,k)
!!$
!!$      end do
!!$    end do
!!$
!!$  end if ! dms , so2
!
   end subroutine chemsox
!
!  *****************************************************************
!  *  calculates reaction rates for the selected mechanism      ****
!  *  Ashraf S. Zakey, 2007                                     ****
!  *  International Center for Theortical Physics (ICTP)        ****
!  *  Trieste -Italy                                            ****
!  *****************************************************************
!
   subroutine chemrate(caircell,temp,rk_com)
     implicit none
     real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: caircell , temp
     real(rkx) , dimension(ici1:ici2,kz,100) , intent(out) ::  rk_com
     real(rkx) :: rk0 , rnn , rki , rmm , te , cair_mlc , alpha
     integer(ik4) :: i , k

     rk_com(:,:,:) = d_zero
     do k = 1 , kz
       do i = ici1 , ici2

          alpha = 0.4_rkx
          te = temp(i,k)

          rk_com(i,k,1) = arr(1.0_rkx,    -234._rkx,te)
          rk_com(i,k,2) = arr(8.46e-10_rkx, 7230._rkx,te)
          rk_com(i,k,3) = arr(2.68e-10_rkx, 7810._rkx,te)
          rk_com(i,k,4) = arr(88.1_rkx,   7460._rkx,te)

          rk_com(i,k,5) = temp(i,k)*rk_com(i,k,1) + rk_com(i,k,2)+ rk_com(i,k,3)
          rk_com(i,k,6) = 1.04e11_rkx*temp(i,k) + rk_com(i,k,4)
          rk_com(i,k,7) = rk_com(i,k,5)/rk_com(i,k,6)
          rk_com(i,k,8) = arr(1.2e-11_rkx,  -260._rkx,te)
          rk_com(i,k,9) = max(0.0_rkx, rk_com(i,k,7) - rk_com(i,k,8))
          rk_com(i,k,10) = rk_com(i,k,8) + rk_com(i,k,9)
          rk_com(i,k,11) = rk_com(i,k,8) + alpha * rk_com(i,k,9)

          cair_mlc = caircell(i,k)
          rk0 = 3.0e-31_rkx
          rnn = 3.3_rkx
          rki = 1.5e-12_rkx
          rmm = 0.0_rkx
          rk_com(i,k,12) = troe(cair_mlc,te,rk0,rnn,rki,rmm)

          rk_com(i,k,13) = abr(1.7e-12_rkx,   600._rkx,te)
          rk_com(i,k,14) = abr(4.9e-32_rkx,  1000._rkx,te)
          rk_com(i,k,15) = abr(1.4e-21_rkx,  2200._rkx,te)
          rk_com(i,k,16) = abr(1.7e-12_rkx,  -200._rkx,te)
          rk_com(i,k,17) = abr(1.0e-12_rkx,   500._rkx,te)

       end do
     end do

     contains

     real(rkx) function troe(cair_mlc,te,rk0,rnn,rki,rmm)
       implicit none
       real(rkx) , intent(in) ::  cair_mlc, te, rnn, rmm
       real(rkx) , intent(inout) :: rk0 , rki
       real(rkx) :: expo
       rk0 = rk0*cair_mlc*(te/300.0_rkx)**(-rnn)
       rki = rki*(te/300.0_rkx)**(-rmm)
       expo= d_one/(d_one + (log10(rk0/rki))**2)
       troe  = (rk0*rki/(rk0+rki))*0.6_rkx**expo
     end function troe

     real(rkx) function arr(aa,bb,te)
       implicit none
       real(rkx), intent(in) :: aa , bb , te
       arr = aa*exp(bb/te)
     end function arr

     real(rkx) function abr(aa,bb,te)
       implicit none
       real(rkx), intent(in) ::  aa , bb , te
       abr = aa*exp(bb*(d_one/te - 0.0033557_rkx))
     end function abr

   end subroutine chemrate

end module mod_che_sox
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
