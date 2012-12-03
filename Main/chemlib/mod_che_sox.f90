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

  private

  public :: chemsox

  contains

    subroutine chemsox(j,wl,fracloud,fracum,rho,ttb)
     implicit none
!
     integer(ik4) , intent(in) :: j
     real(rk8) , dimension(ici1:ici2,kz) , intent(in) :: ttb , wl , rho
     real(rk8) , dimension(ici1:ici2,kz) , intent(in) :: fracloud , fracum
     real(rk8) :: rxs1 , rxs11 , rxs2 , rxs21 , chimol , cldno , &
                 clmin , remcum , oh1int , so2_rate , krembc
     real(rk8) , dimension(ntr) ::  wetrem , wetrem_cvc
     real(rk8) , dimension(ici1:ici2,kz) :: caircell , so2_snk , concmin
     real(rk8) :: rk_com(ici1:ici2,kz,100)
     real(rk8) :: h2o2mol

     integer(ik4) :: i , k

!FAB 
     caircell(:,:) = 1.D-6 * rho(:,:)/amdk * navgdr

     call chemrate(caircell,ttb,rk_com)

     ! clmin = non-precipitating cloud
     ! conversion threshold, clmin=0.01g/m3
     clmin = 0.01D0

     ! remcum= removal rate for cumulus
     ! cloud scavenging (s-1) (in-cloud and not grid level)
     remcum = 1.0D-3

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
           if ( czen(j,i) < 0.001D0 ) then
             oh1int = oh1int * 0.01D0
           else 
             oh1int = oh1int * 1.99D0 
           end if
         else
           oh1int = 30.0D5
           if ( czen(j,i) < 0.001D0 ) then
             oh1int = oh1int * 0.01D0
           else 
             oh1int = oh1int * 1.99D0 
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

         so2_rate = rk_com(i,k,12) * oh1int * d_10
         so2_snk(i,k) = chib(j,i,k,iso2)*(d_one-dexp(-so2_rate*dtche))/dtche

         chiten(j,i,k,iso2) = chiten(j,i,k,iso2) - so2_snk(i,k) * cldno
         chiten(j,i,k,iso4) = chiten(j,i,k,iso4) + 1.5D0*so2_snk(i,k)*cldno 

         !  gazeous conversion diagnostic 

         rxsg(j,i,k,iso2) = rxsg(j,i,k,iso2) + &
                            so2_snk(i,k)*cldno*dtche/d_two
         rxsg(j,i,k,iso4) = rxsg(j,i,k,iso4) + &
                            1.5D0*so2_snk(i,k)*cldno*dtche/d_two

       end do
     end do

     ! AQUEOUS CONVERSION IN CLOUDS AND WET REMOVAL
     ! Aqueous conversion from so2 to so4 : control by h2o2
     do k = 1 , kz
       do i = ici1 , ici2
         chimol = 28.9D0/64.0D0*chib(j,i,k,iso2)/cpsb(j,i) ! kg/kg to mole

         if ( ichaer == 1 ) then 
           h2o2mol =  oxcl(i,k,j,iox_h2o2)
           concmin(i,k) = dmin1(h2o2mol,chimol)*64.0D0/28.9D0*cpsb(j,i)
         else
         ! cb*kg/kg do tests, suppose h2o2 always enough
           concmin(i,k) = chimol*64.D0/28.9D0*cpsb(j,i)     ! cb*kg/kg
         end if
       end do
     end do
     ! conversion in   Large scale clouds
 
     do k = 1 , kz
       do i = ici1 , ici2
         rxs1 = d_zero
         rxs11 = d_zero      ! fraction of conversion, not removed, as SO4 src
         wetrem(iso2) = d_zero
         ! scavenging for SO2, below lsc
         wetrem(iso4) = d_zero
         if ( wl(i,k) > clmin ) then
           ! conversion from so2 to so4
           rxs1 = fracloud(i,k)*chtrsol(iso2)*concmin(i,k) * &
                  (dexp(-wl(i,k)/360.0D0*dtche)-d_one)
           rxs11 = rxs1*1.5D0
           ! SO4 src term and the ratio of molar
           ! mass of SO4 to SO2 is 96/64 = 1.5
 
           ! if removal occurs, a fraction of SO4 src term is also
           ! removed and accounted for in the term  wetrem(iso4)
           ! FAB:
           ! Care , remrat and rembc as calculated in precip are
           ! already grid level average removal rates.
           ! here remart is divided by fracloud
           ! ( large scale cloud fraction) to get the incloud removal rate
           !
! FAB TEST : REMOVE WET DEP AS IT IS CALCULATED IN WETDEPA

!!$           if ( cremrat(j,i,k) > d_zero ) then
!!$             wetrem(iso4) = (fracloud(i,k)*chtrsol(iso4)*chib(j,i,k,iso4) - &
!!$                      rxs11)*(dexp(-cremrat(j,i,k)/fracloud(i,k)*dtche)-d_one)
!!$           end if
 
           ! Below cloud scavenging only for SO2 only stratiform precip !
           ! rembc is in calculated in prec, [mm/hr] and converted to
           ! below cloud scavenging rate for SO2 rate, s^-1)
           !     - Levin & Schwatz
           ! s^-1, it is already a grid scale removal rate!
           krembc = 6.5D0*1.0D-5*crembc(j,i,k)**0.68D0

           if ( crembc(j,i,k) > d_zero ) then
             wetrem(iso2) =  chtrsol(iso2)*concmin(i,k) * &
                             (dexp(-krembc*dtche)-d_one)
           end if
         end if
 
         ! Tendancies large scale cloud
         chiten(j,i,k,iso2) = chiten(j,i,k,iso2) + rxs1/dtche + &
                              wetrem(iso2)/dtche
         chiten(j,i,k,iso4) = chiten(j,i,k,iso4) - rxs11/dtche + &
                              wetrem(iso4)/dtche
 
         ! and wetdep diagnostics
         remlsc(j,i,k,iso2) = remlsc(j,i,k,iso2) - wetrem(iso2)/d_two
         remlsc(j,i,k,iso4) = remlsc(j,i,k,iso4) - wetrem(iso4)/d_two
 
         ! chemical aqueous conversion diagnostic
         rxsaq1(j,i,k,iso2) = rxsaq1(j,i,k,iso2) - rxs1/d_two
         rxsaq1(j,i,k,iso4) = rxsaq1(j,i,k,iso4) - rxs11/d_two
 
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
           wetrem_cvc(iso4) = d_zero
 
           ! conversion from so2 to so4
           rxs2 = fracum(i,k)*chtrsol(iso2)*concmin(i,k) * &
                  (dexp(-d_two/360.0D0*dtche)-d_one)
           rxs21 = rxs2*1.5D0
 
           ! removal (including theremoval on the rxs21 term)
           ! contratily to LS clouds, remcum is already an in cloud removal rate
!!$ FAB TEST DON'T COSIDER REMOVAL HERE
!!$          wetrem_cvc(iso4) = (fracum(i,k)*chtrsol(iso4)*chib(j,i,k,iso4) - &
!!$                              rxs21)*(dexp(-remcum*dtche)-d_one)

           ! tendancies due to convective cloud processes
           chiten(j,i,k,iso2) = chiten(j,i,k,iso2) + rxs2/dtche

           chiten(j,i,k,iso4) = chiten(j,i,k,iso4) + &
                                wetrem_cvc(iso4)/dtche - rxs21/dtche
 
           ! diagnostic of wet deposition
           ! remcvc(j,i,k,1) = remcvc(j,i,k,1) - wetrem_cvc(iso2)/2.
           remcvc(j,i,k,iso4) = remcvc(j,i,k,iso4) - wetrem_cvc(iso4)/d_two
           ! chemical aquesous conversion diagnostic
           rxsaq2(j,i,k,iso2) = rxsaq2(j,i,k,iso2) - rxs2/d_two
           rxsaq2(j,i,k,iso4) = rxsaq2(j,i,k,iso4) - rxs21/d_two
 
         end do
       end if
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
!!$        if ( coszrs(i) < 0.001D0 ) then
!!$          oh1int = oh1int*0.01D0
!!$        end if
!!$
!!$        ratdms_no3 = d_zero
!!$        ratdms_oh  = d_zero
!!$        rattot_dms = d_zero
!!$
!!$        ratdms_no3 = rk_com(i,k,17)  * no3int
!!$        ratdms_oh  = rk_com(i,k,10)  * oh1int
!!$        rattot_dms = ratdms_no3 + ratdms_oh
!!$        ratmsa     = 0.6D0 * rk_com(i,k,9) * oh1int
!!$
!!$        dmsoh_snk(i,k) = 0.4D0*chib(i,k,j,idms) * &
!!$                        (d_one-dexp(-ratdms_oh*dtche))/dtche
!!$        dmsno3_snk(i,k) = chib(i,k,j,idms) * &
!!$                         (d_one-dexp(-ratdms_no3*dtche))/dtche
!!$
!!$        dms_snk(i,k) = chib(i,k,j,idms)*(d_one-dexp(-rattot_dms*dtche))/dtche
!!$        chiten(i,k,j,idms) = chiten(i,k,j,idms) - dms_snk(i,k)
!!$        dms_gas(i,k) = 1.032258D0 * (dmsoh_snk(i,k) + dmsno3_snk(i,k))
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
     real(rk8) , dimension(ici1:ici2,kz) , intent(in) :: caircell , temp
     real(rk8) , dimension(ici1:ici2,kz,100) , intent(out) ::  rk_com
     real(rk8) :: rk0 , rnn , rki , rmm , te , cair_mlc , alpha
     integer(ik4) :: i , k

     rk_com(:,:,:) = d_zero
     do k = 1 , kz
       do i = ici1 , ici2

          alpha = 0.4D0  
          te = temp(i,k)

          rk_com(i,k,1) = arr(1.0D0,    -234.D0,te)
          rk_com(i,k,2) = arr(8.46D-10, 7230.D0,te)
          rk_com(i,k,3) = arr(2.68D-10, 7810.D0,te)
          rk_com(i,k,4) = arr(88.1D0,   7460.D0,te)

          rk_com(i,k,5) = temp(i,k)*rk_com(i,k,1) + rk_com(i,k,2)+ rk_com(i,k,3)
          rk_com(i,k,6) = 1.04D+11*temp(i,k) + rk_com(i,k,4)
          rk_com(i,k,7) = rk_com(i,k,5)/rk_com(i,k,6)
          rk_com(i,k,8) = arr(1.2D-11,  -260.D0,te)
          rk_com(i,k,9) = dmax1(0.0D0, rk_com(i,k,7) - rk_com(i,k,8))
          rk_com(i,k,10) = rk_com(i,k,8) + rk_com(i,k,9)
          rk_com(i,k,11) = rk_com(i,k,8) + alpha * rk_com(i,k,9)

          cair_mlc = caircell(i,k)
          rk0 = 3.0D-31
          rnn = 3.3D0
          rki = 1.5D-12
          rmm = 0.0D0
          rk_com(i,k,12) = troe(cair_mlc,te,rk0,rnn,rki,rmm)

          rk_com(i,k,13) = abr(1.7D-12,   600.D0,te)
          rk_com(i,k,14) = abr(4.9D-32,  1000.D0,te)
          rk_com(i,k,15) = abr(1.4D-21,  2200.D0,te)
          rk_com(i,k,16) = abr(1.7D-12,  -200.D0,te)
          rk_com(i,k,17) = abr(1.0D-12,   500.D0,te)

       end do
     end do

     contains

     real(rk8) function troe(cair_mlc,te,rk0,rnn,rki,rmm)
       implicit none
       real(rk8) , intent(in) ::  cair_mlc, te, rnn, rmm
       real(rk8) , intent(inout) :: rk0 , rki
       real(rk8) :: expo
       rk0 = rk0*cair_mlc*(te/300.0D0)**(-rnn)
       rki = rki*(te/300.0D0)**(-rmm)
       expo= d_one/(d_one + (dlog10(rk0/rki))**2)
       troe  = (rk0*rki/(rk0+rki))*0.6D0**expo
     end function troe
  
     real(rk8) function arr(aa,bb,te)
       implicit none
       real(rk8), intent(in) :: aa , bb , te   
       arr = aa*dexp(bb/te)
     end function arr 

     real(rk8) function abr(aa,bb,te)
       implicit none
       real(rk8), intent(in) ::  aa , bb , te
       abr = aa*dexp(bb*(d_one/te - 0.0033557D0))
     end function abr 

   end subroutine chemrate

end module mod_che_sox
