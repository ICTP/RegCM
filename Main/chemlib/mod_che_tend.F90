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
 
  module mod_che_tend
!
! Tendency and budget for tracer transport and chemicals
!
  use mod_constants
  use mod_realkinds
  use mod_dynparam
  use mod_che_common
  use mod_che_indices
  use mod_che_param
  use mod_che_sox
  use mod_che_drydep
  use mod_che_wetdep
  use mod_che_emission
  use mod_che_dust
  use mod_che_seasalt
  use mod_che_carbonaer
  use mod_che_mppio
  use mod_che_chemistry
  private

  public :: tractend2 , tracbud

  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     This subroutine computes the tendencies for tracer transport and
!     chemistry
!
!     j:             index of j slice in current computation
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine tractend2(ktau,lyear,lmonth,lday,secofday)
      implicit none
!
      integer , intent(in) :: lmonth , lday , lyear
      real(dp), intent(in) :: secofday
      integer(8) , intent(in) :: ktau
!
      real(dp) :: facb , facs , fact , facv , pres10 , qsat10 , &
                 satvp , shu10 , u10 , v10

      real(dp) , dimension(jci1:jci2,iy,kz) :: rho , &
                                    ttb, wl, fracloud, fracum , prec

      integer :: i , j , ibin , k , kk
      integer , dimension(1:jxp,iy) :: ivegcov

      real(dp) , dimension(1:jxp,iy,kz,ntr) :: pdepv
      real(dp) , dimension(1:jxp,iy,ntr) :: ddepa

      real(dp) , dimension(1:jxp,iy) :: psurf , rh10 , soilw , srad ,  &
          temp10 , tsurf , vegfrac , wid10 , zeff , ustar

      real(dp) , dimension(1:jxp,iy,nbin) :: dust_flx
      real(dp), dimension(1:jxp,iy,sbin) :: seasalt_flx
!
      integer(8) :: kchsolv
!
!**************************************************************************
!     A : PRELIMINARY CALCULATIONS
!*************************************************************************
! 
        rho = d_zero
        wl = d_zero
        ttb = d_zero
        prec = d_zero
        ustar =d_zero
        fracloud = d_zero
        fracum = d_zero
        psurf = d_zero  
        dust_flx = d_zero
        seasalt_flx = d_zero
        ivegcov=0
      do j = jci1 , jci2

!       the unit: rho - kg/m3, wl - g/m3
        do k = 1 , kz
          do i = ici1 , ici2
!           rho(j,i,k) = (sfs%psb(i,j)*a(k)+r8pt)* &
!      what the hell   1000./287./atm2%t(i,k,j)*sfs%psb(i,j)
            rho(j,i,k) = crhob3d(j,i,k)
            wl(j,i,k) = cqcb3d(j,i,k)*rho(j,i,k)
            ttb(j,i,k) = ctb3d(j,i,k)
!           precipiation rate is a rquired variable for deposition routines.
!           It is directly taken as rembc (saved in precip routine) in mm/hr !!
            prec(j,i,k) = crembc(j,i,k) / 3600.D0 !passed in mm/s  
          end do
        end do
!       cloud fractionnal cover for wet deposition
!       large scale : fracloud, calculated from fcc coming from pcp.f
!       cumulus scale : fracum, calculated from the total cloud fraction
!       (as defined for the radiation scheme in cldfrac.f routine)
   
        do i = ici1 , ici2
           if ( kcumtop(j,i) > 0 ) then
            do kk = kcumtop(j,i) , kz
              fracum(j,i,kk) = ccldfra(j,i,kk) - cfcc(j,i,kk)
            end do
          end if
        end do
!
!       variables used for natural fluxes and deposition velocities 
! 
        
        do i = ici1 , ici2

! care ocean-lake in veg2d is now back type 14-15 !!

          if ( cveg2d(j,i)== 14 .or. cveg2d(j,i)==15 ) then 
          ivegcov(j,i) = 0
          else
          ivegcov(j,i) = cveg2d(j,i)
          end if 



          psurf(j,i) = cpsb(j,i) * 1.0D3 + ptop
 
!         method based on bats diagnostic in routine interf.
 
          if ( (ivegcov(j,i) /= 0) ) then
            facv = dlog(cza(j,i,kz)/d_10) / &
                   dlog(cza(j,i,kz)/crough(ivegcov(j,i)))
            facb = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zlnd)
            facs = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zsno)
 
!           fact = csfracv2d(j,i)*facv 
            fact = cvegfrac(j,i) * facv + (d_one-cvegfrac(j,i)) * facb
!           FAB REVOIR CETTE partie et definir interface pour sfracs,sfracv
!                   + sfracb2d(j,i)*facb + sfracs2d(j,i)*facs
! 
!           grid level effective roughness lenght (linear averaging for now)
!           zeff(i) = rough(ivegcov(i))*sfracv2d(j,i) + &
!                     zlnd * sfracb2d(j,i) + zsno * sfracs2d(j,i)
            zeff(j,i) = crough(ivegcov(j,i))*cvegfrac(j,i) + &
                      zlnd*(d_one-cvegfrac(j,i))
          else
!           water surface
            fact = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zoce)
            zeff(j,i) = zoce
          end if
!         10 m wind
          u10 = (cubx3d(j,i,kz))*(1-fact)
          v10 = (cvbx3d(j,i,kz))*(1-fact)
          wid10(j,i) = sqrt(u10**2+v10**2)
!         10 m air temperature
          temp10(j,i) = ttb(j,i,kz) - csdeltk2d(j,i)*fact
!         specific  humidity at 10m
          shu10 = cqvb3d(j,i,kz)/ &
                  (d_one+cqvb3d(j,i,kz))-csdelqk2d(j,i)*fact
!         back to mixing ratio
          shu10 = shu10/(1-shu10)
!         saturation mixing ratio at 10m
          if ( temp10(j,i) > tzero ) then
            satvp = svp1*1.0D3*dexp(svp2*(temp10(j,i)-tzero)/(temp10(j,i)-svp3))
          else
            satvp = svp4*1.0D3*dexp(svp5-svp6/temp10(j,i))
          end if
          pres10 = psurf(j,i) - 98.0D0
          qsat10 = ep2*satvp/(pres10-satvp)
!         relative humidity at 10m
          rh10(j,i) = d_zero
          if ( qsat10 > d_zero ) rh10(j,i) = shu10/qsat10
!
!         friction velocity ( from uvdrag so updtaed at  bats or clm frequency : little inconsistenccy here)
!
          ustar(j,i) = sqrt ( cuvdrag(j,i)             *           &
                       sqrt ( (cubx3d(j,i,kz))**2   +  (cvbx3d(j,i,kz))**2 ) /     &
                             crhob3d(j,i,kz) )
 
!           soil wetness
 
          soilw(j,i) = cssw2da(j,i)/cdepuv(idnint(clndcat(j,i)))/(2650.0D0 * &
                (d_one-cxmopor(ciexsol(idnint(clndcat(j,i))))))
!         fraction of vegetation
          vegfrac(j,i) = cvegfrac(j,i)
!         surface temperature
!         over land recalculated from the BATS  deltk air/ surface
!         temperature account for a composite temperature between
!         bare ground and vegetation
          if ( ivegcov(j,i) /= 0 ) then
            tsurf(j,i) = ttb(j,i,kz) - csdeltk2d(j,i)
          else
!           ocean temperature in this case
            tsurf(j,i) = ctg(j,i)
          end if
 
!        aborbed solar radiation (for stb criteria used to calculate
!        aerodynamic resistance)
 
         srad(j,i) = csol2d(j,i)
 
        end do

    end do ! jloop 

!
!       END of preliminary calculations)
!
!*****************************************************************
! B :CALCULATION OF TRACER TENDENCY (except full gas phase chemistry solver)
!*****************************************************************
!
!       SOX CHEMSITRY ( from offline oxidant) 
!
        if ( igaschem == 0 ) then
          if (iso2 > 0 .and. iso4 >0.) then
          do j=jci1,jci2
            call chemsox(j,wl(j,ici1:ici2,:),fracloud(j,ici1:ici2,:), &
                fracum(j,ici1:ici2,:),rho(j,ici1:ici2,:),ttb(j,ici1:ici2,:))
          end do
          end if
        end if


!
!       aging of carboneaceous aerosols
!
        if ( (ibchb > 0 .and. ibchl > 0 ) .or. &
             (iochb > 0 .and. iochl > 0) ) then
          do j=jci1,jci2
          call aging_carb(j)
          end do
        end if

!  before emission and deposition routine set the surfecae netflux used by BL schems to zero
         cchifxuw = d_zero
!

        ! NATURAL EMISSIONS FLUX and tendencies  (dust -sea salt)       
        if ( idust(1) > 0 ) then
        
        do j= jci1,jci2
 
        call sfflux(iy,ici1,ici2,j,ivegcov(j,:),vegfrac(j,:),ustar(j,:), &
                      zeff(j,:),soilw(j,:),wid10(j,:),rho(j,:,kz),dustbsiz,dust_flx(j,:,:))     
        end do 
        end if

        if (  isslt(1) > 0 )then
        do j=jci1,jci2
         call sea_salt(j,wid10(j,ici1:ici2),ivegcov(j,ici1:ici2), &
                       seasalt_flx(j,ici1:ici2,:))
        end do
        end if
!
!       update emission tendencies from inventories


        do j= jci1,jci2
        call emis_tend(ktau,j,lmonth)
        end do
!
!       aerosol settling and drydep 
!       include calculation of dry dep/settling velocities and 
!       updating tendencies
!
        pdepv = d_zero
        ddepa = d_zero
        if ( idust(1) > 0 ) then
        do j=jci1,jci2
          call drydep_aero(j,nbin,idust,rhodust,ivegcov(j,ici1:ici2), &
              ttb(j,ici1:ici2,:),rho(j,ici1:ici2,:),hlev,psurf(j,ici1:ici2), &
              temp10(j,ici1:ici2),tsurf(j,ici1:ici2),srad(j,ici1:ici2),      &
              rh10(j,ici1:ici2),wid10(j,ici1:ici2),zeff(j,ici1:ici2),dustbed, &
              pdepv(j,ici1:ici2,:,:),ddepa(j,ici1:ici2,:))
        end do
        end if

       if ( isslt(1) >0 ) then
       do j=jci1,jci2
         call drydep_aero(j,sbin,isslt,rhosslt,ivegcov(j,ici1:ici2), &
           ttb(j,ici1:ici2,:),rho(j,ici1:ici2,:),hlev,psurf(j,ici1:ici2), &
           temp10(j,ici1:ici2),tsurf(j,ici1:ici2),srad(j,ici1:ici2), &
           rh10(j,ici1:ici2),wid10(j,ici1:ici2),zeff(j,ici1:ici2),ssltbed,  &
           pdepv(j,ici1:ici2,:,:),ddepa(j,ici1:ici2,:))
       end do
       end if 

        if ( icarb(1) > 0 ) then
        ibin = count( icarb > 0 ) 
          do j=jci1,jci2
          call drydep_aero(j,ibin,icarb(1:ibin),rhooc,ivegcov(j,ici1:ici2), &
            ttb(j,ici1:ici2,:),rho(j,ici1:ici2,:),hlev,psurf(j,ici1:ici2), &
            temp10(j,ici1:ici2),tsurf(j,ici1:ici2),srad(j,ici1:ici2), &
            rh10(j,ici1:ici2),wid10(j,ici1:ici2),zeff(j,ici1:ici2),         &
            carbed(1:ibin),pdepv(j,ici1:ici2,:,:),ddepa(j,ici1:ici2,:))
         end do
        end if 
!!$
!       GAS phase dry deposition velocity + tendencies
!       option compatible with BATS and CLM
!!$
        if ( igaschem == 1 ) then
          do j=jci1,jci2
!          call drydep_gas(j,calday, ivegcov(j,ici1:ici2),rh10(j,ici1:ici2), &
!               srad(j,ici1:ici2),tsurf(j,ici1:ici2),prec(j,ici1:ici2,kz),   &
!               temp10(j,ici1:ici2),wid10(j,ici1:ici2),zeff(j,ici1:ici2),    &
!               drydepvg(j,ici1:ici2,:))
          end do
        end if
!!$
!       WET deposition (rainout and washout) for aerosol
!!$
        if ( idust(1) > 0 ) then
         do j=jci1,jci2
          call wetdepa(j,nbin,idust,dustbed,rhodust,ttb(j,ici1:ici2,:), &
                       wl(j,ici1:ici2,:),fracloud(j,ici1:ici2,:),       &
                       fracum(j,ici1:ici2,:),psurf(j,ici1:ici2),        &
                       hlev,rho(j,ici1:ici2,:),prec(j,ici1:ici2,:),     &
                       pdepv(j,ici1:ici2,:,:))  
         end do
        end if

       if ( isslt(1) > 0 )  then   
         do j=jci1,jci2
         call wetdepa(j,sbin,isslt,ssltbed,rhosslt,ttb(j,ici1:ici2,:), &
                      wl(j,ici1:ici2,:),fracloud(j,ici1:ici2,:),       &
                      fracum(j,ici1:ici2,:),psurf(j,ici1:ici2),hlev,   &
                      rho(j,ici1:ici2,:),prec(j,ici1:ici2,:),          &
                      pdepv(j,ici1:ici2,:,:) )  
         end do
       end if
       if ( icarb(1) > 0 )  then   
         ibin = count( icarb > 0 ) 
          do j=jci1,jci2
         call wetdepa(j,ibin,icarb(1:ibin),carbed(1:ibin),rhobchl,  &
                     ttb(j,ici1:ici2,:),wl(j,ici1:ici2,:),          &
                     fracloud(j,ici1:ici2,:),fracum(j,ici1:ici2,:), &
                     psurf(j,ici1:ici2),hlev,rho(j,ici1:ici2,:),    &
                     prec(j,ici1:ici2,:),pdepv(j,ici1:ici2,:,:) ) 
         end do
         end if
!
!!$
!       Wet Deposition for gasphase species 
!!$
        if ( igaschem == 1 ) then
        checum = d_zero !  fix the interface for this variable ! no effect of cumulus scavenging now
        chevap = d_zero
         do j=jci1,jci2        
!          call sethet(j,cza(j,:,:),cht(j,:),ttb(j,:,:),checum(j,:,:),cremrat(j,:,:), &
!                      chevap(j,:,:),dtche,rho(j,:,:),chib(j,:,:,:),iym3,cpsb(j,ici1:ici2))
         end do
        end if
    
!  Gas phase solver 
!  note : solver is called every dtchsolv (900s)- chemten  chemistry raecation tendendy is calculated
!  but chemical tracer tendency is still updated every dtche ( =dt) time step ( insure smoothness)  

      chemten(:,:,:,:) = d_zero
     
      if ( igaschem == 1 .and. ichsolver > 0) then   
      kchsolv = idnint(dtchsolv / dtche)
      kchsolv = 6 ! for the moment
   
      if (mod(ktau+1,kchsolv) == 0 ) then   
       do j = jci1,jci2
        call gas_phase(j,secofday,lyear,lmonth,lday)
       end do
       end if
! add tendency due to chemistry reaction (every dtche)      
       chiten(jci1:jci2,:,:,:) =  chiten(jci1:jci2,:,:,:) +   chemten(jci1:jci2,:,:,:)

      end if

    end subroutine tractend2
!
    subroutine tracbud
      implicit none
      integer :: i , itr , j , k
      dtrace(:,:,:) = d_zero
      wdlsc(:,:,:) = d_zero
      wdcvc(:,:,:) = d_zero
      wxsg(:,:,:) = d_zero
      wxaq(:,:,:) = d_zero
      ddsfc(:,:,:) = d_zero
      ! 
      ! tracers (unit = kg):
      ! 
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              dtrace(j,i,itr) = dtrace(j,i,itr) + chia(j,i,k,itr)*cdsigma(k)
              wdlsc(j,i,itr) = wdlsc(j,i,itr) + remlsc(j,i,k,itr)*cdsigma(k)
              wdcvc(j,i,itr) = wdcvc(j,i,itr) + remcvc(j,i,k,itr)*cdsigma(k)
              wxsg(j,i,itr) = wxsg(j,i,itr) + rxsg(j,i,k,itr)*cdsigma(k)
              ! sum ls and conv contribution
              wxaq(j,i,itr) = wxaq(j,i,itr)                             &
                            & + (rxsaq1(j,i,k,itr)+rxsaq2(j,i,k,itr))   &
                            & *cdsigma(k)
            end do
          end do
        end do
      end do
      do itr = 1 , ntr
        do i = ici1 , ici2
          do j = jci1 , jci2
            ddsfc(j,i,itr) = ddsfc(j,i,itr) + remdrd(j,i,itr)*cdsigma(kz)
            ! Source cumulated diag(care the unit are alredy .m-2)
            cemtrac(j,i,itr) = cemtr(j,i,itr)
          end do
        end do
      end do

      do itr = 1 , ntr
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! unit: mg/m2
            dtrace(j,i,itr) = 1.D6*dtrace(j,i,itr)*d_1000*regrav
            wdlsc(j,i,itr) = 1.D6*wdlsc(j,i,itr)*d_1000*regrav
            wdcvc(j,i,itr) = 1.D6*wdcvc(j,i,itr)*d_1000*regrav
            ddsfc(j,i,itr) = 1.D6*ddsfc(j,i,itr)*d_1000*regrav
            wxsg(j,i,itr) = 1.D6*wxsg(j,i,itr)*d_1000*regrav
            wxaq(j,i,itr) = 1.D6*wxaq(j,i,itr)*d_1000*regrav
            ! cemtrac isbuilt from chsurfem so just need the 1e6*dt/2
            ! factor to to pass im mg/m2
            cemtrac(j,i,itr) = 1.D6*cemtrac(j,i,itr)
          end do
        end do
      end do

    end subroutine tracbud
!
end module mod_che_tend
