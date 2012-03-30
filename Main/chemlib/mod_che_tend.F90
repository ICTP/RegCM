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
    subroutine tractend2(jstart,jend,istart,iend,ktau,lyear,lmonth,lday,calday,secofday)
      implicit none
!
      integer , intent(in) :: jstart , jend , istart , iend , lmonth,lday,lyear
      real(8), intent(in) :: calday,secofday
      integer(8) , intent(in) :: ktau


!
      real(8) :: agct , ak00t , ak0tm , akval , clmin , facb , facs , &
                 fact , facv , pres10 , qsat10 , remcum , satvp ,     &
                 shu10 , u10 , v10 , chias , chibs

      real(8) , dimension(ntr) :: agingtend
      real(8) , dimension(jstart:jend,iy,kz) :: wk, rho , settend , &
                                    ttb, wl, fracloud, fracum , prec

      integer :: i , j , ibin , itr , k , kk , kb , kdwd
      integer , dimension(1:jxp,iy) :: ivegcov

      real(8) , dimension(1:jxp,iy,kz,ntr) :: pdepv
      real(8) , dimension(1:jxp,iy,ntr) :: ddepa


      real(8) , dimension(1:jxp,iy) :: psurf , rh10 , soilw , srad ,  &
          temp10 , tsurf , vegfrac , wid10 , zeff , ustar

      real(dp) , dimension(1:jxp,iy,nbin) :: dust_flx
      real(8), dimension(1:jxp,iy,sbin) :: seasalt_flx
      real(8), dimension(1:jxp,iy,ntr) :: drydepvg
      real(8) , dimension(ntr) :: wetrem , wetrem_cvc
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
      do j = jstart , jend

!       the unit: rho - kg/m3, wl - g/m3
        do k = 1 , kz
          do i = istart , iend
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
   
        do i = istart , iend
           if ( kcumtop(j,i) > 0 ) then
            do kk = kcumtop(j,i) , kz
              fracum(j,i,kk) = ccldfra(j,i,kk) - cfcc(j,i,kk)
            end do
          end if
        end do
!
!       variables used for natural fluxes and deposition velocities 
! 
        
        do i = istart , iend

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
          do j=jstart,jend
            call chemsox(j,wl(j,:,:),fracloud(j,:,:),fracum(j,:,:),rho(j,:,:),ttb(j,:,:))
          end do
          end if
        end if


!
!       aging of carboneaceous aerosols
!
        if ( (ibchb > 0 .and. ibchl > 0 ) .or. &
             (iochb > 0 .and. iochl > 0) ) then
          do j=jstart,jend
          call aging_carb(j)
          end do
        end if

!  before emission and deposition routine set the surfecae netflux used by BL schems to zero
         cchifxuw = d_zero
!

        ! NATURAL EMISSIONS FLUX and tendencies  (dust -sea salt)       
        if ( idust(1) > 0 ) then
        
        do j= jstart,jend
 
        call sfflux(iy,2,iym2,j,ivegcov(j,:),vegfrac(j,:),ustar(j,:), &
                      zeff(j,:),soilw(j,:),wid10(j,:),rho(j,:,kz),dustbsiz,dust_flx(j,:,:))     
        end do 
        end if

        if (  isslt(1) > 0 )then
        do j=jstart,jend
         call sea_salt(j,wid10(j,:),ivegcov(j,:),seasalt_flx(j,:,:))
        end do
        end if
!
!       update emission tendencies from inventories


        do j= jstart,jend
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
        do j=jstart,jend
          call drydep_aero(j,nbin,idust,rhodust,ivegcov(j,:),ttb(j,:,:),rho(j,:,:),hlev,psurf(j,:), &
                            temp10(j,:),tsurf(j,:),srad(j,:),rh10(j,:),wid10(j,:),zeff(j,:),dustbed,      &
                            pdepv(j,:,:,:),ddepa(j,:,:))
        end do
        end if

       if ( isslt(1) >0 ) then
       do j=jstart,jend
         call drydep_aero(j,sbin,isslt,rhosslt,ivegcov(j,:),ttb(j,:,:),rho(j,:,:),hlev,psurf(j,:), &
                          temp10(j,:),tsurf(j,:),srad(j,:),rh10(j,:),wid10(j,:),zeff(j,:),ssltbed,      &
                          pdepv(j,:,:,:),ddepa(j,:,:))
       end do
       end if 

        if ( icarb(1) > 0 ) then
        ibin = count( icarb > 0 ) 
          do j=jstart,jend
          call drydep_aero(j,ibin,icarb(1:ibin),rhooc,ivegcov(j,:),ttb(j,:,:),rho(j,:,:),hlev, &
                           psurf(j,:),temp10(j,:),tsurf(j,:),srad(j,:),rh10(j,:),wid10(j,:),zeff(j,:),         &
                           carbed(1:ibin),pdepv(j,:,:,:),ddepa(j,:,:))
         end do
        end if 
!!$
!       GAS phase dry deposition velocity + tendencies
!       option compatible with BATS and CLM
!!$
        if ( igaschem == 1 ) then
          do j=jstart,jend
!          call drydep_gas(j,calday, ivegcov(j,:),rh10(j,:),srad(j,:),tsurf(j,:),prec(j,:,kz),temp10(j,:),  &
!                          wid10(j,:),zeff(j,:),drydepvg(j,:,:))
          end do
        end if
!!$
!       WET deposition (rainout and washout) for aerosol
!!$
        if ( idust(1) > 0 ) then
         do j=jstart,jend
          call wetdepa(j,nbin,idust,dustbed,rhodust,ttb(j,:,:),wl(j,:,:),fracloud(j,:,:), &
                       fracum(j,:,:),psurf(j,:),hlev,rho(j,:,:),prec(j,:,:),pdepv(j,:,:,:))  
         end do
        end if

       if ( isslt(1) > 0 )  then   
         do j=jstart,jend
         call wetdepa(j,sbin,isslt,ssltbed,rhosslt,ttb(j,:,:),wl(j,:,:),fracloud(j,:,:), &
                      fracum(j,:,:),psurf(j,:),hlev,rho(j,:,:), prec(j,:,:), pdepv(j,:,:,:) )  
         end do
       end if
       if ( icarb(1) > 0 )  then   
         ibin = count( icarb > 0 ) 
          do j=jstart,jend
         call wetdepa(j,ibin,icarb(1:ibin),carbed(1:ibin),rhobchl, &
                     ttb(j,:,:),wl(j,:,:),fracloud(j,:,:), &
                     fracum(j,:,:),psurf(j,:),hlev,rho(j,:,:), prec(j,:,:), pdepv(j,:,:,:) ) 
         end do
         end if
!
!!$
!       Wet Deposition for gasphase species 
!!$
        if ( igaschem == 1 ) then
        checum = d_zero !  fix the interface for this variable ! no effect of cumulus scavenging now
        chevap = d_zero
         do j=jstart,jend        
!          call sethet(j,cza(j,:,:),cht(j,:),ttb(j,:,:),checum(j,:,:),cremrat(j,:,:), &
!                      chevap(j,:,:),dtche,rho(j,:,:),chib(j,:,:,:),iym3,cpsb(j,2:iym2))
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
       do j = jstart,jend
        call gas_phase(j,ktau,secofday,lyear,lmonth,lday)
       end do
       end if
! add tendency due to chemistry reaction (every dtche)      
       chiten(jstart:jend,:,:,:) =  chiten(jstart:jend,:,:,:) +   chemten(jstart:jend,:,:,:)

      end if
 
     

    end subroutine tractend2
!
      subroutine tracbud
      implicit none
!
! Local variables
!
      integer :: i , itr , j , k
!
      do itr = 1 , ntr
        do j = jce1 , jce2
          do i = ice1 , ice2
            dtrace(i,j,itr) = 0.0
            wdlsc(i,j,itr) = 0.0
            wdcvc(i,j,itr) = 0.0
            wxsg(i,j,itr) = 0.0
            wxaq(i,j,itr) = 0.0
            ddsfc(i,j,itr) = 0.0
          end do
        end do
      end do
 
!-----tracers (unit = kg):
      do itr = 1 , ntr
        do j = jce1 , jce2
          do i = ice1 , ice2
            do k = 1 , kz
              dtrace(i,j,itr) = dtrace(i,j,itr) + chia(j,i,k,itr)       &
                              & *cdsigma(k)
              wdlsc(i,j,itr) = wdlsc(i,j,itr) + remlsc(i,k,j,itr)       &
                             & *cdsigma(k)
              wdcvc(i,j,itr) = wdcvc(i,j,itr) + remcvc(i,k,j,itr)       &
                             & *cdsigma(k)
              wxsg(i,j,itr) = wxsg(i,j,itr) + rxsg(i,k,j,itr)*cdsigma(k)
!             sum ls and conv contribution
              wxaq(i,j,itr) = wxaq(i,j,itr)                             &
                            & + (rxsaq1(i,k,j,itr)+rxsaq2(i,k,j,itr))   &
                            & *cdsigma(k)
            end do
            ddsfc(i,j,itr) = ddsfc(i,j,itr) + remdrd(i,j,itr)*cdsigma(kz)
!           Source cumulated diag(care the unit are alredy .m-2)
            cemtrac(i,j,itr) = cemtr(i,j,itr)
          end do
        end do
      end do

      do itr = 1 , ntr
        do j = jce1 , jce2
          do i = ice1 , ice2
            dtrace(i,j,itr) = 1.D6*dtrace(i,j,itr)*d_1000*regrav
                                                        ! unit: mg/m2
            wdlsc(i,j,itr) = 1.D6*wdlsc(i,j,itr)*d_1000*regrav
            wdcvc(i,j,itr) = 1.D6*wdcvc(i,j,itr)*d_1000*regrav
            ddsfc(i,j,itr) = 1.D6*ddsfc(i,j,itr)*d_1000*regrav
            wxsg(i,j,itr) = 1.D6*wxsg(i,j,itr)*d_1000*regrav
            wxaq(i,j,itr) = 1.D6*wxaq(i,j,itr)*d_1000*regrav
!           emtrac isbuilt from chsurfem so just need the 1e6*dt/2
!           factor to to pass im mg/m2
            cemtrac(i,j,itr) = 1.D6*cemtrac(i,j,itr)
          end do
        end do
      end do

      end subroutine tracbud
!
      end module mod_che_tend
