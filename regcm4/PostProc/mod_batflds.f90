!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_batflds

      use mod_regcm_param , only : jxm2 , iym2 , numbat
      use mod_postproc_param , only : nhrbat

      implicit none

      integer , parameter :: nbat2 = numbat + 1

      real(4) , dimension(jxm2,iym2,nbat2,nhrbat) :: b2davg
      real(4) , dimension(jxm2,iym2,nbat2) :: bfld2d

      contains

      subroutine rdsrf(idate,iin,brec,idirect,ierr)
 
      implicit none
!
! Dummy arguments
!
      integer :: brec , idate , idirect , ierr , iin
      intent (in) idirect , iin
      intent (inout) brec , idate , ierr
!
! Local variables
!
      real(4) , dimension(jxm2,iym2) :: b2d
      integer :: i , j , nb
!
      ierr = 0
      if ( idirect/=1 ) then
        read (iin,iostat=ierr) idate
        print * , ' READING SRF (Sequential):' , idate
        if ( ierr/=0 ) return
      else
        print * , ' READING SRF (GrADS):' , idate
      end if
      do nb = 1 , numbat
        if ( idirect==1 ) then
          brec = brec + 1
          read (iin,rec=brec,iostat=ierr) b2d
        else
          read (iin,iostat=ierr) b2d
        end if
        do j = 1 , iym2
          do i = 1 , jxm2
            bfld2d(i,j,nb) = b2d(i,j)
          end do
        end do
      end do
 
      end subroutine rdsrf

      subroutine mmvlubat(vnambat,lnambat,ubat,xmin,xmax,fact,offset)
      use mod_point 
      implicit none
!
! Dummy arguments
!
      real(4) , dimension(nbat2) :: fact , offset , xmax , xmin
      character(64) , dimension(nbat2) :: lnambat
      character(64) , dimension(nbat2) :: ubat
      character(64) , dimension(nbat2) :: vnambat
      intent (out) fact , lnambat , offset , ubat , vnambat
      intent (inout) xmax , xmin
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      lnambat(nux) = 'Anemom Zonal Winds'
      vnambat(nux) = 'UA'
      ubat(nux) = 'm/s'
      xmax(nux) = 50.0
      xmin(nux) = -50.0
 
      lnambat(nvx) = 'Anemom Merid Winds'
      vnambat(nvx) = 'VA'
      ubat(nvx) = 'm/s'
      xmax(nvx) = 50.0
      xmin(nvx) = -50.0
 
      lnambat(ndrag) = 'Surface Drag Stress'
      vnambat(ndrag) = 'DRAG'
      ubat(ndrag) = 'si'
      xmax(ndrag) = 1.0
      xmin(ndrag) = -1.0
 
      vnambat(ntg) = 'TG'
      lnambat(ntg) = 'Ground Temperature'
      ubat(ntg) = 'K'
      xmax(ntg) = 350.0
      xmin(ntg) = 180.0
 
      vnambat(ntf) = 'TF'
      lnambat(ntf) = 'Foliage Temp'
      ubat(ntf) = 'K'
      xmax(ntf) = 350.0
      xmin(ntf) = 180.0
 
      lnambat(ntanm) = 'Anemom Temp'
      vnambat(ntanm) = 'TA'
      ubat(ntanm) = 'K'
      xmax(ntanm) = 350.0
      xmin(ntanm) = 180.0
 
      lnambat(nqanm) = 'Anemom Spec Humidity'
      vnambat(nqanm) = 'QA'
      ubat(nqanm) = 'kg/kg'
      xmax(nqanm) = 0.20
      xmin(nqanm) = -1.0E-5
 
      lnambat(nsmu) = 'Top Layer Soil Moist'
      vnambat(nsmu) = 'SMU'
      ubat(nsmu) = 'mm'
      xmax(nsmu) = 80.0
      xmin(nsmu) = -1.0
 
      lnambat(nsmr) = 'Root Lay Soil Moist'
      vnambat(nsmr) = 'SMR'
      ubat(nsmr) = 'mm'
      xmax(nsmr) = 1200.0
      xmin(nsmr) = -1.0
 
      lnambat(net) = 'Evapotranspiration'
      vnambat(net) = 'ET'
      ubat(net) = 'mm/day'
      xmax(net) = 150.0
      xmin(net) = -5.0
 
      lnambat(nrnfs) = 'Surface Runoff'
      vnambat(nrnfs) = 'RNFS'
      ubat(nrnfs) = 'mm/day'
      xmax(nrnfs) = 2000.0
      xmin(nrnfs) = -200.0
 
      lnambat(nsnow) = 'Snow Depth'
      vnambat(nsnow) = 'SNOW'
      ubat(nsnow) = 'mm H2O'
      xmax(nsnow) = 1000.0
      xmin(nsnow) = -1.0
 
      lnambat(nsh) = 'Sensible Heat'
      vnambat(nsh) = 'SH'
      ubat(nsh) = 'W/m2'
      xmax(nsh) = 1000.0
      xmin(nsh) = -300.0
 
      lnambat(nlwn) = 'Net Longwave'
      vnambat(nlwn) = 'LWN'
      ubat(nlwn) = 'W/m2'
      xmax(nlwn) = 750.0
      xmin(nlwn) = -300.0
 
      lnambat(nlwd) = 'Downward Longwave'
      vnambat(nlwd) = 'LWD'
      ubat(nlwd) = 'W/m2'
      xmax(nlwd) = 750.0
      xmin(nlwd) = -300.0
 
      lnambat(nswn) = 'Net Solar Absorbed'
      vnambat(nswn) = 'SWN'
      ubat(nswn) = 'W/m2'
      xmax(nswn) = 1200.0
      xmin(nswn) = -1.0
 
      lnambat(nswi) = 'Solar Incident'
      vnambat(nswi) = 'SWI'
      ubat(nswi) = 'W/m2'
      xmax(nswi) = 1400.0
      xmin(nswi) = -1.0
 
      lnambat(nprc) = 'Convective Precip'
      vnambat(nprc) = 'RC'
      ubat(nprc) = 'mm/day'
      xmax(nprc) = 1500.0
      xmin(nprc) = -1.0
 
      lnambat(npt) = 'Total Precipitation'
      vnambat(npt) = 'RT'
      ubat(npt) = 'mm/day'
      xmax(npt) = 2500.0
      xmin(npt) = -1.0
 
      lnambat(kxpbl) = 'PBL Height'
      vnambat(kxpbl) = 'ZPBL'
      ubat(kxpbl) = 'm'
      xmax(kxpbl) = 6000.0
      xmin(kxpbl) = -1.0
 
      lnambat(npsrf) = 'Surface Pressure'
      vnambat(npsrf) = 'PSRF'
      ubat(npsrf) = 'hPa'
      xmax(npsrf) = 1500.0
      xmin(npsrf) = 300.0
 
 
      lnambat(nrha) = 'Relative Humidity'
      vnambat(nrha) = 'RHA'
      ubat(nrha) = 'fraction'
      xmax(nrha) = 5.0
      xmin(nrha) = -0.1
 
      lnambat(ntgmax) = 'Max Ground Temp'
      vnambat(ntgmax) = 'TGMAX'
      ubat(ntgmax) = 'K'
      xmax(ntgmax) = 350.0
      xmin(ntgmax) = 200.0
 
      lnambat(ntgmin) = 'Min Ground Temp'
      vnambat(ntgmin) = 'TGMIN'
      ubat(ntgmin) = 'K'
      xmax(ntgmin) = 350.0
      xmin(ntgmin) = 200.0
 
      lnambat(ntamax) = 'Max Anemom Temp'
      vnambat(ntamax) = 'TAMAX'
      ubat(ntamax) = 'K'
      xmax(ntamax) = 350.0
      xmin(ntamax) = 200.0
 
      lnambat(ntamin) = 'Min Anemom Temp'
      vnambat(ntamin) = 'TAMIN'
      ubat(ntamin) = 'K'
      xmax(ntamin) = 350.0
      xmin(ntamin) = 200.0
 
      lnambat(w10max) = 'Max 10m Wind Speed'
      vnambat(w10max) = 'W10MX'
      ubat(w10max) = 'm/s'
      xmax(w10max) = 500.0
      xmin(w10max) = -500.0
 
      lnambat(psmin) = 'Min Surface Pressure'
      vnambat(psmin) = 'PSMIN'
      ubat(psmin) = 'hPa'
      xmax(psmin) = 1500.0
      xmin(psmin) = 300.0
 
      aaa = 2.**16. - 1.
      do l = 1 , nbat2
        fact(l) = (xmax(l)-xmin(l))/aaa
        offset(l) = (xmax(l)+xmin(l))/2.
      end do
 
      end subroutine mmvlubat

      subroutine avgdatabat(ihr,vmisdat)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: ihr
      real(4) :: vmisdat
      intent (in) ihr , vmisdat
!
! Local variables
!
      integer :: i , j , nb
      real(4) :: misdat
! 
      if ( vmisdat>0.0 ) then
        misdat = -1.0*vmisdat
      else
        misdat = vmisdat
      end if
      do nb = 1 , nbat2
!       if ((nb.eq.ntmin.or.nb.eq.ntmax) .and. ihr.lt.nhrbat) then
!       else
        do j = 1 , iym2
          do i = 1 , jxm2
            if ( bfld2d(i,j,nb)>misdat ) then
              b2davg(i,j,nb,ihr) = b2davg(i,j,nb,ihr) + bfld2d(i,j,nb)
            else
              b2davg(i,j,nb,ihr) = vmisdat
            end if
          end do
        end do
!       end if
      end do
      end subroutine avgdatabat
!
      subroutine writebat(vnambat,lnambat,ubat,xmin,xmax,fact,offset,   &
                        & vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,      &
                        & vmisdat,xhr,idout,iotyp,iunt,nrec,u_bat)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr
      real(4) , dimension(nbat2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nbat2) :: lnambat
      character(64) , dimension(nbat2) :: ubat
      integer , dimension(nbat2) :: u_bat
      character(64) , dimension(nbat2) :: vnambat
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , u_bat , xmax , xmin
!
! Local variables
!
      integer :: i , j , nb
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxm2,iym2) :: tmp2d
!
      sig1(1) = 1.
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do nb = 1 , nbat2
        if ( u_bat(nb)==1 ) then
!         if (nb.ne.ntmin .and. nb.ne.ntmax) then
          do i = 1 , jxm2
            do j = 1 , iym2
              tmp2d(i,j) = max(bfld2d(i,j,nb),vmisdat)
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xmin(nb) .or. vmax>xmax(nb) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nb)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nb)
!             stop 999
            end if
            misdat = xmin(nb)
          else if ( iotyp==2 .or. iotyp==3 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idout,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,xhr, &
                        & lnambat(nb),ubat(nb),fact(nb),offset(nb),     &
                        & vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,  &
                        & iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writebat

      subroutine writeavgbat(vmisdat,vnambat,lnambat,ubat,xmin,xmax,    &
                           & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d, &
                           & iadm,ndim,xhr1,nbattime,idbat,iotyp,iunt,  &
                           & nrec,u_bat)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idbat , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nbat2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nbat2) :: lnambat
      integer , dimension(nhrbat) :: nbattime
      character(64) , dimension(nbat2) :: ubat
      integer , dimension(nbat2) :: u_bat
      character(64) , dimension(nbat2) :: vnambat
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) nbattime , ndim , u_bat , xhr1 , xmax , xmin
!
! Local variables
!
      real(4) :: const , misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2,nbat2) :: favgsum
      integer :: i , ihr , j , nb
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(8) :: xhravg
!
      iadm(3) = 1
      sig1(1) = 1.
 
      print * , 'COMPUTING AVERAGE BAT FIELDS:' , nbattime
 
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      call setconst(favgsum,0.0,jxm2,iym2,nbat2,1,1,1,jxm2,1,iym2)
      do ihr = 1 , nhrbat
        do nb = 1 , nbat2
          if ( u_bat(nb)==1 ) then
            const = 1.0
!           if (nb.eq.ntmin .or. nb.eq.ntmax) then
!           xntimes = const/float(nbattime(ihr))
!           else
            xntimes = const/float(nhrbat*nbattime(ihr))
!           end if
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( b2davg(i,j,nb,ihr)>vmisdat ) then
                  favgsum(i,j,nb) = favgsum(i,j,nb) + b2davg(i,j,nb,ihr)&
                                  & *xntimes
                else
                  favgsum(i,j,nb) = vmisdat
                end if
              end do
            end do
          end if
        end do
      end do
      xhravg = xhr1
      do nb = 1 , nbat2
        if ( u_bat(nb)==1 ) then
          if ( xntimes<=0 ) then
            print * , 'NOTHING TO AVERAGE -- nbattime = 0'
            stop 999
          end if
          do j = 1 , iym2
            do i = 1 , jxm2
              tmp2d(i,j) = favgsum(i,j,nb)
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xmin(nb) .or. vmax>xmax(nb) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nb)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nb)
              stop 999
            end if
            misdat = xmin(nb)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idbat,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,     &
                        & xhravg,lnambat(nb),ubat(nb),fact(nb),         &
                        & offset(nb),vvarmin,vvarmax,xlat1d,xlon1d,sig1,&
                        & 0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgbat
!
      subroutine writediurbat(vmisdat,vnambat,lnambat,ubat,xmin,xmax,   &
                            & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,&
                            & iadm,ndim,xhr1,nbattime,idbat,iotyp,iunt, &
                            & nrec,u_bat)
 
      use mod_point
      use mod_postproc_param , only : dtbat
      implicit none
!
! Dummy arguments
!
      integer :: idbat , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nbat2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nbat2) :: lnambat
      integer , dimension(nhrbat) :: nbattime
      character(64) , dimension(nbat2) :: ubat
      integer , dimension(nbat2) :: u_bat
      character(64) , dimension(nbat2) :: vnambat
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) nbattime , ndim , u_bat , xhr1 , xmax , xmin
!
! Local variables
!
      real(4) :: const , misdat , vmax , vmin , xntimes
      integer :: i , ihr , j , nb
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(8) :: xhravg
!
      iadm(3) = 1
      sig1(1) = 1.
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do ihr = 1 , nhrbat
        xhravg = xhr1 + float(ihr-1)*dtbat
        if ( nbattime(ihr)<=0 ) then
          print * , 'NOTHING TO AVERAGE -- nbattime = 0'
          stop 999
        end if
        const = 1.0
        do nb = 1 , nbat2
          if ( u_bat(nb)==1 ) then
!           if (nb.eq.ntmin .or. nb.eq.ntmax) then
!           xntimes = const/float(nbattime(ihr))
!           else
            xntimes = const/float(nbattime(ihr))
!           end if
!           if ((nb.eq.ntmin.or.nb.eq.ntmax) .and. ihr.lt.nhrbat) then
            do j = 1 , iym2
              do i = 1 , jxm2
                tmp2d(i,j) = vmisdat
              end do
            end do
!           else
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( b2davg(i,j,nb,ihr)>vmisdat ) then
                  tmp2d(i,j) = b2davg(i,j,nb,ihr)*xntimes
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
!           end if
            if ( iotyp==1 ) then
              call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
              if ( vmin<xmin(nb) .or. vmax>xmax(nb) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nb)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nb)
                stop 999
              end if
              misdat = xmin(nb)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idbat,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,   &
                          & xhravg,lnambat(nb),ubat(nb),fact(nb),       &
                          & offset(nb),vvarmin,vvarmax,xlat1d,xlon1d,   &
                          & sig1,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
            else
            end if
          end if
        end do
      end do
      end subroutine writediurbat
 
      end module mod_batflds
