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

      module mod_batflds

      use mod_dynparam
      use mod_postproc_param , only : nhrbat

      implicit none

      integer :: nbat2

      real(4) , allocatable , dimension(:,:,:,:) :: b2davg
      real(4) , allocatable , dimension(:,:,:) :: bfld2d
      real(4) , allocatable , dimension(:) :: factbat , offsetbat ,     &
                     & xmaxbat , xminbat
      character(64) , allocatable , dimension(:) :: vnambat
      character(64) , allocatable , dimension(:) :: lnambat
      character(64) , allocatable , dimension(:) :: ubat
      integer , allocatable , dimension(:) :: nbattime
      integer , allocatable , dimension(:) :: u_bat

      contains

      subroutine init_mod_batflds
        implicit none
        nbat2 = numbat + 1
        allocate(b2davg(jxm2,iym2,nbat2,nhrbat))
        allocate(bfld2d(jxm2,iym2,nbat2))
        allocate(factbat(nbat2))
        allocate(offsetbat(nbat2))
        allocate(xmaxbat(nbat2))
        allocate(xminbat(nbat2))
        allocate(vnambat(nbat2))
        allocate(lnambat(nbat2))
        allocate(ubat(nbat2))
        allocate(nbattime(nhrbat))
        allocate(u_bat(nbat2))
      end subroutine init_mod_batflds

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

      subroutine mmvlubat
      use mod_point 
      implicit none
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      lnambat(nux) = 'Anemom Zonal Winds'
      vnambat(nux) = 'UA'
      ubat(nux) = 'm/s'
      xmaxbat(nux) = 50.0
      xminbat(nux) = -50.0
 
      lnambat(nvx) = 'Anemom Merid Winds'
      vnambat(nvx) = 'VA'
      ubat(nvx) = 'm/s'
      xmaxbat(nvx) = 50.0
      xminbat(nvx) = -50.0
 
      lnambat(ndrag) = 'Surface Drag Stress'
      vnambat(ndrag) = 'DRAG'
      ubat(ndrag) = 'si'
      xmaxbat(ndrag) = 1.0
      xminbat(ndrag) = -1.0
 
      vnambat(ntg) = 'TG'
      lnambat(ntg) = 'Ground Temperature'
      ubat(ntg) = 'K'
      xmaxbat(ntg) = 350.0
      xminbat(ntg) = 180.0
 
      vnambat(ntf) = 'TF'
      lnambat(ntf) = 'Foliage Temp'
      ubat(ntf) = 'K'
      xmaxbat(ntf) = 350.0
      xminbat(ntf) = 180.0
 
      lnambat(ntanm) = 'Anemom Temp'
      vnambat(ntanm) = 'TA'
      ubat(ntanm) = 'K'
      xmaxbat(ntanm) = 350.0
      xminbat(ntanm) = 180.0
 
      lnambat(nqanm) = 'Anemom Spec Humidity'
      vnambat(nqanm) = 'QA'
      ubat(nqanm) = 'kg/kg'
      xmaxbat(nqanm) = 0.20
      xminbat(nqanm) = -1.0E-5
 
      lnambat(nsmu) = 'Top Layer Soil Moist'
      vnambat(nsmu) = 'SMU'
      ubat(nsmu) = 'mm'
      xmaxbat(nsmu) = 80.0
      xminbat(nsmu) = -1.0
 
      lnambat(nsmr) = 'Root Lay Soil Moist'
      vnambat(nsmr) = 'SMR'
      ubat(nsmr) = 'mm'
      xmaxbat(nsmr) = 1200.0
      xminbat(nsmr) = -1.0
 
      lnambat(net) = 'Evapotranspiration'
      vnambat(net) = 'ET'
      ubat(net) = 'mm/day'
      xmaxbat(net) = 150.0
      xminbat(net) = -5.0
 
      lnambat(nrnfs) = 'Surface Runoff'
      vnambat(nrnfs) = 'RNFS'
      ubat(nrnfs) = 'mm/day'
      xmaxbat(nrnfs) = 2000.0
      xminbat(nrnfs) = -200.0
 
      lnambat(nsnow) = 'Snow Depth'
      vnambat(nsnow) = 'SNOW'
      ubat(nsnow) = 'mm H2O'
      xmaxbat(nsnow) = 1000.0
      xminbat(nsnow) = -1.0
 
      lnambat(nsh) = 'Sensible Heat'
      vnambat(nsh) = 'SH'
      ubat(nsh) = 'W/m2'
      xmaxbat(nsh) = 1000.0
      xminbat(nsh) = -300.0
 
      lnambat(nlwn) = 'Net Longwave'
      vnambat(nlwn) = 'LWN'
      ubat(nlwn) = 'W/m2'
      xmaxbat(nlwn) = 750.0
      xminbat(nlwn) = -300.0
 
      lnambat(nlwd) = 'Downward Longwave'
      vnambat(nlwd) = 'LWD'
      ubat(nlwd) = 'W/m2'
      xmaxbat(nlwd) = 750.0
      xminbat(nlwd) = -300.0
 
      lnambat(nswn) = 'Net Solar Absorbed'
      vnambat(nswn) = 'SWN'
      ubat(nswn) = 'W/m2'
      xmaxbat(nswn) = 1200.0
      xminbat(nswn) = -1.0
 
      lnambat(nswi) = 'Solar Incident'
      vnambat(nswi) = 'SWI'
      ubat(nswi) = 'W/m2'
      xmaxbat(nswi) = 1400.0
      xminbat(nswi) = -1.0
 
      lnambat(nprc) = 'Convective Precip'
      vnambat(nprc) = 'RC'
      ubat(nprc) = 'mm/day'
      xmaxbat(nprc) = 1500.0
      xminbat(nprc) = -1.0
 
      lnambat(npt) = 'Total Precipitation'
      vnambat(npt) = 'RT'
      ubat(npt) = 'mm/day'
      xmaxbat(npt) = 2500.0
      xminbat(npt) = -1.0
 
      lnambat(kxpbl) = 'PBL Height'
      vnambat(kxpbl) = 'ZPBL'
      ubat(kxpbl) = 'm'
      xmaxbat(kxpbl) = 6000.0
      xminbat(kxpbl) = -1.0
 
      lnambat(npsrf) = 'Surface Pressure'
      vnambat(npsrf) = 'PSRF'
      ubat(npsrf) = 'hPa'
      xmaxbat(npsrf) = 1500.0
      xminbat(npsrf) = 300.0
 
 
      lnambat(nrha) = 'Relative Humidity'
      vnambat(nrha) = 'RHA'
      ubat(nrha) = 'fraction'
      xmaxbat(nrha) = 5.0
      xminbat(nrha) = -0.1
 
      lnambat(ntgmax) = 'Max Ground Temp'
      vnambat(ntgmax) = 'TGMAX'
      ubat(ntgmax) = 'K'
      xmaxbat(ntgmax) = 350.0
      xminbat(ntgmax) = 200.0
 
      lnambat(ntgmin) = 'Min Ground Temp'
      vnambat(ntgmin) = 'TGMIN'
      ubat(ntgmin) = 'K'
      xmaxbat(ntgmin) = 350.0
      xminbat(ntgmin) = 200.0
 
      lnambat(ntamax) = 'Max Anemom Temp'
      vnambat(ntamax) = 'TAMAX'
      ubat(ntamax) = 'K'
      xmaxbat(ntamax) = 350.0
      xminbat(ntamax) = 200.0
 
      lnambat(ntamin) = 'Min Anemom Temp'
      vnambat(ntamin) = 'TAMIN'
      ubat(ntamin) = 'K'
      xmaxbat(ntamin) = 350.0
      xminbat(ntamin) = 200.0
 
      lnambat(w10max) = 'Max 10m Wind Speed'
      vnambat(w10max) = 'W10MX'
      ubat(w10max) = 'm/s'
      xmaxbat(w10max) = 500.0
      xminbat(w10max) = -500.0
 
      lnambat(psmin) = 'Min Surface Pressure'
      vnambat(psmin) = 'PSMIN'
      ubat(psmin) = 'hPa'
      xmaxbat(psmin) = 1500.0
      xminbat(psmin) = 300.0
 
      aaa = 2.**16. - 1.
      do l = 1 , nbat2
        factbat(l) = (xmaxbat(l)-xminbat(l))/aaa
        offsetbat(l) = (xmaxbat(l)+xminbat(l))/2.
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
      subroutine writebat(vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,      &
                        & vmisdat,xhr,idout,iotyp,iunt,nrec)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr
      integer , dimension(ndim) :: iadm
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim
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
            if ( vmin<xminbat(nb) .or. vmax>xmaxbat(nb) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminbat(nb)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxbat(nb)
!             stop 999
            end if
            misdat = xminbat(nb)
          else if ( iotyp==2 .or. iotyp==3 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idout,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,xhr, &
                        & lnambat(nb),ubat(nb),factbat(nb),             &
                        & offsetbat(nb),vvarmin,vvarmax,xlat1d,xlon1d,  &
                        & sig1,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writebat

      subroutine writeavgbat(vmisdat,vvarmin,vvarmax,xlat1d,xlon1d,     &
                           & iadm,ndim,xhr1,idbat,iotyp,iunt,nrec)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idbat , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      integer , dimension(ndim) :: iadm
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , xhr1
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
            if ( vmin<xminbat(nb) .or. vmax>xmaxbat(nb) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminbat(nb)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxbat(nb)
              stop 999
            end if
            misdat = xminbat(nb)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idbat,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,     &
                        & xhravg,lnambat(nb),ubat(nb),factbat(nb),      &
                        & offsetbat(nb),vvarmin,vvarmax,xlat1d,xlon1d,  &
                        & sig1,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgbat
!
      subroutine writediurbat(vmisdat,vvarmin,vvarmax,xlat1d,xlon1d,    &
                            & iadm,ndim,xhr1,idbat,iotyp,iunt,nrec)
 
      use mod_point
      use mod_postproc_param , only : dtbat
      implicit none
!
! Dummy arguments
!
      integer :: idbat , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      integer , dimension(ndim) :: iadm
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , xhr1
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
              if ( vmin<xminbat(nb) .or. vmax>xmaxbat(nb) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambat(nb)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminbat(nb)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxbat(nb)
                stop 999
              end if
              misdat = xminbat(nb)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idbat,vnambat(nb),tmp2d,jxm2,iym2,1,iadm,   &
                          & xhravg,lnambat(nb),ubat(nb),factbat(nb),    &
                          & offsetbat(nb),vvarmin,vvarmax,xlat1d,xlon1d,&
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
