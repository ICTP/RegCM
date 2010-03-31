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

      module mod_bcflds

      use mod_regcm_param , only : jxm2 , iym2 , kz
      use mod_postproc_param , only : nhrbc , npl

      implicit none

      integer , parameter :: ni2d = 2
      integer , parameter :: ni2d2 = 1
      integer , parameter :: ni3d = 4
      integer , parameter :: ni3d2 = 6

      integer , parameter :: nbc2d = ni2d + ni2d2
      integer , parameter :: nbc3d = ni3d + ni3d2
      integer , parameter :: nitot = nbc3d + nbc2d

      real(4) , dimension(jxm2,iym2,nbc2d,nhrbc) :: i2davg
      real(4) , dimension(jxm2,iym2,kz,nbc3d,nhrbc) :: i3davg
      real(4) , dimension(jxm2,iym2,nbc2d) :: ifld2d
      real(4) , dimension(jxm2,iym2,kz,nbc3d) :: ifld3d

      real(4) , dimension(jxm2,iym2,npl,nbc3d,nhrbc) :: i3davg_p
      real(4) , dimension(jxm2,iym2,npl,nbc3d) :: ifld3d_p

      contains

      subroutine rdicbc(idate,iin,irec,ierr)
      use mod_point
 
      implicit none
!
! Dummy arguments
!
      integer :: idate , ierr , iin , irec
      intent (in) iin
      intent (inout) idate , ierr , irec
!
! Local variables
!
      integer :: i , j , k , kk , ni
      real(4) , dimension(jxm2,iym2) :: tmp2d
!
      print * , ''
      ierr = 0
      irec = irec + 1
      read (iin,rec=irec,iostat=ierr) idate
      if ( ierr/=0 ) return
!     print *,'Reading ICBC:  ',idate
      do ni = 1 , ni3d
        do k = 1 , kz
          irec = irec + 1
          read (iin,rec=irec,iostat=ierr) tmp2d
          if ( ierr/=0 ) return
          kk = kz - k + 1
          if ( ni==nui .or. ni==nvi ) then
            do j = 1 , iym2
              do i = 1 , jxm2
                ifld3d(i,j,kk,ni) = 0.25*(tmp2d(i+2,j+2)+tmp2d(i+1,j+1)+&
                                  & tmp2d(i+1,j+2)+tmp2d(i+2,j+1))
              end do
            end do
          else
            do j = 1 , iym2
              do i = 1 , jxm2
                ifld3d(i,j,kk,ni) = tmp2d(i+1,j+1)
!               ifld3d(i,j,k,ni) = tmp2d(i+1,j+1)
              end do
            end do
          end if
        end do
      end do
      do ni = 1 , ni2d
        irec = irec + 1
        read (iin,rec=irec,iostat=ierr) tmp2d
        if ( ierr/=0 ) return
        do j = 1 , iym2
          do i = 1 , jxm2
            if ( ni==npsi ) then
              ifld2d(i,j,ni) = tmp2d(i+1,j+1)*10.
            else
              ifld2d(i,j,ni) = tmp2d(i+1,j+1)
            end if
          end do
        end do
      end do
      print * , 'ICBC DATA READ:' , idate
 
      end subroutine rdicbc

      subroutine mmvlubc(vnambc,lnambc,ubc,xmin,xmax,fact,offset)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      real(4) , dimension(nitot) :: fact , offset , xmax , xmin
      character(64) , dimension(nitot) :: lnambc
      character(64) , dimension(nitot) :: ubc
      character(64) , dimension(nitot) :: vnambc
      intent (out) fact , lnambc , offset , ubc , vnambc
      intent (inout) xmax , xmin
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      vnambc(nui) = 'U'
      vnambc(nvi) = 'V'
      vnambc(nti) = 'TK'
      vnambc(nqvi) = 'QD'
      vnambc(nbc3d+npsi) = 'PS'
      vnambc(nbc3d+ntgi) = 'TGRND'
      vnambc(nbc3d+nslpi) = 'SLP'
      vnambc(nrhi) = 'RH'
      vnambc(nhgti) = 'HGT'
      vnambc(ntdi) = 'TD'
      vnambc(nthi) = 'TH'
      vnambc(nvori) = 'VOR'
      vnambc(ndivi) = 'DIV'
 
      lnambc(nui) = 'Zonal Wind'
      lnambc(nvi) = 'Meridional Wind'
      lnambc(nti) = 'Temperature'
      lnambc(nqvi) = 'Mixing Ratio'
      lnambc(nbc3d+npsi) = 'Surface Pressure'
      lnambc(nbc3d+ntgi) = 'Surface Temperature'
      lnambc(nbc3d+nslpi) = 'Sea Level Temperature'
      lnambc(nrhi) = 'Relative Humidity'
      lnambc(nhgti) = 'Geopotential Height'
      lnambc(ntdi) = 'Dew Point Temperature'
      lnambc(nthi) = 'Potential Temperature'
      lnambc(nvori) = 'Vorticity (Vertical Component)'
      lnambc(nvori) = 'Vorticity (Horizontal Compnent)'
 
      ubc(nui) = 'm/s'
      ubc(nvi) = 'm/s'
      ubc(nti) = 'K'
      ubc(nqvi) = 'kg/kg'
      ubc(nbc3d+npsi) = 'hPa'
      ubc(nbc3d+ntgi) = 'K'
      ubc(nbc3d+nslpi) = 'hPa'
      ubc(nrhi) = 'fraction'
      ubc(nhgti) = 'm'
      ubc(ntdi) = 'K'
      ubc(nthi) = 'K'
      ubc(nvori) = 'm/s'
      ubc(ndivi) = 'm/s'
 
      xmax(nui) = 210.0
      xmax(nvi) = 210.0
      xmax(nti) = 350.0
      xmax(nqvi) = 0.1
      xmax(nbc3d+npsi) = 1200.0
      xmax(nbc3d+ntgi) = 350.0
      xmax(nbc3d+nslpi) = 1200.0
      xmax(nrhi) = 30.0
      xmax(nhgti) = 40000.0
      xmax(ntdi) = 350.0
      xmax(nthi) = 350.0
      xmax(nvori) = 210.0
      xmax(ndivi) = 210.0
 
      xmin(nui) = -210.0
      xmin(nvi) = -210.0
      xmin(nti) = 160.0
      xmin(nqvi) = -0.001
      xmin(nbc3d+npsi) = 300.0
      xmin(nbc3d+ntgi) = 200.0
      xmin(nbc3d+nslpi) = 200.0
      xmin(nrhi) = -0.5
      xmin(nhgti) = -100.0
      xmin(ntdi) = 160.0
      xmin(nthi) = 160.0
      xmin(nvori) = -210.0
      xmin(ndivi) = -210.0
 
 
      aaa = 2.**16. - 1.
      do l = 1 , nitot
        fact(l) = (xmax(l)-xmin(l))/aaa
        offset(l) = (xmax(l)+xmin(l))/2.
      end do
 
      end subroutine mmvlubc

      subroutine writebc(vvarmin,vvarmax,vnambc,lnambc,ubc,xmin,xmax,   &
                       & fact,offset,iadm,ndim,xlat1d,xlon1d,sighrev,   &
                       & vmisdat,idout,xhr,iobctyp,iunt,nrec,lpl,u_bc)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idout , iobctyp , ndim , nrec , iunt
      logical :: lpl
      real(4) :: vmisdat
      real(8) :: xhr
      real(4) , dimension(nitot) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nitot) :: lnambc
      real(4) , dimension(kz) :: sighrev
      character(64) , dimension(nitot) :: ubc
      integer , dimension(nitot) :: u_bc
      character(64) , dimension(nitot) :: vnambc
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , lpl , u_bc , xmax , xmin
!
! Local variables
!
      integer :: i , j , k , ni , nni
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
!
!     **** WRITE OUT 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.lpl ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,ni,vnambc(ni)
            do k = 1 , kz
              do j = 1 , iym2
                do i = 1 , jxm2
                  tmp3d(i,j,k) = max(ifld3d(i,j,k,ni),vmisdat)
                end do
              end do
            end do
            if ( iobctyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                stop 999
              end if
              misdat = xmin(ni)
            else if ( iobctyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iobctyp==1 .or. iobctyp==2 ) then
              call writecdf(idout,vnambc(ni),tmp3d,jxm2,iym2,kz,iadm,   &
                          & xhr,lnambc(ni),ubc(ni),fact(ni),offset(ni), &
                          & vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,    &
                          & misdat,iobctyp)
            else if ( iobctyp==3 ) then
              call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
            else
            end if
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d_p,vmisdat,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,ni,vnambc(ni)
            do k = 1 , npl
              do j = 1 , iym2
                do i = 1 , jxm2
                  tmp3d_p(i,j,k) = max(ifld3d_p(i,j,k,ni),vmisdat)
                end do
              end do
            end do
            if ( iobctyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                stop 999
              end if
              misdat = xmin(ni)
            else if ( iobctyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iobctyp==1 .or. iobctyp==2 ) then
              call writecdf(idout,vnambc(ni),tmp3d_p,jxm2,iym2,npl,iadm,&
                          & xhr,lnambc(ni),ubc(ni),fact(ni),offset(ni), &
                          & vvarmin,vvarmax,xlat1d,xlon1d,plev,0,misdat,&
                          & iobctyp)
            else if ( iobctyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
!     **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do ni = 1 , nbc2d
        nni = ni + nbc3d
        if ( u_bc(nni)==1 ) then
!         print*,ni,nni,vnambc(nni)
          do j = 1 , iym2
            do i = 1 , jxm2
              tmp2d(i,j) = max(ifld2d(i,j,ni),vmisdat)
            end do
          end do
          if ( iobctyp==1 ) then
            misdat = xmin(ni)
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xmin(nni) .or. vmax>xmax(nni) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambc(nni)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nni)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nni)
              stop 999
            end if
          else if ( iobctyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iobctyp==1 .or. iobctyp==2 ) then
            call writecdf(idout,vnambc(nni),tmp2d,jxm2,iym2,1,iadm,xhr, &
                        & lnambc(nni),ubc(nni),fact(nni),offset(nni),   &
                        & vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,      &
                        & misdat,iobctyp)
          else if ( iobctyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
 
      end subroutine writebc

      subroutine writeavgbc(sighrev,vnambc,lnambc,ubc,xmin,xmax,fact,   &
                          & offset,vvarmin,vvarmax,xlat1d,xlon1d,iadm,  &
                          & ndim,xhr1,nbctime,idout,vmisdat,iobctyp,    &
                          & iunt,nrec,plv,u_bc)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idout , iobctyp , ndim , nrec , iunt
      logical :: plv
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nitot) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nitot) :: lnambc
      integer , dimension(nhrbc) :: nbctime
      real(4) , dimension(kz) :: sighrev
      character(64) , dimension(nitot) :: ubc
      integer , dimension(nitot) :: u_bc
      character(64) , dimension(nitot) :: vnambc
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) nbctime , ndim , plv , u_bc , xhr1 , xmax , xmin
!
! Local variables
!
      integer :: i , ihr , j , k , ni , nni
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
      print * , 'COMPUTING AVERAGE ICBC FIELDS:' , nbctime
      xhravg = xhr1
      print * , 'xhravg=' , xhravg
!     **** WRITE ICBC AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,vnambc(ni)
            call setconst(tmp3d,0.0,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrbc
              xntimes = 1./float(nbctime(ihr)*nhrbc)
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( i3davg(i,j,k,ni,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = tmp3d(i,j,k) + i3davg(i,j,k,ni,ihr)&
                                   & *xntimes
                      i3davg(i,j,k,ni,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iobctyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                stop 999
              end if
              misdat = xmin(ni)
            else if ( iobctyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iobctyp==1 .or. iobctyp==2 ) then
              call writecdf(idout,vnambc(ni),tmp3d,jxm2,iym2,kz,iadm,   &
                          & xhravg,lnambc(ni),ubc(ni),fact(ni),         &
                          & offset(ni),vvarmin,vvarmax,xlat1d,xlon1d,   &
                          & sighrev,0,misdat,iobctyp)
            else if ( iobctyp==3 ) then
              call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
            else
            end if
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d_p,vmisdat,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,vnambc(ni)
            call setconst(tmp3d_p,0.0,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrbc
              xntimes = 1./float(nbctime(ihr)*nhrbc)
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( i3davg_p(i,j,k,ni,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = tmp3d_p(i,j,k)                   &
                                     & + i3davg_p(i,j,k,ni,ihr)*xntimes
                      i3davg_p(i,j,k,ni,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iobctyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                stop 999
              end if
              misdat = xmin(ni)
            else if ( iobctyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iobctyp==1 .or. iobctyp==2 ) then
              call writecdf(idout,vnambc(ni),tmp3d_p,jxm2,iym2,npl,iadm,&
                          & xhravg,lnambc(ni),ubc(ni),fact(ni),         &
                          & offset(ni),vvarmin,vvarmax,xlat1d,xlon1d,   &
                          & plev,0,misdat,iobctyp)
            else if ( iobctyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
!     **** WRITE ICBC AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do ni = 1 , nbc2d
        nni = nbc3d + ni
        if ( u_bc(nni)==1 ) then
!         print*,vnambc(nni)
          call setconst(tmp2d,0.0,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
          do ihr = 1 , nhrbc
            xntimes = 1./float(nbctime(ihr)*nhrbc)
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( i2davg(i,j,ni,ihr)>vmisdat ) then
                  tmp2d(i,j) = tmp2d(i,j) + i2davg(i,j,ni,ihr)*xntimes
                  i2davg(i,j,ni,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
          end do
          if ( iobctyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xmin(nni) .or. vmax>xmax(nni) ) then
              print * , 'Values Out of Range:  FIELD=' , vnambc(nni)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nni)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nni)
              stop 999
            end if
            misdat = xmin(ni)
          else if ( iobctyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iobctyp==1 .or. iobctyp==2 ) then
            call writecdf(idout,vnambc(nni),tmp2d,jxm2,iym2,1,iadm,     &
                        & xhravg,lnambc(nni),ubc(nni),fact(nni),        &
                        & offset(nni),vvarmin,vvarmax,xlat1d,xlon1d,    &
                        & sighrev,0,misdat,iobctyp)
          else if ( iobctyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgbc

      subroutine writediurbc(sighrev,vnambc,lnambc,ubc,xmin,xmax,fact,  &
                           & offset,vvarmin,vvarmax,xlat1d,xlon1d,iadm, &
                           & ndim,xhr1,nbctime,idout,vmisdat,iobctyp,   &
                           & iunt,nrec,plv,u_bc)
 
      use mod_point
      use mod_postproc_param , only : plev , dtbc
      implicit none
!
! Dummy arguments
!
      integer :: idout , iobctyp , ndim , nrec , iunt
      logical :: plv
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nitot) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nitot) :: lnambc
      integer , dimension(nhrbc) :: nbctime
      real(4) , dimension(kz) :: sighrev
      character(64) , dimension(nitot) :: ubc
      integer , dimension(nitot) :: u_bc
      character(64) , dimension(nitot) :: vnambc
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) nbctime , ndim , plv , u_bc , xhr1 , xmax , xmin
!
! Local variables
!
      integer :: i , ihr , j , k , ni , nni
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
      print * , 'COMPUTING AVERAGE FIELDS FOR DIURNAL OUTPUT:' , nbctime
 
!     **** WRITE ICBC AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,vnambc(ni)
            do ihr = 1 , nhrbc
              xhravg = xhr1 + float(ihr-1)*dtbc
              xntimes = 1./float(nbctime(ihr))
!             xhravg = float(ihr-1)*dtbc
              if ( nbctime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nbctime = 0'
                stop 999
              end if
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( i3davg(i,j,k,ni,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = i3davg(i,j,k,ni,ihr)*xntimes
                      i3davg(i,j,k,ni,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iobctyp==1 ) then
                call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
                if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                  stop 999
                end if
                misdat = xmin(ni)
              else if ( iobctyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iobctyp==1 .or. iobctyp==2 ) then
                call writecdf(idout,vnambc(ni),tmp3d,jxm2,iym2,kz,iadm, &
                            & xhravg,lnambc(ni),ubc(ni),fact(ni),       &
                            & offset(ni),vvarmin,vvarmax,xlat1d,xlon1d, &
                            & sighrev,0,misdat,iobctyp)
              else if ( iobctyp==3 ) then
                call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
              else
              end if
            end do
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d_p,vmisdat,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
        do ni = 1 , nbc3d
          if ( u_bc(ni)==1 ) then
!           print*,vnambc(ni)
            do ihr = 1 , nhrbc
              xhravg = xhr1 + float(ihr-1)*dtbc
              xntimes = 1./float(nbctime(ihr))
!             xhravg = float(ihr-1)*dtbc
              if ( nbctime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nbctime = 0'
                stop 999
              end if
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( i3davg_p(i,j,k,ni,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = i3davg_p(i,j,k,ni,ihr)*xntimes
                      i3davg_p(i,j,k,ni,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iobctyp==1 ) then
                call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
                if ( vmin<xmin(ni) .or. vmax>xmax(ni) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnambc(ni)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ni)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ni)
                  stop 999
                end if
                misdat = xmin(ni)
              else if ( iobctyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iobctyp==1 .or. iobctyp==2 ) then
                call writecdf(idout,vnambc(ni),tmp3d_p,jxm2,iym2,npl,   &
                            & iadm,xhravg,lnambc(ni),ubc(ni),fact(ni),  &
                            & offset(ni),vvarmin,vvarmax,xlat1d,xlon1d, &
                            & plev,0,misdat,iobctyp)
              else if ( iobctyp==3 ) then
                call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
              else
              end if
 
            end do
          end if
        end do
      end if
 
!     **** WRITE ICBC AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do ni = 1 , nbc2d
        nni = nbc3d + ni
        if ( u_bc(nni)==1 ) then
!         print*,vnambc(nni)
          do ihr = 1 , nhrbc
            xhravg = xhr1 + float(ihr-1)*dtbc
            xntimes = 1./float(nbctime(ihr))
            print * , 'nbctime(ihr)=' , nbctime(ihr) , 'xntimes=' ,     &
                & xntimes , 'ihr=' , ihr , xhravg
            if ( nbctime(ihr)<=0 ) then
              print * , 'NOTHING TO AVERAGE -- nbctime = 0'
              stop 999
            end if
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( i2davg(i,j,ni,ihr)>vmisdat ) then
                  tmp2d(i,j) = i2davg(i,j,ni,ihr)*xntimes
                  i2davg(i,j,ni,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
            if ( iobctyp==1 ) then
              call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
              if ( vmin<xmin(nni) .or. vmax>xmax(nni) ) then
                print * , 'Values Out of Range:  FIELD=' , vnambc(nni)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(nni)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(nni)
                stop 999
              end if
              misdat = xmin(ni)
            else if ( iobctyp==2 ) then
              misdat = vmisdat
            else
            end if
            call writecdf(idout,vnambc(nni),tmp2d,jxm2,iym2,1,iadm,     &
                        & xhravg,lnambc(nni),ubc(nni),fact(nni),        &
                        & offset(nni),vvarmin,vvarmax,xlat1d,xlon1d,    &
                        & sighrev,0,misdat,iobctyp)
          end do
        end if
      end do
 
      end subroutine writediurbc
 
      end module mod_bcflds
