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

      module mod_outflds

      use mod_dynparam
      use mod_postproc_param , only : nhrout , npl

      implicit none

      integer , parameter :: no2d = 5
      integer , parameter :: no2d2 = 1
      integer , parameter :: no3d = 6
      integer , parameter :: no3d2 = 6

      integer , parameter :: nout2d = no2d + no2d2
      integer , parameter :: nout3d = no3d + no3d2
      integer , parameter :: notot = nout3d + nout2d

      real(4) , allocatable , dimension(:,:,:,:) :: o2davg
      real(4) , allocatable , dimension(:,:,:,:,:) :: o3davg
      real(4) , allocatable , dimension(:,:,:) :: ofld2d
      real(4) , allocatable , dimension(:,:,:,:) :: ofld3d

      real(4) , allocatable , dimension(:,:,:,:,:) :: o3davg_p
      real(4) , allocatable , dimension(:,:,:,:) :: ofld3d_p

      real(4) , allocatable , dimension(:) :: factout , offsetout ,     &
                            &                 xmaxout , xminout

      character(64) , allocatable , dimension(:) :: vnamout
      character(64) , allocatable , dimension(:) :: lnamout
      character(64) , allocatable , dimension(:) :: uout

      integer , allocatable , dimension(:) :: u_out
      integer , allocatable , dimension(:) :: nouttime

      contains

      subroutine init_mod_outflds
      implicit none
      allocate(o2davg(jxm2,iym2,nout2d,nhrout))
      allocate(o3davg(jxm2,iym2,kz,nout3d,nhrout))
      allocate(ofld2d(jxm2,iym2,nout2d))
      allocate(ofld3d(jxm2,iym2,kz,nout3d))
      allocate(o3davg_p(jxm2,iym2,npl,nout3d,nhrout))
      allocate(ofld3d_p(jxm2,iym2,npl,nout3d))
      allocate(factout(notot))
      allocate(offsetout(notot))
      allocate(xmaxout(notot))
      allocate(xminout(notot))
      allocate(lnamout(notot))
      allocate(vnamout(notot))
      allocate(uout(notot))
      allocate(u_out(notot))
      allocate(nouttime(nhrout))
      end subroutine init_mod_outflds

      subroutine rdatm(idate,iin,orec,idirect,ierr)
      use mod_point 
      implicit none
!
! Dummy arguments
!
      integer :: idate , idirect , ierr , iin , orec
      intent (in) idirect , iin
      intent (inout) idate , ierr , orec
!
! Local variables
!
      integer :: i , j , k , kk , no
      real(4) , dimension(jxm2,iym2) :: tmp2d
!
      print * , ' '
      ierr = 0
      if ( idirect/=1 ) then
        read (iin,iostat=ierr) idate
        if ( ierr/=0 ) return
      end if
      print * , 'Reading output:  ' , idate
      do no = 1 , no3d
        do k = 1 , kz
          if ( idirect==1 ) then
            orec = orec + 1
            read (iin,rec=orec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
          if ( ierr/=0 ) return
          kk = kz - k + 1
          do j = 1 , iym2
            do i = 1 , jxm2
!             ofld3d(i,j,k,no) = tmp2d(i,j)
              ofld3d(i,j,kk,no) = tmp2d(i,j)
            end do
          end do
        end do
      end do
      do no = 1 , no2d
        if ( idirect==1 ) then
          orec = orec + 1
          read (iin,rec=orec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
        if ( ierr/=0 ) return
        do j = 1 , iym2
          do i = 1 , jxm2
            ofld2d(i,j,no) = tmp2d(i,j)
          end do
        end do
      end do
!     print*,'DONE READING OUTPUT FOR CURRENT TIMESTEP',idate
 
      end subroutine rdatm

      subroutine mmvluout
 
      use mod_point
      implicit none
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      vnamout(nua) = 'U'
      vnamout(nva) = 'V'
      vnamout(nomega) = 'OMEGA'
      vnamout(nta) = 'TK'
      vnamout(nqva) = 'QD'
      vnamout(nqca) = 'QC'
      vnamout(nout3d+npsa) = 'PS'
      vnamout(nout3d+nrt) = 'RT'
      vnamout(nout3d+ntgb) = 'TGRND'
      vnamout(nout3d+nsmt) = 'SMT'
      vnamout(nout3d+nbf) = 'RB'
      vnamout(nout3d+nslp) = 'SLP'
      vnamout(nrh) = 'RH'
      vnamout(nhgt) = 'HGT'
      vnamout(ntda) = 'TD'
      vnamout(ntha) = 'TH'
      vnamout(nvora) = 'VOR'
      vnamout(ndiva) = 'DIV'
 
      lnamout(nua) = 'Zonal Wind'
      lnamout(nva) = 'Meridional Wind'
      lnamout(nomega) = 'Omega'
      lnamout(nta) = 'Temperature'
      lnamout(nqva) = 'Mixing Ratio'
      lnamout(nqca) = 'Cloud Mixing Ratio'
      lnamout(nout3d+npsa) = 'Surface Pressure'
      lnamout(nout3d+ntgb) = 'Ground Temperature'
      lnamout(nout3d+nrt) = 'Total Precip'
      lnamout(nout3d+nsmt) = 'Total Soil Water'
      lnamout(nout3d+nbf) = 'Base Flow'
      lnamout(nout3d+nslp) = 'Sea Level Temperature'
      lnamout(nrh) = 'Relative Humidity'
      lnamout(nhgt) = 'Geopotential Height'
      lnamout(ntda) = 'Dew Point Temperature'
      lnamout(ntha) = 'Potential Temperature'
      lnamout(nvora) = 'Vorticity (Vertical Component)'
      lnamout(nvora) = 'Vorticity (Horizontal Compnent)'
 
      uout(nua) = 'm/s'
      uout(nva) = 'm/s'
      uout(nomega) = 'hPa'
      uout(nta) = 'K'
      uout(nqva) = 'kg/kg'
      uout(nqca) = 'kg/kg'
      uout(nout3d+npsa) = 'hPa'
      uout(nout3d+ntgb) = 'K'
      uout(nout3d+nrt) = 'mm/day'
      uout(nout3d+nsmt) = 'mm'
      uout(nout3d+nbf) = 'mm/day'
      uout(nout3d+nslp) = 'hPa'
      uout(nrh) = 'fraction'
      uout(nhgt) = 'm'
      uout(ntda) = 'K'
      uout(ntha) = 'K'
      uout(nvora) = 'm/s'
      uout(ndiva) = 'm/s'
 
      xmaxout(nua) = 210.0
      xmaxout(nva) = 210.0
      xmaxout(nomega) = 0.1
      xmaxout(nta) = 350.0
      xmaxout(nqva) = 0.1
      xmaxout(nqca) = 0.1
      xmaxout(nout3d+npsa) = 1200.0
      xmaxout(nout3d+ntgb) = 350.0
      xmaxout(nout3d+nrt) = 2500.0
      xmaxout(nout3d+nsmt) = 3000.0
      xmaxout(nout3d+nbf) = 200.0
      xmaxout(nout3d+nslp) = 1200.0
      xmaxout(nrh) = 30.0
      xmaxout(nhgt) = 40000.0
      xmaxout(ntda) = 350.0
      xmaxout(ntha) = 350.0
      xmaxout(nvora) = 210.0
      xmaxout(ndiva) = 210.0
 
      xminout(nua) = -210.0
      xminout(nva) = -210.0
      xminout(nomega) = -0.1
      xminout(nta) = 160.0
      xminout(nqva) = -0.001
      xminout(nqca) = -0.001
      xminout(nout3d+npsa) = 300.0
      xminout(nout3d+ntgb) = 180.0
      xminout(nout3d+nrt) = -10.0
      xminout(nout3d+nsmt) = 0.0
      xminout(nout3d+nbf) = -10.0
      xminout(nout3d+nslp) = 200.0
      xminout(nrh) = -0.5
      xminout(nhgt) = -100.0
      xminout(ntda) = 160.0
      xminout(ntha) = 160.0
      xminout(nvora) = -210.0
      xminout(ndiva) = -210.0
 
      aaa = 2.**16. - 1.
      do l = 1 , notot
        factout(l) = (xmaxout(l)-xminout(l))/aaa
        offsetout(l) = (xmaxout(l)+xminout(l))/2.
      end do
 
      end subroutine mmvluout

      subroutine writeout(vvarmin,vvarmax,iadm,ndim,xlat1d,xlon1d,      &
                        & sighrev,vmisdat,idout,xhr,iotyp,iunt,nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      logical :: plv
      real(4) :: vmisdat
      real(8) :: xhr
      integer , dimension(ndim) :: iadm
      real(4) , dimension(kz) :: sighrev
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , plv
!
! Local variables
!
      integer :: i , j , k , nno , no
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
!
!     **** WRITE OUT 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,no,vnamout(no)
            do k = 1 , kz
              do j = 1 , iym2
                do i = 1 , jxm2
                  tmp3d(i,j,k) = max(ofld3d(i,j,k,no),vmisdat)
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                stop 999
              end if
              misdat = xminout(no)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idout,vnamout(no),tmp3d,jxm2,iym2,kz,iadm,  &
                         & xhr,lnamout(no),uout(no),factout(no),        &
                         & offsetout(no),vvarmin,vvarmax,xlat1d,xlon1d, &
                         & sighrev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
            else
            end if
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d,vmisdat,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,no,vnamout(no)
            do k = 1 , npl
              do j = 1 , iym2
                do i = 1 , jxm2
                  tmp3d_p(i,j,k) = max(ofld3d_p(i,j,k,no),vmisdat)
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                stop 999
              end if
              misdat = xminout(no)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idout,vnamout(no),tmp3d_p,jxm2,iym2,npl,    &
                       & iadm,xhr,lnamout(no),uout(no),factout(no),     &
                       & offsetout(no),vvarmin,vvarmax,xlat1d,xlon1d,   &
                       & plev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
 
!     **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do no = 1 , nout2d
        nno = no + nout3d
        if ( u_out(nno)==1 ) then
!         print*,no,nno,vnamout(nno)
          do j = 1 , iym2
            do i = 1 , jxm2
              tmp2d(i,j) = max(ofld2d(i,j,no),vmisdat)
            end do
          end do
          if ( iotyp==1 ) then
            misdat = xminout(no)
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xminout(nno) .or. vmax>xmaxout(nno) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamout(nno)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(nno)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(nno)
              stop 999
            end if
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idout,vnamout(nno),tmp2d,jxm2,iym2,1,iadm,xhr,&
                        & lnamout(nno),uout(nno),factout(nno),          &
                        & offsetout(nno),vvarmin,vvarmax,xlat1d,xlon1d, &
                        & sighrev,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeout

      subroutine writeavgout(sighrev,vvarmin,vvarmax,xlat1d,xlon1d,     &
                           & iadm,ndim,xhr1,idout,vmisdat,iotyp,iunt,   &
                           & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      logical :: plv
      real(4) :: vmisdat
      real(8) :: xhr1
      integer , dimension(ndim) :: iadm
      real(4) , dimension(kz) :: sighrev
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , plv , xhr1
!
! Local variables
!
      integer :: i , ihr , j , k , nno , no
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
      print * , 'COMPUTING AVERAGE OUT FIELDS:' , nouttime
      xhravg = xhr1
      print * , 'xhravg=' , xhravg
!     **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,vnamout(no)
            call setconst(tmp3d,0.0,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrout
              xntimes = 1./float(nouttime(ihr)*nhrout)
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( o3davg(i,j,k,no,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = tmp3d(i,j,k) + o3davg(i,j,k,no,ihr)&
                                   & *xntimes
                      o3davg(i,j,k,no,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                stop 999
              end if
              misdat = xminout(no)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idout,vnamout(no),tmp3d,jxm2,iym2,kz,iadm,  &
                          & xhravg,lnamout(no),uout(no),factout(no),    &
                          & offsetout(no),vvarmin,vvarmax,xlat1d,xlon1d,&
                          & sighrev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
            else
            end if
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,vnamout(no)
            call setconst(tmp3d,0.0,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrout
              xntimes = 1./float(nouttime(ihr)*nhrout)
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( o3davg_p(i,j,k,no,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = tmp3d_p(i,j,k)                   &
                                     & + o3davg_p(i,j,k,no,ihr)*xntimes
                      o3davg_p(i,j,k,no,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                stop 999
              end if
              misdat = xminout(no)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idout,vnamout(no),tmp3d_p,jxm2,iym2,npl,    &
                       & iadm,xhravg,lnamout(no),uout(no),factout(no),  &
                       & offsetout(no),vvarmin,vvarmax,xlat1d,xlon1d,   &
                       & plev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
 
!     **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do no = 1 , nout2d
        nno = nout3d + no
        if ( u_out(nno)==1 ) then
!         print*,vnamout(nno)
          call setconst(tmp2d,0.0,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
          do ihr = 1 , nhrout
            xntimes = 1./float(nouttime(ihr)*nhrout)
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( o2davg(i,j,no,ihr)>vmisdat ) then
                  tmp2d(i,j) = tmp2d(i,j) + o2davg(i,j,no,ihr)*xntimes
                  o2davg(i,j,no,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xminout(nno) .or. vmax>xmaxout(nno) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamout(nno)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(nno)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(nno)
              stop 999
            end if
            misdat = xminout(no)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idout,vnamout(nno),tmp2d,jxm2,iym2,1,iadm,    &
                        & xhravg,lnamout(nno),uout(nno),factout(nno),   &
                        & offsetout(nno),vvarmin,vvarmax,xlat1d,xlon1d, &
                        & sighrev,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgout

      subroutine writediurout(sighrev,vvarmin,vvarmax,xlat1d,xlon1d,    &
                            & iadm,ndim,xhr1,idout,vmisdat,iotyp,iunt,  &
                            & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev , dtout
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      logical :: plv
      real(4) :: vmisdat
      real(8) :: xhr1
      integer , dimension(ndim) :: iadm
      real(4) , dimension(kz) :: sighrev
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim , plv , xhr1
!
! Local variables
!
      integer :: i , ihr , j , k , nno , no
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
      print * , 'COMPUTING AVERAGE FIELDS FOR DIURNAL OUTPUT:' ,        &
          & nouttime
 
!     **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,vnamout(no)
            do ihr = 1 , nhrout
              xhravg = xhr1 + float(ihr-1)*dtout
              xntimes = 1./float(nouttime(ihr))
!             xhravg = float(ihr-1)*dtout
              if ( nouttime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nouttime = 0'
                stop 999
              end if
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( o3davg(i,j,k,no,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = o3davg(i,j,k,no,ihr)*xntimes
                      o3davg(i,j,k,no,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
                if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                  stop 999
                end if
                misdat = xminout(no)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idout,vnamout(no),tmp3d,jxm2,iym2,kz,iadm,&
                            & xhravg,lnamout(no),uout(no),factout(no),  &
                            & offsetout(no),vvarmin,vvarmax,xlat1d,     &
                            & xlon1d,sighrev,0,misdat,iotyp)
              else if ( iotyp==3 ) then
                call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
              else
              end if
            end do
          end if
        end do
      else
        iadm(3) = npl
        call setconst(tmp3d_p,vmisdat,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
        do no = 1 , nout3d
          if ( u_out(no)==1 ) then
!           print*,vnamout(no)
            do ihr = 1 , nhrout
              xhravg = xhr1 + float(ihr-1)*dtout
              xntimes = 1./float(nouttime(ihr))
!             xhravg = float(ihr-1)*dtout
              if ( nouttime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nouttime = 0'
                stop 999
              end if
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( o3davg_p(i,j,k,no,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = o3davg_p(i,j,k,no,ihr)*xntimes
                      o3davg_p(i,j,k,no,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
                if ( vmin<xminout(no) .or. vmax>xmaxout(no) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamout(no)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(no)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(no)
                  stop 999
                end if
                misdat = xminout(no)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idout,vnamout(no),tmp3d_p,jxm2,iym2,npl,  &
                            & iadm,xhravg,lnamout(no),uout(no),         &
                            & factout(no),offsetout(no),vvarmin,vvarmax,&
                            & xlat1d,xlon1d,plev,0,misdat,iotyp)
              else if ( iotyp==3 ) then
                call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
              else
              end if
            end do
          end if
        end do
      end if
 
!     **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do no = 1 , nout2d
        nno = nout3d + no
        if ( u_out(nno)==1 ) then
!         print*,vnamout(nno)
          do ihr = 1 , nhrout
            xhravg = xhr1 + float(ihr-1)*dtout
            xntimes = 1./float(nouttime(ihr))
            print * , 'nouttime(ihr)=' , nouttime(ihr) , 'xntimes=' ,   &
                & xntimes , 'ihr=' , ihr , xhravg
            if ( nouttime(ihr)<=0 ) then
              print * , 'NOTHING TO AVERAGE -- nouttime = 0'
              stop 999
            end if
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( o2davg(i,j,no,ihr)>vmisdat ) then
                  tmp2d(i,j) = o2davg(i,j,no,ihr)*xntimes
                  o2davg(i,j,no,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
              if ( vmin<xminout(nno) .or. vmax>xmaxout(nno) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamout(nno)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminout(nno)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxout(nno)
                stop 999
              end if
              misdat = xminout(no)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idout,vnamout(nno),tmp2d,jxm2,iym2,1,iadm,  &
                          & xhravg,lnamout(nno),uout(nno),factout(nno), &
                          & offsetout(nno),vvarmin,vvarmax,xlat1d,      &
                          & xlon1d,sighrev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
            else
            end if
          end do
        end if
      end do
 
      end subroutine writediurout

      end module mod_outflds
