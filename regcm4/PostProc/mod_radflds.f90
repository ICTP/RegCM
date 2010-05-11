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

      module mod_radflds

      use mod_dynparam
      use mod_postproc_param , only : nhrrad , npl

      implicit none

      integer , parameter :: nr2d = 9
      integer , parameter :: nr3d = 4
      integer , parameter :: nrtot = nr3d + nr2d

      real(4) , allocatable , dimension(:,:,:,:) :: r2davg
      real(4) , allocatable , dimension(:,:,:,:,:) :: r3davg
      real(4) , allocatable , dimension(:,:,:) :: rfld2d
      real(4) , allocatable , dimension(:,:,:,:) :: rfld3d

      real(4) , allocatable , dimension(:,:,:,:,:) :: r3davg_p
      real(4) , allocatable , dimension(:,:,:,:) :: rfld3d_p

      real(4) , allocatable , dimension(:) :: factrad , offsetrad ,     &
                                  &           xmaxrad , xminrad

      character(64) , allocatable , dimension(:) :: vnamrad
      character(64) , allocatable , dimension(:) :: lnamrad
      character(64) , allocatable , dimension(:) :: urad

      integer , allocatable , dimension(:) :: u_rad
      integer , allocatable , dimension(:) :: nradtime

      contains

      subroutine init_mod_radflds
      implicit none
      allocate(r2davg(jxm2,iym2,nr2d,nhrrad))
      allocate(r3davg(jxm2,iym2,kz,nr3d,nhrrad))
      allocate(rfld2d(jxm2,iym2,nr2d))
      allocate(rfld3d(jxm2,iym2,kz,nr3d))
      allocate(r3davg_p(jxm2,iym2,npl,nr3d,nhrrad))
      allocate(rfld3d_p(jxm2,iym2,npl,nr3d))
      allocate(factrad(nrtot))
      allocate(offsetrad(nrtot))
      allocate(xmaxrad(nrtot))
      allocate(xminrad(nrtot))
      allocate(vnamrad(nrtot))
      allocate(lnamrad(nrtot))
      allocate(urad(nrtot))
      allocate(u_rad(nrtot))
      allocate(nradtime(nhrrad))
      end subroutine init_mod_radflds

      subroutine rdrad(iin,idate,rrec,idirect,ierr)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idate , idirect , ierr , iin , rrec
      intent (in) idirect , iin
      intent (inout) idate , ierr , rrec
!
! Local variables
!
      integer :: i , j , k , nr
      real(4) , dimension(jxm2,iym2) :: tmp2d
! 
      ierr = 0
      if ( idirect/=1 ) then
        read (iin,iostat=ierr) idate
        if ( ierr/=0 ) return
      end if
      print * , 'READING RADIATION DATA:  ' , idate
      do nr = 1 , nr3d
        do k = 1 , kz
          if ( idirect==1 ) then
            rrec = rrec + 1
            read (iin,rec=rrec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
          if ( ierr/=0 ) return
          do j = 1 , iym2
            do i = 1 , jxm2
              rfld3d(i,j,k,nr) = tmp2d(i,j)
            end do
          end do
        end do
      end do
      do nr = 1 , nr2d
        if ( idirect==1 ) then
          rrec = rrec + 1
          read (iin,rec=rrec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
        if ( ierr/=0 ) return
        do j = 1 , iym2
          do i = 1 , jxm2
            rfld2d(i,j,nr) = tmp2d(i,j)
          end do
        end do
      end do
!     Skip added psa not needed here
      rrec = rrec + 1
!     print*,'DONE READING RADIATION FOR CURRENT TIMESTEP',idate
 
      end subroutine rdrad

      subroutine mmvlurad
 
      use mod_point
      implicit none
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      vnamrad(ncld) = 'FC'
      vnamrad(nclwp) = 'CLWP'
      vnamrad(nqrs) = 'QRS'
      vnamrad(nqrl) = 'QRL'
      vnamrad(nr3d+nfsw) = 'FSW'
      vnamrad(nr3d+nflw) = 'FLW'
      vnamrad(nr3d+nclrst) = 'CLRST'
      vnamrad(nr3d+nclrss) = 'CLRSS'
      vnamrad(nr3d+nclrlt) = 'CLRLT'
      vnamrad(nr3d+nclrls) = 'CLRLS'
      vnamrad(nr3d+nsolin) = 'SOLIN'
      vnamrad(nr3d+nsabtp) = 'SABTP'
      vnamrad(nr3d+nfirtp) = 'FIRTP'
 
      lnamrad(ncld) = 'Cloud Fraction'
      lnamrad(nclwp) = 'Cld Liquid H2O Path'
      lnamrad(nqrs) = 'Solar Heating Rate'
      lnamrad(nqrl) = 'LW Cooling Rate'
      lnamrad(nr3d+nfsw) = 'Surface Abs solar'
      lnamrad(nr3d+nflw) = 'LW Cooling of Surf'
      lnamrad(nr3d+nclrst) = 'Clr Sky Col Abs Sol'
      lnamrad(nr3d+nclrss) = 'Clr Sky Surf Abs Sol'
      lnamrad(nr3d+nclrlt) = 'Clr Sky Net Up Flx'
      lnamrad(nr3d+nclrls) = 'Clr Sky LW Surf Cool'
      lnamrad(nr3d+nsolin) = 'Instant Incid Solar'
      lnamrad(nr3d+nsabtp) = 'Column Abs Solar'
      lnamrad(nr3d+nfirtp) = 'Net Up Flux at Top'
 
      urad(ncld) = 'fraction'
      urad(nclwp) = 'g/m2'
      urad(nqrs) = 'K/s'
      urad(nqrl) = 'K/s'
      urad(nr3d+nfsw) = 'W/m2'
      urad(nr3d+nflw) = 'W/m2'
      urad(nr3d+nclrst) = 'W/m2'
      urad(nr3d+nclrss) = 'W/m2'
      urad(nr3d+nclrlt) = 'W/m2'
      urad(nr3d+nclrls) = 'W/m2'
      urad(nr3d+nsolin) = 'W/m2'
      urad(nr3d+nsabtp) = 'W/m2'
      urad(nr3d+nfirtp) = 'W/m2'
 
      xmaxrad(ncld) = 1.1
      xmaxrad(nclwp) = 5000.0
      xmaxrad(nqrs) = 1.0E-2
      xmaxrad(nqrl) = 1.0E-2
      xmaxrad(nr3d+nfsw) = 1200.0
      xmaxrad(nr3d+nflw) = 500.0
      xmaxrad(nr3d+nclrst) = 1500.0
      xmaxrad(nr3d+nclrss) = 1500.0
      xmaxrad(nr3d+nclrlt) = 1500.0
      xmaxrad(nr3d+nclrls) = 500.0
      xmaxrad(nr3d+nsolin) = 1500.0
      xmaxrad(nr3d+nsabtp) = 1500.0
      xmaxrad(nr3d+nfirtp) = 500.0
 
      xminrad(ncld) = -0.1
      xminrad(nclwp) = -10.0
      xminrad(nqrs) = -1.0E-2
      xminrad(nqrl) = -1.0E-2
      xminrad(nr3d+nfsw) = -10.0
      xminrad(nr3d+nflw) = -100.0
      xminrad(nr3d+nclrst) = -10.0
      xminrad(nr3d+nclrss) = -10.0
      xminrad(nr3d+nclrlt) = -10.0
      xminrad(nr3d+nclrls) = -10.0
      xminrad(nr3d+nsolin) = -10.0
      xminrad(nr3d+nsabtp) = -10.0
      xminrad(nr3d+nfirtp) = -10.0
 
      aaa = 2.**16. - 1.
      do l = 1 , nrtot
        factrad(l) = (xmaxrad(l)-xminrad(l))/aaa
        offsetrad(l) = (xmaxrad(l)+xminrad(l))/2.
      end do
 
      end subroutine mmvlurad

      subroutine writerad(vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,      &
                        & sighrev,vmisdat,idout,xhr,iotyp,iunt,nrec)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr
      integer , dimension(ndim) :: iadm
      real(4) , dimension(kz) :: sighrev
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim
!
! Local variables
!
      integer :: i , j , k , nnr , nr
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
!
!     **** WRITE RAD 3-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = kz
      call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
      do nr = 1 , nr3d
!       print*,nr,vnamrad(nr)
        do k = 1 , kz
          do j = 1 , iym2
            do i = 1 , jxm2
              tmp3d(i,j,k) = rfld3d(i,j,k,nr)
            end do
          end do
        end do
        if ( iotyp==1 ) then
          call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
          if ( vmin<xminrad(nr) .or. vmax>xmaxrad(nr) ) then
            print * , 'Values Out of Range:  FIELD=' , vnamrad(nr)
            print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nr)
            print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nr)
            stop 999
          end if
          misdat = xminrad(nr)
        else if ( iotyp==2 ) then
          misdat = vmisdat
        else
        end if
        if ( iotyp==1 .or. iotyp==2 ) then
          call writecdf(idout,vnamrad(nr),tmp3d,jxm2,iym2,kz,iadm,xhr,  &
                      & lnamrad(nr),urad(nr),factrad(nr),offsetrad(nr), &
                      & vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat, &
                      & iotyp)
        else if ( iotyp==3 ) then
          call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
        else
        end if
      end do
 
!     **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do nr = 1 , nr2d
        nnr = nr + nr3d
!       print*,nr,nnr,vnamrad(nnr)
        do j = 1 , iym2
          do i = 2 , jxm2
            tmp2d(i,j) = rfld2d(i,j,nr)
          end do
        end do
        if ( iotyp==1 ) then
          call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
          if ( vmin<xminrad(nnr) .or. vmax>xmaxrad(nnr) ) then
            print * , 'Values Out of Range:  FIELD=' , vnamrad(nnr)
            print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nnr)
            print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nnr)
            stop 999
          end if
          misdat = xminrad(nr)
        else if ( iotyp==2 ) then
          misdat = vmisdat
        else
        end if
        if ( iotyp==1 .or. iotyp==2 ) then
          call writecdf(idout,vnamrad(nnr),tmp2d,jxm2,iym2,1,iadm,xhr,  &
                      & lnamrad(nnr),urad(nnr),factrad(nnr),            &
                      & offsetrad(nnr),vvarmin,vvarmax,xlat1d,xlon1d,   &
                      & sighrev,0,misdat,iotyp)
        else if ( iotyp==3 ) then
          call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
        else
        end if
      end do
      end subroutine writerad
 
      subroutine writeavgrad(xhr1,sighrev,vvarmin,vvarmax,xlat1d,       &
                           & xlon1d,iadm,ndim,vmisdat,idrad,iotyp,iunt, &
                           & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idrad , iotyp , ndim , nrec , iunt
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
      integer :: i , ihr , j , k , nnr , nr
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
      print * , 'COMPUTING AVERAGE RAD FIELDS:' , nradtime
      xhravg = xhr1
      print * , 'nradtime=' , nradtime
      print * , 'xhravg=' , xhravg
 
!     **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do nr = 1 , nr3d
          if ( u_rad(nr)==1 ) then
!           print*,vnamrad(nr)
            call setconst(tmp3d,0.0,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrrad
              xntimes = 1./float(nradtime(ihr)*nhrrad)
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( r3davg(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = tmp3d(i,j,k) + r3davg(i,j,k,nr,ihr)&
                                   & *xntimes
                      r3davg(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xminrad(nr) .or. vmax>xmaxrad(nr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamrad(nr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nr)
                stop 999
              end if
              misdat = xminrad(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idrad,vnamrad(nr),tmp3d,jxm2,iym2,kz,iadm,  &
                          & xhravg,lnamrad(nr),urad(nr),factrad(nr),    &
                          & offsetrad(nr),vvarmin,vvarmax,xlat1d,xlon1d,&
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
        do nr = 1 , nr3d
          if ( u_rad(nr)==1 ) then
!           print*,vnamrad(nr)
            call setconst(tmp3d_p,0.0,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrrad
              xntimes = 1./float(nradtime(ihr)*nhrrad)
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( r3davg_p(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = tmp3d_p(i,j,k)                   &
                                     & + r3davg_p(i,j,k,nr,ihr)*xntimes
                      r3davg_p(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xminrad(nr) .or. vmax>xmaxrad(nr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamrad(nr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nr)
                stop 999
              end if
              misdat = xminrad(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idrad,vnamrad(nr),tmp3d_p,jxm2,iym2,npl,    &
                         & iadm,xhravg,lnamrad(nr),urad(nr),factrad(nr),&
                         & offsetrad(nr),vvarmin,vvarmax,xlat1d,xlon1d, &
                         & plev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
 
!     **** WRITE RAD AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do nr = 1 , nr2d
        nnr = nr3d + nr
        if ( u_rad(nnr)==1 ) then
!         print*,vnamrad(nnr)
          call setconst(tmp2d,0.0,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
          do ihr = 1 , nhrrad
            xntimes = 1./float(nradtime(ihr)*nhrrad)
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( r2davg(i,j,nr,ihr)>vmisdat ) then
                  tmp2d(i,j) = tmp2d(i,j) + r2davg(i,j,nr,ihr)*xntimes
                  r2davg(i,j,nr,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xminrad(nnr) .or. vmax>xmaxrad(nnr) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamrad(nnr)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nnr)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nnr)
              stop 999
            end if
            misdat = xminrad(nr)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idrad,vnamrad(nnr),tmp2d,jxm2,iym2,1,iadm,    &
                        & xhravg,lnamrad(nnr),urad(nnr),factrad(nnr),   &
                        & offsetrad(nnr),vvarmin,vvarmax,xlat1d,xlon1d, &
                        & sighrev,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgrad
 
      subroutine writediurrad(xhr1,sighrev,vvarmin,vvarmax,xlat1d,      &
                            & xlon1d,iadm,ndim,vmisdat,idrad,iotyp,iunt,&
                            & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev , dtrad
      implicit none
!
! Dummy arguments
!
      integer :: idrad , iotyp , ndim , nrec , iunt
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
      integer :: i , ihr , j , k , nnr , nr
      real(4) :: misdat , vmax , vmin , xntimes
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
      real(4) , dimension(jxm2,iym2,npl) :: tmp3d_p
      real(8) :: xhravg
!
!     **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do nr = 1 , nr3d
          if ( u_rad(nr)==1 ) then
!           print*,vnamrad(nr)
            do ihr = 1 , nhrrad
              xhravg = xhr1 + float(ihr-1)*dtrad
              xntimes = 1./float(nradtime(ihr))
!             xhravg = float(ihr-1)*dtrad
              if ( nradtime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nradtime = 0'
                stop 999
              end if
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( r3davg(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = r3davg(i,j,k,nr,ihr)*xntimes
                      r3davg(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
                if ( vmin<xminrad(nr) .or. vmax>xmaxrad(nr) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamrad(nr)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nr)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nr)
                  stop 999
                end if
                misdat = xminrad(nr)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idrad,vnamrad(nr),tmp3d,jxm2,iym2,kz,iadm,&
                            & xhravg,lnamrad(nr),urad(nr),factrad(nr),  &
                            & offsetrad(nr),vvarmin,vvarmax,xlat1d,     &
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
        do nr = 1 , nr3d
          if ( u_rad(nr)==1 ) then
!           print*,vnamrad(nr)
            do ihr = 1 , nhrrad
              xhravg = xhr1 + float(ihr-1)*dtrad
              xntimes = 1./float(nradtime(ihr))
!             xhravg = float(ihr-1)*dtrad
              if ( nradtime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nradtime = 0'
                stop 999
              end if
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( r3davg(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = r3davg_p(i,j,k,nr,ihr)*xntimes
                      r3davg_p(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
                if ( vmin<xminrad(nr) .or. vmax>xmaxrad(nr) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamrad(nr)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nr)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nr)
                  stop 999
                end if
                misdat = xminrad(nr)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idrad,vnamrad(nr),tmp3d_p,jxm2,iym2,npl,  &
                            & iadm,xhravg,lnamrad(nr),urad(nr),         &
                            & factrad(nr),offsetrad(nr),vvarmin,vvarmax,&
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
      do nr = 1 , nr2d
        nnr = nr3d + nr
        if ( u_rad(nnr)==1 ) then
!         print*,vnamrad(nnr)
          do ihr = 1 , nhrrad
            xhravg = xhr1 + float(ihr-1)*dtrad
            xntimes = 1./float(nradtime(ihr))
            if ( nradtime(ihr)<=0 ) then
              print * , 'NOTHING TO AVERAGE -- nradtime = 0'
              stop 999
            end if
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( r2davg(i,j,nr,ihr)>vmisdat ) then
                  tmp2d(i,j) = r2davg(i,j,nr,ihr)*xntimes
                  r2davg(i,j,nr,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
              if ( vmin<xminrad(nnr) .or. vmax>xmaxrad(nnr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamrad(nnr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminrad(nnr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxrad(nnr)
                stop 999
              end if
              misdat = xminrad(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idrad,vnamrad(nnr),tmp2d,jxm2,iym2,1,iadm,  &
                          & xhravg,lnamrad(nnr),urad(nnr),factrad(nnr), &
                          & offsetrad(nnr),vvarmin,vvarmax,xlat1d,      &
                          & xlon1d,sighrev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
            else
            end if
          end do
        end if
      end do
      end subroutine writediurrad

      end module mod_radflds
