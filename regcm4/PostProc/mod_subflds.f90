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

      module mod_subflds
      use mod_regcm_param , only : jxm2 , ixm2 , nsg , numsub
      use mod_postproc_param , only : nhrsub
      implicit none

      integer , parameter :: jxsg = jxm2*nsg
      integer , parameter :: ixsg = ixm2*nsg
      integer , parameter :: nsub2 = numsub+1

      real(4) , dimension(jxsg,ixsg,nsub2,nhrsub) :: s2davg
      real(4) , dimension(jxsg,ixsg,nsub2) :: sfld2d

      contains

      subroutine rdsub(idate,iin,srec,idirect,ierr)
 
      implicit none
!
! Dummy arguments
!
      integer :: idate , idirect , ierr , iin , srec
      intent (in) idirect , iin
      intent (inout) idate , ierr , srec
!
! Local variables
!
      integer :: i , j , ns
      real(4) , dimension(jxsg,ixsg) :: s2d
!
      ierr = 0
      if ( idirect/=1 ) then
        read (iin,iostat=ierr) idate
        if ( ierr/=0 ) return
      end if
      print * , ' READING SUB:' , idate
      do ns = 1 , numsub
        if ( idirect==1 ) then
          srec = srec + 1
          read (iin,rec=srec,iostat=ierr) s2d
        else
          read (iin,iostat=ierr) s2d
        end if
        do j = 1 , ixsg
          do i = 1 , jxsg
            sfld2d(i,j,ns) = s2d(i,j)
          end do
        end do
      end do
      end subroutine rdsub

      subroutine avgdatasub(ihr,vmisdat)
 
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
      do nb = 1 , nsub2
!       if ((nb.eq.nstmin.or.nb.eq.nstmax) .and. ihr.lt.nhrsub) then
!       else
        do j = 1 , ixsg
          do i = 1 , jxsg
            if ( sfld2d(i,j,nb)>misdat ) then
              s2davg(i,j,nb,ihr) = s2davg(i,j,nb,ihr) + sfld2d(i,j,nb)
            else
              s2davg(i,j,nb,ihr) = vmisdat
            end if
          end do
        end do
!       end if
      end do
 
      end subroutine avgdatasub
 
      subroutine mmvlusub(vnamsub,lnamsub,usub,xmin,xmax,fact,offset)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      real(4) , dimension(nsub2) :: fact , offset , xmax , xmin
      character(64) , dimension(nsub2) :: lnamsub
      character(64) , dimension(nsub2) :: usub
      character(64) , dimension(nsub2) :: vnamsub
      intent (out) fact , lnamsub , offset , usub , vnamsub
      intent (inout) xmax , xmin
!
! Local variables
!
      real(4) :: aaa
      integer :: l
!
      lnamsub(nsux) = 'Anemom Zonal Winds'
      vnamsub(nsux) = 'UA'
      usub(nsux) = 'm/s'
      xmax(nsux) = 50.0
      xmin(nsux) = -50.0
 
      lnamsub(nsvx) = 'Anemom Merid Winds'
      vnamsub(nsvx) = 'VA'
      usub(nsvx) = 'm/s'
      xmax(nsvx) = 50.0
      xmin(nsvx) = -50.0
 
      lnamsub(nsdrag) = 'Surface Drag Stress'
      vnamsub(nsdrag) = 'DRAG'
      usub(nsdrag) = 'si'
      xmax(nsdrag) = 1.0
      xmin(nsdrag) = -1.0
 
      vnamsub(nstg) = 'TG'
      lnamsub(nstg) = 'Ground Temperature'
      usub(nstg) = 'K'
      xmax(nstg) = 350.0
      xmin(nstg) = 180.0
 
      vnamsub(nstf) = 'TF'
      lnamsub(nstf) = 'Foliage Temp'
      usub(nstf) = 'K'
      xmax(nstf) = 350.0
      xmin(nstf) = 180.0
 
      lnamsub(nstanm) = 'Anemom Temp'
      vnamsub(nstanm) = 'TA'
      usub(nstanm) = 'K'
      xmax(nstanm) = 350.0
      xmin(nstanm) = 180.0
 
      lnamsub(nsqanm) = 'Anemom Spec Humidity'
      vnamsub(nsqanm) = 'QA'
      usub(nsqanm) = 'kg/kg'
      xmax(nsqanm) = 0.20
      xmin(nsqanm) = -1.0E-5
 
      lnamsub(nssmu) = 'Top Layer Soil Moist'
      vnamsub(nssmu) = 'SMU'
      usub(nssmu) = 'mm'
      xmax(nssmu) = 80.0
      xmin(nssmu) = -1.0
 
      lnamsub(nssmr) = 'Root Lay Soil Moist'
      vnamsub(nssmr) = 'SMR'
      usub(nssmr) = 'mm'
      xmax(nssmr) = 1200.0
      xmin(nssmr) = -1.0
 
      lnamsub(nset) = 'Evapotranspiration'
      vnamsub(nset) = 'ET'
      usub(nset) = 'mm/day'
      xmax(nset) = 150.0
      xmin(nset) = -5.0
 
      lnamsub(nsrnfs) = 'Surface Runoff'
      vnamsub(nsrnfs) = 'RNFS'
      usub(nsrnfs) = 'mm/day'
      xmax(nsrnfs) = 2000.0
      xmin(nsrnfs) = -200.0
 
      lnamsub(nssnow) = 'Snow Depth'
      vnamsub(nssnow) = 'SNOW'
      usub(nssnow) = 'mm H2O'
      xmax(nssnow) = 1000.0
      xmin(nssnow) = -1.0
 
      lnamsub(nssh) = 'Sensible Heat'
      vnamsub(nssh) = 'SH'
      usub(nssh) = 'W/m2'
      xmax(nssh) = 1000.0
      xmin(nssh) = -300.0
 
      lnamsub(nsprc) = 'Convective Precip'
      vnamsub(nsprc) = 'RC'
      usub(nsprc) = 'mm/day'
      xmax(nsprc) = 1500.0
      xmin(nsprc) = -1.0
 
      lnamsub(nspt) = 'Total Precipitation'
      vnamsub(nspt) = 'RT'
      usub(nspt) = 'mm/day'
      xmax(nspt) = 2500.0
      xmin(nspt) = -1.0
 
      lnamsub(nspsrf) = 'Surface Pressure'
      vnamsub(nspsrf) = 'PSRF'
      usub(nspsrf) = 'hPa'
      xmax(nspsrf) = 1500.0
      xmin(nspsrf) = 300.0
 
      lnamsub(nsrha) = 'Relative Humidity'
      vnamsub(nsrha) = 'RHA'
      usub(nsrha) = 'fraction'
      xmax(nsrha) = 5.0
      xmin(nsrha) = -0.1
 
      aaa = 2.**16. - 1.
      do l = 1 , nsub2
        fact(l) = (xmax(l)-xmin(l))/aaa
        offset(l) = (xmax(l)+xmin(l))/2.
      end do
 
      end subroutine mmvlusub

      subroutine writesub(vnamsub,lnamsub,usub,xmin,xmax,fact,offset,   &
                        & vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,      &
                        & vmisdat,xhr,idout,iotyp,iunt,nrec)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr
      real(4) , dimension(nsub2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nsub2) :: lnamsub
      character(64) , dimension(nsub2) :: usub
      character(64) , dimension(nsub2) :: vnamsub
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(ixsg) :: xlat1d
      real(4) , dimension(jxsg) :: xlon1d
      intent (in) ndim , xmax , xmin
!
! Local variables
!
      integer :: i , j , ns
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxsg,ixsg) :: tmp2d
!
      sig1(1) = 1.0
      call setconst(tmp2d,vmisdat,jxsg,ixsg,1,1,1,1,jxsg,1,ixsg)
      do ns = 1 , nsub2
!       if (ns.ne.nstmin .and. ns.ne.nstmax) then
        do i = 1 , jxsg
          do j = 1 , ixsg
            tmp2d(i,j) = max(sfld2d(i,j,ns),vmisdat)
          end do
        end do
        if ( iotyp==1 ) then
          call getminmax(tmp2d,jxsg,ixsg,1,vmin,vmax,vmisdat)
          if ( vmin<xmin(ns) .or. vmax>xmax(ns) ) then
            print * , 'Values Out of Range:  FIELD=' , vnamsub(ns)
            print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ns)
            print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ns)
!           stop 999
          end if
          misdat = xmin(ns)
        else if ( iotyp==2 .or. iotyp==3 ) then
          misdat = vmisdat
        else
        end if
!       print*,vnamsub(ns),nrec+1
        if ( iotyp==1 .or. iotyp==2 ) then
          call writecdf(idout,vnamsub(ns),tmp2d,jxsg,ixsg,1,iadm,xhr,   &
                      & lnamsub(ns),usub(ns),fact(ns),offset(ns),       &
                      & vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,    &
                      & iotyp)
        else if ( iotyp==3 ) then
          call writegrads(iunt,tmp2d,jxsg,ixsg,1,nrec)
        else
        end if
!       end if
      end do
      end subroutine writesub
 
      subroutine writeavgsub(vmisdat,vnamsub,lnamsub,usub,xmin,xmax,    &
                           & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d, &
                           & iadm,ndim,xhr1,nsubtime,idsub,iotyp,iunt,  &
                           & nrec)
 
      use mod_point
      implicit none
!
! Dummy arguments
!
      integer :: idsub , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nsub2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nsub2) :: lnamsub
      integer , dimension(nhrsub) :: nsubtime
      character(64) , dimension(nsub2) :: usub
      character(64) , dimension(nsub2) :: vnamsub
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(ixsg) :: xlat1d
      real(4) , dimension(jxsg) :: xlon1d
      intent (in) ndim , nsubtime , xhr1 , xmax , xmin
!
! Local variables
!
      real(4) :: const , misdat , vmax , vmin , xntimes
      real(4) , dimension(jxsg,ixsg,nsub2) :: favgsum
      integer :: i , ihr , j , ns
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxsg,ixsg) :: tmp2d
      real(8) :: xhravg
!
      iadm(3) = 1
      sig1(1) = 1.
 
      print * , 'COMPUTING AVERAGE SUB FIELDS:' , nsubtime
 
      call setconst(tmp2d,vmisdat,jxsg,ixsg,1,1,1,1,jxsg,1,ixsg)
      call setconst(favgsum,0.0,jxsg,ixsg,nsub2,1,1,1,jxsg,1,ixsg)
      do ihr = 1 , nhrsub
        do ns = 1 , nsub2
          const = 1.0
!         if (ns.eq.nstmin .or. ns.eq.nstmax) then
!         xntimes = const/float(nsubtime(ihr))
!         else
          xntimes = const/float(nhrsub*nsubtime(ihr))
!         end if
          do j = 1 , ixsg
            do i = 1 , jxsg
              if ( s2davg(i,j,ns,ihr)>vmisdat ) then
                favgsum(i,j,ns) = favgsum(i,j,ns) + s2davg(i,j,ns,ihr)  &
                                & *xntimes
              else
                favgsum(i,j,ns) = vmisdat
              end if
            end do
          end do
        end do
      end do
      xhravg = xhr1
      do ns = 1 , nsub2
        if ( xntimes<=0 ) then
          print * , 'NOTHING TO AVERAGE -- nsubtime = 0'
          stop 999
        end if
        do j = 1 , ixsg
          do i = 1 , jxsg
            tmp2d(i,j) = favgsum(i,j,ns)
          end do
        end do
        if ( iotyp==1 ) then
          call getminmax(tmp2d,jxsg,ixsg,1,vmin,vmax,vmisdat)
          if ( vmin<xmin(ns) .or. vmax>xmax(ns) ) then
            print * , 'Values Out of Range:  FIELD=' , vnamsub(ns)
            print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ns)
            print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ns)
            stop 999
          end if
          misdat = xmin(ns)
        else if ( iotyp==2 ) then
          misdat = vmisdat
        else
        end if
        if ( iotyp==1 .or. iotyp==2 ) then
          call writecdf(idsub,vnamsub(ns),tmp2d,jxsg,ixsg,1,iadm,xhravg,&
                      & lnamsub(ns),usub(ns),fact(ns),offset(ns),       &
                      & vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,    &
                      & iotyp)
        else if ( iotyp==3 ) then
          call writegrads(iunt,tmp2d,jxsg,ixsg,1,nrec)
        else
        end if
      end do
      end subroutine writeavgsub

      subroutine writediursub(vmisdat,vnamsub,lnamsub,usub,xmin,xmax,   &
                            & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,&
                            & iadm,ndim,xhr1,nsubtime,idsub,iotyp,iunt, &
                            & nrec)
 
      use mod_point
      use mod_postproc_param , only : dtsub , nhrbat

      implicit none
!
! Dummy arguments
!
      integer :: idsub , iotyp , ndim , nrec , iunt
      real(4) :: vmisdat
      real(8) :: xhr1
      real(4) , dimension(nsub2) :: fact , offset , xmax , xmin
      integer , dimension(ndim) :: iadm
      character(64) , dimension(nsub2) :: lnamsub
      integer , dimension(nhrbat) :: nsubtime
      character(64) , dimension(nsub2) :: usub
      character(64) , dimension(nsub2) :: vnamsub
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(ixsg) :: xlat1d
      real(4) , dimension(jxsg) :: xlon1d
      intent (in) ndim , nsubtime , xhr1 , xmax , xmin
!
! Local variables
!
      real(4) :: const , misdat , vmax , vmin , xntimes
      integer :: i , ihr , j , ns
      real(4) , dimension(1) :: sig1
      real(4) , dimension(jxsg,ixsg) :: tmp2d
      real(8) :: xhravg
!
      iadm(3) = 1
      sig1(1) = 1.
      call setconst(tmp2d,vmisdat,jxsg,ixsg,1,1,1,1,jxsg,1,ixsg)
      do ihr = 1 , nhrsub
        xhravg = xhr1 + float(ihr-1)*dtsub
        if ( nsubtime(ihr)<=0 ) then
          print * , 'NOTHING TO AVERAGE -- nsubtime = 0'
          stop 999
        end if
        const = 1.0
        do ns = 1 , nsub2
!         if (ns.eq.nstmin .or. ns.eq.nstmax) then
!         xntimes = const/float(nsubtime(ihr))
!         else
          xntimes = const/float(nsubtime(ihr))
!         end if
!         if ((ns.eq.nstmin.or.ns.eq.nstmax) .and. ihr.lt.nhrsub) then
          do j = 1 , ixsg
            do i = 1 , jxsg
              tmp2d(i,j) = vmisdat
            end do
          end do
!         else
          do j = 1 , ixsg
            do i = 1 , jxsg
              if ( s2davg(i,j,ns,ihr)>vmisdat ) then
                tmp2d(i,j) = s2davg(i,j,ns,ihr)*xntimes
              else
                tmp2d(i,j) = vmisdat
              end if
            end do
          end do
!         end if
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxsg,ixsg,1,vmin,vmax,vmisdat)
            if ( vmin<xmin(ns) .or. vmax>xmax(ns) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamsub(ns)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xmin(ns)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmax(ns)
              stop 999
            end if
            misdat = xmin(ns)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idsub,vnamsub(ns),tmp2d,jxsg,ixsg,1,iadm,     &
                        & xhravg,lnamsub(ns),usub(ns),fact(ns),         &
                        & offset(ns),vvarmin,vvarmax,xlat1d,xlon1d,sig1,&
                        & 0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxsg,ixsg,1,nrec)
          else
          end if
        end do
      end do
      end subroutine writediursub
 
      end module mod_subflds
