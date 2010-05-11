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

      module mod_cheflds

      use mod_dynparam
      use mod_postproc_param , only : nhrche , npl

      implicit none

      integer :: nc3d
      integer :: nc2d
      integer :: nctot

      real(4) , allocatable , dimension(:,:,:,:) :: c2davg
      real(4) , allocatable , dimension(:,:,:,:,:) :: c3davg
      real(4) , allocatable , dimension(:,:,:) :: cfld2d
      real(4) , allocatable , dimension(:,:,:,:) :: cfld3d

      real(4) , allocatable , dimension(:,:,:,:,:) :: c3davg_p
      real(4) , allocatable , dimension(:,:,:,:) :: cfld3d_p

      real(4) , allocatable , dimension(:) :: factche , offsetche ,     &
                               &              xmaxche , xminche

      character(64) , allocatable , dimension(:) :: vnamche
      character(64) , allocatable , dimension(:) :: lnamche
      character(64) , allocatable , dimension(:) :: uche

      integer , allocatable , dimension(:) :: nchetime
      integer , allocatable , dimension(:) :: u_che

      contains

      subroutine init_mod_cheflds
      implicit none
      nc3d  = ntr*1 + 3
      nc2d  = ntr*7 + 2
      nctot = nc3d + nc2d
      allocate(c2davg(jxm2,iym2,nc2d,nhrche))
      allocate(c3davg(jxm2,iym2,kz,nc3d,nhrche))
      allocate(cfld2d(jxm2,iym2,nc2d))
      allocate(cfld3d(jxm2,iym2,kz,nc3d))
      allocate(c3davg_p(jxm2,iym2,npl,nc3d,nhrche))
      allocate(cfld3d_p(jxm2,iym2,npl,nc3d))
      allocate(factche(nctot))
      allocate(offsetche(nctot))
      allocate(xmaxche(nctot))
      allocate(xminche(nctot))
      allocate(vnamche(nctot))
      allocate(lnamche(nctot))
      allocate(uche(nctot))
      allocate(u_che(nctot))
      allocate(nchetime(nhrche))
      end subroutine init_mod_cheflds

      subroutine rdche(iin,idate,crec,idirect,ierr)
      implicit none
!
! Dummy arguments
!
      integer :: crec , idate , idirect , ierr , iin
      intent (in) idirect , iin
      intent (inout) crec , idate , ierr
!
! Local variables
!
      integer :: i , j , k , nc
      real(4) , dimension(jxm2,iym2) :: tmp2d
!
      if ( idirect/=1 ) read (iin,iostat=ierr) idate
      if ( ierr/=0 ) return
      print * , 'READING CHEM-TRACER DATA:  ' , idate
 
      do nc = 1 , nc3d
        do k = 1 , kz
          if ( idirect==1 ) then
            crec = crec + 1
            read (iin,rec=crec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
          if ( ierr/=0 ) return
          do j = 1 , iym2
            do i = 1 , jxm2
              cfld3d(i,j,k,nc) = tmp2d(i,j)
            end do
          end do
        end do
      end do
 
      do nc = 1 , nc2d
        if ( idirect==1 ) then
          crec = crec + 1
          read (iin,rec=crec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
        if ( ierr/=0 ) return
        do j = 1 , iym2
          do i = 1 , jxm2
            cfld2d(i,j,nc) = tmp2d(i,j)
          end do
        end do
      end do
!     Skip added psa not needed here
      crec = crec + 1
      end subroutine rdche
 
      subroutine mmvluche
 
      implicit none
!
! Local variables
!
      real(4) :: aaa
      integer :: l , n , r
      character(10) , dimension(ntr) :: tracname
!
      tracname(1) = 'TR1'
      tracname(2) = 'TR2'
      tracname(3) = 'TR3'
      tracname(4) = 'TR4'
      tracname(5) = 'TR5'
      tracname(6) = 'TR6'
      tracname(7) = 'TR7'
      tracname(8) = 'TR8'
      tracname(9) = 'TR9'
      tracname(10) = 'TR10'
 
      r = (nc3d-3)/ntr
      do n = 1 , ntr
        vnamche(r*(n-1)+1) = tracname(n)
        lnamche(r*(n-1)+1) = 'MMR_'//tracname(n)
        uche(r*(n-1)+1) = 'micro-g/Kg'
      end do
      print * , 'r is ' , r
      lnamche(7) = 'aer mix. ext. coef'
      vnamche(7) = 'aext8'
      uche(7) = 'na'
      xmaxche(7) = 0.5
      xminche(7) = 0.
 
      lnamche(8) = 'aer mix. scat. alb'
      vnamche(8) = 'assa9'
      uche(8) = 'na'
      xmaxche(8) = 0.5
      xminche(8) = 0.
 
      lnamche(9) = 'aer mix. scat. alb'
      vnamche(9) = 'agfu8'
      uche(9) = 'na'
      xmaxche(9) = 0.5
      xminche(9) = 0.
 
      r = (nc2d-2)/ntr
 
      do n = 1 , ntr
        vnamche(nc3d+r*(n-1)+1) = 'BURD'//tracname(n)
        vnamche(nc3d+r*(n-1)+2) = 'WDLS'//tracname(n)
        vnamche(nc3d+r*(n-1)+3) = 'WDCV'//tracname(n)
        vnamche(nc3d+r*(n-1)+4) = 'DRDP'//tracname(n)
        vnamche(nc3d+r*(n-1)+5) = 'XGAZ'//tracname(n)
        vnamche(nc3d+r*(n-1)+6) = 'XSAQ'//tracname(n)
        vnamche(nc3d+r*(n-1)+7) = 'EMRA'//tracname(n)
 
        lnamche(nc3d+r*(n-1)+1) = 'Int column '//tracname(n)
        lnamche(nc3d+r*(n-1)+2) = 'Wetdep lsc '//tracname(n)
        lnamche(nc3d+r*(n-1)+3) = 'Wetdep cvc '//tracname(n)
        lnamche(nc3d+r*(n-1)+4) = 'Drydep surf'//tracname(n)
        lnamche(nc3d+r*(n-1)+5) = 'gaz conv'//tracname(n)
        lnamche(nc3d+r*(n-1)+6) = 'aq conv'//tracname(n)
        lnamche(nc3d+r*(n-1)+7) = 'Emission rate'//tracname(n)
 
        uche(nc3d+r*(n-1)+1) = 'mg/m2'
        uche(nc3d+r*(n-1)+2) = 'mg/m2'
        uche(nc3d+r*(n-1)+3) = 'mg/m2'
        uche(nc3d+r*(n-1)+4) = 'mg/m2'
        uche(nc3d+r*(n-1)+5) = 'mg/m2'
        uche(nc3d+r*(n-1)+6) = 'mg/m2'
        uche(nc3d+r*(n-1)+7) = 'micro-g/m2.s'
 
 
!       print*,'1 :',vnamche(1)
!       print*,'2 :',vnamche(2)
!       print*,'3 :',vnamche(3)
!       print*,'4 :',vnamche(4)
!       print*,'5 :',vnamche(5)
!       print*,'6 :',vnamche(6)
!       print*,'7 :',vnamche(7)
!       print*,'8 :',vnamche(8)
 
      end do
 
      do n = 1 , nc3d - 3
        xmaxche(n) = 10.E6
        xminche(n) = 0.
      end do
 
!     xmaxche(5) = 10.E9
!     xmaxche(6) = 10.E9
      do n = nc3d + 1 , nc3d + nc2d - 2
        xmaxche(n) = 10.E10
        xminche(n) = 0.
      end do
 
!     xmaxche(ncld)        = 1.1
!     xmaxche(nclwp)       = 5000.0
!     xmaxche(nqrs)        = 1.0e-3
!     xmaxche(nqrl)        = 1.0e-3
!     xmaxche(nr3d+nfsw)   = 1200.0
!     xmaxche(nr3d+nflw)   = 500.0
!     xmaxche(nr3d+nclrst) = 1500.0
!     xmaxche(nr3d+nclrss) = 1500.0
!     xmaxche(nr3d+nclrlt) = 1500.0
!     xmaxche(nr3d+nclrls) = 500.0
!     xmaxche(nr3d+nsolin) = 1500.0
!     xmaxche(nr3d+nsabtp) = 1500.0
!     xmaxche(nr3d+nfirtp) = 500.0
 
!     xminche(ncld)        = -0.1
!     xminche(nclwp)       = -10.0
!     xminche(nqrs)        = -1.0e-3
!     xminche(nqrl)        = -1.0e-3
!     xminche(nr3d+nfsw)   = -10.0
!     xminche(nr3d+nflw)   = -100.0
!     xminche(nr3d+nclrst) = -10.0
!     xminche(nr3d+nclrss) = -10.0
!     xminche(nr3d+nclrlt) = -10.0
!     xminche(nr3d+nclrls) = -10.0
!     xminche(nr3d+nsolin) = -10.0
!     xminche(nr3d+nsabtp) = -10.0
!     xminche(nr3d+nfirtp) = -10.0
 
      lnamche(52) = 'TOArad forcing'
      vnamche(52) = 'acsto'
      uche(52) = 'W/m^2'
      xmaxche(52) = 200.0
      xminche(52) = -200.0
 
      lnamche(53) = 'TOArad forcing'
      vnamche(53) = 'acsto'
      uche(53) = 'W/m^2'
      xmaxche(53) = 200.0
      xminche(53) = -200.0
 
      aaa = 2.**16. - 1.
      do l = 1 , nctot
        factche(l) = (xmaxche(l)-xminche(l))/aaa
        offsetche(l) = (xmaxche(l)+xminche(l))/2.
      end do
 
      end subroutine mmvluche

      subroutine writeche(vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,      &
                        & sighrev,vmisdat,idout,xhr,iotyp,iunt,nrec)
 
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
      integer :: i , j , k , nc , nnc
      real(4) :: misdat , vmax , vmin
      real(4) , dimension(jxm2,iym2) :: tmp2d
      real(4) , dimension(jxm2,iym2,kz) :: tmp3d
!
!     **** WRITE RAD 3-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = kz
      call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
 
      do nc = 1 , nc3d
        if ( u_che(nc)==1 ) then
!         print*,nr,vnamrad(nr)
          do k = 1 , kz
            do j = 1 , iym2
              do i = 1 , jxm2
                tmp3d(i,j,k) = cfld3d(i,j,k,nc)*1.E9
              end do
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
            if ( vmin<xminche(nc) .or. vmax>xmaxche(nc) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamche(nc)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nc)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nc)
              stop 999
            end if
            misdat = xminche(nc)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
 
            call writecdf(idout,vnamche(nc),tmp3d,jxm2,iym2,kz,iadm,xhr,&
                        & lnamche(nc),uche(nc),factche(nc),             &
                        & offsetche(nc),vvarmin,vvarmax,xlat1d,xlon1d,  &
                        & sighrev,0,misdat,iotyp)
 
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp3d,jxm2,iym2,kz,nrec)
          else
          end if
        end if
      end do
 
!     **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      call setconst(tmp2d,vmisdat,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
      do nc = 1 , nc2d
        nnc = nc + nc3d
        if ( u_che(nnc)==1 ) then
!         print*,nr,nnr,vnamrad(nnr)
          do j = 1 , iym2
            do i = 2 , jxm2
              tmp2d(i,j) = cfld2d(i,j,nc)
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xminche(nnc) .or. vmax>xmaxche(nnc) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamche(nnc)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nnc)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nnc)
              stop 999
            end if
            misdat = xminche(nc)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
 
            call writecdf(idout,vnamche(nnc),tmp2d,jxm2,iym2,1,iadm,xhr,&
                        & lnamche(nnc),uche(nnc),factche(nnc),          &
                        & offsetche(nnc),vvarmin,vvarmax,xlat1d,xlon1d, &
                        & sighrev,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeche

      subroutine writeavgche(xhr1,sighrev,vvarmin,vvarmax,xlat1d,       &
                           & xlon1d,iadm,ndim,vmisdat,idche,iotyp,iunt, &
                           & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : plev
      implicit none
!
! Dummy arguments
!
      integer :: idche , iotyp , ndim , nrec , iunt
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
      print * , 'COMPUTING AVERAGE CHE FIELDS:' , nchetime
      xhravg = xhr1
      print * , 'nchetime=' , nchetime
      print * , 'xhravg=' , xhravg
 
!     **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if ( .not.plv ) then
        iadm(3) = kz
        call setconst(tmp3d,vmisdat,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
        do nr = 1 , nc3d
          if ( u_che(nr)==1 ) then
!           print*,vnamrad(nr)
            call setconst(tmp3d,0.0,jxm2,iym2,kz,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrche
              xntimes = 1./float(nchetime(ihr)*nhrche)
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( c3davg(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = tmp3d(i,j,k) + c3davg(i,j,k,nr,ihr)&
                                   & *xntimes
                      c3davg(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
              if ( vmin<xminche(nr) .or. vmax>xmaxche(nr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamche(nr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nr)
                stop 999
              end if
              misdat = xminche(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idche,vnamche(nr),tmp3d,jxm2,iym2,kz,iadm,  &
                          & xhravg,lnamche(nr),uche(nr),factche(nr),    &
                          & offsetche(nr),vvarmin,vvarmax,xlat1d,xlon1d,&
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
        do nr = 1 , nc3d
          if ( u_che(nr)==1 ) then
!           print*,vnamrad(nr)
            call setconst(tmp3d_p,0.0,jxm2,iym2,npl,1,1,1,jxm2,1,iym2)
            do ihr = 1 , nhrche
              xntimes = 1./float(nchetime(ihr)*nhrche)
              do k = 1 , npl
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( c3davg_p(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = tmp3d_p(i,j,k)                   &
                                     & + c3davg_p(i,j,k,nr,ihr)*xntimes
                      c3davg_p(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
              if ( vmin<xminche(nr) .or. vmax>xmaxche(nr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamche(nr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nr)
                stop 999
              end if
              misdat = xminche(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idche,vnamche(nr),tmp3d_p,jxm2,iym2,npl,    &
                          & iadm,xhravg,lnamche(nr),uche(nr),           &
                          & factche(nr),offsetche(nr),vvarmin,vvarmax,  &
                          & xlat1d,xlon1d,plev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp3d_p,jxm2,iym2,npl,nrec)
            else
            end if
          end if
        end do
      end if
      print * , 'repere1'
!     **** WRITE RAD AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      iadm(3) = 1
      do nr = 1 , nc2d
        nnr = nc3d + nr
        if ( u_che(nnr)==1 ) then
!         print*,vnamrad(nnr)
          call setconst(tmp2d,0.0,jxm2,iym2,1,1,1,1,jxm2,1,iym2)
          do ihr = 1 , nhrche
            xntimes = 1./float(nchetime(ihr)*nhrche)
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( c2davg(i,j,nr,ihr)>vmisdat ) then
                  tmp2d(i,j) = tmp2d(i,j) + c2davg(i,j,nr,ihr)*xntimes
                  c2davg(i,j,nr,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
          end do
          if ( iotyp==1 ) then
            call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
            if ( vmin<xminche(nnr) .or. vmax>xmaxche(nnr) ) then
              print * , 'Values Out of Range:  FIELD=' , vnamche(nnr)
              print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nnr)
              print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nnr)
              stop 999
            end if
            misdat = xminche(nr)
          else if ( iotyp==2 ) then
            misdat = vmisdat
          else
          end if
          if ( iotyp==1 .or. iotyp==2 ) then
            call writecdf(idche,vnamche(nnr),tmp2d,jxm2,iym2,1,iadm,    &
                        & xhravg,lnamche(nnr),uche(nnr),factche(nnr),   &
                        & offsetche(nnr),vvarmin,vvarmax,xlat1d,xlon1d, &
                        & sighrev,0,misdat,iotyp)
          else if ( iotyp==3 ) then
            call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
          else
          end if
        end if
      end do
      end subroutine writeavgche
 
      subroutine writediurche(xhr1,sighrev,vvarmin,vvarmax,xlat1d,      &
                            & xlon1d,iadm,ndim,vmisdat,idche,iotyp,iunt,&
                            & nrec,plv)
 
      use mod_point
      use mod_postproc_param , only : dtche , plev
      implicit none
!
! Dummy arguments
!
      integer :: idche , iotyp , ndim , nrec , iunt
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
        do nr = 1 , nc3d
          if ( u_che(nr)==1 ) then
            print * , vnamche(nr)
            do ihr = 1 , nhrche
              xhravg = xhr1 + float(ihr-1)*dtche
              xntimes = 1./float(nchetime(ihr))
!             xhravg = float(ihr-1)*dtche
              if ( nchetime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nchetime = 0'
                stop 999
              end if
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( c3davg(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d(i,j,k) = c3davg(i,j,k,nr,ihr)*xntimes
                      c3davg(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d,jxm2,iym2,kz,vmin,vmax,vmisdat)
                if ( vmin<xminche(nr) .or. vmax>xmaxche(nr) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamche(nr)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nr)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nr)
                  stop 999
                end if
                misdat = xminche(nr)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idche,vnamche(nr),tmp3d,jxm2,iym2,kz,iadm,&
                            & xhravg,lnamche(nr),uche(nr),factche(nr),  &
                            & offsetche(nr),vvarmin,vvarmax,xlat1d,     &
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
        do nr = 1 , nc3d
          if ( u_che(nr)==1 ) then
            print * , vnamche(nr)
            do ihr = 1 , nhrche
              xhravg = xhr1 + float(ihr-1)*dtche
              xntimes = 1./float(nchetime(ihr))
!             xhravg = float(ihr-1)*dtche
              if ( nchetime(ihr)<=0 ) then
                print * , 'NOTHING TO AVERAGE -- nchetime = 0'
                stop 999
              end if
              do k = 1 , kz
                do j = 1 , iym2
                  do i = 1 , jxm2
                    if ( c3davg_p(i,j,k,nr,ihr)>vmisdat ) then
                      tmp3d_p(i,j,k) = c3davg_p(i,j,k,nr,ihr)*xntimes
                      c3davg_p(i,j,k,nr,ihr) = 0.0
                    else
                      tmp3d_p(i,j,k) = vmisdat
                    end if
                  end do
                end do
              end do
              if ( iotyp==1 ) then
                call getminmax(tmp3d_p,jxm2,iym2,npl,vmin,vmax,vmisdat)
                if ( vmin<xminche(nr) .or. vmax>xmaxche(nr) ) then
                  print * , 'Values Out of Range:  FIELD=' , vnamche(nr)
                  print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nr)
                  print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nr)
                  stop 999
                end if
                misdat = xminche(nr)
              else if ( iotyp==2 ) then
                misdat = vmisdat
              else
              end if
              if ( iotyp==1 .or. iotyp==2 ) then
                call writecdf(idche,vnamche(nr),tmp3d_p,jxm2,iym2,npl,  &
                            & iadm,xhravg,lnamche(nr),uche(nr),         &
                            & factche(nr),offsetche(nr),vvarmin,vvarmax,&
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
      do nr = 1 , nc2d
        nnr = nc3d + nr
        if ( u_che(nnr)==1 ) then
!         print*,vnamche(nnr)
          do ihr = 1 , nhrche
            xhravg = xhr1 + float(ihr-1)*dtche
            xntimes = 1./float(nchetime(ihr))
            if ( nchetime(ihr)<=0 ) then
              print * , 'NOTHING TO AVERAGE -- nchetime = 0'
              stop 999
            end if
            do j = 1 , iym2
              do i = 1 , jxm2
                if ( c2davg(i,j,nr,ihr)>vmisdat ) then
                  tmp2d(i,j) = c2davg(i,j,nr,ihr)*xntimes
                  c2davg(i,j,nr,ihr) = 0.0
                else
                  tmp2d(i,j) = vmisdat
                end if
              end do
            end do
            if ( iotyp==1 ) then
              call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
              if ( vmin<xminche(nnr) .or. vmax>xmaxche(nnr) ) then
                print * , 'Values Out of Range:  FIELD=' , vnamche(nnr)
                print * , 'MINVAL=' , vmin , 'XMIN=' , xminche(nnr)
                print * , 'MAXVAL=' , vmax , 'XMAX=' , xmaxche(nnr)
                stop 999
              end if
              misdat = xminche(nr)
            else if ( iotyp==2 ) then
              misdat = vmisdat
            else
            end if
            if ( iotyp==1 .or. iotyp==2 ) then
              call writecdf(idche,vnamche(nnr),tmp2d,jxm2,iym2,1,iadm,  &
                          & xhravg,lnamche(nnr),uche(nnr),factche(nnr), &
                          & offsetche(nnr),vvarmin,vvarmax,xlat1d,      &
                          & xlon1d,sighrev,0,misdat,iotyp)
            else if ( iotyp==3 ) then
              call writegrads(iunt,tmp2d,jxm2,iym2,1,nrec)
            else
            end if
          end do
        end if
      end do
      end subroutine writediurche

      end module mod_cheflds
