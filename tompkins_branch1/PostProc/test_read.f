      implicit none
      integer mdate,mdate0,ibltyp,icup,imoist,iboudy,igrads,ibigend
     &      , iin,ni,nj,nk,nx,ny,nz,i,j,k,kk,ierr,idirect,ibyte,ierrsum
      character proj*6, inhead*70
      logical there

      parameter(nx=100,ny=80,nz=23,ibyte=1,inhead='OUT_HEAD_K')
c     parameter(nx=122,ny=113,nz=18,ibyte=1,inhead='OUT_HEAD_J')

      real v2d(nx-2,ny-2), sigf(nz), ds, clat, clon, pt, xplat, xplon
     &   , dto, dtb, dtr, dtc

      iin=10
      print*,'OPENING: ',inhead
      print*,'  nx=',nx,'ny=',ny,'nk=',nk
      print*,'  ibyte=',ibyte
      inquire(file=inhead,exist=there)
      if (there) then
        print*,'INHEAD FILE EXISTS'
      else
        stop 'INHEAD FILE DOES NOT EXIST'
      end if

      open (iin,file=inhead,status='old',form='unformatted'
     &         ,recl=(nx-2)*(ny-2)*ibyte,access='direct')
      print*,'OUTHEAD HAS BEEN OPENED'
      read(iin,rec=1,iostat=ierr) mdate0,ibltyp,icup,imoist,iboudy
     &        ,ni,nj,nk,(sigf(k),k=nz+1,1,-1),ds,pt,clat,clon,xplat
     &        ,xplon,proj,dto,dtb,dtr,dtc,idirect

      print*,'mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds='
      print*,mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds
      print*,'sigf='
      print*,sigf
      print*,'pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,dtc='
      print*,pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,dtc
      print*,'ibyte= ',ibyte

      if (ni.ne.nyf .or. nj.ne.nxf .or. nz.ne.nk) then
        print*,'Grid Dimensions DO NOT MATCH'
        print*,'  nyf=',nyf,'nxf=',nxf,'nz=',nz
        print*,'  ni=',ni,'nj=',nj,'nk=',nk
        print*,'  Also check ibyte in postproc.param: ibyte= ',ibyte
        stop 'BAD DIMENSIONS'
      end if

      if (ierr.ne.0) then
        print*,'ERROR READING FILE'
        print*,'  Check ibyte in postproc.param: ibyte= ',ibyte
        stop 'EOF'
      end if

      ierrsum=0
      read (iin,rec=2,iostat=ierr) v2d
      print*,'ZS:     ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=3,iostat=ierr) v2d
      print*,'ZSSD:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=4,iostat=ierr) v2d
      print*,'LS:     ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=5,iostat=ierr) v2d
      print*,'SATBRT: ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=6,iostat=ierr) v2d
      print*,'XLAT:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=7,iostat=ierr) v2d
      print*,'XLON:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=8,iostat=ierr) v2d
      print*,'XMAP:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=9,iostat=ierr) v2d
      print*,'DMAP:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=10,iostat=ierr) v2d
      print*,'F:      ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1
      read (iin,rec=11,iostat=ierr) v2d
      print*,'MASK:   ierr=',ierr
      if (ierr.ne.0) ierrsum=ierrsum+1

      print*,'DONE READING FILE'
      print*,'  There were ',ierrsum,'read errors'


      end
