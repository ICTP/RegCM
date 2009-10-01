C CRU05 1901-1998 0.5 degree Latitude/Longitude Monthly Climate Data
C 
C Variables                    Directory     Units         Filename
C Precipitation                pre           mm            glo.pre.YYYY
C Wet day frequency            rd0           days          glo.rd0.YYYY
C Mean temperature             tmp           degC*10       glo.tmp.YYYY
C Diurnal temperature range    dtr           degC*10       glo.dtr.YYYY
C Vapour pressure              vap           hPa*10        glo.vap.YYYY
C Cloud cover                  cld           oktas*10      glo.cld.YYYY
C 
C Format (ascii)
C One file per year per variable.
C The first and second lines of the file contain information on the grid size,
C  for example, for a global field:
C grd_sz    xmin    ymin    xmax    ymax  n_cols  n_rows n_months 
C  0.5    -179.75  -89.75  189.75  89.75    720     360     12
C 
C This is followed by 12 monthly grids that are n_cols by n_rows in size.
C Each record contains n_cols longitude values, format n_colsi5, 
C starting at xmin and ending at xmax.  The first record starts at January,
C ymin and the last record is the row for December, ymax.
C Missing values (Antarctica, oceans and major inland water bodies) are 
C assigned values of -9999.
C To read in one year of data the following program should work:

      program rd_ascii
      implicit none
      include '../../Main/regcm.param'
c f90 program to read in an integer ascii grid with
c  variable dimensions into a global grid (720x360x12)
c
      integer n_cols,n_rows
      parameter (n_cols=720,n_rows=360)  
      integer*2 grid(n_cols,n_rows,12)
      real*4  prec(n_cols,n_rows)
      integer imx,jmx
      parameter (imx=jx-2,jmx=ix-2)
      real*4  xlon(jmx,imx),xlat(jmx,imx),precOUT(jmx,imx)
      character infl*72,fmt*20
      real*4  gs,xmin,ymin,xmax,ymax
      integer ncols,nrows,nmonths
      integer im,lat,lon,i,j,ix0,iy,ix1
      real*4  x,y,add,sum
C
      open(49,file='OUT_HEAD',form='unformatted'
     &         ,access='direct',recl=jmx*imx*4)
      read(49,rec=6) xlat
      read(49,rec=7) xlon
      write(*,*) 'Enter ascii grid file name'
      read(*,'(a72)')infl
      open(1,file=infl,status='old')
      read(1,*)
      read(1,*)gs,xmin,ymin,xmax,ymax,ncols,nrows,nmonths
      do im=1,nmonths
         do lat=1,n_rows
            do lon=1,n_cols
               grid(lon,lat,im)=-9999
            enddo
         enddo
      enddo
      fmt='(   i5)'
      write(fmt(2:4),'(i3)')n_cols
      open(10,file='CRU_DAT',form='unformatted'
     &       ,access='direct',recl=jmx*imx*4)
      do im=1,nmonths
         do lat=1,n_rows
            read(1,fmt)(grid(lon,lat,im),lon=1,n_cols)
            do lon=1,n_cols
               prec(lon,lat)=grid(lon,lat,im)
            enddo
         enddo
         do j=1,imx
         do i=1,jmx
            x=(xlon(i,j)+179.75+0.50000001)*2.
            y=(xlat(i,j)+ 89.75+0.50000001)*2.
            ix0=x
            ix1=ix0+1
            if(ix0.eq.720) ix1=1
            iy=y
            if(x-ix0.le.0.00001) then
               if(prec(ix0,iy).lt.-9998.) then
                  if(prec(ix0,iy+1).lt.-9998.) then
                     precOUT(i,j) = -9999.
                  else
                     precOUT(i,j) = prec(ix0,iy+1)
                  endif
               else 
                  if(prec(ix0,iy+1).lt.-9998.) then
                     precOUT(i,j) = prec(ix0,iy)
                  else
                     precOUT(i,j) = prec(ix0,iy)*(iy+1-y)
     &                            + prec(ix0,iy+1)*(y-iy)
                  endif
               endif
            else if(y-iy.le.0.00001) then
               if(prec(ix0,iy).lt.-9998.) then
                  if(prec(ix1,iy).lt.-9998.) then
                     precOUT(i,j) = -9999.
                  else
                     precOUT(i,j) = prec(ix1,iy)
                  endif
               else
                  if(prec(ix1,iy).lt.-9998.) then
                     precOUT(i,j) = prec(ix0,iy)
                  else
                     precOUT(i,j) = prec(ix0,iy)*(ix0+1-x)
     &                            + prec(ix1,iy)*(x-ix0)
                  endif
               endif
            else
               if(prec(ix0,iy).lt.-9998.0.and.
     &            prec(ix1,iy).lt.-9998.0.and.
     &            prec(ix0,iy+1).lt.-9998.0.and.
     &            prec(ix1,iy+1).lt.-9998.0) then
                  precOUT(i,j) = -9999.
               else
                  sum=0.0
                  add=0.0
                  if(prec(ix0,iy).gt.-9998.0) then
                     sum=sum+prec(ix0,iy)*(ix0+1-x)*(iy+1-y)
                     add=add+(ix0+1-x)*(iy+1-y)
                  endif
                  if(prec(ix1,iy).gt.-9998.0) then
                     sum=sum+prec(ix1,iy)*(x-ix0)*(iy+1-y)
                     add=add+(x-ix0)*(iy+1-y)
                  endif
                  if(prec(ix0,iy+1).gt.-9998.0) then
                     sum=sum+prec(ix0,iy+1)*(ix0+1-x)*(y-iy)
                     add=add+(ix0+1-x)*(y-iy)
                  endif
                  if(prec(ix1,iy+1).gt.-9998.0) then
                     sum=sum+prec(ix1,iy+1)*(x-ix0)*(y-iy)
                     add=add+(x-ix0)*(y-iy)
                  endif
                  precOUT(i,j)=sum/add
                  if(add.lt.0.25) precOUT(i,j)=-9999.0
               endif
            endif
         enddo
         enddo
         write(10,rec=im) precOUT
      enddo
      stop
      end

