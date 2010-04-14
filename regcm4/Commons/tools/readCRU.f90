! CRU05 1901-1998 0.5 degree Latitude/Longitude Monthly Climate Data
!
! Variables                    Directory     Units         Filename
! Precipitation                pre           mm            glo.pre.YYYY
! Wet day frequency            rd0           days          glo.rd0.YYYY
! Mean temperature             tmp           degC*10       glo.tmp.YYYY
! Diurnal temperature range    dtr           degC*10       glo.dtr.YYYY
! Vapour pressure              vap           hPa*10        glo.vap.YYYY
! Cloud cover                  cld           oktas*10      glo.cld.YYYY
!
! Format (ascii)
! One file per year per variable.
! The first and second lines of the file contain information on the grid size,
!  for example, for a global field:
! grd_sz    xmin    ymin    xmax    ymax  n_cols  n_rows n_months
!  0.5    -179.75  -89.75  189.75  89.75    720     360     12
!
! This is followed by 12 monthly grids that are n_cols by n_rows in size.
! Each record contains n_cols longitude values, format n_colsi5,
! starting at xmin and ending at xmax.  The first record starts at January,
! ymin and the last record is the row for December, ymax.
! Missing values (Antarctica, oceans and major inland water bodies) are
! assigned values of -9999.
! To read in one year of data the following program should work:
 
      program readCRU
      use mod_regcm_param , only : ix => iy
      use mod_regcm_param , only : jx => jx
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: n_cols = 720 , n_rows = 360 ,              &
                          &  imx = jx - 2 , jmx = ix - 2
!
! Local variables
!
      real(4) :: add , gs , x , xmax , xmin , xsum , y , ymax , ymin
      character(20) :: cfmt
      integer(2) , dimension(n_cols,n_rows,12) :: grid
      integer :: i , im , ix0 , ix1 , iy , j , lat , lon , ncols ,      &
               & nmonths , nrows
      character(72) :: infl
      real(4) , dimension(n_cols,n_rows) :: prec
      real(4) , dimension(jmx,imx) :: precout , xlat , xlon
!
!     f90 program to read in an integer ascii grid with
!     variable dimensions into a global grid (720x360x12)
!
      open (49,file='OUT_HEAD',form='unformatted',access='direct',      &
          & recl=jmx*imx*4)
      read (49,rec=6) xlat
      read (49,rec=7) xlon
      write (*,*) 'Enter ascii grid file name'
      read (*,'(a72)') infl
      open (1,file=infl,status='old')
      read (1,*)
      read (1,*) gs , xmin , ymin , xmax , ymax , ncols , nrows ,       &
               & nmonths
      do im = 1 , nmonths
        do lat = 1 , n_rows
          do lon = 1 , n_cols
            grid(lon,lat,im) = -9999
          end do
        end do
      end do
      cfmt = '(   i5)'
      write (cfmt(2:4),'(i3)') n_cols
      open (10,file='CRU_DAT',form='unformatted',access='direct',       &
          & recl=jmx*imx*4)
      do im = 1 , nmonths
        do lat = 1 , n_rows
          read (1,cfmt) (grid(lon,lat,im),lon=1,n_cols)
          do lon = 1 , n_cols
            prec(lon,lat) = grid(lon,lat,im)
          end do
        end do
        do j = 1 , imx
          do i = 1 , jmx
            x = (xlon(i,j)+179.75+0.50000001)*2.
            y = (xlat(i,j)+89.75+0.50000001)*2.
            ix0 = x
            ix1 = ix0 + 1
            if ( ix0==720 ) ix1 = 1
            iy = y
            if ( x-ix0<=0.00001 ) then
              if ( prec(ix0,iy)<-9998. ) then
                if ( prec(ix0,iy+1)<-9998. ) then
                  precout(i,j) = -9999.
                else
                  precout(i,j) = prec(ix0,iy+1)
                end if
              else if ( prec(ix0,iy+1)<-9998. ) then
                precout(i,j) = prec(ix0,iy)
              else
                precout(i,j) = prec(ix0,iy)*(iy+1-y) + prec(ix0,iy+1)   &
                             & *(y-iy)
              end if
            else if ( y-iy<=0.00001 ) then
              if ( prec(ix0,iy)<-9998. ) then
                if ( prec(ix1,iy)<-9998. ) then
                  precout(i,j) = -9999.
                else
                  precout(i,j) = prec(ix1,iy)
                end if
              else if ( prec(ix1,iy)<-9998. ) then
                precout(i,j) = prec(ix0,iy)
              else
                precout(i,j) = prec(ix0,iy)*(ix0+1-x) + prec(ix1,iy)    &
                             & *(x-ix0)
              end if
            else if ( prec(ix0,iy)<-9998.0 .and. prec(ix1,iy)           &
                    & <-9998.0 .and. prec(ix0,iy+1)<-9998.0 .and.       &
                    & prec(ix1,iy+1)<-9998.0 ) then
              precout(i,j) = -9999.
            else
              xsum = 0.0
              add = 0.0
              if ( prec(ix0,iy)>-9998.0 ) then
                xsum = xsum + prec(ix0,iy)*(ix0+1-x)*(iy+1-y)
                add = add + (ix0+1-x)*(iy+1-y)
              end if
              if ( prec(ix1,iy)>-9998.0 ) then
                xsum = xsum + prec(ix1,iy)*(x-ix0)*(iy+1-y)
                add = add + (x-ix0)*(iy+1-y)
              end if
              if ( prec(ix0,iy+1)>-9998.0 ) then
                xsum = xsum + prec(ix0,iy+1)*(ix0+1-x)*(y-iy)
                add = add + (ix0+1-x)*(y-iy)
              end if
              if ( prec(ix1,iy+1)>-9998.0 ) then
                xsum = xsum + prec(ix1,iy+1)*(x-ix0)*(y-iy)
                add = add + (x-ix0)*(y-iy)
              end if
              precout(i,j) = xsum/add
              if ( add<0.25 ) precout(i,j) = -9999.0
            end if
          end do
        end do
        write (10,rec=im) precout
      end do
      end program readCRU
