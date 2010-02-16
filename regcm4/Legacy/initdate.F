C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine initdate
      implicit none
c
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/message.cb"
      integer mdatez
      common /datenum/mdatez(289276)
      integer mbase,nbase,nrec,nyear,mon,nday,i,m
c
      if (ibdyfrq.eq.6) then
        nrec=0
        do nyear=1948,2145
          mbase = nyear*1000000
          do mon=1,12
            mbase = mbase+10000
            if (mon.eq.1.or.mon.eq.3.or.mon.eq.5.or.mon.eq.7
     &                 .or.mon.eq.8.or.mon.eq.10.or.mon.eq.12) then
              nday=31
            else if (mon.eq.4.or.mon.eq.6.or.mon.eq.9.or.mon.eq.11)then
              nday=30
            else
              nday=28
              if (mod(nyear,4)  .eq.0) nday=nday+1
              if (mod(nyear,100).eq.0) nday=nday-1
              if (mod(nyear,400).eq.0) nday=nday+1
            end if
            nbase = mbase
            do i=1,nday
              nbase = nbase+100
              do m=1,4
                nrec=nrec+1
                if (m.eq.1) then
                  mdatez(nrec)=nbase
                else if (m.eq.2) then
                  mdatez(nrec)=nbase+6
                else if (m.eq.3) then
                  mdatez(nrec)=nbase+12
                else
                  mdatez(nrec)=nbase+18
                end if
              end do
            end do
          end do
        end do
      else if (ibdyfrq.eq.3) then
        nrec=0
        do nyear=1948,2046
          mbase = nyear*1000000
          do mon=1,12
            mbase = mbase+10000
            if (mon.eq.1.or.mon.eq.3.or.mon.eq.5.or.mon.eq.7
     &                  .or.mon.eq.8.or.mon.eq.10.or.mon.eq.12) then
              nday=31
            else if (mon.eq.4.or.mon.eq.6.or.mon.eq.9.or.mon.eq.11)then
              nday=30
            else
              nday=28
              if (mod(nyear,4)  .eq.0) nday=nday+1
              if (mod(nyear,100).eq.0) nday=nday-1
              if (mod(nyear,400).eq.0) nday=nday+1
            end if
            nbase = mbase
            do i=1,nday
              nbase = nbase+100
              do m=1,8
                nrec=nrec+1
                if (m.eq.1) then
                  mdatez(nrec)=nbase
                else if (m.eq.2) then
                  mdatez(nrec)=nbase+3
                else if (m.eq.3) then
                  mdatez(nrec)=nbase+6
                else if (m.eq.4) then
                  mdatez(nrec)=nbase+9
                else if (m.eq.5) then
                  mdatez(nrec)=nbase+12
                else if (m.eq.6) then
                  mdatez(nrec)=nbase+15
                else if (m.eq.7) then
                  mdatez(nrec)=nbase+18
                else
                  mdatez(nrec)=nbase+21
                end if
              end do
            end do
          end do
        end do
      else if (ibdyfrq.eq.12) then
        nrec=0
        do nyear=1948,2145
          mbase = nyear*1000000
          do mon=1,12
            mbase = mbase+10000
            if (mon.eq.1.or.mon.eq.3.or.mon.eq.5.or.mon.eq.7
     &                 .or.mon.eq.8.or.mon.eq.10.or.mon.eq.12) then
              nday=31
            else if (mon.eq.4.or.mon.eq.6.or.mon.eq.9.or.mon.eq.11)then
              nday=30
            else
              nday=28
              if (mod(nyear,4)  .eq.0) nday=nday+1
              if (mod(nyear,100).eq.0) nday=nday-1
              if (mod(nyear,400).eq.0) nday=nday+1
            end if
            nbase = mbase
            do i=1,nday
              nbase = nbase+100
              do m=1,2
                nrec=nrec+1
                if (m.eq.1) then
                  mdatez(nrec)=nbase
                else if (m.eq.2) then
                  mdatez(nrec)=nbase+12
                end if
              end do
            end do
          end do
        end do
      else
        write(aline,*) 'please double check your ibdyfrq'
        call say
        write(aline,*) 'we only support ibdyfrq = 3, 6, or 12'
        call say
        write(aline,*) 'and we use ibdyfrq = 6 as the default'
        call say
        call fatal(__FILE__,__LINE__,'UNSUPPORTED BOUNDARY FREQUENCY')
      end if
      write(aline,*) 'nrec = ',nrec
      call say
      return
      end subroutine initdate
