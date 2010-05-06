!!>
!!c  ROUTINE: HEADER
!!c  ACTION: print out initial information
!!c          about machine/time,compilation stamp 
!!c          and compilation flags 
!!c  
!!c  
!!************************************************************
!!<

SUBROUTINE header(myid)

  IMPLICIT NONE 
  !! local variables:
  INTEGER ihost,idir
  INTEGER hostnm
  INTEGER getcwd
  CHARACTER (len=24) :: data='?'
  CHARACTER (len=32) :: hostname='?' 
  CHARACTER (len=30) :: user='?' 
  CHARACTER (len=100) :: directory='?'
  INTEGER :: nrite=6, myid
  
  if (myid.eq.0)  then 

     !!open here the output file:
!     OPEN ( 20,file='OUTPUT')
     WRITE (nrite,"(/,2x,'This is RegCM version 4 ')")
     WRITE (nrite,100)  SVN_REV, __DATE__ , __TIME__   
100  FORMAT(2x,' SVN Revision: ',A,' compiled at: data : ',A,'  time: ',A,/)

#ifdef IBM
     hostname='ibm platform '
     user= 'Unknown'
     call fdate_(data)
#else
     Ihost = hostnm(hostname)
     call getlog(user)
     CALL fdate(data)
#endif 

     Idir=GETCWD(directory)


     WRITE(nrite,*) ": this run start at    : ",data
     WRITE(nrite,*) ": it is submitted by   : ",trim(user)
     WRITE(nrite,*) ": it is running on     : ",trim(hostname)
     WRITE(nrite,*) ": in directory         : ",trim(directory)
     WRITE(nrite,*) "                      " 
     end if 
  RETURN 
END SUBROUTINE header
