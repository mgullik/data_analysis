!-----------------------------------------------------------------------
      subroutine printerror(status)
!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      integer status
      character errtext*30,errmessage*80
!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return
!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext
!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
       end do
       stop 
    end subroutine printerror
!-----------------------------------------------------------------------
