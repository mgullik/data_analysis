!---------------------------------------------------------------------!
      subroutine deletefile(filename,status)
!  A simple little routine to delete a FITS file
      implicit none
      integer status,unit,blocksize
      character*(*) filename
!  Simply return if status is greater than zero
      if (status .gt. 0)return

!  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
          call ftdelt(unit,status)
          !write(*,*)"Deleted a file"
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          !write(*,*)"Didn't delete a file"
          call ftcmsg
      else
!         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
!-----------------------------------------------------------------------
