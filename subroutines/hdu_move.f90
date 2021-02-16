!---------------------------------------------------------------------!
    subroutine hdu_move(unit,hdu_name,status)
!PURPOSE:  Move to the HDU based on the name. 
      implicit none
      integer :: unit,hdutype,hdu_extver,status
      character (len=30) hdu_name,error_description

! The hdutype parameter may have a value of IMAGE HDU(0), ASCII TBL (1),
! BINARY TBL (2), or ANY HDU (-1) where ANY HDU means that
! only the extname and extver values will be used to locate the correct extension.
! If the input value of extver is 0 then the EXTVER keyword is ignored and the first HDU with a matching
! EXTNAME (or HDUNAME) keyword will be found. If no matching HDU is found in the file
! then the current HDU will remain unchanged and a status = BAD HDU NUM (301) will be returned.     
      hdutype = -1
      hdu_extver = 0
      call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      call ftgerr(status,error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
      
! if the extension is not RATE the code asks for a new name until it gets a right one      
      do while(status .eq. 301)
         write(*,*) "   No extension called: ",hdu_name
         write(*,*) "   Please specify the correct name or type 'no' to quit (in the case of GTI extension you can choose the GTI manually)"
         read(*,*) hdu_name
         if (hdu_name .eq. 'no') goto 111
         hdutype = -1
         hdu_extver = 0
         status = 0
         call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      enddo
111   continue 
    end subroutine  hdu_move
!---------------------------------------------------------------------!
