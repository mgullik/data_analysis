!---------------------------------------------------------------------!
    subroutine extract_lc(filename)
      use dyn_lc
      implicit none

! CFITSIO variable
      character (len=500) :: filename

      character (len=30)  :: error_description, hdu_name &
                             ,col_name_temp, keyword
      logical             :: anynul, yes_no
      integer             :: i, status, readwrite, blocksize, unit, chdu  &
                             ,nrow, nrow_GTI, ncol, ncol_GTI, colnum & 
                             ,felem, nelem, datacode, repeat, width, frow
      integer             :: colnum_f
      double precision    :: get_keyword_double, dt_temp
      real                :: nullval

      double precision    :: first_time_bin, first_time_bin_GTI

      real  , dimension(:), allocatable :: time_e, rate_e, bkg_e, start_GTI_e, end_GTI_e, err_rate_e
      double precision, dimension(:), allocatable :: time_d, rate_d, bkg_d, start_GTI_d, end_GTI_d, err_rate_d


! Cross-spectum variables
!      integer :: i,count,count2,dim
      

      status = 0
      readwrite = 0  !file access mode: readonly = 0 , read and write = 1
! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
      ! write(*,*) status 
      
! Open the FITS file
      call ftopen(unit,filename,readwrite,blocksize,status)
     ! call ftnopn(unit,filename,readwrite,status)  !thisi is used with the extension in the name (e.g. namefile.fits+2)
      ! write(*,*) status 
      
      call ftgerr(status,error_description)
      if (status .ne. 0) then
         write(*,*) '!! ATTENTION !! ', error_description
         write(*,*) '    Exit...'
         stop
      endif

      write(*,*) '------------------------------------------'
      write(*,*) '   Extraction of the light curve: ', trim(filename)
      write(*,*) 
      write(*,*) '   Check TIME, RATE, and GTI'
      
      
!Let's move in the correct HDU (with the name in hdu_name)      
      hdu_name = 'RATE'
      call hdu_move(unit,hdu_name,status)
      call ftghdn(unit,chdu)
      if (status .eq. 301) stop 
      status = 0 
!      write(*,*) '   Current hdu',chdu

! Get the number of rows and columns in the CHDU
      call ftgnrw(unit,nrow,status)
      call ftgncl(unit,ncol,status)
!      write(*,*) "number of rows", nrow

! extract the dt from fits file
      keyword = 'TIMEDEL'
      dt_temp = dt
      dt = get_keyword_double(unit,keyword,status)
!      write(*,*) 'dt', dt

! Check if the dt is equal to the previous one, besides the first time the dt is calculated (starting dt=-1)
      if(dt_temp .ne. -1) then 
         if (dt .ne. dt_temp) then
            write(*,*) '   dt of this new light curve is different from the previous one'
            write(*,*) '   to make the cross spectrum you need to have the same dt. EXIT...'
            stop
         endif
      endif


! extract the telescope starting time from fits
      ! keyword = 'TSTART'
      ! starting_telescope_time = get_keyword_double(unit,keyword,status)
!      write(*,*) 'dt', dt
      starting_telescope_time = 0.d0

      
!**************************************************************!      
! Let's get the column number for a particolar name with colnum_f (function)
      col_name_temp = '*TIME*'
      ! write(*,*) 'Look at the *TIME* column'
      colnum = colnum_f(unit,col_name_temp,status)
      if (colnum .gt. 0.0) then 
         write(*,*) '     TIME column found'

! Get the datatype of a column and read that column
!   repeat tells you if there is more than one element in very space of the table  

! Get the datatype
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         if( status .gt. 0 )call printerror(status)

! Get the the column values      
         frow = 1 !starting row
         nelem = nrow ! last row to read
         felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
         nullval = -1
         if(datacode .eq. 42) then
            write(*,*) '     real type'
            if(.not.allocated(time_e)) allocate(time_e(nrow))
            call ftgcve(unit,colnum,frow,felem,nelem,nullval,time_e,anynul,status)         
         else if (datacode .eq. 82) then
            write(*,*) '     double type'
            if(.not.allocated(time_d)) allocate(time_d(nrow))
            call ftgcvd(unit,colnum,frow,felem,nelem,nullval,time_d,anynul,status)
         else
            write(*,*) "   TIME column is niether a real nor a double"
            write(*,*) "   Modify the source code"
         end if
      else
         write(*,*) '    TIME NOT found, exit...'
         write(*,*) '    Check the fits file for exact name of the column '
         stop 
      endif
!**************************************************************!      
     
!**************************************************************!      
      col_name_temp = '*RATE*'
      ! write(*,*) 'Look at the *RATE* column'
      colnum = colnum_f(unit,col_name_temp,status)
      if (colnum .gt. 0.0) then 
         write(*,*) '     RATE column found'
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         if( status .gt. 0 )call printerror(status)

         frow = 1 !starting row
         nelem = nrow ! last row to read
         felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
         nullval = -1
         if(datacode .eq. 42) then
            write(*,*) '     real type'
            if(.not.allocated(rate_e)) allocate(rate_e(nrow))
            call ftgcve(unit,colnum,frow,felem,nelem,nullval,rate_e,anynul,status)         
         else if (datacode .eq. 82) then
            write(*,*) '     double type'
            if(.not.allocated(rate_d)) allocate(rate_d(nrow))
            call ftgcvd(unit,colnum,frow,felem,nelem,nullval,rate_d,anynul,status)         
         else
            write(*,*) "   RATE column is niether an real nor a double"
            write(*,*) "   Modify the source code"
         end if
      else
         write(*,*) '    RATE NOT found, exit...'
         write(*,*) '    Check the fits file for exact name of the column '
         stop 
      end if
!**************************************************************!      

! !**************************************************************!      
      col_name_temp = '*ERROR*'
      ! write(*,*) 'Look at the *ERROR* column'
      colnum = colnum_f(unit, col_name_temp, status)
      if (colnum .gt. 0.0) then 
         write(*,*) '     ERROR column found'
         call ftgtcl(unit, colnum, datacode, repeat, width, status)
         if( status .gt. 0 )call printerror(status)

         frow    = 1    !starting row
         nelem   = nrow ! last row to read
         felem   = 1    ! first pixel of the element vector (ignored for ASCII tables)
         nullval = -1
         if(datacode .eq. 42) then
            write(*,*) '     real type'
            if(.not.allocated(err_rate_e)) allocate(err_rate_e(nrow))
            call ftgcve(unit, colnum, frow, felem, nelem, nullval, err_rate_e, anynul, status)
         else if (datacode .eq. 82) then
            write(*,*) '     double type'
            if(.not.allocated(err_rate_d)) allocate(err_rate_d(nrow))
            call ftgcvd(unit, colnum, frow, felem, nelem, nullval, err_rate_d, anynul, status) 
         else
            write(*,*) "   ERROR column is niether an real nor a double"
            write(*,*) "   Modify the source code "
         end if
      else
         write(*,*) '    ERROR NOT found'
         write(*,*)
      endif
! !**************************************************************!      


!**************************************************************!      
      ! col_name_temp = '*BACKV*'
      ! write(*,*) 'Look at the *BACKV* column'
      ! colnum        = colnum_f(unit, col_name_temp, status)
      ! if (colnum .gt. 0.0) then 
      !    write(*,*) 'BACKV found'
      !    call ftgtcl(unit, colnum, datacode, repeat, width, status)
      !    if( status .gt. 0 )call printerror(status)

      !    frow    = 1    !starting row
      !    nelem   = nrow ! last row to read
      !    felem   = 1    ! first pixel of the element vector (ignored for ASCII tables)
      !    nullval = -1
      !    if(datacode .eq. 42) then
      !       if(.not.allocated(bkg_e)) allocate(bkg_e(nrow))
      !       call ftgcve(unit, colnum, frow, felem, nelem, nullval, bkg_e, anynul, status)
      !    else if (datacode .eq. 82) then
      !       if(.not.allocated(bkg_d)) allocate(bkg_d(nrow))
      !       call ftgcvd(unit, colnum, frow, felem, nelem, nullval, bkg_d, anynul, status) 
      !    else
      !       write(*,*) "   BACKV column is niether an real nor a double"
      !       write(*,*) "   Modify the source code "
      !    end if
      ! endif
!**************************************************************!      
      
!Allocation of the general array
      if(.not. allocated(lc) ) then 
         allocate(lc (nrow))
      else 
         if (nrow .ne. dim_lc) then 
            write(*,*) '   ATTENTION!! The light curves do not have the same length!'
            stop 
         endif
      endif
      
      ! if(.not. allocated(bkg)  ) allocate(bkg  (nrow))

!Set the dimension of the lc and write the lc, time, and bkg in the common arrays
      dim_lc = nrow

!CONVERSION FROM REAL TO DOUBLE PRECISION OF TIME AND 
      if (allocated(time_e)) then
         if(.not. allocated(time)  ) allocate(time(dim_lc))
         ! if(.not. allocated(time_tel)) allocate(time_tel(dim_lc))
         first_time_bin = dble(time_e(1))
         starting_telescope_time = first_time_bin
         do i = 1, dim_lc
            time(i) = dble(time_e(i)) - first_time_bin
            ! time_tel(i) = dble(time_e(i))
         enddo
         deallocate(time_e)
      else
         if(.not. allocated(time)  ) allocate(time(dim_lc))
         ! if(.not. allocated(time_tel)) allocate(time_tel(nrow))
         first_time_bin = time_d(1)
         starting_telescope_time = first_time_bin
         do i = 1, nrow
            time(i) = time_d(i) - first_time_bin
            ! time_tel(i) = time_d(i)
         enddo         
         deallocate(time_d)
      endif

!CONVERSION FROM REAL TO DOUBLE PRECISION OF RATE_E AND STORE INTO THE LC ARRAY
      if (allocated(rate_e)) then
         if(.not. allocated(lc)) allocate(lc(dim_lc))
         do i = 1, dim_lc
            lc(i) = dble(rate_e(i))
         enddo
         deallocate(rate_e)
      else
         if(.not. allocated(lc)) allocate(lc(dim_lc))
         do i = 1, dim_lc
            lc(i) = rate_d(i)
         enddo         
         deallocate(rate_d)
      endif

      if (allocated(err_rate_e)) then
         if(.not. allocated(err_rate)) allocate(err_rate(dim_lc))
         do i = 1, dim_lc
            err_rate(i) = dble(err_rate_e(i))
         enddo
         deallocate(err_rate_e)
      else
         if(.not. allocated(err_rate)) allocate(err_rate(dim_lc))
         do i = 1, dim_lc
            err_rate(i) = err_rate_d(i)
         enddo
      endif

            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      ! ! GET THE GTI ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!Let's move in the correct HDU (with the name in hdu_name)      
      hdu_name = 'SRC_GTIS'
      call hdu_move(unit, hdu_name,status)
      if (status .eq. 301) then
         if (yes_no('   Do you want to consider the full light curve length as GTI')) then
            nrow_GTI = 1
            if(.not.allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI)) 
            if(.not.allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
            start_GTI_d = 0.d0
            end_GTI_d   = time_d(nrow) + 1.d0
            status = 0
         else
            if (yes_no('   Do you want to manually write the GTI?')) then
               write(*,*) '     How many intervals in the GTI? (integer)'
               read(*,*) nrow_GTI
               if(.not.allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI)) 
               if(.not.allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
               do i = 1, nrow_GTI
                  write(*,'(A,I2)') '     Interval number:', i
                  write(*,*) '     Insert start GTI time'
                  read (*,*) start_GTI_d(i) 
                  write(*,*) '     Insert end GTI time'
                  read (*,*) end_GTI_d(i) 
               enddo
            else
               write(*,*)'    Options are finished...exit'               
               stop
            endif
         endif

      else

         call ftghdn(unit, chdu)
         if( status .gt. 0 ) call printerror(status)
         status = 0 
! Get the number of rows and columns in the CHDU
! and check if the number of columns is equal to 2      
         call ftgnrw(unit, nrow_GTI, status)
         call ftgncl(unit, ncol_GTI, status)
         if (ncol_GTI .ne. 2) then
            write(*,*) '    The selected GTI extension has not 2 columns! Exit...'
            stop
         endif
      
!**************************************************************!      
         col_name_temp = '*START*'
         colnum = colnum_f(unit, col_name_temp, status)
         if (colnum .gt. 0.0) then 
            write(*,*) '     GTI START found'
            call ftgtcl(unit, colnum, datacode, repeat, width, status)
            if( status .gt. 0 )call printerror(status)

            frow    = 1 !starting row
            nelem   = nrow_GTI ! last row to read
            felem   = 1 ! first pixel of the element vector (ignored for ASCII tables)
            nullval = -1
            if(datacode .eq. 42) then
               write(*,*) '     real type'
               if(.not.allocated(start_GTI_e)) allocate(start_GTI_e(nrow_GTI))
               call ftgcve(unit, colnum, frow, felem, nelem, nullval, start_GTI_e, anynul, status)         
            else if (datacode .eq. 82) then
               write(*,*) '     double type'
               if(.not.allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI))
               call ftgcvd(unit, colnum, frow, felem, nelem, nullval, start_GTI_d, anynul, status)         
            else
               write(*,*) "   START GTI column is niether an real nor a double"
               write(*,*) "   Modify the source code"
            endif
         else
            write(*,*) '    GTI START NOT found, exit...'
            write(*,*) '    Check the fits file for exact name of the column '
            stop 

         end if
!**************************************************************!      

!**************************************************************!      
         col_name_temp = '*STOP*'
         colnum = colnum_f(unit, col_name_temp, status)
         if (colnum .gt. 0.0) then 
            write(*,*) '     GTI STOP found'

            call ftgtcl(unit, colnum, datacode, repeat, width, status)
            if( status .gt. 0 )call printerror(status)

            frow    = 1        !starting row
            nelem   = nrow_GTI ! last row to read
            felem   = 1        ! first pixel of the element vector (ignored for ASCII tables)
            nullval = -1
            if(datacode .eq. 42) then
               write(*,*) '     real type'
               if(.not.allocated(end_GTI_e)) allocate(end_GTI_e(nrow_GTI))
               call ftgcve(unit, colnum, frow, felem, nelem, nullval, end_GTI_e, anynul, status)         
            else if (datacode .eq. 82) then
               write(*,*) '     double type'
               if(.not.allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
               call ftgcvd(unit, colnum, frow, felem, nelem, nullval, end_GTI_d, anynul, status)         
            else
               write(*,*) "   END GTI column is niether an real nor a double"
               write(*,*) "   Modify the source code"
            end if
         else
            write(*,*) '    GTI STOP NOT found, exit...'
            write(*,*) '    Check the fits file for exact name of the column '
            stop
         endif
         write(*,*) 
         write(*,*) '   Extraction ended'
         write(*,*) '------------------------------------------'
 !**************************************************************!


         !CONVERSION FROM REAL TO DOUBLE PRECISION
         if (first_time_bin .lt. 1.d0) then ! this means that the starting time of the light curve is 0.0 then we need to set the GTI starting from zero 
            
            if (allocated(start_GTI_e)) then
               if(.not. allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI))
               ! if(.not. allocated(start_GTI_tel_tim)) allocate(start_GTI_tel_time(nrow_GTI))
               first_time_bin_GTI = dble(start_GTI_e(1)) 
               do i = 1, nrow_GTI
                  start_GTI_d(i) = dble(start_GTI_e(i)) - first_time_bin_GTI
                  ! start_GTI_tel_time(i) =  dble(start_GTI_e(i))  
               enddo

               deallocate(start_GTI_e)
            else
               ! if(.not. allocated(start_GTI_tel_tim)) allocate(start_GTI_tel_time(nrow_GTI))
               first_time_bin_GTI = start_GTI_d(1) 
               do i = 1, nrow_GTI
                  start_GTI_d(i) = start_GTI_d(i) - first_time_bin_GTI
                  ! start_GTI_tel_time(i) =  start_GTI_d(i)  
               enddo

            endif

            if (allocated(end_GTI_e)) then
               if(.not. allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
               ! if(.not. allocated(end_GTI_tel_time)) allocate(end_GTI_tel_time(nrow_GTI))
               do i = 1, nrow_GTI
                  end_GTI_d(i) = dble(end_GTI_e(i)) - first_time_bin_GTI
                  ! end_GTI_tel_time(i) = dble(end_GTI_e(i))
               enddo
               deallocate(end_GTI_e)
            else
               ! if(.not. allocated(end_GTI_tel_time)) allocate(end_GTI_tel_time(nrow_GTI))
               do i = 1, nrow_GTI
                  end_GTI_d(i) = end_GTI_d(i)  - first_time_bin_GTI
                  ! end_GTI_tel_time(i) = end_GTI_d(i) 
               enddo
            endif

            
! The starting time is the telescope time so we need to subtract the GTI from the starting starting time of the light curve
         else 

            if (allocated(start_GTI_e)) then
               if(.not. allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI))
               ! if(.not. allocated(start_GTI_tel_tim)) allocate(start_GTI_tel_time(nrow_GTI))
               do i = 1, nrow_GTI
                  start_GTI_d(i) = dble(start_GTI_e(i)) - first_time_bin
                  ! start_GTI_tel_time(i) =  dble(start_GTI_e(i))  
               enddo

               deallocate(start_GTI_e)
            else
               ! if(.not. allocated(start_GTI_tel_tim)) allocate(start_GTI_tel_time(nrow_GTI))
               first_time_bin_GTI = start_GTI_d(1) 
               do i = 1, nrow_GTI
                  start_GTI_d(i) = start_GTI_d(i) - first_time_bin
                  ! start_GTI_tel_time(i) =  start_GTI_d(i)  
               enddo

            endif

            if (allocated(end_GTI_e)) then
               if(.not. allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
               ! if(.not. allocated(end_GTI_tel_time)) allocate(end_GTI_tel_time(nrow_GTI))
               do i = 1, nrow_GTI
                  end_GTI_d(i) = dble(end_GTI_e(i)) - first_time_bin
                  ! end_GTI_tel_time(i) = dble(end_GTI_e(i))
               enddo
               deallocate(end_GTI_e)
            else
               ! if(.not. allocated(end_GTI_tel_time)) allocate(end_GTI_tel_time(nrow_GTI))
               do i = 1, nrow_GTI
                  end_GTI_d(i) = end_GTI_d(i)  - first_time_bin
                  ! end_GTI_tel_time(i) = end_GTI_d(i) 
               enddo
            endif

         endif
      endif

      
!HERE it is possible to print the GTI and the gaps
      ! write(*,*) '-----------------------------------------'
      ! write(*,*) 'GTI: ', nrow_GTI
      ! do i = 1, nrow_GTI - 1
      !    write(*,*) start_GTI_d(i), end_GTI_d(i)
      !    write(*,*) 'ok time ',  end_GTI_d(i) -  start_GTI_d(i)
      !    write(*,*) 'gap ', start_GTI_d(i + 1) - end_GTI_d(i)
      ! enddo
      ! write(*,*) start_GTI_d(nrow_GTI), end_GTI_d(nrow_GTI)
      ! write(*,*) '-----------------------------------------'

!HERE it is possible to print the light curve 

      ! do i=1, nrow
      !    write(99,*) time_d(i), rate_d(i)
      ! enddo

      ! do i=1, nrow_GTI
      !    write(98,*) start_GTI_d(i),end_GTI_d(i)
      ! enddo

      ! write(*,*) 'First time bin of the GTI', first_time_bin
      ! write(*,*)

      
!Fill the GTI 
!It is complicated because we have to distinguish between the first call and the others
! and check if the GTIs are the same in all the calls (checking if the GTIs are identical to the previous call)      
      if(.not. allocated(start_GTI)) then  
         allocate(start_GTI(nrow_GTI))
         allocate(end_GTI(nrow_GTI))
!Set the dimension of the GTI and write them in the common arrays      
         dim_GTI = nrow_GTI
!Fill the actual GTI         
         do i = 1, dim_GTI
            start_GTI(i) = start_GTI_d(i)
            end_GTI  (i) = end_GTI_d(i)
            ! write(*,*) 'GTI', start_GTI_d(i), end_GTI_d(i)
         enddo
         
      else
         ! write(*, *) ' check_gap_num, nrow_GTI, dim_GTI', check_gap_num, nrow_GTI, dim_GTI
         if (check_gap_num .ne. nrow_GTI) then 
            write(*,*) '  ATTENTION!! The GTIs are not the same in all the light curves'
            stop 
         else 
!deallocate the GTIs because they have been modified by the interpolation 
            deallocate(start_GTI)
            deallocate(end_GTI)
!allocate the new GTIs (they will be modify by interpolation)
            allocate(start_GTI(nrow_GTI))
            allocate(end_GTI(nrow_GTI))
            dim_GTI = nrow_GTI
         do i = 1, dim_GTI
            start_GTI(i) = start_GTI_d(i)
            end_GTI  (i) = end_GTI_d(i)
         enddo
         endif
      endif

      
! Close and show if there are errors      
      call ftclos(unit,status)
      call ftgerr(status,error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description

      if(allocated(time_e)) deallocate(time_e)
      if(allocated(time_d)) deallocate(time_d)
      if(allocated(rate_e)) deallocate(rate_e)
      if(allocated(rate_d)) deallocate(rate_d)
      if(allocated(bkg_e) ) deallocate(bkg_e )
      if(allocated(bkg_d) ) deallocate(bkg_d )
      if(allocated(err_rate_d) ) deallocate(err_rate_d )
      if(allocated(err_rate_e) ) deallocate(err_rate_e )
      
      if(allocated(start_GTI_d)) deallocate(start_GTI_d)
      if(allocated(end_GTI_d)) deallocate(end_GTI_d)
      if(allocated(start_GTI_e)) deallocate(start_GTI_e)
      if(allocated(end_GTI_e)) deallocate(end_GTI_e)

    end subroutine extract_lc
!---------------------------------------------------------------------!
