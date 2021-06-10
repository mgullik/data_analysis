!---------------------------------------------------------------------!
    subroutine split_lc()
      use dyn_lc
      implicit none

      integer               :: i, new_dim_GTI
      double precision      :: tot_time_lc
      logical               :: power2, check_interval, check, yes_no 
      double precision, allocatable  :: temp_array(:), temp_GTI1(:), temp_GTI2(:)

! The first time (when int_number is -1 and split_ind is not allocated) we call the lc_split_first to work out split_index array
! Then we call only lc_split_index

!At this point there is the possibility to print the complete light curve without the interpolated gaps            
      ! do i=1, nrow
      !    write(98,*) time(i), lc(i)
      ! enddo


      if (.not. allocated(split_ind)) then
         write(*,*) 
         write(*,*) '   TOTAL length of the light curve (sec) considering the gaps (final_time - starting time): ', time(dim_lc),' - ', time(1) , ' = ',  time(dim_lc) - time(1) 
         write(*,*) 
         write(*,*) '   Number of gaps in the light curve: ', dim_GTI - 1
         write(*,*) 

         check_gap_num = dim_GTI


! Check if you can interpolate with the current routine: the light curve needs no jumps in time 
         check = .true.
         do i = 2 , dim_lc
!         write(*,*) 'time index', i
            ! write(*,*) time(i) - time(i - 1)
            if ((time(i) - time(i - 1)) .gt. 2.d0 * dt  ) then
               check = .false.
            endif
         enddo

         write(*,*)
         write(*,*) '   GTI situation of the observation'
         do i = 1, dim_GTI
            write(*,'(A,I4,A,F12.2,A)') '   GTI number ', i, ' length ', end_GTI(i) - start_GTI(i), ' sec'
         enddo
         write(*,*)
            
         if (check) then
            if (verbose_merge) then !*** TO RUN AUTOMATICALLY
            else
               goto 101
            endif
            if (yes_no('    Do you want to interpolate the light curve gaps?')) then
                  
101            continue
                  
! This is to save the rate before the interpolation 
               if(.not. allocated(temp_array)) allocate(temp_array(dim_lc))
               if(.not. allocated(temp_GTI1) ) allocate(temp_GTI1 (dim_GTI))
               if(.not. allocated(temp_GTI2) ) allocate(temp_GTI2 (dim_GTI))
               temp_array = lc 
               temp_GTI1  = start_GTI
               temp_GTI2  = end_GTI
            
! This is a while loop because after the interpolation the user might want to interpolate differently
               do             
               ! write(*,*) 'start do loop'
! CALL THE INTERPOLATION SUBROUTINE TO FILL THE SMALL GAPS      
                  call interpol_split_silent(new_dim_GTI)
                  write(*,*)
                  write(*,*) 'After the interpolation the number of gaps is: ', new_dim_GTI - 1
                  write(*,*)
                  
                  if (verbose_merge) then !*** TO RUN AUTOMATICALLY
                  else
                     exit
                  endif
                  
                  if (yes_no('   Do you want to interpolate differently?')) then
                     lc      = temp_array
                     start_GTI = temp_GTI1
                     end_GTI   = temp_GTI2
                     gap = -1
                  else
                     exit
                  endif

               ! write(*,*) 'end do loop'
               enddo
            
               if(allocated(temp_array)) deallocate(temp_array)
               if(allocated(temp_GTI1) ) deallocate(temp_GTI1)
               if(allocated(temp_GTI2) ) deallocate(temp_GTI2)

               write(*,*)
               write(*,*) '   GTI situation of the observation aftert the interpolation'
               do i = 1, dim_GTI
                  write(*,'(A,I4,A,F12.2,A)') '   GTI number ', i, ' length ', end_GTI(i) - start_GTI(i), ' sec'
               enddo
               write(*,*)
               
! RE-SET THE GTI INTERVALS BASED ON THE INTERPOLATION
               allocate(temp_array(new_dim_GTI))
               do i=1, new_dim_GTI
                  temp_array(i) = start_GTI(i)
               enddo
               deallocate(start_GTI)      
               allocate(start_GTI(new_dim_GTI))
               do i = 1, new_dim_GTI
                  start_GTI(i) = temp_array(i)  ! write the new stat_GTI
                  temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
               enddo
               deallocate(end_GTI)
               allocate(end_GTI(new_dim_GTI))
               do i = 1, new_dim_GTI
                  end_GTI(i) = temp_array(i)  ! write the new end_GTI
               enddo
               deallocate(temp_array)

               dim_GTI = new_dim_GTI

            else
               check_gap_num = dim_GTI
            endif

         else
            write(*,*) 'No interpolation of the light curve is possible with the current routine'
            check_gap_num = dim_GTI
         endif

         write(*,*)
         tot_time_lc = 0.d0
         do i = 1, dim_GTI
            tot_time_lc = tot_time_lc +  end_GTI(i) - start_GTI(i)
         enddo
         write(*,'(A, F8.1, A)') '   Total exposure time of the light curve ', tot_time_lc, ' s'
         write(*,*)

         do 
! !length of the intervals
            write(*,*) '   Enter the length of the interval in steps: '
            write(*,'(A, F10.4, A)') '      !! Remember that dt is ', dt, ' s'
            if (verbose_merge) then 
               read(*,*) int_len_dim
            else
               int_len_dim = int(tot_time_lc / dt)
               write(*,*) '    length of the interval is: ',  int_len_dim
            endif
            write(*,*)

!Compute the split in the light curve and create the split_ind array
            call lc_split_first(check_interval)
            if (check_interval) then 
               write(*,*) 'Number of segments in the light curve: ', int_number
               write(*,*)
               if(.not. allocated(lc_int)  ) allocate(lc_int  (int_number, int_len_dim) )
               if(.not. allocated(time_int)) allocate(time_int(int_number, int_len_dim) )

            !Fill the time_int, lc_int and bkg_int which are the light curves separated in intervals
               call lc_split_index_time()

! Check if int_len_dim is a power of 2
               if (power2(int_len_dim) .eqv. .false.)then
                  write(*,*) '   Precedure without the FFT (it is not a power of 2)'
                  check_power2 = .false.
               else
               write(*, *) '   Procedure uses FFT'
                  check_power2 = .true.
               endif
               
               exit
            endif
            write(*,*) '  '
            write(*,*) '  Try again...'
         enddo
         
      else
!THIS is called if it's not the first time

!Call interpol_split_silent because the new lc needs interpolation. The number of gaps is saved from the first time          
!There is no need to adjust the GTI arrays because the split is already done. 
! The IF statement before calling the subroutine is because if it is false the first time it should be false every time (no gap interpolation!!)
         write(*,*) '    Automatic interpolation of gaps shorter than (number of bins): ', gap

         if (gap .ne. -1) then 
            call interpol_split_silent(new_dim_GTI)
         ! write(*,*) 'ciao 1'
            allocate(temp_array(new_dim_GTI))
         ! write(*,*) 'ciao 1.1'
            do i=1, new_dim_GTI
         ! write(*,*) 'ciao loop', i 
               temp_array(i) = start_GTI(i)
            enddo
         ! write(*,*) 'ciao 2'
            deallocate(start_GTI)      
            allocate(start_GTI(new_dim_GTI))
            do i = 1, new_dim_GTI
               start_GTI(i) = temp_array(i)  ! write the new stat_GTI
               temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
            enddo
            ! write(*,*) 'ciao 3'
            deallocate(end_GTI)
            allocate(end_GTI(new_dim_GTI))
            do i = 1, new_dim_GTI
               end_GTI(i) = temp_array(i)  ! write the new end_GTI
            enddo
            deallocate(temp_array)
            ! write(*,*) 'ciao 4'

            dim_GTI = new_dim_GTI
         
         endif
         print *, '   Number of gaps in the total light curve: : ', dim_GTI - 1

         call lc_split_index_time()
      endif

      write(*,*) '  Split completed'
     
    end subroutine split_lc
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
  subroutine interpol_split_silent(new_dim_GTI)
    !This subroutine differs from interpol_split only because when it asks to  
    !set the maximum number of bins that you want to interpolate, it has the 
    !option to skip this if it has been already asked. This is why is called 
    !SILENT, it has the option to be silent. 

    !This subroutine fills the gaps in the light curve according to the maximum number
    ! of bins that the user chooses.
    !Then the GTI arrays (start and end) are modified and a new dimension is set.
    !Basically they become larger.
    !NOTE: after this function the GTI arrays are re-set (deallocate and reallocate)
    !      with the correct number of bins
    
    use dyn_lc
    ! use rand
    ! use rand_fun_real
    implicit none
    integer, intent(OUT)   :: new_dim_GTI

    integer                :: j, i, k, w, ee,  max_gap
    double precision       :: m, q, lc_interpol, m_b, max_gap_sec, q_b, bkg_interpol
    double precision, allocatable   :: new_start_GTI(:), new_end_GTI(:)
    integer :: poisson
    integer :: funPoissonSingle

    ! do j = 1, dim_GTI
    !    write(*,*) start_GTI(j), end_GTI(j)
    ! enddo
    

    allocate(new_start_GTI(dim_GTI))
    allocate(new_end_GTI(dim_GTI))
    
    w = 1 
    new_start_GTI(w) = start_GTI(w)

    if (gap .eq. -1) then 

       write(*,*)
       write(*,*) '  The gap situation is: '
       do j = 1, dim_GTI - 1
          ee = int( ( start_GTI(j + 1) - end_GTI(j) ) / dt)
          write(*,'(A, I3, A, F8.2, A, I4, A)') 'Gap number: ', j,' ', real(ee) * real(dt), ' sec,   ', ee, ' bins'
       enddo

       write(*,*)
       write(*,*) '   Set the maximum length that you want to interpolate (in sec).'
       write(*,'(A,F6.3)') '   Remember that dt is ', real(dt)

       if (verbose_merge) then 
          read(*,*) max_gap_sec
       else
          max_gap_sec = max_gap_sec_init
       endif      
       ! max_gap_sec = 2000 !*** TO RUN AUTOMATICALLY                  
       gap = int(max_gap_sec / dt)
       max_gap = gap
    else
       max_gap = gap
    endif 

    
    do j = 1, dim_GTI - 1
       ee = int( ( start_GTI(j + 1) - end_GTI(j) ) / dt)

       if (ee .ge. 1 .and. ee .le. max_gap ) then
          i = 1
          do while(time(i) .lt. end_GTI(j))
             i = i + 1
             if (i .gt. dim_lc) then
                exit
             endif
          enddo

          ! write(*,*)  i, ee + i,  time(ee+i), lc(ee+i), time(i-1), lc(i-1)

          m = (lc(ee + i) - lc(i - 1)) / (time(ee + i) - time(i - 1))
          q = lc(i - 1) - m * time(i - 1)
          if (allocated(bkg)) then
             ! write(*,*) 'bkg active'
             m_b = (bkg(ee + i) - bkg(i - 1)) / (time(ee + i) - time(i - 1))
             q_b = bkg(i - 1) - m_b * time(i - 1)
          endif
          
          do k = i, i + ee - 1
             lc_interpol = m * time(k) + q
             lc_interpol = lc_interpol * dt
             ! write(*,*) 'time interpol: ', time(k), lc_interpol
             
             if (allocated(bkg)) then 
                bkg_interpol  = m_b * time(k) + q_b
                bkg_interpol = bkg_interpol * dt
             endif
!Extracting the interpolation value from a Poisson distribution with the same mean
             poisson = funPoissonSingle(real(lc_interpol))
             lc(k) = dble(poisson) / dt
             ! write(*,*) 'extracted lc values: ', lc(k)
             if (lc(k) .lt. 0.0) lc(k) = 0.0

             lc(k) = m * time(k) + q
             if (allocated(bkg)) then
                poisson = funPoissonSingle(real(lc_interpol))
                bkg(k)  = dble(poisson) / dt
                write(*,*) 'bkg: ', bkg(k)
                if (bkg(k) .lt. 0.0) bkg(k) = 0.0
             endif
             
          enddo
          new_end_GTI(w) = end_GTI(j+1)
       else
          new_end_GTI(w) = end_GTI(j)
          w = w + 1
          new_start_GTI(w) = start_GTI(j+1)
       endif
       
    enddo
    new_end_GTI(w) = end_GTI(j)

    write(*,*) "   Number of GTI after interpolation " ,w 
    do j = 1, dim_GTI
       start_GTI(j) = new_start_GTI(j)
       end_GTI(j) = new_end_GTI(j)
       ! write(*,*)   new_start_GTI(j),new_end_GTI(j)
    enddo
    
    new_dim_GTI = w

    
  end subroutine interpol_split_silent
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
    subroutine lc_split_first(check_interval)
!! Working out the intervals in the light curve
! The best way to do that is to check with the GTI
      use dyn_lc
      implicit none
      logical, intent(OUT) :: check_interval
      
      integer :: i, j, count, count2
      logical :: check, check_jump


      if(.not. allocated(split_ind)) allocate(split_ind(dim_lc))
      split_ind = -1

      
! Counter: this first part work out how many intervals there are in the complete light curve
!      the calculation is based on the time difference
      count  = 0   !count how many elements to form an interval  
      count2 = 1  !count how many intervals
      split_ind(1) = 1
      
!First step
      ! write(*,*) 'First time of the light curve ', time(1)
      ! write(*,*) 'First time of the first GTI', start_GTI(1)
      if( (time(1) .ge. start_GTI(1))  .and. (time(1) .le. end_GTI(1)) ) then
         count = count + 1
      endif

      
      do i = 2 , dim_lc
         ! write(*,*) 'time index', i, time(i)
         check_jump = .false. 
         check = .false.
         if ((time(i) - time(i - 1)) .gt. 2.d0 * dt  ) then
            ! write(*,*) 'jump in time of the light curve'
            ! write(*,*) 'from ', time(i - 1), 'to ', time(i)
            do j = 1, dim_GTI
               if( (time(i) .ge. start_GTI(j))  .and. (time(i) .le. end_GTI(j)) ) then
                  check_jump = .true. 
                  exit
               endif
            enddo
         else 
            
            do j = 1, dim_GTI
               ! write(*,*) 'ciao ', time(i), start_GTI(j), end_GTI(j) 
               if( (time(i) .ge. start_GTI(j))  .and. (time(i) .le. end_GTI(j)) ) then
                  check = .true.
                  ! write(*,*) 'exit', check
                  exit
               endif
            enddo
         endif
            
        if (check) then
           count = count + 1
           ! write(*,*) check, count
        else 
           count = 0
           ! write(*,*) check, time(i), count
           if (check_jump) then
              count = 1
              split_ind((count2 * 2) - 1) = i
           else
              split_ind((count2 * 2) - 1) = i + 1
           endif
              ! split_ind((count2 * 2) - 1) = i + 1
        end if
        
        if (count .eq. int_len_dim) then
           split_ind( count2 * 2 ) = i  
           count2 = count2 + 1
           count = 0
           split_ind((count2 * 2) - 1) = i + 1
        end if        
     enddo

     ! write(*,*) split_ind
     
     int_number = count2 - 1 ! There is always one more than what we need 
! Check if there are at least one interval in the light curve (it can happen that the light curve has too many holes or the int_len_dim is too big)
     check_interval = .true.
     if (int_number .lt. 1 ) then
        write(*,*) '   The interval chosen is too long for the light curve.'
        if (allocated(split_ind)) deallocate(split_ind)
        check_interval = .false.
     endif

     write(*,*)  'End lc_split_first'
     
   end subroutine lc_split_first
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
    subroutine lc_split_index_time()
! This suroutine writes in lc matrix all the intervals based on the split_ind array      
      use dyn_lc
      implicit none
      integer   :: i, j, count

     do i = 1, int_number
     count = 1 
        do j = split_ind(2 * i - 1), split_ind(2 * i)
           lc_int  (i, count) = lc(j)
           time_int(i, count) = time(j)
           ! if(allocated(bkg)) bkg_int (i, count) = bkg(j)
           ! if (lc(j) .lt. 0.0) lc(j) = 0.0
           ! if (bkg(j) .lt. 0.0) bkg(j) = 0.0

           ! if (lc(j) .lt. 0.0) write(*,*) 'ATTENTION!! (in lc_split_index_time) negative lc (time, lc)', time(j), lc(j)
           ! if (bkg(j) .lt. 0.0) write(*,*) 'ATTENTION!! (in lc_split_index_time) negative bkg (time, bkg)', time(j), bkg(j)
           count = count + 1
        enddo
     enddo
   end subroutine lc_split_index_time
!---------------------------------------------------------------------!
