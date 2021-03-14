subroutine load_lc_lag_freq(filename1, filename2)
  use dyn_lc
  use rand_fun
  implicit none 

  integer               :: i, j
  character (len=500)   :: filename1, filename2

  
  ! call set_seed(seed)

!******************************************************************
    
      ! call execute_command_line('ls')

  ave_rate_freq1 = 0.0
  ave_rate_freq2 = 0.0
      
    write(*,*) 'First light curve upload'
    
    call extract_lc(filename1)
    if (allocated(lc)) then 
       call split_lc()
    else
       write(*,*) '    !!! NO LIGHT CURVE !!! '
       write(*,*) '        STOP HERE   '
       stop
    endif
    
    ! if (.not. allocated(time_lc)) allocate(time_lc(int_number, int_len_dim))
    if (.not. allocated(lc_freq1)) allocate(lc_freq1(int_number, int_len_dim))

    do j = 1, int_len_dim
       do i = 1, int_number
          ! time_lc  (i, j) = time_int(i, j)
          lc_freq1 (i, j) = lc_int(i, j)
          ! bkg_freq1(i, j) = bkg_int(i, j)
          ave_rate_freq1 = ave_rate_freq1 + lc_int(i, j)
          ! ave_bkg_freq1  = ave_bkg_freq1  + bkg_int(i, j)
       enddo
    enddo
    
    ave_rate_freq1 = ave_rate_freq1  / real(int_len_dim * int_number)
    ! ave_bkg_freq1  = ave_bkg_freq1   / real(int_len_dim * int_number)
    
    write(*,*) 
    write(*,*) '********************************************'
    write(*,*) 'Second light curve upload'
    write(*,*) 
    
    call extract_lc(filename2)
    call split_lc()

    
    if (.not. allocated(lc_freq2)) allocate(lc_freq2(int_number, int_len_dim))

    do j = 1, int_len_dim
       do i = 1, int_number
          lc_freq2 (i, j) = lc_int(i, j)
          ! bkg_freq2(i, j) = bkg_int(i, j)
          ave_rate_freq2 = ave_rate_freq2 + lc_int(i, j)
          ! ave_bkg_freq2  = ave_bkg_freq2  + bkg_int(i, j)
       enddo
    enddo
    
    ave_rate_freq2 = ave_rate_freq2  / real(int_len_dim * int_number)
    ! ave_bkg_freq2  = ave_bkg_freq2   / real(int_len_dim * int_number)
  
         write(*,*) 
         write(*,*) '*************************************'
         write(*,*) '   Light curve NAMEs:    '
         write(*,*) trim(filename1), trim(filename2)
         write(*,*) '   Number of intervals '          , int_number
         write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
         write(*,*) '   Exposure (sec) '        , int_len_dim * dt * int_number
         write(*,*) '   Average count rate (count/s) in both lc:'
         write(*,*)  ave_rate_freq1, ave_rate_freq2
         write(*,*) '*************************************'
         write(*,*)

    if(allocated(lc)   ) deallocate(lc  )
    if(allocated(time) ) deallocate(time)
    if(allocated(bkg)  ) deallocate(bkg )

    if(allocated(start_GTI)) deallocate(start_GTI)
    if(allocated(end_GTI)  ) deallocate(end_GTI  )
    if(allocated(split_ind)) deallocate(split_ind)

    ! if(allocated(lc_int)   ) deallocate(lc_int   )
    ! if(allocated(time_int) ) deallocate(time_int )
    ! if(allocated(bkg_int)  ) deallocate(bkg_int  )
    ! gap = -1 
    ! check_gap_num = -1

end subroutine load_lc_lag_freq
