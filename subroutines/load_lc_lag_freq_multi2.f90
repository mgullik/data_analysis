subroutine load_lc_lag_freq_multi2()
  use dyn_lc
  use rand_fun
  implicit none 

  integer               :: i, j, o
  character (len=500)   :: filename1, filename2, name_path, name, name_lc

  
  ! call set_seed(seed)

!******************************************************************
    
      ! call execute_command_line('ls')

  ! int_len_dim  = 2048

  ! if (.not. allocated(int_len_dim_obs) ) allocate(int_len_dim_obs(obs_num))
  if (.not. allocated(int_number_obs) ) allocate(int_number_obs(obs_num))



    ! if (.not. allocated(ave_rate_freq1_obs)) allocate(ave_rate_freq1_obs(obs_num))
    ! if (.not. allocated(ave_rate_freq2_obs)) allocate(ave_rate_freq2_obs(obs_num))
    ! ave_rate_freq1_obs = 0.0
    ! ave_rate_freq2_obs = 0.0
    
    ! write(*,*) 'First light curve upload'

    do o = 1, obs_num
       ! print *, ' '
       ! print *, '   Write the name of the path (end with "/"):'
       ! read(*,'(A)') name_path
       ! write(*,*) 'Path read ', trim(name_path)

!Aternatively, comment the 3 lines above and use this one below         
         write(name_path, '(A,I1,A)') '/Users/gullo/Work/AGN_project/ark564/0670130', o + 1, '01/lc_lag_freq/'
         write(*,*) 'Path is: (enter to continue)'
         write(*,*) trim(name_path)
         read(*,*)

         
         ! write(*,*) '   Write the name of the light curve:'
         ! read(*,'(A)') name
         ! write(*,*)  trim(name)
!Aternatively, comment the 3 lines above and use this one below       
         name = 'PN_lccorr_en300-1000_lag_freq.lc'


         write(filename1, '(A, A)') trim(name_path), trim(name)
         write(*,*) trim(filename1)
         write(*,*)
         
         call extract_lc(filename1)
         ! check_gap_num = dim_GTI
         if (allocated(lc)) then 
            call split_lc()
         else
            write(*,*) '    !!! NO LIGHT CURVE !!! '
            write(*,*) '        STOP HERE   '
            stop
         endif

         write(*,'(A, I2, A)') '   Observation number ', o, ' first light cureve extracted and splitted'

         ! int_number = 1

         int_number_obs(o)  = int_number

         if (int_number_obs(o) .gt. int_number_max) then
            write(*,*) '   ERROR! The number of the intervals is larger than the limit (modify the limit in the code)'
            write(*,*) '          max length for the interval is ', int_number_max, ' bins'
            write(*,*) '          the error is for obs number: ', o, ' with number of intervals ', int_number_obs(o)  
            write(*,*) '   Exit...'
            stop
         endif
         
         ! if (.not. allocated(time_lc)) allocate(time_lc(int_number_max, int_len_dim, obs_num))
         if (.not. allocated(lc_freq1_obs)) allocate(lc_freq1_obs(int_number_max, int_len_dim, obs_num))
  
         if (.not. allocated(lc_freq2_obs)) allocate(lc_freq2_obs(int_number_max, int_len_dim, obs_num))

         write(name_lc, '(A,I1,A)') 'lc_first',o, '.lc'
         open(81, file = name_lc)
         write(81, *) 'skip on'
         do i = 1, int_number_obs(o)
            do j = 1, int_len_dim
               write(81,*) time_int(i, j), lc_int(i, j)
            enddo
            write(81,*) 'no no'
         enddo
         close(81)

         
         do j = 1, int_len_dim
            do i = 1, int_number_obs(o)
               ! time_lc  (i, j) = time_int(i, j, o)
               lc_freq1_obs(i, j, o) = lc_int(i, j)
               ! ave_rate_freq1_obs(o) = ave_rate_freq1_obs(o) + lc_int(i, j)
               ! ave_bkg_freq1  = ave_bkg_freq1  + bkg_int(i, j)

            enddo
         enddo
         ! ave_rate_freq1_obs(o) = ave_rate_freq1_obs(o)  / real(int_len_dim * int_number_obs(o))
       ! ave_bkg_freq1  = ave_bkg_freq1   / real(int_len_dim * int_number)

         write(*,*) 
         write(*,*) '********************************************'
         write(*,*) 'Second light curve upload'
         write(*,*) 

         name = 'PN_lccorr_en1200-4000_lag_freq.lc'

         write(filename2, '(A, A)') trim(name_path), trim(name)
         write(*,*) trim(filename2)
         write(*,*)

         call extract_lc(filename2)
         ! write(*,*) 'ciao'
         call split_lc()

         write(*,'(A, I2, A)') '   Observation number ', o, ' second light cureve extracted and splitted'

         write(name_lc, '(A,I1,A)') 'lc_second',o, '.lc'
         open(82, file = name_lc)
         write(82, *) 'skip on'
         do i = 1, int_number_obs(o)
            do j = 1, int_len_dim
               write(82,*) time_int(i, j), lc_int(i, j)
            enddo
            write(82,*) 'no no'
         enddo
         close(82)
         
         do j = 1, int_len_dim
            do i = 1, int_number_obs(o)
               lc_freq2_obs (i, j, o) = lc_int(i, j)
               ! bkg_freq2(i, j) = bkg_int(i, j)
               ! ave_rate_freq2_obs(o) = ave_rate_freq2_obs(o) + lc_int(i, j)
               ! ave_bkg_freq2  = ave_bkg_freq2  + bkg_int(i, j)
            enddo
         enddo
         
         ! ave_rate_freq2_obs(o) = ave_rate_freq2_obs(o)  / real(int_len_dim * int_number_obs(o))
         ! ave_bkg_freq2  = ave_bkg_freq2   / real(int_len_dim * int_number)

         write(*,*) 
         write(*,*) '*************************************'
         write(*,*) '   Light curve NAMEs:    '
         write(*,*) trim(filename1), trim(filename2)
         write(*,*) '   Number of intervals '          , int_number_obs(o)
         write(*,*) '   Exposure (sec) '        , int_len_dim * dt * int_number_obs(o)
         ! write(*,*) '   Average count rate (count/s) in both lc:'
         ! write(*,*)  ave_rate_freq1_obs(o), ave_rate_freq2_obs(o)
         write(*,*) '*************************************'
         write(*,*)

         if(allocated(lc      )) deallocate(lc  )
         if(allocated(time    )) deallocate(time)
         if(allocated(bkg     )) deallocate(bkg )
         if(allocated(lc_int  )) deallocate(lc_int)
         if(allocated(time_int)) deallocate(time_int)
         if(allocated(err_rate)) deallocate(err_rate)

         if(allocated(start_GTI)) deallocate(start_GTI)
         if(allocated(end_GTI)  ) deallocate(end_GTI  )
         if(allocated(split_ind)) deallocate(split_ind)

         gap = -1
         check_gap_num = -1
         int_number = -1

      enddo

    ! if(allocated(lc_int)   ) deallocate(lc_int   )
    ! if(allocated(time_int) ) deallocate(time_int )
    ! if(allocated(bkg_int)  ) deallocate(bkg_int  )
    ! gap = -1 
    ! check_gap_num = -1

  end subroutine load_lc_lag_freq_multi2
