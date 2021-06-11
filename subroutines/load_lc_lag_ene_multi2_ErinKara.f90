!this subroutine loads lc from different observations and keep them separated, so the arrays have one dimention more than usual

subroutine load_lc_lag_ene_obs2_ek()
  use dyn_lc
  implicit none

  integer              :: i, j, k
  double precision     :: en_units
  character (len = 500)  :: name_path,  name_base, filename, &
       filename_en_bin, name_extension, filename_ref, prefix_name
   logical              :: yes_no

  integer              :: file_line_num
  integer, allocatable :: l_bin(:), r_bin(:)

  integer                       :: gti_dim_obs, o
  double precision, allocatable :: ave_rate_en_obs(:,:)
  double precision, allocatable :: time_obs(:,:,:)
  double precision              :: start_GTI_obs(max_GTI_dim), end_GTI_obs(max_GTI_dim)

  check_merge = .true.
  
  gti_dim_obs = 1
  
  obs_num = 8
      print *, ' '

      dt = 10.d0
      en_num = 35

      if (.not. allocated(l_bin)    ) allocate(l_bin   (en_num))
      if (.not. allocated(r_bin)    ) allocate(r_bin   (en_num))
      if (.not. allocated(l_en_bin) ) allocate(l_en_bin(en_num))
      if (.not. allocated(r_en_bin) ) allocate(r_en_bin(en_num))

      
      l_en_bin = (/0.3, 0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1, &
              1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0/)
      r_en_bin = (/     0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1, &
              1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0/)

      
      print *, '   Number of energies in the spectrum: ', en_num
      print *, ' '

      if (.not. allocated(ave_rate_en_obs) ) allocate(ave_rate_en_obs(en_num, obs_num))
      ! if (.not. allocated(ave_bkg_en)  ) allocate(ave_bkg_en (en_num))
      ave_rate_en_obs = 0.d0

      if (.not. allocated(int_len_dim_obs) ) allocate(int_len_dim_obs(obs_num))
      if (.not. allocated(int_number_obs ) ) allocate(int_number_obs (obs_num))
            
      if (.not. allocated(time_obs)   ) allocate(time_obs   (1, int_len_dim_max, obs_num))
      if (.not. allocated(lc_en_obs)  ) allocate(lc_en_obs   (1, int_len_dim_max, en_num, obs_num))
      
!GET ALL THE LIGHT CURVES IN ALL THE OBSERVATIONS
      print *, ' '
      ! print *, '   Write the name of the path (end with "/"):'
      ! read(*,'(A)') name_path
      name_path = '/Users/gullo/Dropbox/forGullo/'

      write(*,*) 'Path read ', trim(name_path)

      do o = 1, obs_num
!Aternatively, comment the 3 lines above and use this one below         
         write(filename, '(A,A,I1,A)') trim(name_path),'lc_obs', o, '.dat'
         write(*,*) 'Extracting light curve: (enter to continue)'
         write(*,*) trim(filename)
         read(*,*)

         
         call extract_split_erin(filename, o)

            
         int_len_dim_obs(o) = int_len_dim
         int_number_obs(o)  = int_number

         
         if (int_number_obs(o) .gt. int_number_max) then
            write(*,*) '   ERROR! The number of the intervals is larger than the limit (modify the limit in the code)'
            write(*,*) '          max length for the interval is ', int_number_max, ' bins'
            write(*,*) '          the error is for obs number: ', o, ' with number of intervals ', int_number_obs(o)  
            write(*,*) '   Exit...'
            stop
         endif
         if (int_len_dim_obs(o) .gt. int_len_dim_max) then
            write(*,*) '   ERROR! The length of the interval is larger than the limit (modify the limit in the code)'
            write(*,*) '          max length for the interval is ', int_len_dim_max, ' bins'
            write(*,*) '          the error is for obs number: ', o, ' with bin length ', int_len_dim_obs(o)  
            write(*,*) '   Exit...'
            stop
         endif

         ! if (.not. allocated(bkg_en)  ) allocate(bkg_en  (int_number, int_len_dim, en_num))

         
         
         do j = 1, int_len_dim_obs(o)
            do i = 1, int_number_obs(o)
               ! lc_en_obs (i, j, k, o) = lc_int(i, j)
               ! bkg_en(i, j, k) = bkg_int(i, j)
               ave_rate_en_obs(k, o) = ave_rate_en_obs(k, o) + lc_en_obs (i, j, k, o)
               ! ave_bkg_en(k)  = ave_bkg_en(k)  + bkg_int(i, j)
            enddo
         enddo
         ave_rate_en_obs(k, o) = ave_rate_en_obs(k, o)  / real(int_len_dim * int_number)
         ! ave_bkg_en (k) = ave_bkg_en (k)  / real(int_len_dim * int_number)

         
         write(*,*) '*************************************'
         write(*,*) '   Light curve NAME     '               , trim(filename)
         write(*,*) '   Average count rate (count/s) '       , ave_rate_en_obs(k, o)            
         write(*,*) '*************************************'
         write(*,*) 

         write(*,*) '   Number of intervals in light curves '  , int_number_obs(o)
         write(*,*) '   Number of bins per interval         '  , int_len_dim_obs(o) 
         write(*,*) 

         ! do i = 1, int_number_obs(o)
         !    do j = 1, int_len_dim_obs(o)
         !       time_obs(i, j, o) = time_int(i, j) + starting_telescope_time
         !    enddo
         ! enddo

         ! do j = 1, dim_GTI

         !    start_GTI_obs(gti_dim_obs) = start_GTI(j) + starting_telescope_time
         !    end_GTI_obs(gti_dim_obs) = end_GTI(j) + starting_telescope_time
         !    gti_dim_obs = gti_dim_obs + 1 
         ! enddo

         !deallocation of the usuful array for the next observations 
         ! if(allocated(split_ind)) deallocate(split_ind)
         ! if(allocated(time_int) ) deallocate(time_int)
         ! if(allocated(lc_int)   ) deallocate(lc_int)
         ! if(allocated(time)     ) deallocate(time)
         ! if(allocated(lc)       ) deallocate(lc)
         ! if(allocated(err_rate) ) deallocate(err_rate)
         ! if(allocated(start_GTI)) deallocate(start_GTI)
         ! if(allocated(end_GTI)  ) deallocate(end_GTI  )

         ! gap = -1
         ! check_gap_num = -1
         ! int_number = -1


      enddo

      !reference band
      !so far it's calculated summing all the energy bands and subtracting the subject band 

      ! if (yes_no('   Do you want to use your reference band for the cross-spectrum?  If no, the code is going to calculate the reference summing all the single energy band light curves, apart from the subject band.')) then      
      !    print *, '   Enter the name of the reference band light cureve (with the full path if it is not in this folder). '

      !    read (*,'(A)') filename_ref
      !    print *, ' '

      !    call extract_lc(filename_ref)
      !    call split_lc()

      !    if(.not. allocated(lc_ref)) allocate(lc_ref(int_number, int_len_dim))
      !    ave_rate_ref = 0.d0
      !    do j = 1, int_len_dim
      !       do i = 1, int_number
      !          lc_ref(i, j) = lc_int(i, j)
      !          ave_rate_ref = ave_rate_ref + lc_int(i, j)
      !       enddo
      !    enddo
      !    ave_rate_ref = ave_rate_ref / real(int_number * int_len_dim)
      
      !    print *, '   The light curves have been divided in segments. '
      !    print *, ' '     
      !    write(*,*) '   Number of intervals '          , int_number
      !    write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
      !    write(*,*) '   Actual exposure (sec) '        , int_len_dim * dt * int_number
      !    write(*,*) '   Average count rate of the reference light curve (count/s)'  , ave_rate_ref
      !    print *, ' '
      !    print *, '   Type enter to continue'
      !    read(*,*)
      !    print *, ' '

      !    ! call make_pds_ref()

      ! else
      !    print *, ' '
      ! endif

      
    end subroutine load_lc_lag_ene_obs2_ek
