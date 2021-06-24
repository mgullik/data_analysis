!this subroutine loads lc from different observations and keep them separated, so the arrays have one dimention more than usual

subroutine load_lc_lag_ene_obs2()
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

  
  gti_dim_obs = 1
  
      print *, ' '
      ! print *, '   Write the name of the channel/energy bin'
      ! read(*,'(A)') filename_en_bin
      ! write(*,*) trim(filename_en_bin)
      ! print *, ' '
!Aternatively, comment the 4 lines above and use this one below         
      filename_en_bin = '/Users/gullo/Work/AGN_project/ark564/binning_en.txt'
      
      if (yes_no('   Are the energies expressed in eV? ')) then
         en_units = 1000.d0
      else
         write(*,*) '   Write the number the energies has to be divided to obtain keV'
         read(*,*) en_units
      end if

      print *, ' '

      en_num = file_line_num(filename_en_bin)

      if (.not. allocated(l_bin)    ) allocate(l_bin   (en_num))
      if (.not. allocated(r_bin)    ) allocate(r_bin   (en_num))
      if (.not. allocated(l_en_bin) ) allocate(l_en_bin(en_num))
      if (.not. allocated(r_en_bin) ) allocate(r_en_bin(en_num))

      open(55, file = filename_en_bin)
         do k = 1, en_num
            read(55, *) l_bin(k), r_bin(k)
            l_en_bin(k) = real(l_bin(k)) / en_units
            r_en_bin(k) = real(r_bin(k)) / en_units
            ! write(*,*) l_bin_ch(k), r_bin_ch(k), l_bin(k), r_bin(k)
         enddo
      close(55)
      print *, '   Number of energies in the spectrum: ', en_num
      print *, ' '
      ! write(*, '(A, F6.3, A, F6.3)') '   Energy boundaries in keV: ', (l_en_bin(1) + r_en_bin(1)) * 0.5d0 ,' - ', (l_en_bin(en_num) + r_en_bin(en_num)) * 0.5d0  
      write(*, '(A, F6.3, A, F6.3)') '   Energy boundaries in keV: ', real(l_en_bin(1)),' - ',  real(r_en_bin(en_num))  
      print *, ' '

      if (.not. allocated(ave_rate_en_obs) ) allocate(ave_rate_en_obs(en_num, obs_num))
      ! if (.not. allocated(ave_bkg_en)  ) allocate(ave_bkg_en (en_num))
      ave_rate_en_obs = 0.d0

      if (.not. allocated(int_len_dim_obs) ) allocate(int_len_dim_obs(obs_num))
      if (.not. allocated(int_number_obs ) ) allocate(int_number_obs (obs_num))
            
      if (.not. allocated(time_obs)   ) allocate(time_obs   (int_number_max, int_len_dim_max, obs_num))
      if (.not. allocated(lc_en_obs)  ) allocate(lc_en_obs   (int_number_max, int_len_dim_max, en_num, obs_num))
      
!GET ALL THE LIGHT CURVES IN ALL THE OBSERVATIONS

      do o = 1, obs_num
         print *, ' '
         ! print *, '   Write the name of the path (end with "/"):'
         ! read(*,'(A)') name_path
         ! write(*,*) 'Path read ', trim(name_path)
!Aternatively, comment the 3 lines above and use this one below         
         write(name_path, '(A,I1,A)') '/Users/gullo/Work/AGN_project/ark564/0670130', o + 1, '01/lc_lag_en_pileup_nofil/'
         write(*,*) 'Path is: (enter to continue)'
         write(*,*) trim(name_path)
         read(*,*)

         
         ! write(*,*) '   Write the prefix of the light curve filename:'
         ! read(*,'(A)') prefix_name
         ! write(*,*)  trim(prefix_name)
!Aternatively, comment the 3 lines above and use this one below       
         prefix_name = 'PN_lccorr_en'

         write(name_base, '(A, A)')  trim(name_path), trim(prefix_name)

         ! write(*,*) '   Write the extension of the light curve filename:'
         ! read(*,'(A)') name_extension
         ! write(*,*)  trim(name_extension)
!Aternatively, comment the 3 lines above and use this one below       
         name_extension = '_lag_en_pileup_nofil.lc'
         
         ! en_num = 1
         do k = 1, en_num
         
!Work out the name of the light curve 
            write(filename, '(A, I0, A, I0, A)') trim(name_base), l_bin(k), '-', r_bin(k), trim(name_extension)
            write(*,*) trim(filename)
            write(*,*)
            call extract_lc(filename)

            if (allocated(lc)) then 
               call split_lc()
            else
               write(*,*) '    !!! NO LIGHT CURVE !!! '
               write(*,*) '        STOP HERE   '
               stop
            endif

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
                  lc_en_obs (i, j, k, o) = lc_int(i, j)
                  ! bkg_en(i, j, k) = bkg_int(i, j)
                  ave_rate_en_obs(k, o) = ave_rate_en_obs(k, o) + lc_int(i, j)
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
         enddo
         
         write(*,*) '   Number of intervals in light curves '  , int_number_obs(o)
         write(*,*) '   Number of bins per interval         '  , int_len_dim_obs(o) 
         write(*,*) 
         ! starting_telescope_time = 0.d0
         do i = 1, int_number_obs(o)
            do j = 1, int_len_dim_obs(o)
               time_obs(i, j, o) = time_int(i, j) + starting_telescope_time
            enddo
         enddo
         
         do j = 1, dim_GTI

            start_GTI_obs(gti_dim_obs) = start_GTI(j) + starting_telescope_time
            end_GTI_obs(gti_dim_obs) = end_GTI(j) + starting_telescope_time
            gti_dim_obs = gti_dim_obs + 1 
         enddo

         ! write(*,*) ' ********************* STARTING TIME *************************', starting_telescope_time
         
!deallocation of the usuful array for the next observations 
         if(allocated(split_ind)) deallocate(split_ind)
         if(allocated(time_int) ) deallocate(time_int)
         if(allocated(lc_int)   ) deallocate(lc_int)
         if(allocated(time)     ) deallocate(time)
         if(allocated(lc)       ) deallocate(lc)
         if(allocated(err_rate) ) deallocate(err_rate)
         if(allocated(start_GTI)) deallocate(start_GTI)
         if(allocated(end_GTI)  ) deallocate(end_GTI  )

         gap = -1
         check_gap_num = -1
         int_number = -1
         
      enddo

      call print_lc_multi2(time_obs)

!reference band

      if (yes_no('   Do you want to use your reference band for the cross-spectrum?  If no, the code is going to calculate the reference summing all the single energy band light curves, apart from the subject band.')) then      
         print *, '   Enter the name of the reference band light cureve (with the full path if it is not in this folder). '

         read (*,'(A)') filename_ref
         print *, ' '

         call extract_lc(filename_ref)
         call split_lc()

         if(.not. allocated(lc_ref)) allocate(lc_ref(int_number, int_len_dim))
         ave_rate_ref = 0.d0
         do j = 1, int_len_dim
            do i = 1, int_number
               lc_ref(i, j) = lc_int(i, j)
               ave_rate_ref = ave_rate_ref + lc_int(i, j)
            enddo
         enddo
         ave_rate_ref = ave_rate_ref / real(int_number * int_len_dim)
      
         print *, '   The light curves have been divided in segments. '
         print *, ' '     
         write(*,*) '   Number of intervals '          , int_number
         write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
         write(*,*) '   Actual exposure (sec) '        , int_len_dim * dt * int_number
         write(*,*) '   Average count rate of the reference light curve (count/s)'  , ave_rate_ref
         print *, ' '
         print *, '   Type enter to continue'
         read(*,*)
         print *, ' '

         ! call make_pds_ref()

      else
         print *, ' '
      endif      


      
    end subroutine load_lc_lag_ene_obs2
