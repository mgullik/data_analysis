subroutine make_freq_intervals()
  use dyn_lc
  implicit none

  integer              :: jj
  double precision     :: min_freq, max_freq
  double precision , allocatable :: fq_min_arr(:), fq_max_arr(:), fq(:)
  character (len=200)  :: filename 

  filename = 'products/ener_freq_info.txt'

  open(12, file = filename)
   
  write(*,*)
  write(*,*) '   Calculate the frequency intervals for the lag spectrum'
  write(*,*)
  write(*,*) '   Enter the minimum frequency in Hz '
  read (*,*)   min_freq
  write(*,*) '   Enter the maximum frequency in Hz '
  read (*,*)   max_freq
   ! min_freq = 2.9
   ! max_freq = 30.0
  write(*,*) '   Enter the number of frequency intervals (they will be log spaced) '
  read (*,*)    freq_num
   
   if(.not. allocated(upper_fq  )) allocate(upper_fq(freq_num))
   if(.not. allocated(lower_fq  )) allocate(lower_fq(freq_num))
      
   if(.not. allocated(fq_min_arr)) allocate(fq_min_arr(freq_num))
   if(.not. allocated(fq_max_arr)) allocate(fq_max_arr(freq_num))
   if(.not. allocated(fq)        ) allocate(fq        (freq_num + 1))

!log shape for the frequency ranges 
   call make_array_logspace(min_freq, max_freq, freq_num + 1, fq)
   
   write(12,*) '   FREQUNECY INTERVALS INFORMATION for the lag vs energy spectra'
   write(12,*) '   User input:'   
   write(12,*) '   Minimum frequency: ', min_freq
   write(12,*) '   Maximum frequency: ', max_freq
   write(12,*) '   Number of frequency intervals: ', freq_num
   write(12,*)   
   write(12, *) '   Frequency intervals calculated from the boundaries give by the user'
   write(12,*)
   do jj = 1, freq_num
      fq_min_arr(jj) = fq(jj) 
      fq_max_arr(jj) = fq(jj + 1)
      write(12,*) jj, fq_min_arr(jj), fq_max_arr(jj)
   enddo
   write(12,*)

    
   do jj = 1, freq_num 
      lower_fq(jj) = ceiling(fq_min_arr(jj) * real(int_len_dim) * dt)
      upper_fq(jj) = floor  (fq_max_arr(jj) * real(int_len_dim) * dt)
      write(12,*) '----------------------------------------------'
      write(12,*) "   Frequency interval: ", jj
      write(12,*) "   Min frequency bin       = ", lower_fq(jj)
      write(12,*) "   Max frequency bin       = ", upper_fq(jj)
      write(12,*) "   Min true frequency (Hz) = ", lower_fq(jj) / (real(int_len_dim) * dt) 
      write(12,*) "   Max true frequency (Hz) = ", upper_fq(jj) / (real(int_len_dim) * dt)
      write(12,*) '   Number of frequency bins in this interval ', upper_fq(jj) - lower_fq(jj) + 1
      write(12,*) '----------------------------------------------'
      write(12,*)         
   enddo

   if (.not. allocated(relevant_freq)) allocate(relevant_freq(freq_num))

   do jj = 1, freq_num       
      relevant_freq(jj) = (fq_max_arr(jj) + fq_min_arr(jj)) * 0.5d0 
   enddo

   write(*,*) '   Frequnecy intervals have been calculated'
   write(*,*) '   Informations are in the file: ', trim(filename)
   write(*,*) 
   write(*,*) '   Press Enter to continue'
   read (*,*) 
   write(*,*) 
        
   deallocate(fq_min_arr)
   deallocate(fq_max_arr)
   deallocate(fq        )

 end subroutine make_freq_intervals





 subroutine make_freq_intervals_obs()
  use dyn_lc
  implicit none

  integer              :: jj, o
  double precision     :: min_freq, max_freq
  double precision , allocatable :: fq_min_arr(:), fq_max_arr(:), fq(:)
  character (len=200)  :: filename 

  filename = 'ener_freq_obs_info.txt'

  open(12, file = filename)
   
  write(*,*)
  write(*,*) '   Calculate the frequency intervals for the lag spectrum'
  write(*,*)
  write(*,*) '   Enter the minimum frequency in Hz '
  read (*,*)   min_freq
  write(*,*) '   Enter the maximum frequency in Hz '
  read (*,*)   max_freq
   ! min_freq = 2.9
   ! max_freq = 30.0
  write(*,*) '   Enter the number of frequency intervals (they will be log spaced) '
  read (*,*)    freq_num
   
   if(.not. allocated(upper_fq_obs)) allocate(upper_fq_obs(freq_num, obs_num))
   if(.not. allocated(lower_fq_obs)) allocate(lower_fq_obs(freq_num, obs_num))
      
   if(.not. allocated(fq_min_arr)  ) allocate(fq_min_arr(freq_num))
   if(.not. allocated(fq_max_arr)  ) allocate(fq_max_arr(freq_num))
   if(.not. allocated(fq)          ) allocate(fq        (freq_num + 1))

!log shape for the frequency ranges 
   call make_array_logspace(min_freq, max_freq, freq_num + 1, fq)
   
   write(12,*) '   FREQUNECY INTERVALS INFORMATION for the lag vs energy spectra'
   write(12,*) '   User input:'   
   write(12,*) '   Minimum frequency: ', min_freq
   write(12,*) '   Maximum frequency: ', max_freq
   write(12,*) '   Number of frequency intervals: ', freq_num
   write(12,*)   
   write(12, *) '   Frequency intervals calculated from the boundaries give by the user'
   write(12,*)
   do jj = 1, freq_num
      fq_min_arr(jj) = fq(jj) 
      fq_max_arr(jj) = fq(jj + 1)
      write(12,*) jj, fq_min_arr(jj), fq_max_arr(jj)
   enddo
   write(12,*)

   do o = 1, obs_num 
   
      do jj = 1, freq_num 
         lower_fq_obs(jj, o) = ceiling(fq_min_arr(jj) * real(int_len_dim_obs(o)) * dt)
         upper_fq_obs(jj, o) = floor  (fq_max_arr(jj) * real(int_len_dim_obs(o)) * dt)
         write(12,*) '----------------------------------------------'
         write(12,*) "   Frequency interval: ", jj
         write(12,*) "   Min frequency bin       = ", lower_fq_obs(jj,o)
         write(12,*) "   Max frequency bin       = ", upper_fq_obs(jj,o)
         write(12,*) "   Min true frequency (Hz) = ", lower_fq_obs(jj,o) / (real(int_len_dim_obs(o)) * dt) 
         write(12,*) "   Max true frequency (Hz) = ", upper_fq_obs(jj,o) / (real(int_len_dim_obs(o)) * dt)
         write(12,*) '   Number of frequency bins in this interval ', upper_fq_obs(jj,o) - lower_fq_obs(jj,o) + 1
         write(12,*) '----------------------------------------------'
         write(12,*)         
      enddo


   enddo

   if (.not. allocated(relevant_freq)) allocate(relevant_freq(freq_num))
   do jj = 1, freq_num       
      relevant_freq(jj) = (fq_max_arr(jj) + fq_min_arr(jj)) * 0.5d0 
   enddo
   
   write(*,*) '   Frequnecy intervals have been calculated for all the observations'
   write(*,*) '   Informations are in the file: ', trim(filename)
   write(*,*) 
   write(*,*) '   Press Enter to continue'
   read (*,*) 
   write(*,*) 
        
   deallocate(fq_min_arr)
   deallocate(fq_max_arr)
   deallocate(fq        )
   

 end subroutine make_freq_intervals_obs
