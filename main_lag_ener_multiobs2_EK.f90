include 'header.h'

!This main is specific to analyse Erin's light curves 

program cross
  use dyn_lc
  use rand_fun
  implicit none 

  double precision, parameter   :: pi = acos(-1.0)
  integer              :: i, j, k, kk, jj, o, en_ind
  integer, allocatable :: num_freq_bins(:)
  double precision, allocatable :: rc_en_obs(:,:,:,:), ic_en_obs(:,:,:,:), &
                          pw_en_obs(:,:,:,:), pw_ref_en_obs(:,:,:,:), &
                          lc_ref_cross_obs(:,:,:,:)

  double precision, allocatable :: rc_avefq_en(:,:),ic_avefq_en(:,:),&
                          rc2_avefq_en(:,:), ic2_avefq_en(:,:), &
                          rc_ic_avefq_en(:,:), pw_avefq_en(:,:), &
                          pw_avefq_en_ref(:,:)
  double precision, allocatable :: var_rc(:,:), var_ic(:,:), deriv_rc(:,:), &
                          deriv_ic(:,:), covariance(:,:), &
                          err_std_rc(:,:), err_std_ic(:,:), &
                          coher2_en(:,:), err_Aform_rc_ic(:,:), &
                          err_Aform_lag(:,:), err_prop_lag(:,:), &
                          err_cohe_lag(:,:), lag_en(:,:)
!P_noise
  double precision              :: freq_limit, bias2
  double precision              :: get_p_noise_pw  !function
  double precision, allocatable :: P_noise_ext_obs(:,:), P_noise_ext_ref_obs(:,:), &
                          P_noise_ext(:), P_noise_ext_ref(:), &
                          pw_aveint(:), pw_aveint_ref(:)

!NOT power of 2 segments
  double precision, allocatable :: re_lc(:), im_lc(:), re_ref(:), im_ref(:)
  
!******************************************************************
    
!GET LIGHT CURVEs 

  ! call load_lc_lag_ene_obs2()
  call load_lc_lag_ene_obs2_ek()

!Arrays to save the fourier trasform analysis for every energy       
      if(.not. allocated(rc_en_obs)    ) allocate(rc_en_obs    (int_number_max, int_len_dim_max / 2, en_num, obs_num))
      if(.not. allocated(ic_en_obs)    ) allocate(ic_en_obs    (int_number_max, int_len_dim_max / 2, en_num, obs_num))
      if(.not. allocated(pw_en_obs)    ) allocate(pw_en_obs    (int_number_max, int_len_dim_max / 2, en_num, obs_num))
      if(.not. allocated(pw_ref_en_obs)) allocate(pw_ref_en_obs(int_number_max, int_len_dim_max / 2, en_num, obs_num))

      print *, ' -----------------------------------------------'
      print *, '   Starting the cross spectrum analysis'


      if (.not. allocated(lc_ref_cross_obs)) allocate(lc_ref_cross_obs(int_number_max, int_len_dim_max, en_num, obs_num))
      lc_ref_cross_obs = 0.d0

      do o = 1, obs_num
      
!Calculate the reference band for every light curve 
         ! if (.not. allocated(lc_ref)) then          
!In this case the reference light curve is the sum of all the light curve but the subject band 
         write(*,*)
         write(*,*) '   Reference light curve done by adding all the lc but subject one'
         do k = 1, en_num             
            do i = 1, int_number_obs(o)
               do kk = 1, en_num 
                  do j = 1, int_len_dim_obs(o)
                     if (kk .ne. k) then 
                        lc_ref_cross_obs(i, j, k, o) = lc_ref_cross_obs(i, j, k, o) + lc_en_obs(i, j, kk, o)
                     endif
                  enddo
               enddo
            enddo
         enddo
         ! else
!In this case the reference light curve is the loaded reference light curve subtracted by the subject light curve
            
            ! write(*,*) '   Reference light curve done by subtracting the subject lc to the loaded one'
            ! do k = 1, en_num             
            !    do i = 1, int_number_obs(o)
            !       do j = 1, int_len_dim_obs(o)
            !          lc_ref_cross(i, j, k, o) =  lc_ref(i, j) - lc_en_obs(i, j, k, o)
            !       enddo
            !    enddo
            ! enddo
            
         ! endif

      
         if (.not. check_power2) then
            allocate(re_lc (int_len_dim_obs(o) / 2))
            allocate(im_lc (int_len_dim_obs(o) / 2))
            allocate(re_ref(int_len_dim_obs(o) / 2))
            allocate(im_ref(int_len_dim_obs(o) / 2))
         endif

         do k = 1, en_num             
            do i = 1, int_number_obs(o)
               if (check_power2) then 
                  
                  ! call ncperiodogram(lc_en_obs(i, :, k, o), lc_ref_cross_obs(i, :, k, o), rc_en_obs(i, :, k, o), ic_en_obs(i, :, k, o), int_len_dim_obs(o)) 
                  ! call periodogram(lc_en_obs(i, :, k, o), pw_en_obs(i, :, k, o), int_len_dim_obs(o))
                  ! call periodogram(lc_ref_cross_obs(i, :, k, o), pw_ref_en_obs(i, :, k, o), int_len_dim_obs(o))
                  
                  ! call ncperiodogram_frac_rms(lc_en_obs(i, :, k, o), lc_ref_cross_obs(i, :, k, o), rc_en_obs(i, :, k, o), ic_en_obs(i, :, k, o), int_len_dim_obs(o)) 
                  ! call periodogram_frac_rms(lc_en_obs(i, :, k, o), pw_en_obs(i, :, k, o), int_len_dim_obs(o))
                  ! call periodogram_frac_rms(lc_ref_cross_obs(i, :, k, o), pw_ref_en_obs(i, :, k, o), int_len_dim_obs(o))
                  
                  call ncperiodogram_no_norm(lc_en_obs(i, :, k, o), lc_ref_cross_obs(i, :, k, o), rc_en_obs(i, :, k, o), ic_en_obs(i, :, k, o), int_len_dim_obs(o)) 
                  call periodogram_no_norm(lc_en_obs(i, :, k, o), pw_en_obs(i, :, k, o), int_len_dim_obs(o))
                  call periodogram_no_norm(lc_ref_cross_obs(i, :, k, o), pw_ref_en_obs(i, :, k, o), int_len_dim_obs(o))

                  ! write(*,*) 'FFT done', i
               else
                  ! write(*,*) '   The length of the intervals is not a power of 2 ', 

                  !Calculate the Fourier transform of the light curves and reference bands
                  call FT_not_fast(lc_en_obs(i, :, k, o), re_lc, im_lc, int_len_dim_obs(o))
                  call FT_not_fast(lc_ref_cross_obs(i, :, k, o), re_ref, im_ref, int_len_dim_obs(o))

                  !Calculate Power spectra and cross spectrum of light curves and ref band
                  do jj = 1, int_len_dim_obs(o) / 2    
                     rc_en_obs    (i, jj, k, o) = (re_lc(jj) * re_ref(jj)) + (im_lc(jj) * im_ref(jj))
                     ic_en_obs    (i, jj, k, o) = (im_lc(jj) * re_ref(jj)) - (re_lc(jj) * im_ref(jj))
                     pw_en_obs    (i, jj, k, o) = (re_lc(jj) * re_lc(jj)) + (im_lc(jj) * im_lc(jj)) 
                     pw_ref_en_obs(i, jj, k, o) = (re_ref(jj) * re_ref(jj)) + (im_ref(jj) * im_ref(jj))  
                  enddo
                  write(*,'(A,I4,A,I8)') '   Slow Fourier transform done for interval ', i, ' with interval length ', int_len_dim_obs(o) 

               endif

            enddo ! loop over intervarls (int_number)

            write(*,'(A, F6.2,A,F6.2,A)') '   Finished calculating the cross spectrum between the reference band and the energy band between: ', l_en_bin(k), ' - ', r_en_bin(k), ' keV' 
            write(*,*)
         enddo !loop over enery bins (en_num)

         if(allocated(re_lc )) deallocate(re_lc )
         if(allocated(im_lc )) deallocate(im_lc )
         if(allocated(re_ref)) deallocate(re_ref)
         if(allocated(im_ref)) deallocate(im_ref)

      enddo

      if (allocated(lc_en_obs)       ) deallocate(lc_en_obs)
      if (allocated(lc_ref_cross_obs)) deallocate(lc_ref_cross_obs)

      write(*,*) '   Press Enter to continue'
      read(*,*)

      write(*,*)
      write(*,*) '   Start lag calculation'
      write(*,*)
      
!Compute the poisson noise for every pw of each energy band and ref band  
      if( .not. allocated(P_noise_ext_obs    )) allocate(P_noise_ext_obs    (en_num, obs_num))
      if( .not. allocated(P_noise_ext_ref_obs)) allocate(P_noise_ext_ref_obs(en_num, obs_num))
      if( .not. allocated(P_noise_ext    ))     allocate(P_noise_ext        (en_num))
      if( .not. allocated(P_noise_ext_ref))     allocate(P_noise_ext_ref    (en_num))
      P_noise_ext     = 0.d0
      P_noise_ext_ref = 0.d0
      
      call make_freq_array_obs()
      write(*,*) '   Above which frequency do you want to calculate the poisson noise?'
      read (*,*) freq_limit
      write(*,*)
      ! write(*,'(A,F5.1)') '   Calculating the Poisson noise for every energy PDS. Above frequency  ', freq_limit
      ! write(*,*)

      ! write(21,*) 'skip on'
 
      do o = 1, obs_num
         if( .not. allocated(pw_aveint    )) allocate(pw_aveint    (int_len_dim_obs(o)))
         if( .not. allocated(pw_aveint_ref)) allocate(pw_aveint_ref(int_len_dim_obs(o)))

         do k = 1, en_num
            pw_aveint     = 0.d0
            pw_aveint_ref = 0.d0
            do i = 1, int_number_obs(o)
               do j = 1, int_len_dim_obs(o)
                  ! write(21,*) freq_obs(j,o), pw_en_obs(i, j, k, o)
                  pw_aveint    (j) = pw_aveint    (j) + pw_en_obs    (i, j, k, o)
                  pw_aveint_ref(j) = pw_aveint_ref(j) + pw_ref_en_obs(i, j, k, o)
               enddo
               ! write(21,*) 'no no'
            enddo
            pw_aveint     = pw_aveint     / dble(int_number_obs(o))
            pw_aveint_ref = pw_aveint_ref / dble(int_number_obs(o))

            P_noise_ext_obs    (k, o) = get_p_noise_pw(pw_aveint    , int_len_dim_obs(o) / 2, freq_obs(:,o), freq_limit )
            P_noise_ext_ref_obs(k, o) = get_p_noise_pw(pw_aveint_ref, int_len_dim_obs(o) / 2, freq_obs(:,o), freq_limit )
            write(*,'(A, I4, A, I4)') '   Poisson noise level from obs ', o, ' and energy bin ', k
            write(*,*) P_noise_ext_obs(k, o), P_noise_ext_ref_obs(k, o)
            P_noise_ext    (k) = P_noise_ext    (k) + P_noise_ext_obs    (k, o)
            P_noise_ext_ref(k) = P_noise_ext_ref(k) + P_noise_ext_ref_obs(k, o)
         enddo
         deallocate(pw_aveint    )
         deallocate(pw_aveint_ref)
         write(*,*)
         ! P_noise_ext     = 0.d0
         ! P_noise_ext_ref = 0.d0

      enddo

      write(*,*)

      call make_freq_intervals_obs()

!-----------------------------------------------------------------------------------------------------
!Energy rebinning (only factor of 2) and the last 3 energy bins together
!-----------------------------------------------------------------------------------------------------
      en_num = (en_num - 1) / 2  ! half of the number of energy bins 

!change the energy bins
      do k = 1, en_num - 1
         en_ind = (2 * k) - 1
         l_en_bin(k) = l_en_bin(en_ind)
         r_en_bin(k) = r_en_bin(en_ind + 1)
      enddo
!Last energy bin
      k = en_num
      en_ind = (2 * k) - 1
      l_en_bin(k) = l_en_bin(en_ind)
      r_en_bin(k) = r_en_bin(en_ind + 2)

      
      if(.not. allocated(rc_avefq_en    )) allocate(rc_avefq_en    (freq_num, en_num))
      if(.not. allocated(ic_avefq_en    )) allocate(ic_avefq_en    (freq_num, en_num))
      if(.not. allocated(rc2_avefq_en   )) allocate(rc2_avefq_en   (freq_num, en_num))
      if(.not. allocated(ic2_avefq_en   )) allocate(ic2_avefq_en   (freq_num, en_num))
      if(.not. allocated(rc_ic_avefq_en )) allocate(rc_ic_avefq_en (freq_num, en_num))
      if(.not. allocated(pw_avefq_en    )) allocate(pw_avefq_en    (freq_num, en_num))
      if(.not. allocated(pw_avefq_en_ref)) allocate(pw_avefq_en_ref(freq_num, en_num))
      rc_avefq_en     = 0.d0
      ic_avefq_en     = 0.d0
      rc2_avefq_en    = 0.d0
      ic2_avefq_en    = 0.d0
      rc_ic_avefq_en  = 0.d0
      pw_avefq_en     = 0.d0 
      pw_avefq_en_ref = 0.d0

      
      do o = 1, obs_num
         
         do k = 1, en_num - 1
            do jj = 1, freq_num
               do j = lower_fq_obs(jj, o), upper_fq_obs(jj, o)
                  do i = 1, int_number_obs(o)
                     en_ind = (2 * k) - 1
                     rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) + rc_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) 
                     ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) + ic_en_obs(i, j, en_ind, o) + ic_en_obs(i, j, en_ind + 1, o)
                     rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) + rc_en_obs(i, j, en_ind, o) * rc_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) * rc_en_obs(i, j, en_ind + 1, o)
                     ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) + ic_en_obs(i, j, en_ind, o) * ic_en_obs(i, j, en_ind, o) + ic_en_obs(i, j, en_ind + 1, o) * ic_en_obs(i, j, en_ind + 1, o)
                     rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) + rc_en_obs(i, j, en_ind, o) * ic_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) * ic_en_obs(i, j, en_ind + 1, o)
                     pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) + pw_en_obs(i, j, en_ind, o) + pw_en_obs(i, j, en_ind + 1, o)
                     pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) + pw_ref_en_obs(i, j, en_ind, o) + pw_ref_en_obs(i, j, en_ind + 1, o)
                  enddo
               enddo
            enddo
         enddo
         
!last 3 energy bins averaged together
!----------------------------------------------------------------------------------------------------
         k = en_num
         do jj = 1, freq_num
            do j = lower_fq_obs(jj, o), upper_fq_obs(jj, o)
               do i = 1, int_number_obs(o)
                  en_ind = (2 * k) - 1
                  rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) + rc_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) + rc_en_obs(i, j, en_ind + 2, o) 
                  ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) + ic_en_obs(i, j, en_ind, o) + ic_en_obs(i, j, en_ind + 1, o) + ic_en_obs(i, j, en_ind + 2, o)
                  rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) + rc_en_obs(i, j, en_ind, o) * rc_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) * rc_en_obs(i, j, en_ind + 1, o) + rc_en_obs(i, j, en_ind + 2, o) * rc_en_obs(i, j, en_ind + 2, o)
                  ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) + ic_en_obs(i, j, en_ind, o) * ic_en_obs(i, j, en_ind, o) + ic_en_obs(i, j, en_ind + 1, o) * ic_en_obs(i, j, en_ind + 1, o) + ic_en_obs(i, j, en_ind + 2, o) * ic_en_obs(i, j, en_ind + 2, o)
                  rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) + rc_en_obs(i, j, en_ind, o) * ic_en_obs(i, j, en_ind, o) + rc_en_obs(i, j, en_ind + 1, o) * ic_en_obs(i, j, en_ind + 1, o) + rc_en_obs(i, j, en_ind + 2, o) * ic_en_obs(i, j, en_ind + 2, o)
                  pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) + pw_en_obs(i, j, en_ind, o) + pw_en_obs(i, j, en_ind + 1, o) + pw_en_obs(i, j, en_ind + 2, o)
                  pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) + pw_ref_en_obs(i, j, en_ind, o) + pw_ref_en_obs(i, j, en_ind + 1, o) + pw_ref_en_obs(i, j, en_ind + 2, o)
               enddo
            enddo
         enddo
         
         
      enddo
   
      if(.not. allocated(num_freq_bins)) allocate(num_freq_bins(freq_num))
      num_freq_bins = 0
      do o = 1, obs_num
         do jj = 1, freq_num
            num_freq_bins(jj) = num_freq_bins(jj) + ((upper_fq_obs(jj, o) - lower_fq_obs(jj, o) + 1) * int_number_obs(o) * 2)
         enddo
      enddo
!----------------------------------------------------------------------------------------------------

      
!-----------------------------------------------------------------------------------------------------


   if (allocated(rc_en_obs))        deallocate(rc_en_obs)
   if (allocated(ic_en_obs))        deallocate(ic_en_obs) 
   if (allocated(pw_en_obs)    )    deallocate(pw_en_obs    )
   if (allocated(pw_ref_en_obs))    deallocate(pw_ref_en_obs)


   
      do k = 1, en_num - 1
         do jj = 1, freq_num
            ! write(10,*) 'Energy ', k, ' frequency range ', jj
            rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) / dble(num_freq_bins(jj))
            ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) / dble(num_freq_bins(jj))
            rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) / dble(num_freq_bins(jj))
            ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) / dble(num_freq_bins(jj))
            rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) / dble(num_freq_bins(jj))
            pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) / dble(num_freq_bins(jj))
            pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) / dble(num_freq_bins(jj))
            ! write(*,'(A,I8,A)') '   Number of frequency bins for the averaged cross-spectrum: ', num_freq_bins(jj)
            ! write(10,*) rc_avefq_en(jj, k), ic_avefq_en(jj, k)
            ! write(10,*) pw_avefq_en(jj, k), pw_avefq_en_ref(jj, k)
         enddo
      enddo

      k = en_num
! this is for the last energy bin that averages 3 energy bins and NOT 2 
!----------------------------------------------------------------------------------------------------
      do jj = 1, freq_num
         ! write(10,*) 'Energy ', k, ' frequency range ', jj
         rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) / dble(num_freq_bins(jj) / 2 * 3)
         pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) / dble(num_freq_bins(jj) / 2 * 3)
      enddo
!----------------------------------------------------------------------------------------------------
      

!!! Error calculation !!!
   
   if(.not. allocated(var_rc         )) allocate(var_rc         (freq_num, en_num)) 
   if(.not. allocated(var_ic         )) allocate(var_ic         (freq_num, en_num)) 
   if(.not. allocated(covariance     )) allocate(covariance     (freq_num, en_num))
   if(.not. allocated(deriv_rc       )) allocate(deriv_rc       (freq_num, en_num))
   if(.not. allocated(deriv_ic       )) allocate(deriv_ic       (freq_num, en_num))
   if(.not. allocated(coher2_en      )) allocate(coher2_en      (freq_num, en_num))
   if(.not. allocated(err_Aform_rc_ic)) allocate(err_Aform_rc_ic(freq_num, en_num))
   if(.not. allocated(err_std_rc     )) allocate(err_std_rc     (freq_num, en_num))
   if(.not. allocated(err_std_ic     )) allocate(err_std_ic     (freq_num, en_num))
   if(.not. allocated(lag_en         )) allocate(lag_en         (freq_num, en_num)) 
   if(.not. allocated(err_prop_lag   )) allocate(err_prop_lag   (freq_num, en_num))
   if(.not. allocated(err_Aform_lag  )) allocate(err_Aform_lag  (freq_num, en_num))
   if(.not. allocated(err_cohe_lag   )) allocate(err_cohe_lag   (freq_num, en_num))
   
   !Error std of real and im part, propagation for the lag, Adam's formula for rc and ic and lag  
   do k = 1, en_num - 1
      do jj = 1, freq_num 

!propagation error 
      var_rc(jj, k) = rc2_avefq_en(jj, k) - rc_avefq_en(jj, k)**2
      var_ic(jj, k) = ic2_avefq_en(jj, k) - ic_avefq_en(jj, k)**2

      covariance(jj, k) = rc_ic_avefq_en(jj, k) - (rc_avefq_en(jj, k) * ic_avefq_en(jj, k))

      deriv_rc(jj, k) =  (-1.d0 *  ic_avefq_en(jj, k) / rc_avefq_en(jj, k)**2) / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )
      deriv_ic(jj, k) =                        (1.d0  / rc_avefq_en(jj, k))    / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )

!bias term
      ! if ((num_freq_bins(jj) * int_number_obs(o)) .gt. 500 ) then 
      !    bias2 = ((pw_avefq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(num_freq_bins(jj) * int_number_obs(o)))
      ! else
      !    bias2 = 0.d0 
      ! endif

      bias2 = 0.d0 

!Adam's formula real and imaginary part
      err_Aform_rc_ic(jj, k) = sqrt (pw_avefq_en_ref(jj, k) * (pw_avefq_en(jj, k) - ( (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) ) ) / (2.d0 * real(num_freq_bins(jj)) ) )

! Standat error on the real and imaginary part 
      err_std_rc(jj, k) = sqrt(var_rc(jj, k) / dble(num_freq_bins(jj)))
      err_std_ic(jj, k) = sqrt(var_ic(jj, k) / dble(num_freq_bins(jj)))
      
!lag calculation 
      lag_en(jj, k) = atan2(ic_avefq_en(jj, k), rc_avefq_en(jj, k)) / (2.d0 * pi * relevant_freq(jj))

!error lag with propagation formula
      err_prop_lag(jj, k) = sqrt( (deriv_rc(jj, k)**2 * var_rc(jj, k) + deriv_ic(jj, k)**2 * var_ic(jj, k) + 2.d0 * deriv_rc(jj, k) * deriv_ic(jj, k) * covariance(jj, k)) / dble(num_freq_bins(jj)) ) / (2.d0 * pi * relevant_freq(jj))

!coherence calculations      
      coher2_en(jj, k) = (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en(jj, k) * pw_avefq_en_ref(jj, k))

!Coherence error lag 
         err_cohe_lag(jj, k) = sqrt( (1.0 - coher2_en(jj, k)) / (2.0 * coher2_en(jj, k) * dble(num_freq_bins(jj) ) ) ) / (2 * pi * relevant_freq(jj)) 

!Adam's formula lag
         err_Aform_lag(jj, k) = sqrt( pw_avefq_en_ref(jj, k) * ( (pw_avefq_en(jj, k) / (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) ) - ( 1.d0 / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  dble(num_freq_bins(jj))) ) / (2.d0 * pi * relevant_freq(jj)) 


         write(*,*) 'Energy num, Freq range num,  relevant freq,  num_freq_average,   coherence,    cross^2,     bias^2 ' 
         write(*,'(I12, I12, E15.4, I12, E15.4, E15.5, E15.5)') k, jj,  relevant_freq(jj), num_freq_bins(jj), sqrt(coher2_en(jj, k)), rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2, bias2 

      enddo
   enddo

! this is for the last energy bin that averages 3 energy    
!----------------------------------------------------------------------------------------------------
   k = en_num 
   
   do jj = 1, freq_num 

      !propagation error 
      var_rc(jj, k) = rc2_avefq_en(jj, k) - rc_avefq_en(jj, k)**2
      var_ic(jj, k) = ic2_avefq_en(jj, k) - ic_avefq_en(jj, k)**2

      covariance(jj, k) = rc_ic_avefq_en(jj, k) - (rc_avefq_en(jj, k) * ic_avefq_en(jj, k))

      deriv_rc(jj, k) =  (-1.d0 *  ic_avefq_en(jj, k) / rc_avefq_en(jj, k)**2) / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )
      deriv_ic(jj, k) =                        (1.d0  / rc_avefq_en(jj, k))    / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )

      !bias term
      ! if ((num_freq_bins(jj) * int_number_obs(o)) .gt. 500 ) then 
      !    bias2 = ((pw_avefq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(num_freq_bins(jj) / 2 * 3 * int_number_obs(o)))
      ! else
      !    bias2 = 0.d0 
      ! endif

      bias2 = 0.d0 

      !Adam's formula real and imaginary part
      err_Aform_rc_ic(jj, k) = sqrt (pw_avefq_en_ref(jj, k) * (pw_avefq_en(jj, k) - ( (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) ) ) / (2.d0 * real(num_freq_bins(jj) / 2 * 3) ) )

      ! Standat error on the real and imaginary part 
      err_std_rc(jj, k) = sqrt(var_rc(jj, k) / dble(num_freq_bins(jj) / 2 * 3))
      err_std_ic(jj, k) = sqrt(var_ic(jj, k) / dble(num_freq_bins(jj) / 2 * 3))

      !lag calculation 
      lag_en(jj, k) = atan2(ic_avefq_en(jj, k), rc_avefq_en(jj, k)) / (2.d0 * pi * relevant_freq(jj))

      !error lag with propagation formula
      err_prop_lag(jj, k) = sqrt( (deriv_rc(jj, k)**2 * var_rc(jj, k) + deriv_ic(jj, k)**2 * var_ic(jj, k) + 2.d0 * deriv_rc(jj, k) * deriv_ic(jj, k) * covariance(jj, k)) / dble(num_freq_bins(jj) / 2 * 3) ) / (2.d0 * pi * relevant_freq(jj))

      !coherence calculations      
      coher2_en(jj, k) = (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en(jj, k) * pw_avefq_en_ref(jj, k))

      !Coherence error lag 
      err_cohe_lag(jj, k) = sqrt( (1.0 - coher2_en(jj, k)) / (2.0 * coher2_en(jj, k) * dble(num_freq_bins(jj)  / 2 * 3 ) ) ) / (2 * pi * relevant_freq(jj)) 

      !Adam's formula lag
      err_Aform_lag(jj, k) = sqrt( pw_avefq_en_ref(jj, k) * ( (pw_avefq_en(jj, k) / (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) ) - ( 1.d0 / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  dble(num_freq_bins(jj) / 2 * 3)) ) / (2.d0 * pi * relevant_freq(jj)) 


      write(*,*) 'Energy num, Freq range num,  relevant freq,  num_freq_average,   coherence,    cross^2,     bias^2 ' 
      write(*,'(I12, I12, E15.4, I12, E15.4, E15.5, E15.5)') k, jj,  relevant_freq(jj), num_freq_bins(jj), sqrt(coher2_en(jj, k)), rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2, bias2 

   enddo
!----------------------------------------------------------------------------------------------------



   
   call print_lag_ener(rc_avefq_en, ic_avefq_en, err_std_rc, &
        err_std_ic, err_Aform_rc_ic, lag_en, err_Aform_lag, &
        err_cohe_lag, err_prop_lag)


   deallocate(rc_avefq_en    )
   deallocate(ic_avefq_en    )
   deallocate(rc2_avefq_en   )
   deallocate(ic2_avefq_en   )
   deallocate(rc_ic_avefq_en )
   deallocate(pw_avefq_en    )
   deallocate(pw_avefq_en_ref)
   deallocate(num_freq_bins  )
   deallocate(P_noise_ext    )
   deallocate(P_noise_ext_ref)
   deallocate(var_rc         )
   deallocate(var_ic         )
   deallocate(covariance     )
   deallocate(deriv_rc       )
   deallocate(deriv_ic       )
   deallocate(coher2_en      )
   deallocate(err_Aform_rc_ic)
   deallocate(err_std_rc     )
   deallocate(err_std_ic     )
   deallocate(lag_en         )
   deallocate(err_prop_lag   )
   deallocate(err_Aform_lag  )
   deallocate(err_cohe_lag   )
   
  end program cross
