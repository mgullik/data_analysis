include 'header.h'

program cross
  use dyn_lc
  use rand_fun
  implicit none 

  double precision, parameter   :: pi = acos(-1.0)
  integer              :: i, j, k, kk, jj
  integer, allocatable :: num_freq_bins(:)
  double precision , allocatable :: rc_en(:,:,:), ic_en(:,:,:), &
                          pw_en(:,:,:), pw_ref_en(:,:,:), &
                          lc_ref_cross(:,:,:)

  double precision , allocatable :: rc_avefq_en(:,:),ic_avefq_en(:,:),&
                          rc2_avefq_en(:,:), ic2_avefq_en(:,:), &
                          rc_ic_avefq_en(:,:), pw_avefq_en(:,:), &
                          pw_avefq_en_ref(:,:)
  double precision , allocatable :: var_rc(:,:), var_ic(:,:), deriv_rc(:,:), &
                          deriv_ic(:,:), covariance(:,:), &
                          err_std_rc(:,:), err_std_ic(:,:), &
                          coher2_en(:,:), err_Aform_rc_ic(:,:), &
                          err_Aform_lag(:,:), err_prop_lag(:,:), &
                          err_cohe_lag(:,:), lag_en(:,:)
!P_noise
  double precision               :: freq_limit, bias2
  double precision               :: get_p_noise_pw  !function
  double precision , allocatable :: P_noise_ext(:), P_noise_ext_ref(:), &
                          pw_aveint(:), pw_aveint_ref(:)

  
!******************************************************************
    
!GET LIGHT CURVEs 

  call load_lc_lag_ene()

!Arrays to save the fourier trasform analysis for every energy       
      if(.not. allocated(rc_en)) allocate(rc_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(ic_en)) allocate(ic_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(pw_en)) allocate(pw_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(pw_ref_en)) allocate(pw_ref_en(int_number, int_len_dim / 2, en_num))

      print *, ' -----------------------------------------------'
      print *, '   Starting the cross spectrum analysis'


      if (.not. allocated(lc_ref_cross)) allocate(lc_ref_cross(int_number, int_len_dim, en_num))
      lc_ref_cross = 0.d0
      
!Calculate the reference band for every light curve 
      if (.not. allocated(lc_ref)) then          
!In this case the reference light curve is the sum of all the light curve but the subject band 
         write(*,*) '   Reference light curve done by adding all the lc but subject one'
         do k = 1, en_num             
            do i = 1, int_number
               do kk = 1, en_num 
                  do j = 1, int_len_dim
                     if (kk .ne. k) then 
                        lc_ref_cross(i, j, k) = lc_ref_cross(i, j, k) + lc_en(i, j, kk)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else
!In this case the reference light curve is the loaded reference light curve subtracted by the subject light curve
         write(*,*) '   Reference light curve done by subtracting the subject lc to the loaded one'
          do k = 1, en_num             
            do i = 1, int_number
                  do j = 1, int_len_dim
                     lc_ref_cross(i, j, k) =  lc_ref(i, j) - lc_en(i, j, k)
                  enddo
               enddo
            enddo

         endif
         
            
         do k = 1, en_num             
            do i = 1, int_number
               if (check_power2) then 
                  ! call ncperiodogram(lc_en(i, :, k), lc_ref_cross(i, :, k), rc_en(i, :, k), ic_en(i, :, k), dt, int_len_dim) 
                  ! call periodogram(lc_en(i, :, k), pw_en(i, :, k), dt, int_len_dim)
                  ! call periodogram(lc_ref_cross(i, :, k), pw_ref_en(i, :, k), dt, int_len_dim)
                  
                  ! call ncperiodogram_frac_rms(lc_en(i, :, k), lc_ref_cross(i, :, k), rc_en(i, :, k), ic_en(i, :, k), dt, int_len_dim) 
                  ! call periodogram_frac_rms(lc_en(i, :, k), pw_en(i, :, k), dt, int_len_dim)
                  ! call periodogram_frac_rms(lc_ref_cross(i, :, k), pw_ref_en(i, :, k), dt, int_len_dim)

                  call ncperiodogram_no_norm(lc_en(i, :, k), lc_ref_cross(i, :, k), rc_en(i, :, k), ic_en(i, :, k), dt, int_len_dim) 
                  call periodogram_no_norm(lc_en(i, :, k), pw_en(i, :, k), dt, int_len_dim)
                  call periodogram_no_norm(lc_ref_cross(i, :, k), pw_ref_en(i, :, k), dt, int_len_dim)
                  
               ! write(*,*) 'FFT done', i
            else
               write(*,*) '   The length of the intervals is not a power of 2 ', int_len_dim
               stop 
               
! !Calculate the Fourier transform of the light curves and reference bands
!                call FT_not_fast(lc_en(i, :, k), re_lc, im_lc, int_len_dim)
!                call FT_not_fast(lc_ref_cross, re_ref, im_ref, int_len_dim)

! !Calculate Power spectra and cross spectrum of light curves and ref band
!                do jj = 1, int_len_dim / 2    
!                   rc_en    (i, jj, k) = (re_lc(jj) * re_ref(jj)) + (im_lc(jj) * im_ref(jj))
!                   ic_en    (i, jj, k) = (im_lc(jj) * re_ref(jj)) - (re_lc(jj) * im_ref(jj))
!                   pw_en    (i, jj, k) = (re_lc(jj) * re_lc(jj)) + (im_lc(jj) * im_lc(jj)) 
!                   pw_ref_en(i, jj, k) = (re_ref(jj) * re_ref(jj)) + (im_ref(jj) * im_ref(jj))  
!                enddo               
!                write(*,*) 'Slow Fourier transform done', k, i  

            endif

         enddo ! loop over intervarls (int_number)
         
         write(*,'(A, I3)') '   Finished to calculate the cross spectrum between the reference band and the energy bin number ', k
      enddo !loop over enery bins (en_num)

      write(*,*) '   Start error calculation'
      write(*,*)
      ! write(*,*) '   Press Enter to continue'
      ! read(*,*)

!Compute the poisson noise for every pw of each energy band and ref band  
      if( .not. allocated(P_noise_ext    )) allocate(P_noise_ext    (en_num))
      if( .not. allocated(P_noise_ext_ref)) allocate(P_noise_ext_ref(en_num))
      if( .not. allocated(pw_aveint    )) allocate(pw_aveint    (int_len_dim))
      if( .not. allocated(pw_aveint_ref)) allocate(pw_aveint_ref(int_len_dim))

      call make_freq_array()
      write(*,*) '   Above which frequency do you want to calculate the poisson noise?'
      read (*,*) freq_limit
      write(*,*)
      ! write(*,'(A,F5.1)') '   Calculating the Poisson noise for every energy PDS. Above frequency  ', freq_limit
      ! write(*,*)
      do k = 1, en_num
         pw_aveint     = 0.d0
         pw_aveint_ref = 0.d0
         do i = 1, int_number
            do j = 1, int_len_dim
               ! write(*,*) 'ciao 123 ', pw_en(i, j, k), pw_ref_en(i, j, k)
               pw_aveint    (j) = pw_aveint    (j) + pw_en    (i, j, k)
               pw_aveint_ref(j) = pw_aveint_ref(j) + pw_ref_en(i, j, k)
            enddo
         enddo
         pw_aveint     = pw_aveint     / real(int_number)
         pw_aveint_ref = pw_aveint_ref / real(int_number)

         P_noise_ext    (k) = get_p_noise_pw(pw_aveint    , int_len_dim / 2, freq, freq_limit )
         P_noise_ext_ref(k) = get_p_noise_pw(pw_aveint_ref, int_len_dim / 2, freq, freq_limit )         
         write(*,*) P_noise_ext(k), P_noise_ext_ref(k)
      enddo
      deallocate(pw_aveint    )
      deallocate(pw_aveint_ref)
      write(*,*)
      P_noise_ext     = 0.d0
      P_noise_ext_ref = 0.d0         

      call make_freq_intervals()

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

      do k = 1, en_num 
         do jj = 1, freq_num
            do j = lower_fq(jj), upper_fq(jj)
               do i = 1, int_number             
                  rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) + rc_en(i, j, k)
                  ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) + ic_en(i, j, k)
                  rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) + rc_en(i, j, k) * rc_en(i, j, k)
                  ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) + ic_en(i, j, k) * ic_en(i, j, k)
                  rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) + rc_en(i, j, k) * ic_en(i, j, k)
                  pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) + pw_en(i, j, k)
                  pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) + pw_ref_en(i, j, k)
               enddo
            enddo
         enddo
      enddo

   if (allocated(rc_en)) deallocate(rc_en)
   if (allocated(ic_en)) deallocate(ic_en) 

   if(.not. allocated(num_freq_bins)) allocate(num_freq_bins(freq_num))
   do jj = 1, freq_num
      num_freq_bins(jj) = upper_fq(jj) - lower_fq(jj) + 1
   enddo
   do k = 1, en_num 
      do jj = 1, freq_num
         rc_avefq_en    (jj, k) = rc_avefq_en    (jj, k) / real(num_freq_bins(jj) * int_number)
         ic_avefq_en    (jj, k) = ic_avefq_en    (jj, k) / real(num_freq_bins(jj) * int_number)
         rc2_avefq_en   (jj, k) = rc2_avefq_en   (jj, k) / real(num_freq_bins(jj) * int_number)
         ic2_avefq_en   (jj, k) = ic2_avefq_en   (jj, k) / real(num_freq_bins(jj) * int_number)
         rc_ic_avefq_en (jj, k) = rc_ic_avefq_en (jj, k) / real(num_freq_bins(jj) * int_number)
         pw_avefq_en    (jj, k) = pw_avefq_en    (jj, k) / real(num_freq_bins(jj) * int_number)
         pw_avefq_en_ref(jj, k) = pw_avefq_en_ref(jj, k) / real(num_freq_bins(jj) * int_number)
      enddo
   enddo


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
   do k = 1, en_num 
      do jj = 1, freq_num 

!propagation error 
      var_rc(jj, k) = rc2_avefq_en(jj, k) - rc_avefq_en(jj, k)**2
      var_ic(jj, k) = ic2_avefq_en(jj, k) - ic_avefq_en(jj, k)**2

      covariance(jj, k) = rc_ic_avefq_en(jj, k) - (rc_avefq_en(jj, k) * ic_avefq_en(jj, k))

      deriv_rc(jj, k) =  (-1.d0 *  ic_avefq_en(jj, k) / rc_avefq_en(jj, k)**2) / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )
      deriv_ic(jj, k) =                        (1.d0  / rc_avefq_en(jj, k))    / (1.d0 + (ic_avefq_en(jj, k) / rc_avefq_en(jj, k))**2 )

!bias term
      if ((num_freq_bins(jj) * int_number) .gt. 500 ) then 
         bias2 = ((pw_avefq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(num_freq_bins(jj) * int_number))
      else
         bias2 = 0.d0 
      endif

      ! bias2 = 0.d0 

!Adam's formula real and imaginary part
      err_Aform_rc_ic(jj, k) = sqrt (pw_avefq_en_ref(jj, k) * (pw_avefq_en(jj, k) - ( (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k)) ) ) / (2.d0 * real(num_freq_bins(jj) * int_number) ) )

! Standat error on the real and imaginary part 
      err_std_rc(jj, k) = sqrt(var_rc(jj, k) / real(num_freq_bins(jj) * int_number))
      err_std_ic(jj, k) = sqrt(var_ic(jj, k) / real(num_freq_bins(jj) * int_number))
      
!lag calculation 
      lag_en(jj, k) = atan2(ic_avefq_en(jj, k), rc_avefq_en(jj, k)) / (2.d0 * pi * relevant_freq(jj))

!error lag with propagation formula
      err_prop_lag(jj, k) = sqrt( (deriv_rc(jj, k)**2 * var_rc(jj, k) + deriv_ic(jj, k)**2 * var_ic(jj, k) + 2.d0 * deriv_rc(jj, k) * deriv_ic(jj, k) * covariance(jj, k)) / real(num_freq_bins(jj) * int_number) ) / (2.d0 * pi * relevant_freq(jj))

!coherence calculations      
      coher2_en(jj, k) = (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) / (pw_avefq_en(jj, k) * pw_avefq_en_ref(jj, k))

!Coherence error lag 
         err_cohe_lag(jj, k) = sqrt( (1.0 - coher2_en(jj, k)) / (2.0 * coher2_en(jj, k) * real(num_freq_bins(jj) * int_number) ) ) / (2 * pi * relevant_freq(jj)) 

!Adam's formula lag
         err_Aform_lag(jj, k) = sqrt( pw_avefq_en_ref(jj, k) * ( (pw_avefq_en(jj, k) / (rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2 - bias2) ) - ( 1.d0 / (pw_avefq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  real(num_freq_bins(jj) * int_number)) ) / (2.d0 * pi * relevant_freq(jj)) 


         write(*,*) '          k,          jj,    relevant freq,    coherence,       cross^2,              bias^2,            num_freq_average ' 
         write(*,*) k, jj,  relevant_freq(jj), sqrt(coher2_en(jj, k)), rc_avefq_en(jj, k)**2 + ic_avefq_en(jj, k)**2, bias2, num_freq_bins(jj)

      enddo
   enddo


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
