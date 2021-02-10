subroutine make_PDS_ref()
  use dyn_lc
  use rebin
  implicit none

  integer              :: i, j, reb_dim
  double precision, allocatable :: pw_ref_int(:,:), pw2_ref(:)
  ! real   , allocatable :: reb_freq(:), reb_pds_ref(:), reb_err_pds_ref(:)
  character (len = 200):: name_base_2
  logical              :: yes_no
  
  if (.not. allocated(pw_ref_int)) allocate(pw_ref_int(int_number, int_len_dim / 2))
  do i = 1, int_number
     call periodogram_frac_rms(lc_ref(i,:), pw_ref_int(i, :), int_len_dim)
  enddo

!Declared in the modulae file
  if(.not. allocated(pw_ref)    ) allocate(pw_ref    (int_len_dim/2))
  if(.not. allocated(err_pw_ref)) allocate(err_pw_ref(int_len_dim/2))

  if(.not. allocated(pw2_ref)   ) allocate(pw2_ref   (int_len_dim/2))
  pw_ref  = 0.0
  pw2_ref = 0.0
  do j = 1, int_len_dim / 2 
     do i = 1, int_number
        pw_ref  (j) = pw_ref  (j) + pw_ref_int(i, j)
        pw2_ref (j) = pw2_ref (j) + pw_ref_int(i, j)**2
     enddo
  enddo
  pw_ref  = pw_ref  / real(int_number)
  pw2_ref = pw2_ref / real(int_number)

  do j = 1, int_len_dim / 2 - 1
     err_pw_ref(j) = sqrt((pw2_ref(j)-pw_ref(j)**2)/real(int_number))
  enddo

  if (.not. allocated(freq)) then
     call make_freq_array()
  endif
  

!------------------- PRINT and REBIN ---------------------
  
  if (yes_no('   Do you want to save the PDS of the reference band light curve?')) then
     
      name_base_2 = 'ener_PDS_ref_band.dat'  
      write(*,*) '   The name the reference band PDS is ', trim(name_base_2)
      open(71, file = trim(name_base_2))
      write(71, *) 'skip on'
      write(71, *) 'read serr 1 2'
      do j = 1, int_len_dim / 2 - 1
         write(71, *) (freq(j + 1) + freq(j)) * 0.5 , df, pw_ref(j), err_pw_ref(j)  
      enddo
      write(71, *) 'no no'
      write(71, *) 'log x y on'
      close(71)

! -------------- REBIN --------------      
      if (yes_no('   Do you want to rebin the PDS?')) then
         ! If(.not. allocated(Reb_Freq) ) allocate(Reb_freq (int_len_dim /2))
         call rebin_PDS(freq, pw_ref_int, int_number, int_len_dim, reb_dim)
         name_base_2 = 'ener_PDS_ref_band_reb.dat'  
         write(*,*) '   The name the rebinned PDS is ', trim(name_base_2)
         open(71, file = trim(name_base_2))
         write(71, *) 'skip on'
         write(71, *) 'read serr 1 2'
         do j = 1, reb_dim - 1
            
            write(71, *) (reb_freq(j + 1) + reb_freq(j)) * 0.5 ,(reb_freq(j + 1) - reb_freq(j)) * 0.5, reb_pds(j), reb_err_pds(j)  
         enddo
         write(71, *) 'no no'
         write(71, *) 'log x y on'
         close(71)
         
      endif


!--------------- Print final statements ---------------
      write(*,*)
      write(*,*) '   The normalisation of the PDS is rms-squared (see Uttley et al 2014 Eq.3)'
      write(*,*)
      write(*,*) '   Press Enter to continue'
      read(*,*)
      write(*,*)

   endif  ! end of the print and rebin stuff

   if (allocated(pw2_ref        )) deallocate(pw2_ref        )

  
 end subroutine make_PDS_ref
