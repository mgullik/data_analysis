subroutine make_PDS(freq, lc_int, int_number, int_len_dim, pds, err_pds)
  use rebin
  implicit none

  integer, intent(IN)  :: int_number, int_len_dim
  real   , intent(IN)  :: lc_int(int_number, int_len_dim), freq(int_len_dim /2)
  real   , intent(OUT) :: pds(int_len_dim /2), err_pds(int_len_dim /2)

  integer              :: i, j, reb_dim
  real   , allocatable :: pds_int(:,:), pds2(:)
  character (len = 200):: name_base_2
  logical              :: yes_no
  
  if(.not. allocated(pds_int)) allocate(pds_int(int_number, int_len_dim/2))
  do i = 1, int_number
     call periodogram_frac_rms(lc_int(i,:), pds_int(i, :), int_len_dim)
  enddo

  if(.not. allocated(pds2)) allocate(pds2(int_len_dim/2))
  pds  = 0.0
  pds2 = 0.0
  do j = 1, int_len_dim / 2 
     do i = 1, int_number
        pds  (j) = pds  (j) + pds_int(i, j)
        pds2 (j) = pds2 (j) + pds_int(i, j) * pds_int(i, j)
     enddo
  enddo
  pds  = pds  / real(int_number)
  pds2 = pds2 / real(int_number)

  do j = 1, int_len_dim / 2 - 1
     err_pds(j) = sqrt((pds2(j) - pds(j)**2) / real(int_number))
  enddo
  

!------------------- PRINT and REBIN ---------------------
  
  if (yes_no('   Do you want to save the PDS of the reference band light curve?')) then
     
      name_base_2 = 'PDS.dat'  
      write(*,*) '   The name the reference band PDS is ', trim(name_base_2)
      open(71, file = trim(name_base_2))
      write(71, *) 'skip on'
      write(71, *) 'read serr 1 2'
      do j = 1, int_len_dim / 2 - 1
         write(71, *) (freq(j + 1) + freq(j)) * 0.5 , (freq(j + 1) - freq(j)) * 0.5, pds(j), err_pds(j)  
      enddo
      write(71, *) 'no no'
      write(71, *) 'log x y on'
      close(71)

! -------------- REBIN --------------      
      if (yes_no('   Do you want to rebin the PDS?')) then
         ! If(.not. allocated(Reb_Freq) ) allocate(Reb_freq (int_len_dim /2))
         call rebin_PDS(freq, pds_int, int_number, int_len_dim, reb_dim)
         name_base_2 = 'PDS_reb.dat'  
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

   if (allocated(pds_int)) deallocate(pds_int)
   if (allocated(pds2   )) deallocate(pds2   )

  
 end subroutine make_PDS
