subroutine make_PDS(freq, lc_int, int_number, dt, int_len_dim, pds, err_pds)
  use rebin
  implicit none

  integer, intent(IN)  :: int_number, int_len_dim
  double precision, intent(IN)  :: dt, lc_int(int_number, int_len_dim), freq(int_len_dim /2)
  double precision, intent(OUT) :: pds(int_len_dim /2), err_pds(int_len_dim /2)

  integer              :: i, j, reb_dim, norm_type
  double precision, allocatable :: pds_int(:,:), pds2(:)
  character (len = 200):: name_base_2
  logical              :: yes_no


  write(*,*) '    Choose the PDS normalisation:'
  write(*,*) '    (1) rms-squared  '
  write(*,*) '    (2) absolute rms '
  write(*,*) '    (3) Leahy  '
  write(*,*)
  read (*,*) norm_type 
  

  if(.not. allocated(pds_int)) allocate(pds_int(int_number, int_len_dim/2))

  select case(norm_type)
     case (1)
        ! write(*,*) '   rms-squared'
        do i = 1, int_number
           call periodogram_frac_rms(lc_int(i,:), pds_int(i, :), dt, int_len_dim)
        enddo
     case (2)
        ! write(*,*) '   Absolute rms'
        do i = 1, int_number
           call periodogram(lc_int(i,:), pds_int(i, :), dt, int_len_dim)
        enddo
     case (3)
        ! write(*,*) '   Leahy'
        do i = 1, int_number
           call   periodogram_leahy(lc_int(i,:), pds_int(i, :), dt, int_len_dim)
        enddo
     case default
        write(*,*) '   No valid normalisation. Exit...'
        stop
     end select
  

  
  if(.not. allocated(pds2)) allocate(pds2(int_len_dim/2))
  pds  = 0.0
  pds2 = 0.0
  do j = 1, int_len_dim / 2 
     do i = 1, int_number
        pds  (j) = pds  (j) + pds_int(i, j)
        pds2 (j) = pds2 (j) + pds_int(i, j) * pds_int(i, j)
     enddo
  enddo
  pds  = pds  / dble(int_number)
  pds2 = pds2 / dble(int_number)

  do j = 1, int_len_dim / 2 - 1
     err_pds(j) = sqrt((pds2(j) - pds(j)**2) / dble(int_number))
  enddo
  

!------------------- PRINT and REBIN ---------------------
  
  if (yes_no('   Do you want to save the PDS of the light curve?')) then
     
      select case(norm_type)
      case (1)
         name_base_2 = 'PDS_rms2.qdp'  
      case (2) 
         name_base_2 = 'PDS_abs_rms.qdp'  
      case (3)
         name_base_2 = 'PDS_leahy.qdp'  
      case default
         write(*,*) '   No valid normalisation.'
      end select
      write(*,*) '   The name the PDS file is ', trim(name_base_2)
      open(71, file = trim(name_base_2))
      write(71, *) 'skip on'
      write(71, *) 'read serr 1 2'
      do j = 1, int_len_dim / 2 - 1
         write(71, *) (freq(j + 1) + freq(j)) * 0.5 , (freq(j + 1) - freq(j)) * 0.5, pds(j), err_pds(j)  
      enddo
      write(71, *) 'no no'
      write(71, *) 'log x y on'
      write(71, *) 'scr white'
      write(71, *) 'lw 5'
      write(71, *) 'la x Frequency [Hz]'
      write(71, *) 't off'
      select case(norm_type)
      case (1)
         write(71, *) 'la y Power [rms\u2\d]'
      case (2) 
         write(71, *) 'la y Power [abs rms]'
      case (3)
         write(71, *) 'la y Power [Leahy]'
      case default
         write(*,*) '   No valid normalisation.'
      end select

      close(71)

! -------------- REBIN --------------      
666   continue 
      if (yes_no('   Do you want to rebin the PDS?')) then
         ! If(.not. allocated(Reb_Freq) ) allocate(Reb_freq (int_len_dim /2))
         call rebin_PDS(freq, pds_int, int_number, int_len_dim, reb_dim)
         select case(norm_type)
         case (1)
            name_base_2 = 'PDS_reb_rms2.qdp'  
         case (2) 
            name_base_2 = 'PDS_reb_abs_rms.qdp'  
         case (3)
            name_base_2 = 'PDS_reb_leahy.qdp'  
         case default
            write(*,*) '   No valid normalisation.'
         end select
         write(*,*) '   The name the rebinned PDS is ', trim(name_base_2)
         open(71, file = trim(name_base_2))
         write(71, *) 'skip on'
         write(71, *) 'read serr 1 2'
         do j = 1, reb_dim - 1            
            write(71, *) (reb_freq(j + 1) + reb_freq(j)) * 0.5 ,(reb_freq(j + 1) - reb_freq(j)) * 0.5, reb_pds(j), reb_err_pds(j)  
         enddo
         write(71, *) 'no no'
         write(71, *) 'log x y on'
         write(71, *) 'scr white'
         write(71, *) 'lw 5'
         write(71, *) 'la x Frequency [Hz]'
         write(71, *) 't off'
         select case(norm_type)
         case (1)
            write(71, *) 'la y Power [rms\u2\d]'
         case (2) 
            write(71, *) 'la y Power [abs rms]'
         case (3)
            write(71, *) 'la y Power [Leahy]'
         case default
            write(*,*) '   No valid normalisation.'
         end select

         close(71)

      else
         go to 111
      endif

      if (yes_no('   Are you satified with the rebin?')) then
      else
         goto 666
      endif

      if (yes_no('   Do you want to create xspec file of the rebin PDS?')) then
         select case(norm_type)
         case (1)
            name_base_2 = 'PDS_reb_rms2_xspec.dat'  
         case (2) 
            name_base_2 = 'PDS_reb_abs_rms_xspec.dat'  
         case (3)
            name_base_2 = 'PDS_reb_leahy_xspec.dat'  
         case default
            write(*,*) '   No valid normalisation.'
         end select

         open(11, file = trim(name_base_2))
         ! write(11, *) '0.1 0.3 0.01 0.01 '
         do j = 1, reb_dim - 1            
            write(11, *) reb_freq(j), reb_freq(j + 1), (reb_freq(j + 1) - reb_freq(j))* reb_pds(j), (reb_freq(j + 1) - reb_freq(j)) * reb_err_pds(j)  
         enddo
         ! write(11, *)  l_bin(k), r_bin(k),  (r_bin(k) - l_bin(k)) * lag_freq(jj, k), (r_bin(k) - l_bin(k)) *  error_Aform_lag(jj, k)
      ! write(11, *) '10.0 20. 0.01 0.01 '
         close(11)
         write(*,*) '   The name the xspec PDS is ', trim(name_base_2)
         write(*,*) '   Use flx2xsp to create the pha files from the PDS file'

         
      endif
      
!--------------- Print final statements ---------------
      ! write(*,*)
      ! write(*,*) '   The normalisation of the PDS is rms-squared (see Uttley et al 2014 Eq.3)'
111   continue
      write(*,*)
  ! select case(norm_type)
  !    case (1)
  !       write(*,*) '   The normalisation is rms-squared'
  !    case (2)
  !       write(*,*) '   The normalisation is absolute rms'
  !    case (3)
  !       write(*,*) '   The normalisation is leahy'
  !    case default
  !       write(*,*) '   No valid normalisation.'
  !    end select

      write(*,*) '   Press Enter to continue'
      read(*,*)
      write(*,*)

   endif  ! end of the print and rebin stuff

   if (allocated(pds_int)) deallocate(pds_int)
   if (allocated(pds2   )) deallocate(pds2   )

  
 end subroutine make_PDS
