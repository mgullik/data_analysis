subroutine load_lc_lag_ene()
  use dyn_lc
  implicit none

  integer              :: i, j, k
  double precision     :: en_units
  character (len = 500)  :: name_path,  name_base, filename, &
       filename_en_bin, name_extension, filename_ref
  logical              :: yes_no

  integer              :: file_line_num
  integer, allocatable :: l_bin(:), r_bin(:)
  
      print *, ' '
      print *, '   Name of the path!'
      name_path = '/Users/gullo/Work/BHB_project/CygX1_Ole/energy_resolved_lcs'
      write(*,*) trim(name_path)
      
      print *, ' '
      print *, '   Name of the channel/energy bin'
      filename_en_bin = '/Users/gullo/Work/BHB_project/CygX1_Ole/bins.txt'
      write(*,*) trim(filename_en_bin)
      print *, ' '
      
      if (yes_no('   Are the energies in the bins.txt file expressed in eV? ')) then
         en_units = 1000.d0
      else
         write(*,*) '   Write the number the energies has to be divided to obtain keV'
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
      write(*, '(A, F6.3, A, F6.3)') '   Energy boundaries in keV: ', (l_en_bin(1) + r_en_bin(1)) * 0.5d0 ,' - ', (l_en_bin(en_num) + r_en_bin(en_num)) * 0.5d0  
      print *, ' '

      if (.not. allocated(ave_rate_en) ) allocate(ave_rate_en(en_num))
      ! if (.not. allocated(ave_bkg_en)  ) allocate(ave_bkg_en (en_num))
      ave_rate_en = 0.d0

!GET ALL THE LIGHT CURVES  

      write(name_base, '(A, A)')  trim(name_path), '/ni2636010101_0.01s_'  
      ! write(*,*)  name_base
      write(name_extension, '(A)') 'eV.lc'

      do k = 1, en_num
         
!Work out the name of the light curve 
         write(filename, '(A, I0, A, I0, A)') trim(name_base), l_bin(k), 'eV-', r_bin(k), trim(name_extension)
         write(*,*)
         call extract_lc(filename)

         if (allocated(lc)) then 
            call split_lc()
         else
            write(*,*) '    !!! NO LIGHT CURVE !!! '
            write(*,*) '        STOP HERE   '
            stop
         endif

         if (.not. allocated(lc_en)   ) allocate(lc_en   (int_number, int_len_dim, en_num))
         ! if (.not. allocated(bkg_en)  ) allocate(bkg_en  (int_number, int_len_dim, en_num))

         do j = 1, int_len_dim
            do i = 1, int_number
               lc_en (i, j, k) = lc_int(i, j)
               ! bkg_en(i, j, k) = bkg_int(i, j)
               ave_rate_en(k) = ave_rate_en(k) + lc_int(i, j)
               ! ave_bkg_en(k)  = ave_bkg_en(k)  + bkg_int(i, j)
            enddo
         enddo
         ave_rate_en(k) = ave_rate_en(k)  / real(int_len_dim * int_number)
         ! ave_bkg_en (k) = ave_bkg_en (k)  / real(int_len_dim * int_number)

         write(*,*) '*************************************'
         write(*,*) '   Light curve NAME    '          , trim(filename)
         write(*,*) '   Average count rate (count/s)'  , ave_rate_en(k)
         write(*,*) '*************************************'
         write(*,*) 
      enddo

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

      
end subroutine load_lc_lag_ene
