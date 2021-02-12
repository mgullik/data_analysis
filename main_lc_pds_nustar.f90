include 'header.h'
program ls_to_pds_nustar
  use dyn_lc
  implicit none

  integer                       :: i, j
  double precision              :: average_rate_A, average_rate_B 
  double precision, allocatable :: lc_FPMA(:,:), lc_FPMB(:,:)  
  double precision, allocatable :: pdsAB(:), err_pdsAB(:)
  double precision, allocatable :: pdsA(:), err_pdsA(:), &
                                   pdsB(:), err_pdsB(:)
  character (len = 200) :: filename, lc_nameA, lc_nameB 
  logical               :: yes_no
  
  filename = '/Users/gullo/Work/BHB_project/Gx339/2021_outburst/spec/nu90702303003_srcA_jg_sr.lc'
  
  call load_single_lc(filename)

  allocate(lc_FPMA(int_number, int_len_dim))
  allocate(lc_FPMB(int_number, int_len_dim))
  
  average_rate_A = 0.d0
  do i = 1, int_number
     do j = 1, int_len_dim
        lc_FPMA(i,j) =  lc_int(i, j)
        average_rate_A = average_rate_A + lc_int(i, j)
     enddo
  enddo
  average_rate_A = average_rate_A / dble(int_number * int_len_dim)

  write(*,'(A, F8.3)') '   Average rate of FPMA: ', average_rate_A 


  filename = '/Users/gullo/Work/BHB_project/Gx339/2021_outburst/spec/nu90702303003_srcB_jg_sr.lc'
  call load_single_lc(filename)

  average_rate_B = 0.d0
  do i = 1, int_number
     do j = 1, int_len_dim
        lc_FPMB(i,j) =  lc_int(i, j)
        average_rate_B = average_rate_B + lc_int(i, j)
     enddo
  enddo
  average_rate_B = average_rate_B / dble(int_number * int_len_dim)

  write(*,'(A, F8.3)') '   Average rate of FPMA: ', average_rate_B 
  
  write(*,*)
  if (yes_no('   Do you want to print the lcs?')) then
     lc_nameA = 'lc_FPMA_split.qdp'
     lc_nameB = 'lc_FPMB_split.qdp'
     open(20, file = lc_nameA)
     open(21, file = lc_nameB)
     write(20,*) 'skip on'
     write(21,*) 'skip on'
     do i = 1, int_number
        do j = 1, int_len_dim
           write(20,*) time_int(i, j), lc_FPMA(i, j)
           write(21,*) time_int(i, j), lc_FPMB(i, j)
        enddo
        write(20,*) 'no no '
        write(21,*) 'no no '
     enddo
     write(20,*) 'scr white'
     write(21,*) 'scr white'
     write(*,*) '   The file name of the light curves are ', trim(lc_nameA), trim(lc_nameB)
  endif
  write(*,*)

  call make_freq_array()

  allocate(pdsAB    (int_len_dim/2))
  allocate(err_pdsAB(int_len_dim/2))
  call make_pds_FPMA_B(freq, lc_FPMA, lc_FPMB, int_number, int_len_dim, pdsAB, err_pdsAB)


  if(.not. allocated(pdsA    )) allocate(pdsA    (int_len_dim/2))
  if(.not. allocated(pdsB    )) allocate(pdsB    (int_len_dim/2))
  if(.not. allocated(err_pdsA)) allocate(err_pdsA(int_len_dim/2))
  if(.not. allocated(err_pdsB)) allocate(err_pdsB(int_len_dim/2))
  
  call make_pds(freq, lc_FPMA, int_number, int_len_dim, pdsA, err_pdsA)

  write(*,*) 
  call execute_command_line('mv PDS.qdp PDS_FPMA.qdp')
  call execute_command_line('mv PDS_reb.qdp PDS_FPMA_reb.qdp')
  write(*,*) '   PDS.qdp changed in PDS_FPMA.qdp and PDS_reb.qdp in PDS_FPMA_reb.qdp'
  write(*,*) 
  write(*,*) '   Press enter to continue'
  read (*,*)
  
  call make_pds(freq, lc_FPMB, int_number, int_len_dim, pdsB, err_pdsB)
  write(*,*) 
  call execute_command_line('mv PDS.qdp PDS_FPMB.qdp')
  call execute_command_line('mv PDS_reb.qdp PDS_FPMB_reb.qdp')
  write(*,*) '   PDS.qdp changed in PDS_FPMB.qdp and PDS_reb.qdp in PDS_FPMB_reb.qdp'
  write(*,*) 
  write(*,*) '   Press enter to continue'
  read (*,*)

  
end program ls_to_pds_nustar
  
  
