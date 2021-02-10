subroutine make_freq_array()
  use dyn_lc
  implicit none
  integer :: j

! !Frequency array
    df = 0.5d0 / (dt * int_len_dim)
    if (.not. allocated(freq)) allocate(freq(int_len_dim / 2))
    ! print *, dt, int_len_dim

    do j = 1, int_len_dim / 2 
       freq(j) = real(j) / (dt * int_len_dim)
       ! write(*,*) freq(j)
    enddo

  
  end subroutine make_freq_array
