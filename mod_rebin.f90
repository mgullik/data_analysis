!------------------------------------------------------------------
module rebin 

  double precision, allocatable :: reb_freq(:), reb_pds(:), reb_err_pds(:)
  
contains 
  subroutine rebin_PDS(freq, pds_int, int_number, int_len_dim, reb_dim)
    implicit none

    integer         , intent(IN)  :: int_number, int_len_dim
    double precision, intent(IN)  :: freq(int_len_dim / 2), pds_int(int_number, int_len_dim / 2)

    integer              :: i, j, k, reb_dim, start, end_reb
    double precision     :: reb_fac
    integer, allocatable :: reb_array(:)
    double precision, allocatable :: reb_pds2(:)    


    if(allocated(reb_array  )) deallocate(reb_array)
    if(allocated(reb_freq   )) deallocate(reb_freq)
    if(allocated(reb_pds    )) deallocate(reb_pds )
    if(allocated(reb_err_pds)) deallocate(reb_err_pds)
    if(allocated(reb_pds2   )) deallocate(reb_pds2 )


    write(*,*) 'Choose exponensial rebinning factor'
    read(*,*) reb_fac
    print *, ' '

    allocate(reb_freq (int_len_dim /2))
    allocate(reb_array(int_len_dim /2))

!Finds out the dimension of the rebinned array and the rebin indexes
    call rebin_log(freq, reb_freq, reb_array, int_len_dim / 2, reb_dim, reb_fac)

    write(*,*) 'ciao ', reb_dim
    
    allocate(reb_pds    (reb_dim))
    allocate(reb_pds2   (reb_dim))
    allocate(reb_err_pds(reb_dim))

    reb_pds  = 0.d0
    reb_pds2 = 0.d0

    write(*,*) 'ciao ciao ', reb_array

       
       
    do i = 1, int_number
       start   = 0
       end_reb = 0
       do j = 1, reb_dim
          start   = end_reb + 1 
          end_reb = end_reb + reb_array(j)

          do k = start, end_reb
             reb_pds (j) = reb_pds (j) + pds_int(i, k)
             reb_pds2(j) = reb_pds2(j) + pds_int(i, k) * pds_int(i, k)
          enddo
       enddo
    enddo

    do j = 1, reb_dim
       reb_pds (j) = reb_pds (j) / real(int_number * reb_array(j))
       reb_pds2(j) = reb_pds2(j) / real(int_number * reb_array(j))
       reb_err_pds(j) = sqrt((reb_pds2(j) - reb_pds(j)**2 ) / real(int_number * reb_array(j)))
       
    enddo

    if (allocated(reb_array)) deallocate(reb_array)
    if (allocated(reb_pds2 )) deallocate(reb_pds2)

  end subroutine rebin_PDS

!-----------------------------------------------------------------------
!     PURPOUSE:  Logarithmic rebin of a vector. Find how many bins (nf) 
!     the logarithmic rebinned vector has. c is the binning factor, 
!     dim is the dimension of the lin vector. It takes the linear x_axis 
!        
!     INPUTS:  vector -- linear vector
!              dim -- # bin of the linear vector
!              c -- rebinning factor
!      
!     OUTPUTS: rebin_vector -- logarithmic rebinned vector
!              rebin_array -- how many numbers for each rebin
!
!     ROUTINE CALLED: none
!     AUTHOR: Mastroserio 
!     DATE WRITTEN: Feb 2019
!     LAST CHANGE: 
!     NOTE: 
!           
!----------------------------------------------------------------
!----------------------------------------------------------------
      subroutine rebin_log(vector, rebin_vector, rebin_array, dim, nf, c)
        implicit none

        integer         , intent(IN)    :: dim
        integer         , intent(INOUT) :: nf
        integer         , intent(OUT)   :: rebin_array(dim)
        double precision, intent(IN)    :: vector(dim), c 
        double precision, intent(OUT)   :: rebin_vector(dim)

        integer                :: i, j, nnp, remain, cont
        integer, allocatable   :: iar(:), np(:)
 

        i = 0
        j = 0
        cont = 0 
        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i
           if (floor( c**j ).eq.remain) cont = 1
           nnp  = min( floor( c**j ) , remain )
           nnp  = max( 1 , nnp )
           i   = i + nnp
        end do
        if (cont.eq.1) nf = j
        if (cont.ne.1) nf = j - 1

        allocate(iar(0:nf))
        allocate(np(nf))

        np = 0
! First calculate binning scheme
        i = 0
        j = 0
        iar(0) = 0

        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i

           if (floor( c**j ).gt.remain) then
              np(nf)  = np(nf) + remain
              iar(nf) = dim
              i      = i + np(nf)
           else
              np(j)  = min( floor( c**j ) , remain )
              np(j)  = max( 1 , np(j) )
              i      = i + np(j)
              iar(j) = i
           endif

        end do
        
! Now do the binning
        do j = 0, nf - 1
           rebin_vector(j + 1)  = 0.d0
           cont = 0
           do i = iar(j) + 1,iar(j + 1)
              cont = cont + 1
              rebin_vector(j + 1)  = rebin_vector(j + 1) + vector(i)
           end do
           rebin_array (j + 1) = cont
           rebin_vector(j + 1) = rebin_vector(j + 1) / float( np(j + 1) )
        end do

        if(allocated(iar)) deallocate(iar)
        if(allocated(np)) deallocate(np)
        
        return
      end subroutine rebin_log
 
  
end module rebin
!------------------------------------------------------------------
