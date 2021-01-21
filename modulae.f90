MODULE dyn_lc
!------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!-----------------------------------------------------------------
! int_len_dim: length of the interval (in units of element) -> this is decided by the user (must be a power of 2)
! int_number: number of interval in the light curve -> this is calculated automatically
  implicit none 
  integer     :: dim_lc, dim_GTI, int_number = -1
  integer     :: int_len_dim, gap = -1, check_gap_num = -1, en_num
  real                 :: dt = -1
  logical              :: check_power2

  real   , allocatable :: lc(:), time(:), err_rate(:), bkg(:)
  real   , allocatable :: start_GTI(:), end_GTI(:)

!split lc variables   
  real   , allocatable :: lc_int(:,:), time_int(:,:), bkg_int(:,:)

  real                 :: df
  real   , allocatable :: freq(:)
  
!lag vs freq variables
  real                 :: ave_rate_freq1, ave_rate_freq2
  real   , allocatable :: lc_freq1(:,:), lc_freq2(:,:)
  ! real   , allocatable :: pw_freq1_ave(:), pw_freq2_ave(:), &
  !                         pw2_freq1_ave(:), pw2_freq2_ave(:)
  real   , allocatable :: lc_en(:,:,:), bkg_en(:,:,:)
  ! real   , allocatable :: time_int_o(:,:,:), lc_en_o(:,:,:,:), bkg_en_o(:,:,:,:)
  integer, allocatable :: split_ind(:), int_len_dim_o(:)
END MODULE dyn_lc
!------------------------------------------------------------------

!------------------------------------------------------------------
!------------------------------------------------------------------
module rand_fun

  implicit none 
  private
  integer  :: idum 

  public   :: set_seed, poidev, ran1

  contains 
    
    subroutine set_seed(seed)
      implicit none 
      integer :: seed
      idum = seed
    end subroutine set_seed

!------------------------------------------------------------------
   function poidev(xm)
     implicit none
     real poidev, xm, pi
     parameter (pi=3.141592654)
     real alxm, em, g, oldm, sq, t, y!, gammln, ran1
     save alxm, g, oldm, sq
      data oldm /-1./

      if (xm .lt. 12) then
        if(xm .ne. oldm) then
          oldm = xm
          g    = exp(-xm)
        end if
        em = -1.
        t  = 1.
 2      em = em + 1.
        t  = t * ran1()
        if( t .gt. g) goto 2
      else
        if(xm .ne. oldm)then
          oldm = xm
          sq   = sqrt(2. * xm)
          alxm = log(xm)
          g    = xm * alxm - gammln( xm + 1.)
        end if
 1      y = tan(pi * ran1())
        em = sq * y + xm
        if(em .lt. 0.) goto 1
        em = int(em)
        t = 0.9*(1. + y**2.) * exp(em * alxm - gammln(em + 1.) - g)
        if(ran1() .gt. t) goto 1
      end if
      poidev = em
      return
    end function poidev
!-----------------------------------------------------------------

!------------------------------------------------------------------
      FUNCTION gammln(xx)
        implicit none 
        REAL gammln, xx
        INTEGER j
        DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
        SAVE cof,stp
        DATA cof, stp / 76.18009172947146d0, -86.50532032941677d0, 24.01409824083091d0, &
             -1.231739572450155d0, 0.1208650973866179d-2, -0.5395239384953d-5, 2.5066282746310005d0/
        x   = xx
        y   = x
        tmp = x + 5.5d0
        tmp = (x + 0.5d0) * log(tmp) - tmp
        ser = 1.000000000190015d0
        do 11 j = 1, 6
           y   = y + 1.d0
           ser = ser +cof(j) / y
11         continue
           gammln = real(tmp + log(stp * ser / x))
           return
        END function gammln
!-----------------------------------------------------------------

!------------------------------------------------------------------
        FUNCTION ran1()
          implicit none 
          INTEGER IA,IM,IQ,IR,NTAB,NDIV
          REAL ran1,AM,EPS,RNMX
          PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
          NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
          INTEGER j,k,iv(NTAB),iy
          SAVE iv,iy
          DATA iv /NTAB*0/, iy /0/
          if (idum.le.0.or.iy.eq.0) then
             idum=max(-idum,1)
             do 11 j=NTAB+8,1,-1
                k=idum/IQ
                idum=IA*(idum-k*IQ)-IR*k
                if (idum.lt.0) idum=idum+IM
                if (j.le.NTAB) iv(j)=idum
11              continue
                iy=iv(1)
             endif
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum.lt.0) idum=idum+IM
             j=1+iy/NDIV
             iy=iv(j)
             iv(j)=idum
             ran1=min(AM*iy,RNMX)
             return
          END function ran1
!----------------------------------------------------------------

  end module rand_fun
!----------------------------------------------------------------

