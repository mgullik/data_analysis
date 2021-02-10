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

   function poidev(xm)
     implicit none
     double precision poidev, xm, pi
     parameter (pi=3.141592654)
     double precision alxm, em, g, oldm, sq, t, y!, gammln, ran1
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
          sq   = sqrt(2.d0 * xm)
          alxm = log(xm)
          g    = xm * alxm - gammln( xm + 1.)
        end if
 1      y = tan(pi * ran1())
        em = sq * y + xm
        if(em .lt. 0.d0) goto 1
        em = int(em)
        t = 0.9d0*(1. + y**2.) * exp(em * alxm - gammln(em + 1.d0) - g)
        if(ran1() .gt. t) goto 1
      end if
      poidev = em
      return
    end function poidev

    
    FUNCTION gammln(xx)
        implicit none 
        DOUBLE PRECISION gammln, xx
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

        FUNCTION ran1()
          implicit none 
          INTEGER IA,IM,IQ,IR,NTAB,NDIV
          DOUBLE PRECISION ran1,AM,EPS,RNMX
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

  end module rand_fun
!----------------------------------------------------------------





! !------------------------------------------------------------------
! module rand_fun_real

!   implicit none 
!   private
!   integer  :: idum 

!   public   :: set_seed_real, poidev_real, ran1_real

!   contains 
    
!     subroutine set_seed_real(seed)
!       implicit none 
!       integer :: seed
!       idum = seed
!     end subroutine set_seed_real
    
!    function poidev_real(xm)
!      implicit none
!      real poidev_real, xm, pi
!      parameter (pi=3.141592654)
!      real alxm, em, g, oldm, sq, t, y!, gammln, ran1
!      save alxm, g, oldm, sq
!       data oldm /-1./

!       if (xm .lt. 12) then
!         if(xm .ne. oldm) then
!           oldm = xm
!           g    = exp(-xm)
!         end if
!         em = -1.
!         t  = 1.
!  2      em = em + 1.
!         t  = t * ran1_real()
!         if( t .gt. g) goto 2
!       else
!         if(xm .ne. oldm)then
!           oldm = xm
!           sq   = sqrt(2. * xm)
!           alxm = log(xm)
!           g    = xm * alxm - gammln_real( xm + 1.)
!         end if
!  1      y = tan(pi * ran1_real())
!         em = sq * y + xm
!         if(em .lt. 0.) goto 1
!         em = int(em)
!         t = 0.9*(1. + y**2.) * exp(em * alxm - gammln_real(em + 1.) - g)
!         if(ran1_real() .gt. t) goto 1
!       end if
!       poidev_real = em
!       return
!     end function poidev_real

    
!     FUNCTION gammln_real(xx)
!         implicit none 
!         REAL gammln_real, xx
!         INTEGER j
!         DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
!         SAVE cof,stp
!         DATA cof, stp / 76.18009172947146d0, -86.50532032941677d0, 24.01409824083091d0, &
!              -1.231739572450155d0, 0.1208650973866179d-2, -0.5395239384953d-5, 2.5066282746310005d0/
!         x   = xx
!         y   = x
!         tmp = x + 5.5d0
!         tmp = (x + 0.5d0) * log(tmp) - tmp
!         ser = 1.000000000190015d0
!         do 11 j = 1, 6
!            y   = y + 1.d0
!            ser = ser +cof(j) / y
! 11         continue
!            gammln_real = real(tmp + log(stp * ser / x))
!            return
!         END function gammln_real

!         FUNCTION ran1_real()
!           implicit none 
!           INTEGER IA,IM,IQ,IR,NTAB,NDIV
!           REAL ran1_real,AM,EPS,RNMX
!           PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
!           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!           INTEGER j,k,iv(NTAB),iy
!           SAVE iv,iy
!           DATA iv /NTAB*0/, iy /0/
!           if (idum.le.0.or.iy.eq.0) then
!              idum=max(-idum,1)
!              do 11 j=NTAB+8,1,-1
!                 k=idum/IQ
!                 idum=IA*(idum-k*IQ)-IR*k
!                 if (idum.lt.0) idum=idum+IM
!                 if (j.le.NTAB) iv(j)=idum
! 11              continue
!                 iy=iv(1)
!              endif
!              k=idum/IQ
!              idum=IA*(idum-k*IQ)-IR*k
!              if (idum.lt.0) idum=idum+IM
!              j=1+iy/NDIV
!              iy=iv(j)
!              iv(j)=idum
!              ran1_real=min(AM*iy,RNMX)
!              return
!           END function ran1_real

!   end module rand_fun_real
! !----------------------------------------------------------------

