!-------------------------------------------------------------------
      subroutine make_array_logspace(a_start, a_end, dim, out)
      implicit none
      integer :: i,dim
      real :: a_start,a_end,out(dim)
      real :: du,u1,u2,exp
      if (dim.gt.1) then
         u1=log10(a_start)
         u2=log10(a_end)      
         du=(u2-u1)/float(dim-1)
         exp=log10(a_start)
         do i=1, dim
            out(i)=10**exp
            exp=exp+du
         enddo
      else if(dim.eq.1) then
       out(1)=a_start
      else if(dim.lt.1) then
       write(*,'(A)',advance="no") 'Error! It is impossible create' 
       write(*,*) ' the log scaled array: dim is less than 1'
       stop
      endif
      return
    end subroutine make_array_logspace
!-----------------------------------------------------------------
