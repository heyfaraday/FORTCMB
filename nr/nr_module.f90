module nr_module

  implicit none

  interface
    subroutine four1(data,nn,isign)
      integer, intent(in) :: isign, nn
      real, dimension(2*nn), intent(in) :: data
    end subroutine four1

    real function gasdev(idum)
      integer, intent(in) :: idum
    end function gasdev

    real function ran1(idum)
      integer, intent(in) :: idum
    end function ran1

  end interface

  contains

end module nr_module
