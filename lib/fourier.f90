module fourier

  implicit none
  double precision, parameter :: my_PI = 4.d0 * datan(1.d0)

  contains

    subroutine p_lm_gen(theta, p_lm_max, p_lm_array)

      implicit none
      double precision, intent(in) :: theta
      integer, intent(in) :: p_lm_max
      double precision, dimension(0:p_lm_max, 0:p_lm_max), &
                                  intent(out) :: p_lm_array

      integer :: m, l

      p_lm_array(0, 0) = 1.d0 / sqrt(4.d0 * my_PI)

      do m = 0, p_lm_max - 1, 1
        p_lm_array(m + 1, m + 1) = -p_lm_array(m, m)*sin(theta) &
        *sqrt(2.d0*m+3.d0)/sqrt(2.d0*m+2.d0)
      end do
      do m = 0, p_lm_max - 1, 1
        p_lm_array(m, m + 1) = p_lm_array(m, m)*cos(theta)*sqrt(2.d0*m+3.d0)
      end do
      do m = 0, p_lm_max - 2, 1
        do l = m + 2, p_lm_max, 1
          p_lm_array(m, l) = ((2.d0*l-1.d0)*dsqrt((l-m)*(2.d0*l+1.d0)) &
          /sqrt((l+m)*(2.d0*l-1.d0))*p_lm_array(m, l-1)*dcos(theta)-(l+m-1.d0) &
          *sqrt((l-m)*(l-1.d0-m)*(2.d0*l+1.d0)) /sqrt((l+m)*(l-1.d0+m) &
          *(2.d0*l-3.d0))*p_lm_array(m, l-2))/(l-m)
        end do
      end do
      do m = 1, p_lm_max, 1
        do l = m, p_lm_max, 1
          p_lm_array(m, l) = p_lm_array(m, l) * sqrt(2.d0)
        end do
      end do

    end subroutine p_lm_gen

end module fourier
