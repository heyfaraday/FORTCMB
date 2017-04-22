module map_coef

  ! Healpix module
  use healpix_types

  implicit none

  contains

    subroutine a_lm_gasdev(p_lm_max, coef, seed, mean, std)

      ! p_lm_max - maximum of l and m
      ! coef - array for writing
      ! seed - number for nr ran. num. generator
      ! mean, std - mean and std values for normal distribution

      ! Numerical recipes modules
      use nr, only: gasdev
      use ran_state, only: ran_seed

      implicit none
      integer(kind=i8b), intent(in) :: p_lm_max
      complex(kind=dpc), dimension(0:p_lm_max, 0:p_lm_max), &
                                  intent(out) :: coef
      integer(kind=i4b), intent(in) :: seed
      real(kind=dp), intent(in) :: mean, std

      integer(kind=i8b) :: m, l ! Vars for iterating
      real(kind=sp) :: rand1_sp, rand2_sp ! Vars for complex ran. num.

      ! Initialization
      coef = (0.0_dp, 0.0_dp)

      call ran_seed(seed)

      do l = 0, p_lm_max, 1
        call gasdev(rand1_sp)
        call gasdev(rand2_sp)
        coef(0, l) = dcmplx(std * rand1_sp + mean, 0.0_dp + mean)
        coef(l, l) = dcmplx(std * rand2_sp + mean, 0.0_dp + mean)
      end do

      do m = 1, p_lm_max, 1
        do l = m + 1, p_lm_max, 1
          call gasdev(rand1_sp)
          call gasdev(rand2_sp)
          coef(m, l) = dcmplx(std * rand1_sp + mean, std * rand2_sp + mean)
        end do
      end do

      ! 0- and 1- mode disabled
      coef(0, 0) = (0.0_dp, 0.0_dp)
      coef(0, 1) = (0.0_dp, 0.0_dp)
      coef(1, 1) = (0.0_dp, 0.0_dp)

    end subroutine a_lm_gasdev

end module map_coef
