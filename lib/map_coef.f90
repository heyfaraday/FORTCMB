module map_coef

  implicit none

  contains

    subroutine a_lm_gasdev(p_lm_max, coef, seed, mean, std)

      ! p_lm_max = maximum of l and m
      ! coef = array for writing
      ! seed = number for nr ran. num. generator
      ! mean, std = mean and std values for normal distribution

      ! Numerical recipes modules
      use nr, only: gasdev
      use ran_state, only: ran_seed

      implicit none
      integer, intent(in) :: p_lm_max
      complex, dimension(0:p_lm_max, 0:p_lm_max), &
                                  intent(out) :: coef
      integer, intent(in) :: seed
      double precision, intent(in) :: mean, std

      integer :: m, l ! Var for iterating
      real :: rand1, rand2 ! Vars for complex ran. num.

      ! Initialization
      coef = (0.d0, 0.d0)

      call ran_seed(seed)

      do l = 0, p_lm_max, 1
        call gasdev(rand1)
        coef(0, l) = cmplx(std * rand1 + mean, 0.d0 + mean)
      end do

      do m = 1, p_lm_max, 1
        do l = m, p_lm_max, 1
          call gasdev(rand1)
          call gasdev(rand2)
          coef(m, l) = cmplx(std * rand1 + mean, std * rand2 + mean)
        end do
      end do

      do l = 0, p_lm_max, 1
        coef(l, l) = coef(l, l) - dcmplx(0.d0, aimag(coef(l, l)))
      end do

      ! 0- and 1- mode disabled
      coef(0, 0) = (0.d0, 0.d0)
      coef(0, 1) = (0.d0, 0.d0)
      coef(1, 1) = (0.d0, 0.d0)

    end subroutine a_lm_gasdev

end module map_coef
