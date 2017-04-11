module fourier

  implicit none
  double precision, parameter :: fourier_PI = 4.d0 * datan(1.d0)

  contains

    subroutine p_lm_gen(theta, p_lm_max, p_lm_array)

      ! theta = theta as parameter for polynoms
      ! p_lm_max = maximum of l and m
      ! p_lm_array = array for writing

      implicit none
      double precision, intent(in) :: theta
      integer, intent(in) :: p_lm_max
      double precision, dimension(0:p_lm_max, 0:p_lm_max), &
                                  intent(out) :: p_lm_array

      integer :: m, l ! Var for iterating

      ! Initialization
      p_lm_array = 0.d0

      p_lm_array(0, 0) = 1.d0 / dsqrt(4.d0 * fourier_PI)

      do m = 0, p_lm_max - 1, 1
        p_lm_array(m + 1, m + 1) = - p_lm_array(m, m) * dsin(theta) &
        * dsqrt(2.d0 * m + 3.d0)/sqrt(2.d0 * m + 2.d0)
      end do

      do m = 0, p_lm_max - 1, 1
        p_lm_array(m, m + 1) = p_lm_array(m, m) * dcos(theta) &
        * dsqrt(2.d0 * m + 3.d0)
      end do

      do m = 0, p_lm_max - 2, 1
        do l = m + 2, p_lm_max, 1
          p_lm_array(m, l) = ((2.d0 * l - 1.d0) * dsqrt((l - m) &
          * (2.d0 * l + 1.d0)) / dsqrt((l + m) * (2.d0 * l - 1.d0)) &
          * p_lm_array(m, l - 1) * dcos(theta) - (l + m - 1.d0) &
          * dsqrt((l - m) * (l - 1.d0 - m) * (2.d0 * l + 1.d0)) &
          / dsqrt((l + m) * (l - 1.d0 + m) * (2.d0 * l - 3.d0)) &
          * p_lm_array(m, l - 2)) / (l - m)
        end do
      end do

      do m = 1, p_lm_max, 1
        do l = m, p_lm_max, 1
          p_lm_array(m, l) = p_lm_array(m, l) * dsqrt(2.d0)
        end do
      end do

    end subroutine p_lm_gen


    subroutine direct_fourier(n_max, map, p_lm_max, coef)

      ! n_max = number of pixels for different phi
      ! map = array for writing map
      ! p_lm_max = maximum of l and m
      ! coef = array with a_lm

      ! C module for fftw
      use, intrinsic :: iso_c_binding

      implicit none
      include 'fftw3.f03' ! Header for fftw
      integer, intent(in) :: n_max
      double precision, dimension(1:n_max+1, 1:n_max/2+1), intent(out) :: map
      integer, intent(in) :: p_lm_max
      complex, dimension(0:p_lm_max, 0:p_lm_max), intent(in) :: coef

      integer :: j ! Var for iterating in the map
      double precision :: theta ! Theta as parameter for polynoms and map
      ! Array for p_lm_gen
      double precision, dimension(0:p_lm_max, 0:p_lm_max) :: p_lm
      integer :: m, l ! Var for iterating
      ! Arrays for our fftw on sphere method
      double precision, dimension(:), allocatable :: p1, p2
      integer :: err_p1, err_p2 = 0! Error flags memory allocating


      type(C_PTR) :: plan ! Pointer for fftw plan
      ! Arrays for fftw
      complex(C_DOUBLE_COMPLEX), dimension(1:n_max) :: in, out

      allocate(p1(1:n_max), stat=err_p1)
      if (err_p1 /= 0) print *, "p1: Allocation request denied"

      allocate(p2(1:n_max), stat=err_p2)
      if (err_p2 /= 0) print *, "p2: Allocation request denied"

      ! Initialization
      map = 0.d0

      do j = 2, n_max / 2, 1

        p1 = 0.d0
        p2 = 0.d0

        theta = 2.d0 * fourier_PI * (j - 1) / n_max

        call p_lm_gen(theta, p_lm_max, p_lm)

        do m = 0, p_lm_max, 1
          do l = m, p_lm_max, 1
            p1(m + 1) = p1(m + 1) + real(coef(m, l)) * p_lm(m, l)
            p2(m + 1) = p2(m + 1) + aimag(coef(m, l)) * p_lm(m, l)
          end do
        end do

        in = dcmplx(p1, 0.d0)
        out = (0.d0, 0.d0)

        plan = fftw_plan_dft_1d(n_max, in, out, FFTW_FORWARD, FFTW_MEASURE)

        call fftw_execute_dft(plan, in, out)

        call fftw_destroy_plan(plan)

        map(1:n_max, j) = real(out(1:n_max))
        map(n_max+1, j) = real(out(1))

        in = dcmplx(p2, 0.d0)
        out = (0.d0, 0.d0)

        plan = fftw_plan_dft_1d(n_max, in, out, FFTW_FORWARD, FFTW_MEASURE)

        call fftw_execute_dft(plan, in, out)

        call fftw_destroy_plan(plan)

        map(1:n_max, j) = map(1:n_max, j) + aimag(out(1:n_max))
        map(n_max+1, j) = map(n_max+1, j) + aimag(out(1))

      end do

      if (allocated(p1)) deallocate(p1, stat=err_p1)
      if (err_p1 /= 0) print *, "p1: Deallocation request denied"

      if (allocated(p2)) deallocate(p2, stat=err_p2)
      if (err_p2 /= 0) print *, "p2: Deallocation request denied"

    end subroutine direct_fourier

end module fourier
