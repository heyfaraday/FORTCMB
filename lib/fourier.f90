module fourier

  ! Healpix module
  use healpix_types

  implicit none
  real(kind=dp), parameter :: fourier_PI = 4.0_dp * atan(1.0_dp)

  contains

    subroutine p_lm_gen(theta, p_lm_max, p_lm_array)

      ! theta = theta as parameter for polynoms
      ! p_lm_max = maximum of l and m
      ! p_lm_array = array for writing

      implicit none
      real(kind=dp), intent(in) :: theta
      integer(kind=i8b), intent(in) :: p_lm_max
      real(kind=dp), dimension(0:p_lm_max, 0:p_lm_max), &
                                  intent(out) :: p_lm_array

      integer(kind=i8b) :: m, l ! Var for iterating

      ! Initialization
      p_lm_array = 0.0_dp

      p_lm_array(0, 0) = 1.0_dp / sqrt(4.0_dp * fourier_PI)

      do m = 0, p_lm_max - 1, 1
        p_lm_array(m + 1, m + 1) = - p_lm_array(m, m) * sin(theta) &
        * sqrt(2.0_dp * m + 3.0_dp)/sqrt(2.0_dp * m + 2.0_dp)
      end do

      do m = 0, p_lm_max - 1, 1
        p_lm_array(m, m + 1) = p_lm_array(m, m) * cos(theta) &
        * sqrt(2.0_dp * m + 3.0_dp)
      end do

      do m = 0, p_lm_max - 2, 1
        do l = m + 2, p_lm_max, 1
          p_lm_array(m, l) = ((2.0_dp * l - 1.0_dp) * sqrt((l - m) &
          * (2.0_dp * l + 1.0_dp)) / sqrt((l + m) * (2.0_dp * l - 1.0_dp)) &
          * p_lm_array(m, l - 1) * cos(theta) - (l + m - 1.0_dp) &
          * dsqrt((l - m) * (l - 1.0_dp - m) * (2.0_dp * l + 1.0_dp)) &
          / dsqrt((l + m) * (l - 1.0_dp + m) * (2.0_dp * l - 3.0_dp)) &
          * p_lm_array(m, l - 2)) / (l - m)
        end do
      end do

      do m = 1, p_lm_max, 1
        do l = m, p_lm_max, 1
          p_lm_array(m, l) = p_lm_array(m, l) * dsqrt(2.0_dp)
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
      integer(kind=i4b), intent(in) :: n_max ! i4b
      real(kind=dp), dimension(1:n_max+1, 1:n_max/2+1), intent(out) :: map
      integer(kind=i8b), intent(in) :: p_lm_max
      complex(kind=dpc), dimension(0:p_lm_max, 0:p_lm_max), intent(in) :: coef

      integer(kind=i8b) :: j ! Var for iterating in the map
      real(kind=dp) :: theta ! Theta as parameter for polynoms and map
      ! Array for p_lm_gen
      real(kind=dp), dimension(0:p_lm_max, 0:p_lm_max) :: p_lm
      integer(kind=i8b) :: m, l ! Var for iterating
      ! Arrays for our fftw on sphere method
      real(kind=dp), dimension(:), allocatable :: p1, p2
      integer(kind=i8b) :: err_p1, err_p2 = 0! Error flags memory allocating


      type(C_PTR) :: plan ! Pointer for fftw plan
      ! Arrays for fftw
      complex(C_DOUBLE_COMPLEX), dimension(1:n_max) :: in, out

      allocate(p1(1:n_max), stat=err_p1)
      if (err_p1 /= 0) print *, "p1: Allocation request denied"

      allocate(p2(1:n_max), stat=err_p2)
      if (err_p2 /= 0) print *, "p2: Allocation request denied"

      ! Nyquist warning
      if ( p_lm_max > n_max / 2 + 1 ) then
        print *, "Nyquist frequency warning!"
      end if

      ! Initialization
      map = 0.0_dp

      do j = 2, n_max / 2, 1

        p1 = 0.0_dp
        p2 = 0.0_dp

        theta = 2.0_dp * fourier_PI * (j - 1) / n_max

        call p_lm_gen(theta, p_lm_max, p_lm)

        do m = 0, p_lm_max, 1
          do l = m, p_lm_max, 1
            p1(m + 1) = p1(m + 1) + real(coef(m, l)) * p_lm(m, l)
            p2(m + 1) = p2(m + 1) + aimag(coef(m, l)) * p_lm(m, l)
          end do
        end do

        in = cmplx(p1, 0.0_dp)
        out = (0.0_dp, 0.0_dp)

        plan = fftw_plan_dft_1d(n_max, in, out, FFTW_FORWARD, FFTW_ESTIMATE)

        call fftw_execute_dft(plan, in, out)

        call fftw_destroy_plan(plan)

        map(1:n_max, j) = real(out(1:n_max))
        map(n_max+1, j) = real(out(1))

        in = cmplx(p2, 0.0_dp)
        out = (0.0_dp, 0.0_dp)

        plan = fftw_plan_dft_1d(n_max, in, out, FFTW_FORWARD, FFTW_ESTIMATE)

        call fftw_execute_dft(plan, in, out)

        call fftw_destroy_plan(plan)

        map(1:n_max, j) = map(1:n_max, j) - aimag(out(1:n_max))
        map(n_max+1, j) = map(n_max+1, j) - aimag(out(1))

      end do

      if (allocated(p1)) deallocate(p1, stat=err_p1)
      if (err_p1 /= 0) print *, "p1: Deallocation request denied"

      if (allocated(p2)) deallocate(p2, stat=err_p2)
      if (err_p2 /= 0) print *, "p2: Deallocation request denied"

    end subroutine direct_fourier


    subroutine inverse_fourier(n_max, map, p_lm_max, coef)

      ! n_max = number of pixels for different phi
      ! map = array for reading map
      ! p_lm_max = maximum of l and m
      ! coef = array for writing

      ! C module for fftw
      use, intrinsic :: iso_c_binding

      implicit none
      include 'fftw3.f03' ! Header for fftw
      integer(kind=i4b), intent(in) :: n_max
      real(kind=dp), dimension(1:n_max+1, 1:n_max/2+1), intent(in) :: map
      integer(kind=i8b), intent(in) :: p_lm_max
      complex(kind=dpc), dimension(0:p_lm_max, 0:p_lm_max), intent(out) :: coef

      integer(kind=i8b) :: i, j ! Vars for iterating in the map
      real(kind=dp) :: theta ! Theta as parameter for polynoms and map
      ! Array for p_lm_gen
      real(kind=dp), dimension(0:p_lm_max, 0:p_lm_max) :: p_lm
      integer(kind=i8b) :: m, l ! Var for iterating
      ! Normalization for backward fourier transform
      real(kind=dp) :: norm = 0.0_dp

      type(C_PTR) :: plan ! Pointer for fftw plan
      ! Arrays for fftw
      complex(C_DOUBLE_COMPLEX), dimension(1:n_max) :: in, out

      ! Nyquist warning
      if ( p_lm_max > n_max / 2 + 1 ) then
        print *, "Nyquist frequency warning!"
      end if

      ! Initialization
      coef = (0.0_dp, 0.0_dp)

      do j = 2, n_max / 2, 1

        in = (0.0_dp, 0.0_dp)
        out = (0.0_dp, 0.0_dp)

        theta = 2.0_dp * fourier_PI * (j - 1) / n_max

        call p_lm_gen(theta, p_lm_max, p_lm)

        do i = 1, n_max, 1
            in(i) = cmplx(map(i,j), 0.0_dp)
        end do

        plan = fftw_plan_dft_1d(n_max, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)

        call fftw_execute_dft(plan, in, out)

        call fftw_destroy_plan(plan)

        do i = 2, n_max/2, 1
            out(i) = cmplx(real(out(i)) + real(out(n_max + 2 - i)), &
                            aimag(out(i)) - aimag(out(n_max + 2 - i)))
        end do

        out = out / n_max
        out(n_max/2+2:n_max) = (0.0_dp, 0.0_dp)

        norm = norm + sin(theta)

        do m = 1, p_lm_max, 1
          do l = m, p_lm_max, 1
            coef(m, l) = coef(m, l) &
            + cmplx(real(out(m + 1)) * p_lm(m, l), 0.0_dp) * sin(theta) &
            * 4 * fourier_PI / 2.0_dp
            coef(m, l) = coef(m, l) &
            + cmplx(0.0_dp, aimag(out(m+1)) * p_lm(m, l)) * sin(theta) &
            * 4 * fourier_PI / 2.0_dp
          end do
        end do

        do m = 0, 0, 1
          do l = m, p_lm_max, 1
            coef(m, l) = coef(m, l) &
            + cmplx(real(out(m + 1)), 0.0_dp) * p_lm(m, l) * sin(theta) &
            * 4 * fourier_PI
            coef(m, l) = coef(m, l) &
            + cmplx(0.0_dp, aimag(out(m+1))) * p_lm(m, l) * sin(theta) &
            * 4 * fourier_PI
          end do
        end do

      end do

      coef = coef / norm

    end subroutine inverse_fourier

end module fourier
