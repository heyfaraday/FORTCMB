module routines

  use healpix_types

  implicit none

  contains

    subroutine QuadraticEquationSolver(a, b, c, root1, root2)

      implicit none
      real(kind=dp), intent(in) :: a, b, c
      real(kind=dp), intent(out) :: root1, root2

      real(kind=dp) :: d

      d = b * b - 4.0_dp * a * c

      if ( d >= 0.0_dp ) then

        d = dsqrt(d)

        root1 = (-b + d) / (2.0_dp * a)
        root2 = (-b - d) / (2.0_dp * a)

      else

        write(*, *) "ERROR: There is no real roots!"
        write(*, *) "Discriminant = ", d

      end if

    end subroutine QuadraticEquationSolver


end module routines
