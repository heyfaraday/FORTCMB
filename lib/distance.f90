module distance

  ! Healpix module
  use healpix_types

  implicit none

  contains

    real(kind=dp) function s2(phi1, phi2, theta1, theta2)

      real(kind=dp), intent(in) :: phi1, phi2, theta1, theta2

      s2 = dacos(dcos(theta1) * dcos(theta2) + dsin(theta1) * dsin(theta2) &
      * dcos(phi1 - phi2))

    end function s2


    real(kind=dp) function r2(x1, x2, y1, y2)

      real(kind=dp), intent(in) :: x1, x2, y1, y2

      r2 = dsqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))

    end function r2


    subroutine cross(phi1a, theta1a, phi1b, &
      theta1b, phi2a, theta2a, phi2b, theta2b, phi, theta)

      real(kind=dp), intent(in) :: phi1a, theta1a, phi1b, &
      theta1b, phi2a, theta2a, phi2b, theta2b
      real(kind=dp), intent(out) :: phi, theta

      real(kind=dp), dimension(0:2) :: a1, b1, c1, a2, b2, c2, a
      real(kind=dp) :: mod

      a1 = (/dcos(theta1a) * dcos(phi1a), &
      dcos(theta1a) * dsin(phi1a), dsin(theta1a)/)
      b1 = (/dcos(theta1b) * dcos(phi1b), &
      dcos(theta1b) * dsin(phi1b), dsin(theta1b)/)
      c1 = (/a1(1) * b1(2) - a1(2) * b1(1), &
      a1(2) * b1(0) - a1(0) * b1(2), a1(0) * b1(1) - a1(1) * b1(0)/)

      a2 = (/dcos(theta2a) * dcos(phi2a), &
      dcos(theta2a) * dsin(phi2a), dsin(theta2a)/)
      b2 = (/dcos(theta2b) * dcos(phi2b), &
      dcos(theta2b) * dsin(phi2b), dsin(theta2b)/)
      c2 = (/a2(1) * b2(2) - a2(2) * b2(1), &
      a2(2) * b2(0) - a2(0) * b2(2), a2(0) * b2(1) - a2(1) * b2(0)/)

      a = (/c1(1) * c2(2) - c1(2) * c2(1), &
      c1(2) * c2(0) - c1(0) * c2(2), c1(0) * c2(1) - c1(1) * c2(0)/)
      mod = dsqrt(a(0) * a(0) + a(1) * a(1) + a(2) * a(2))

      phi = datan(a(1) / a(0))
      theta = dasin(a(2) / mod)

    end subroutine

end module distance
