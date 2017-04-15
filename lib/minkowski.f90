module minkowski

  use healpix_types

  implicit none
  real(kind=dp), parameter :: minkowski_PI = 4.0_dp * datan(1.0_dp)

  contains

    real(kind=dp) function area(n_max, map, level)

      ! area = area for parameter 'level'
      ! n_max = number of pixels for different phi
      ! map = array for reading map
      ! level = F_1(level) - minkowski parameter

      implicit none
      integer(kind=i4b), intent(in) :: n_max
      real(kind=dp), dimension(1:n_max+1, 1:n_max/2+1), intent(in) :: map
      real(kind=dp), intent(in) :: level

      integer(kind=i8b) :: i, j ! Vars for iterating in the map
      real(kind=dp) :: theta ! Theta as parameter for polynoms and map
      real(kind=dp) :: mean ! mean value in pixel
      ! Area without normalization and normalization factor
      real(kind=dp) :: a, na

      do j = 2, n_max / 2 , 1

        theta = 2.0_dp * minkowski_PI * (j - 1) / n_max

        do i = 1, n_max, 1

          mean = (map(i, j))
          ! mean = (map(i, j) + map(i, j + 1) &
          ! + map(i + 1, j) + map(i + 1, j + 1)) / 4.0_dp

          if ( mean > level ) then
            a = a + dsin(theta)
          end if

          na = na + dsin(theta)

        end do
      end do

      area = a / na

    end function area


    real(kind=dp) function length(n_max, map, level)

      use distance, only: s2

      ! length = length for parameter 'level'
      ! n_max = number of pixels for different phi
      ! map = array for reading map
      ! level = F_2(level) - minkowski parameter

      implicit none
      integer(kind=i8b), intent(in) :: n_max
      real(kind=dp), dimension(1:n_max+1, 1:n_max/2+1), &
      intent(out) :: map
      real(kind=dp), intent(in) :: level

      integer(kind=i8b) i, j ! Vars for iterating in the map
      real(kind=dp) :: h_phi, h_theta !
      real(kind=dp) :: l = 0.0_dp! Lattice period
      real(kind=dp) :: phi1, theta1, phi2, theta2 ! Intersection coordinates
      real(kind=dp) :: phi, phi_1, theta, theta_1 ! Pixels coordinates

      map = map - level

      do j = 2, n_max/2 - 1, 1

        h_theta = 2.0_dp * minkowski_PI / n_max
        theta = (j - 1) * h_theta
        theta_1 = (j) * h_theta

        do i = 1, n_max, 1

          h_phi = 2.0_dp * minkowski_PI / n_max
          phi = (i - 1) * h_phi
          phi_1 = (i) * h_phi

          if ( map(i, j) == 0.0_dp .and. &
          (map(i, j + 1) == 0.0_dp .or. map(i + 1, j) == 0.0_dp)) then
            l = l + 0.0_dp
          end if

          if ( map(i, j) * map(i, j + 1) < 0.0_dp) then

            if ( map(i, j) * map(i + 1, j) < 0.0_dp) then

              phi1 = phi
              theta1 = theta + h_theta * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i, j + 1)))

              phi2 = phi + h_phi * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i + 1, j)))
              theta2 = theta

              l = l + s2(phi1, phi2, theta1, theta2)

            elseif ( map(i + 1, j) * map(i + 1, j + 1) < 0.0_dp) then

              phi1 = phi
              theta1 = theta + h_theta * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i, j + 1)))

              phi2 = phi_1
              theta2 = theta + h_theta * abs(map(i + 1, j)) &
              / (abs(map(i + 1, j)) + abs(map(i + 1, j + 1)))

              l = l + s2(phi1, phi2, theta1, theta2)

            else if ( map(i, j + 1) * map(i + 1, j + 1) < 0.0_dp) then

              phi1 = phi
              theta1 = theta + h_theta * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i, j + 1)))

              phi2 = phi + h_phi * abs(map(i, j + 1)) &
              / (abs(map(i, j + 1)) + abs(map(i + 1, j + 1)))
              theta2 = theta_1

              l = l + s2(phi1, phi2, theta1, theta2)

            end if

          else if ( map(i, j) * map(i + 1, j) < 0.0_dp) then

            if ( map(i + 1, j) * map(i + 1, j + 1) < 0.0_dp) then

              phi1 = phi + h_phi * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i + 1, j)))
              theta1 = theta

              phi2 = phi_1
              theta2 = theta + h_theta * abs(map(i + 1, j)) &
              / (abs(map(i + 1, j)) + abs(map(i + 1, j + 1)))

              l = l + s2(phi1, phi2, theta1, theta2)

            else if ( map(i, j + 1) * map(i + 1, j + 1) < 0.0_dp) then

              phi1 = phi + h_phi * abs(map(i, j)) &
              / (abs(map(i, j)) + abs(map(i + 1, j)))
              theta1 = theta

              phi2 = phi + h_phi * abs(map(i, j + 1)) &
              / (abs(map(i, j + 1)) + abs(map(i + 1, j + 1)))
              theta2 = theta_1

              l = l + s2(phi1, phi2, theta1, theta2)

            end if

          else if ( map(i + 1, j) * map(i + 1, j + 1) < 0.0_dp) then

            if ( map(i, j + 1) * map(i + 1, j + 1) < 0.0_dp) then

              phi1 = phi_1
              theta1 = theta + h_theta * abs(map(i + 1, j)) &
              / (abs(map(i + 1, j)) + abs(map(i + 1, j + 1)))

              phi2 = phi + h_phi * abs(map(i, j + 1)) &
              / (abs(map(i, j + 1)) + abs(map(i + 1, j + 1)))
              theta2 = theta_1

              l = l + s2(phi1, phi2, theta1, theta2)

            end if

          end if

        end do
      end do

    map = map + level
    length = l / 4.0_dp / minkowski_PI

    end function length


    integer(kind=i8b) function condition_1(xx, yy, xy)

      implicit none
      real(kind=dp), intent(in) :: xx, yy, xy

      if (( xx * yy - xy * xy >= 0.0_dp .and. xx >= 0.0_dp ) &
      .or. (xx * yy - xy * xy >= 0.0_dp .and. yy >= 0.0_dp )) then
        condition_1 = 0
      else
        if (( xx * yy - xy * xy >= 0.0_dp .and. 0.0_dp > xx ) &
        .or. (xx * yy - xy * xy >= 0.0_dp .and. 0.0_dp > yy )) then
          condition_1 = 2
        else
          condition_1 = 1
        end if
      end if

    end function condition_1

    integer(kind=i8b) function condition_2(xx, yy, xy)

      use routines, only: QuadraticEquationSolver

      implicit none
      real(kind=dp), intent(in) :: xx, yy, xy

      if (( xx * yy - xy * xy >= 0.0_dp .and. xx >= 0.0_dp ) &
      .or. (xx * yy - xy * xy >= 0.0_dp .and. yy >= 0.0_dp )) then
          condition_2 = 0
      else
        if (( xx * yy - xy * xy >= 0.0_dp .and. 0.0_dp > xx ) &
        .or. (xx * yy - xy * xy >= 0.0_dp .and. 0.0_dp > yy )) then
          condition_2 = 2
          else
            condition_2 = 1
        end if
      end if

    end function condition_2

end module minkowski
