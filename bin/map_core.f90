program map_core

  use write_in_file, only: write_iteration
  use nrtype

  implicit none
  double precision, dimension(:, :), allocatable :: map
  double precision, dimension(:, :), allocatable :: polynom
  integer :: err_map, err_polynom
  integer :: i, j
  integer :: pix_num = 8
  integer :: legendre_num = 4
  double precision :: theta
  integer :: j_map
  double precision :: my_PI = 4.d0 * datan(1.d0)
  integer :: m, l
  real :: ran_num

  ! Allocation
  allocate(map(1:pix_num, 1:pix_num/2), stat=err_map)
  if (err_map /= 0) print *, "map: Allocation request denied"

  allocate(polynom(0:legendre_num, 0:legendre_num), stat=err_polynom)
  if (err_polynom /= 0) print *, "polynom: Allocation request denied"
  ! -

  ! Initialization
  do i = 1, pix_num, 1
    do j = 1, pix_num/2, 1
      map(i, j) = 0.d0
    end do
  end do

  do i = 0, legendre_num, 1
    do j = 0, legendre_num, 1
      polynom(i, j) = 0.d0
    end do
  end do
  ! -

  j_map = 2
  theta = 2.d0 * my_PI * (j_map - 1) / pix_num

  polynom(0, 0) = 1.d0 / sqrt(4.d0 * my_PI)

  do m = 0, legendre_num - 1, 1
    polynom(m + 1, m + 1) = -polynom(m, m)*dsin(theta)*sqrt(2.d0*m+3.d0)/sqrt(2.d0*m+2.d0)
  end do
  do m = 0, legendre_num - 1, 1
    polynom(m, m + 1) = polynom(m, m)*dcos(theta)*sqrt(2.d0*m+3.d0)
  end do
  do m = 0, legendre_num - 2, 1
    do l = m + 2, legendre_num, 1
      polynom(m, l) = ((2.d0*l-1.d0)*dsqrt((l-m)*(2.d0*l+1.d0)) &
      /sqrt((l+m)*(2.d0*l-1.d0))*polynom(m, l-1)*dcos(theta)-(l+m-1.d0) &
      *sqrt((l-m)*(l-1.d0-m)*(2.d0*l+1.d0)) /sqrt((l+m)*(l-1.d0+m) &
      *(2.d0*l-3.d0))*polynom(m, l-2))/(l-m)
    end do
  end do
  do m = 1, legendre_num, 1
    do l = m, legendre_num, 1
      polynom(m, l) = polynom(m, l) * sqrt(2.d0)
    end do
  end do

  ! write(*, *) polynom(:, 2)

  ! Deallocation
  if (allocated(polynom)) deallocate(polynom, stat=err_polynom)
  if (err_polynom /= 0) print *, "polynom: Deallocation request denied"

  if (allocated(map)) deallocate(map, stat=err_map)
  if (err_map /= 0) print *, "map: Deallocation request denied"
  ! -

  do i = 1, 10, 1
    call ran1_s(ran_num)
    write(*, *) ran_num
  end do

end program map_core
