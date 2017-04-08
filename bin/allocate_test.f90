program alocate_test

  use write_in_file, only: write_iteration
  use nr_module, only: four1, gasdev, ran1

  implicit none
  real, dimension(:, :), allocatable :: array
  real, dimension(4) :: a = (/0.0, 0.0, 0.0, 0.0/)
  real :: random_num
  integer :: err, idum_num

  allocate(array(0:4, 0:2), stat=err)
  if (err /= 0) print *, "array: Allocation request denied"

  call four1(a(:), 2, 1)

  idum_num = 0
  write(*, "(4F)") a(:)
  random_num = gasdev(idum_num)
  write(*, *) random_num

  if (allocated(array)) deallocate(array, stat=err)
  if (err /= 0) print *, "array: Deallocation request denied"

end program alocate_test
