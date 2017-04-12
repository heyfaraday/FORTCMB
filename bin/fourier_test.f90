program fourier_test

  use, intrinsic :: iso_c_binding

  implicit none
  include 'fftw3.f03'

  type(C_PTR) :: plan
  complex(C_DOUBLE_COMPLEX), dimension(1:8) :: in, out

  plan = fftw_plan_dft_1d(8, in, out, FFTW_FORWARD, FFTW_ESTIMATE)

  in(1) = (0.0, 0.d0)
  in(2) = (0.0, 1.d0)
  in(3) = (0.0, 0.d0)
  in(4) = (0.d0, 0.d0)

  write(unit=*, fmt=*) 'in', in
  write(unit=*, fmt=*) 'out', out
  write(unit=*, fmt=*)

  call fftw_execute_dft(plan, in, out)

  write(unit=*, fmt=*) 'in', in
  write(unit=*, fmt=*) 'out', out
  write(unit=*, fmt=*)

  call fftw_destroy_plan(plan)

  plan = fftw_plan_dft_1d(8, in, out, FFTW_BACKWARD, FFTW_ESTIMATE) ! FFTW_MEASURE

  call fftw_execute_dft(plan, out, in)

  write(unit=*, fmt=*) 'in', in
  write(unit=*, fmt=*) 'out', out
  write(unit=*, fmt=*)

  call fftw_destroy_plan(plan)

end program fourier_test
