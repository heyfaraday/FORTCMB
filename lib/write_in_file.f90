module write_in_file

  implicit none

  contains

    subroutine write_iteration(out_num, t, r, v)

      implicit none
      integer, intent(in) :: out_num
      double precision, intent(in) :: t
      double precision, dimension(1:3), intent(in) :: r, v

      write(unit=out_num, fmt="(F)", advance='no') t
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') r(1)
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') r(2)
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') r(3)
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') v(1)
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') v(2)
      write(unit=out_num, fmt="(A)", advance='no') ','
      write(unit=out_num, fmt="(F)", advance='no') v(3)
      write(unit=out_num, fmt="(A)", advance='no') NEW_LINE('F')

    end subroutine write_iteration

end module write_in_file
