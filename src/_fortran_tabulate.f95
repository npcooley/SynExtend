module tabulate_mod
  implicit none
contains
  pure subroutine tabulate_double(v, l, out_val, out_count, ctr) bind(C, name="tabulate_double_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    integer(c_int), intent(in) :: l
    real(c_double), intent(in) :: v(l)
    integer(c_int), intent(out) :: out_count(l), ctr
    real(c_double), intent(out) :: out_val(l)

    integer(c_int) :: j, k

    ctr = 1
    out_val(:) = -1.0
    out_count(:) = 0
    out_val(1) = v(1)
    out_count(1) = 1

    outer: do j=2, l
      inner: do k=1, ctr
        if (v(j) == out_val(k)) then
          out_count(k) = out_count(k) + 1
          cycle outer
        end if
      end do inner
      ctr = ctr+1
      out_val(ctr) = v(j)
      out_count(ctr) = 1
    end do outer
  end subroutine tabulate_double

  pure subroutine tabulate_int(v, l, out_val, out_count, ctr) bind(C, name="tabulate_int_")
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none

    integer(c_int), intent(in) :: l
    integer(c_int), intent(in) :: v(l)
    integer(c_int), intent(out) :: out_count(l), out_val(l), ctr

    integer(c_int) :: j, k
    ctr = 1
    out_val(:) = -1.0
    out_count(:) = 0
    out_val(1) = v(1)
    out_count(1) = 1

    outer: do j=2, l
      inner: do k=1, ctr
        if (v(j) == out_val(k)) then
          out_count(k) = out_count(k) + 1
          cycle outer
        end if
      end do inner
      ctr = ctr+1
      out_val(ctr) = v(j)
      out_count(ctr) = 1
    end do outer
  end subroutine tabulate_int

end module tabulate_mod
