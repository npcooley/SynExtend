module utilfuncs
  use iso_c_binding
  interface
    function unifRand() result(x) bind(C, name="unif_rand")
      use iso_c_binding, only: c_double
      real(c_double) :: x
    end function unifRand

    function normRand() result(x) bind(C, name="norm_rand")
      use iso_c_binding, only: c_double
      real(c_double) :: x
    end function normRand

    subroutine getRNGseed() bind(C, name="GetRNGstate")
    end subroutine getRNGseed

    subroutine putRNGseed() bind(C, name="PutRNGstate")
    end subroutine putRNGseed
  end interface

contains
  subroutine shuffle_vec(v, l)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    integer(c_int), intent(in) :: l
    integer(c_int), intent(out) :: v(l)

    integer :: i, j
    integer(c_int) :: tmp
    v = [(i, i=1, l)]

    call getRNGseed()
    do i=l, 1, -1
      j = floor((unifRand() * i) + 1.0)
      tmp = v(j)
      v(j) = v(i)
      v(i) = tmp
    end do

    call putRNGseed()
  end subroutine shuffle_vec
end module utilfuncs
