module cart_methods
  implicit none
  private
  ! any exported subroutines should be included here
  ! public ...
  public reorder_matrix

contains
  ! make sure to purify all your subroutines before finalizing
  pure subroutine reorder_matrix(mat, r, c, split_col, split_val, split_point) bind(C, name="f_reorder_matrix")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: split_col
    integer(c_int), intent(in) :: r, c
    real(c_double), intent(in) :: split_val
    real(c_double), intent(inout) :: mat(r,c)
    integer(c_int), intent(out) :: split_point

    logical :: vec(r)
    integer(c_int) :: i

    vec = mat(:,split_col) > split_val

    !dir$ ivdep
    do i=1,c
      mat(:,i) = [pack(mat(:,i), .not. vec), pack(mat(:,i), vec)]
    end do
    split_point = count(vec)
  end subroutine reorder_matrix

  ! subroutine find_split(mat, r, c, num_to_check, out_column, out_value)
  !   ! I think I'm going to do this directly from C
  !   use, intrinsic :: iso_c_binding, only: c_int, c_double
  !   implicit none

  !   integer(c_int), intent(in) :: r, c, num_to_check
  !   real(c_double), intent(in) :: mat(r,c)
  !   integer(c_int), intent(out) :: out_column
  !   real(c_double), intent(out) :: out_value

  !   integer(c_int) :: to_check(c-1), i
  !   real(c_double) :: scores(num_to_check)

  !   ! determine which columns we'll check
  !   call shuffle_vec(to_check, c-1)

  !   !do i=1,num_to_check
  !   !  call find_gini_split(mat(:,to_check(i)), mat(:,c), r, scores, i, c-1)
  !   !end do
  ! end subroutine find_split


  subroutine gini_imp(classes, l, nclass, o_v)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: l, nclass
    integer(c_int), intent(in) :: classes(l)
    real(c_double), intent(out) :: o_v

    real(c_double) :: class_counts(nclass), sum_amt
    integer(c_int) :: i, j
    if(l == 0) then
      o_v = 1.0
      return
    end if

    ! i have a feeling this is slow
    ! do i=1, nclass
    !   class_counts(i) = 0.0+count(classes==i) ! cast to double for later
    ! end do

    ! total is always just l
    ! total = sum(class_counts)
    class_counts(:) = 0.0
    sum_amt = 1.0/l
    do i=1, l
      j = classes(i)
      class_counts(j) = class_counts(j) + sum_amt
    end do
    o_v = 1.0-sum(class_counts**2)
  end subroutine gini_imp

  subroutine double_gini_imp(classes, l, nclass, mask, mask_count, o_v)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: l, nclass, mask_count
    integer(c_int), intent(in) :: classes(l)
    logical, intent(in) :: mask(l)
    real(c_double), intent(out) :: o_v

    real(c_double) :: cc_left(nclass), cc_right(nclass), sum_amt, divisor
    integer(c_int) :: i, j

    cc_left(:) = 0.0
    cc_right(:) = 0.0
    sum_amt = 1.0
    do i=1, l
      j = classes(i)
      if(mask(i)) then
        cc_left(j) = cc_left(j) + 1.0
      else
        cc_right(j) = cc_right(j) + 1.0
      end if
    end do

    if(mask_count .ne. 0) cc_left = cc_left / mask_count
    if(mask_count .ne. l) cc_right = cc_right / (l-mask_count)

    ! return weighted gini impurity
    divisor = real(mask_count) / l
    o_v = divisor*(1.0-sum(cc_left**2)) + (1-divisor)*(1.0-sum(cc_right**2))
  end subroutine double_gini_imp

  subroutine find_gini_sim_anneal(v, response, l, nclass, tempmax, o_v, o_gini_score)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use :: utilfuncs
    implicit none
    integer(c_int), intent(in) :: l, nclass, tempmax
    integer(c_int), intent(in) :: response(l)
    real(c_double), intent(in) :: v(l)
    real(c_double), intent(out) :: o_gini_score, o_v


    ! local variables
    real(c_double) :: cur_score, cur_thresh, roll, temp, new_thresh, new_score, scale_fac, vmax, vmin
    integer(c_int) :: i
    logical :: shouldSwap, tmpmask(l)

    ! If we have too many data points, use simulated annealing
    ! This works really well regardless of the number of iterations -- beating RF with only 10
    vmax = maxval(v)
    vmin = minval(v)
    cur_thresh = (vmax+vmin)/2.0
    tmpmask(:) = v <= cur_thresh
    call double_gini_imp(response, l, nclass, tmpmask, count(tmpmask), cur_score)

    ! minimize the weighted gini importance, corresponds to maximizing gini gain
    scale_fac = 0.33 * (vmax-vmin)
    call getRNGseed()
    do i=1, tempmax
      shouldSwap = .false.
      temp = 1 - ((i-1) / tempmax)
      ! 1. generate a random new candidate
      !roll = normRand() * scale_fac
      roll = unifRand() * scale_fac
      new_thresh = cur_thresh + roll

      ! if the roll is out of bounds, I'm just going to reflect it back in range
      if(new_thresh > vmax) new_thresh = 2*vmax - new_thresh
      if(new_thresh < vmin) new_thresh = 2*vmin - new_thresh
      tmpmask(:) = v <= new_thresh

      if(new_thresh == cur_thresh) cycle
      call double_gini_imp(response, l, nclass, tmpmask, count(tmpmask), new_score)

      ! 2. Probabalistically move to the new state based on temperature
      if(new_score < cur_score) then
        shouldSwap = .true.
      else
        temp = exp((cur_score-new_score) / temp)
        roll = unifRand()
        if(roll <= temp) shouldSwap = .true.
      end if

      if(shouldSwap) then
        cur_score = new_score
        cur_thresh = new_thresh
      end if
    end do
    call putRNGseed()

    o_v = cur_thresh
    o_gini_score = cur_score
  end subroutine find_gini_sim_anneal

  subroutine find_gini_split(v, response, l, nclass, o_v, o_gini_score) bind(C, name="find_gini_split_")
    ! Here I'm going to assume that scores are INTEGERS on scale 1:n
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    integer(c_int), intent(in) :: l, nclass
    integer(c_int), intent(in) :: response(l)
    real(c_double), intent(in) :: v(l)
    real(c_double), intent(out) :: o_gini_score, o_v

    ! local variables
    integer(c_int) :: i, j, mloc, tempmax
    real(c_double) :: total_gini, gains(l), tmpscore
    logical :: tmpmask(l)

    tempmax = 10

    ! calculate the base gini impurity
    call gini_imp(response, l, nclass, total_gini)

    ! if the data are small, just calculate gain for every possible split point
    if(l <= tempmax) then
      gains(:) = total_gini
      do i=1, l
        tmpmask(:) = v <= v(i)
        j = count(tmpmask)
        if(j == l) then
          gains(i) = -1.0
        else
          call double_gini_imp(response, l, nclass, tmpmask, j, total_gini)
          gains(i) = gains(i) - total_gini
        end if
      end do

      mloc = maxloc(gains, dim=1)
      o_v = v(mloc)
      o_gini_score = gains(mloc)
    else
      ! otherwise, use simulated annealing
      call find_gini_sim_anneal(v, response, l, nclass, tempmax, o_v, tmpscore)
      o_gini_score = total_gini - tmpscore
    end if
  end subroutine find_gini_split

  subroutine find_sse_split(v, response, l, o_v, o_sse_score) bind(C, name="find_sse_split_")
    ! scores here are DOUBLES
    ! SSE calculation is much easier than gini
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    integer(c_int), intent(in) :: l
    real(c_double), intent(in) :: v(l), response(l)
    real(c_double), intent(out) :: o_sse_score, o_v

    ! local variables
    integer(c_int) :: i, j, mloc, tempmax
    real(c_double) :: total_sse, gains(l), tmpscore, lmean, rmean
    logical :: tmpmask(l)

    tempmax = 10

    ! original sse is just the mean squared error
    total_sse = sum(((response - sum(response)/l))**2)

    ! eventually I should do something smarter here
    ! is it even that slow though to just check everything?

    ! small number, just check all splits
    gains(:) = total_sse
    do i=1, l
      tmpmask(:) = v <= v(i)
      j = count(tmpmask)
      if(j == l) then
        gains(i) = -1.0
      else ! sse for both sides
        ! first get the mean for each side
        lmean = sum(response, mask=tmpmask) / j
        rmean = sum(response, mask=(.not. tmpmask)) / (l-j)
        ! then get the sse
        lmean = sum((response-lmean)**2, mask=tmpmask)
        rmean = sum((response-rmean)**2, mask=(.not. tmpmask))
        tmpscore = lmean + rmean
        gains(i) = total_sse - tmpscore
      end if
    end do



    ! get result
    mloc = maxloc(gains, dim=1)
    o_v = v(mloc)
    o_sse_score = gains(mloc)
  end subroutine find_sse_split

end module cart_methods
