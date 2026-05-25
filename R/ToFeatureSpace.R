###### -- space shift from genomic to feature space ---------------------------

ToFeatureSpace <- function(hit_blocks,
                           qranges,
                           qstrands,
                           sranges,
                           sstrands) {
  res <- vector(mode = "list",
                length = length(hit_blocks))
  for (d1 in seq_along(hit_blocks)) {
    # space shift for query
    if (qstrands[d1] == 1) {
      # forward strand
      if (all(hit_blocks[[d1]][, "QuerySubKey"] == 1)) {
        # no funny business
        curr_q_start <- hit_blocks[[d1]][, "QLeftPos"] - start(qranges[[d1]][1]) + 1L
        curr_q_end <- hit_blocks[[d1]][, "QRightPos"] - start(qranges[[d1]][1]) + 1L
      } else {
        # need to find the interruption space offsets
        curr_q_start <- curr_q_end <- vector(mode = "integer",
                                             length = nrow(hit_blocks[[d1]]))
        ranges_select <- hit_blocks[[d1]][, "QuerySubKey"]
        # don't actually need to subset anything here because I'll never actually
        # index to that position
        prior_widths <- c(0, width(qranges[[d1]]))
        prior_widths <- cumsum(prior_widths)
        for (d2 in seq_along(curr_q_start)) {
          curr_q_start[d2] <- hit_blocks[[d1]][d2, "QLeftPos"] - start(qranges[[d1]][ranges_select[d2]]) + 1L
          curr_q_end[d2] <- hit_blocks[[d1]][d2, "QRightPos"] - start(qranges[[d1]][ranges_select[d2]]) + 1L
        }
        curr_q_start <- curr_q_start + prior_widths[ranges_select]
        curr_q_end <- curr_q_end + prior_widths[ranges_select]
      } # end subkey check
    } else {
      # reverse strand
      # using the subkey here is a bit funny because adjustments need to happen
      # in descending order now
      range_len <- length(qranges[[d1]])
      if (all(hit_blocks[[d1]][, "QuerySubKey"] == range_len)) {
        # no funny business
        curr_q_start <- end(qranges[[d1]][range_len]) - hit_blocks[[d1]][, "QRightPos"] + 1L
        # curr_q_start <- rev(curr_q_start)
        curr_q_end <- end(qranges[[d1]][range_len]) - hit_blocks[[d1]][, "QLeftPos"] + 1L
        # curr_q_end <- rev(curr_q_end)
      } else {
        # all the funny business
        curr_q_start <- curr_q_end <- vector(mode = "integer",
                                             length = nrow(hit_blocks[[d1]]))
        ranges_select <- hit_blocks[[d1]][, "QuerySubKey"]
        # don't actually need to subset out the last position here because I'll
        # never actually index to that position
        prior_widths <- c(width(qranges[[d1]]), 0)
        # reverse to grab get the cumulative sum correctly, and then reverse back
        prior_widths <- rev(cumsum(rev(prior_widths)))
        for (d2 in seq_along(curr_q_start)) {
          curr_q_start[d2] <- end(qranges[[d1]][ranges_select[d2]]) - hit_blocks[[d1]][d2, "QRightPos"] + 1L
          curr_q_end[d2] <- end(qranges[[d1]][ranges_select[d2]]) - hit_blocks[[d1]][d2, "QLeftPos"] + 1L
        }
        curr_q_start <- curr_q_start + prior_widths[ranges_select]
        # curr_q_start <- rev(curr_q_start)
        curr_q_end <- curr_q_end + prior_widths[ranges_select]
        # curr_q_end <- rev(curr_q_end)
      }
      
    } # end qstrand check
    
    # space shift for subject
    if (sstrands[d1] == 1) {
      # forward strand
      if (all(hit_blocks[[d1]][, "SubjectSubKey"] == 1)) {
        # no funny business
        curr_s_start <- hit_blocks[[d1]][, "SLeftPos"] - start(sranges[[d1]][1]) + 1L
        curr_s_end <- hit_blocks[[d1]][, "SRightPos"] - start(sranges[[d1]][1]) + 1L
      } else {
        # need to find the interruption space offsets
        curr_s_start <- curr_s_end <- vector(mode = "integer",
                                             length = nrow(hit_blocks[[d1]]))
        ranges_select <- hit_blocks[[d1]][, "SubjectSubKey"]
        # don't actually need to subset anything here because I'll never actually
        # index to that position
        prior_widths <- c(0, width(sranges[[d1]]))
        prior_widths <- cumsum(prior_widths)
        for (d2 in seq_along(curr_s_start)) {
          curr_s_start[d2] <- hit_blocks[[d1]][d2, "SLeftPos"] - start(sranges[[d1]][ranges_select[d2]]) + 1L
          curr_s_end[d2] <- hit_blocks[[d1]][d2, "SRightPos"] - start(sranges[[d1]][ranges_select[d2]]) + 1L
        }
        curr_s_start <- curr_s_start + prior_widths[ranges_select]
        curr_s_end <- curr_s_end + prior_widths[ranges_select]
      } # end subkey check
    } else {
      # reverse strand
      # using the subkey here is a bit funny because adjustments need to happen
      # in descending order now
      range_len <- length(sranges[[d1]])
      if (all(hit_blocks[[d1]][, "SubjectSubKey"] == range_len)) {
        # no funny business
        curr_s_start <- end(sranges[[d1]][range_len]) - hit_blocks[[d1]][, "SRightPos"] + 1L
        # curr_s_start <- rev(curr_s_start)
        curr_s_end <- end(sranges[[d1]][range_len]) - hit_blocks[[d1]][, "SLeftPos"] + 1L
        # curr_s_end <- rev(curr_s_end)
      } else {
        # all the funny business
        curr_s_start <- curr_s_end <- vector(mode = "integer",
                                             length = nrow(hit_blocks[[d1]]))
        ranges_select <- hit_blocks[[d1]][, "SubjectSubKey"]
        # don't actually need to subset out the last position here because I'll
        # never actually index to that position
        prior_widths <- c(width(sranges[[d1]]), 0)
        # reverse to grab get the cumulative sum correctly, and then reverse back
        prior_widths <- rev(cumsum(rev(prior_widths)))
        for (d2 in seq_along(curr_s_start)) {
          curr_s_start[d2] <- end(sranges[[d1]][ranges_select[d2]]) - hit_blocks[[d1]][d2, "SRightPos"] + 1L
          curr_s_end[d2] <- end(sranges[[d1]][ranges_select[d2]]) - hit_blocks[[d1]][d2, "SLeftPos"] + 1L
        }
        curr_s_start <- curr_s_start + prior_widths[ranges_select]
        # curr_s_start <- rev(curr_s_start)
        curr_s_end <- curr_s_end + prior_widths[ranges_select]
        # curr_s_end <- rev(curr_s_end)
      }
    } # end strand check
    # if hits are marching through the features in opposite directions -- in
    # feature space then they are non-sensical for alignment
    # there could be scenarios where blocks of sensical hits can be salvaged
    # here but i'm just not in the mood to figure that out right now
    o1 <- order(curr_q_start)
    o2 <- order(curr_s_start)
    if (any(o1 != o2)) {
      res[[d1]] <- matrix(data = integer(),
                          nrow = 4)
    } else {
      ph <- rbind(curr_q_start,
                  curr_q_end,
                  curr_s_start,
                  curr_s_end)
      rownames(ph) <- NULL
      ph <- ph[, o1, drop = FALSE]
      res[[d1]] <- ph
    }
    
  } # end d1 loop
  return(res)
}

