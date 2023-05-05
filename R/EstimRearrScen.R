##### -- Function to estimate genome rearrangement scenarios -------------------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

##########
#Function to generate most likely rearrangement scenario
#Expects as input a synteny object
#Returns a list containing the rearrangement steps in order

##########
EstimRearrScen <- function(SyntenyObject,
                                           NumRuns = -1,
                                           Mean = FALSE,
                                           MinBlockLength = -1,
                                           Verbose = TRUE) {
  if (!is(object = SyntenyObject,
          class2 = "Synteny")) {
    stop ("Expected class of type 'Synteny'.")
  }
  starttime <- Sys.time()
  
  rearrange_chromosome <- function(synblocks, syn, NumRuns=-1,
                                   Verbose = FALSE, Mean = FALSE, #main user parameters
                                   MinBlockLength = -1, #extra user parameters
                                  first_gen=1, second_gen=2, chrom=-1, EarlyOutput=FALSE){ #internal parameters  
    # ERROR CHECKING  --------------------------------------------------------------
    #set.seed(13)
    stopifnot("Error: Incorrect genome indices provided. Make sure indices are not the same."=
                second_gen!=first_gen)
    if ( second_gen < first_gen ){
      temp <- second_gen
      second_gen <- first_gen
      first_gen <- temp
    }
    
    #case where both genomes are identical or have a single synteny block causes problems
    #this calculates the number of unique genome sequences and then returns before the main
    #function to avoid any potential errors
    if(nrow(syn[[second_gen,first_gen]]) == 1){
      len_1 <- syn[[1,1]]
      len_2 <- syn[[2,2]]
      start_1 <- SyntenyObject[[2,1]][,5]
      start_2 <- SyntenyObject[[2,1]][,6]
      end_1 <- SyntenyObject[[2,1]][,7]
      end_2 <- SyntenyObject[[2,1]][,8]
      
      ins_del_trans <- 0
      if(start_1 != 1)
        ins_del_trans <- 1
      if(start_2 != 1)
        ins_del_trans <- ins_del_trans + 1
      if(end_1 != len_1)
        ins_del_trans <- ins_del_trans + 1
      if(end_2 != len_2)
        ins_del_trans <- ins_del_trans + 1
      
      if(EarlyOutput)
        return(c(0, 0, ins_del_trans))
      rearrangements <- c('Original: 1')
      block_key <- matrix(c(1,1,1,SyntenyObject[[1,1]], 1), nrow=1)
      colnames(block_key) <- c('start1', 'start2', 'length', 'rel_direction_on_2', 'index1')
      
      return(list('counts'=c(0, 0, ins_del_trans),
                  'scenario'=rearrangements,
                  'block_key'=block_key))
    }
    
    # END ERROR CHECKING------------------------------------------------------------
    
    # HELPER FUNCTIONS -------------------------------------------------------------
    
    ## Adds connection between two vertices of graph
    add_connection <- function(mat, i1, i2, col){
      mat[i1, col] <- i2
      mat[i2, col] <- i1
      return(mat)
    }
    
    ## Double cut and join operation
    DCJ <- function(graph, v1, v2){
      ##v1 and v2 are joined by a gray line
      v1b <- graph[v1, 1] #black connection of v1
      v2b <- graph[v2, 1] #black connection of v2
      
      graph <- add_connection(graph, v1, v2, 1)
      graph <- add_connection(graph, v1b, v2b, 1)
      
      return(graph)
    }
    
    ## Convert adjacency graph to genome (as vector)
    graph_to_genome <- function(graph){
      start <- nrow(graph)
      cur <- graph[start, 1]
      genome <- c()
      
      while (cur != (start - 1)){
        genome <- c(genome, (cur + (cur%%2)) / 2 * graph[cur,3])
        cur <- graph[cur + graph[cur,3], 1]
      }
      
      return(genome)
    }
    
    condense_genome <- function(genome){
      i <- 1
      while(i <= length(genome)){
        start <- i
        while(genome[i] + 1 == genome[i+1]){
          i <- i+1
          genome <- genome[-(i)]
        }
        if (start == i)
          i <- i+1
      }
      return(genome)
    }
    
    #faster and more memory efficient way to append to a list
    lappend <- function(lst, val, increase=10){
      
      #find first empty value
      first_empty <- length(lst) + 1
      for(i in seq_len(length(lst))){
        if(isEmpty(lst[[i]])){
          first_empty <- i
          break()
        }
      }
      
      #get current length of list
      len <- length(lst)
      tmp <- len
      
      #find new length to be able to add new value(s)
      #should overestimate slightly, this is intentional to reduce calls to this function
      while((tmp - first_empty) < length(val)){
        tmp <- tmp + increase
      }
      
      #if we increased the size of the list, copy the old list into the new one  
      if(tmp != len){
        new_lst <- vector("list", length = tmp)
        new_lst[seq_len(len)] <- lst
        lst <- new_lst
      }
      
      #add new value to list
      lst[first_empty:(first_empty+length(val)-1)] <- val 
      return(lst)
    }
    
    update_direc <- function(graph){
      start <- nrow(graph) #the last entry is the begin point for the mixed up genome
      cur <- graph[start, 1] #travel to the vertex connected by a black line
      
      #note here that the genome is represented as a sequence of permutations with a direction
      #ex. < -3 2 -4 5 1 > becomes < +3 -3 -2 +2 +4 -4 -5 +5 -1 +1 >, where genes are read from negative to positive
      #in the graph, instead of representing the ends of block n as -n and +n, we use 2n-1 and 2n
      #so the program represents the above as < 6 5 3 4 7 8 9 10 1 2 >
      #thus is a gene is even -> odd it's reversed, and if it's odd -> even it's in the correct orientation
      
      
      while (cur != (start-1)){
        if(cur%%2 == 1){ #if value is odd, gene is in the correct orientation (+)
          #we store a direction of 1 in both edges of the gene, and set cur to the other side of the block
          graph[cur, 3] <- 1 
          cur <- cur + 1
          graph[cur, 3] <- 1
        }
        else{ #if value is even, gene is in the reversed orientation (-)
          #same as above, except setting direction to -1
          graph[cur, 3] <- -1
          cur <- cur-1
          graph[cur, 3] <- -1
        }
        
        #follow the black path to the next gene in the mixed up genome
        cur <- graph[cur, 1]
      }
      
      
      #This may change the end caps' direction, which is not allowed
      #this fixes it
      graph[nrow(graph), 3] <- 1
      graph[nrow(graph)-1, 3] <- 1
      return(graph)
    }
    
    
    generate_rearrangements <- function(graph, bi_skip_prob=0, inv_skip_prob=0){
      rearrangements <- vector("list", length=nrow(graph))
      rearrangements[[1]] <- paste("Original: ", paste(graph_to_genome(graph), collapse=" "))
      block_count <- 0
      invert_count <- 0
      gen <- graph_to_genome(graph)
      t <- 0
      len <- length(gen) + 1
      gen <- c(0, gen, len)
      done <- FALSE
      while(!identical(as.integer(gen), 0:(length(gen)-1))){
        t <- t + 1
        if (t == (3*len)) {
          #The number of operations is at most double the number of blocks
          #(as we could sort any vector into place by transposing each block one at a time and then inverting it)
          #Thus if we ever get to triple the number of blocks, we can be sure it's in an infinite loop
          fname <- "./ERRONEOUS_SYNT_OBJ.RData"
          save(SyntenyObject, file=fname)
          err_msg <- paste("Infinite Loop, bug in code. Error-producing synteny object saved to ", fname, sep="")
          stop(err_msg)
        }
        
        gen <- graph_to_genome(graph)
        len <- length(gen) + 1
        gen <- c(0, gen, len)
        done <- TRUE
        order <- sample(seq_len(len))
        for(i in order){
          sign1 <- sign(gen[i])
          sign2 <- sign(gen[i+1])
          if (sign1==0) sign1 <- 1
          
          #block interchange====
          for(p in c(1,-1)){
            if (done == FALSE)
              next()
            if (p == 1)
              subsec <- gen[i:len]
            else
              subsec <- gen[seq_len(i-1)]
            
            bi_skip_roll <- runif(1) # naive way to correct for overestimation of BI
            if((bi_skip_roll > bi_skip_prob) && ((gen[i]+p) %in% subsec) && (gen[i] + p) != gen[i+p]){
              #no block interchanges to connect inverted elements
              #getting length of interchange
              j <- i+p
              len_block <- 0
              cont <- FALSE
              while(gen[j] != (gen[i]+p)){
                len_block <- len_block + 1
                if(gen[j] == gen[length(gen)]){
                  cont <- TRUE
                  break()
                }
                j <- j+p
              }
              if(cont) next()
              #ex. 0 1 4 3 2 5, block interchange required
              #first cut creates the circular intermediate
              if((gen[i] < 0 && p==1) || (gen[i] > 0 && p==-1))
                v1 <-  abs(gen[i])*2 - 1
              else
                v1 <- abs(gen[i])*2
              
              new_graph <- DCJ(graph, v1, graph[v1,2])
              circular_intermediate <- sort(gen[(!gen%in%(graph_to_genome(new_graph)))])
              #we now have 0 1 2 5, with a CI of (4 3)
              #take genome as a vector and sort it, adding 0 to the front
              cut_vertex <- -1
              
              for(vert in circular_intermediate){
                if(vert == 0 || vert==max(circular_intermediate))
                  next()
                if((vert-1)%in%gen && !(vert-1)%in%circular_intermediate){
                  if (vert < 0){
                    cut_vertex <- (abs(vert))*2
                    break()
                  }
                  else{
                    cut_vertex <- (vert)*2-1
                    break()
                  }
                }
                if((vert + 1)%in%gen && !(vert+1)%in%circular_intermediate){
                  if (vert < 0){
                    cut_vertex <- (abs(vert))*2 -1
                    break()
                  }
                  else{
                    cut_vertex <- (vert)*2 
                    break()
                  }
                }
                
                if((vert-1) %in% gen && !(vert-1)%in%circular_intermediate){
                  if(vert < 0){
                    cut_vertex <- abs(vert)*2
                    break()
                  }
                  else{
                    cut_vertex <- (vert)*2 - 1
                    break()
                  }
                }
              }
              
              if(cut_vertex != -1){
                graph <- new_graph
                #then we just do a DCJ on the left vertex of this block and its connected gray vertex
                graph <- DCJ(graph, cut_vertex, graph[cut_vertex, 2])
                rearrangements[[block_count + invert_count +2]] <- paste("transposition: ", paste(graph_to_genome(graph), collapse=" "), "{", len_block, "}")
                block_count <- block_count + 1
                done <- FALSE
                break()   
              }
              else {
                gen <- graph_to_genome(graph)
                gen <- c(0, gen, length(gen) + 1)
              }
            }
            #end block interchange====
            
            inv_skip_roll <- runif(1)
            
            #inversion====
            if ((inv_skip_prob < inv_skip_roll) && !(p==-1 && i==1) && gen[i] < 0 && ((-1*(gen[i]+p)) %in% gen)){
              v1 <- -1
              if(gen[i]*p > 0)
                v1 <- (abs(gen[i]*2)) 
              else 
                v1 <- (abs(gen[i])*2-1)
              
              #ex. 0 1 -3 -2 4
              i1 <- i
              i2 <- match((-1*(gen[i]+p)), gen)
              len_block <- abs(i1 - i2)
              
              if(v1 != -1){
                graph <- DCJ(graph, v1, graph[v1, 2])
                graph <- update_direc(graph)
                rearrangements[[invert_count + block_count + 2]] <- paste("inversion: ", paste(graph_to_genome(graph), collapse=" "), "{", len_block, "}")
                invert_count <- invert_count + 1
                done <- FALSE
                break()
              }
            }
          }
          #end inversion====
        }
      }
      total_ops <- block_count + invert_count
      
      rearrangements <- rearrangements[-((total_ops+2):length(rearrangements))]
      return(list("rearrangements"=rearrangements, "graph"=graph, "invert_count"=invert_count, "block_count"=block_count))
    }
    
    
    # END HELPER FUNCTIONS----------------------------------------------------------
    
    
    # FUNCTION BODY # 
    
    # SYNTENY OBJECT TO VECTOR OF PERMUTATIONS 
    
    # ==== Synteny Object to Blocks ====
    block_matrix <- matrix(nrow=nrow(synblocks), ncol=4) #setting up the matrix to eventually return
    
    for(rowindex in seq_len(nrow(synblocks))){
      row <- synblocks[rowindex,] #grab information on the block
      
      start1 <- row[5] #start of block on the first genome
      blocklength <- row[7] - row[5] #length of block
      
      start2 <- row[6]
      #get direction of the block
      if(row[3] == 0){ 
        dir <- 1
      }
      else{ 
        dir <- -1
      }
      if (blocklength > MinBlockLength){ #drop blocks that are SNPs
        block_matrix[rowindex, 1] <- start1
        block_matrix[rowindex, 2] <- start2
        block_matrix[rowindex, 3] <- blocklength
        block_matrix[rowindex, 4] <- dir
      }
    }
    
    block_matrix <- block_matrix[complete.cases(block_matrix),]
    
    # ==== end synteny object to blocks ====
    
    # ==== Blocks to Permutation Matrix ====
    sorted_mat <- block_matrix
    
    if (!is(sorted_mat, 'array')){
      return(list('counts'=c(0, 0, 0),
                  'scenario'='none',
                  'block_key'=as.matrix(sorted_mat)))
    }
    
    sorted_mat <- sorted_mat[order(sorted_mat[,1]),] #sort based on the first genome
    indices <- seq_len(nrow(block_matrix)) #simple vector from 1 to {n}, to order the blocks
    sorted_mat <- cbind(sorted_mat, indices) #attach the indices, now sorted_mat[,5] represents permutation order for genome1
    
    sorted_mat <- sorted_mat[order(sorted_mat[,2]),] #sort based on start coord on genome2
    sorted_mat <- cbind(sorted_mat, indices) #append indices, now sorted_mat[,6] represents permutation order for genome2
    sorted_mat[,6] <- sorted_mat[,4] * sorted_mat[,6] #multiply by direction to get pos/neg permutations for genome2
    
    sorted_mat <- sorted_mat[order(sorted_mat[,5]),] #sort by permutation order for genome1
    
    permutation_order <- cbind(sorted_mat[,5], sorted_mat[,6]) #now columns 5 and 6 represent permutation order for genomes 1 and 2, resp.
    genome <- permutation_order[,2]
    block_key <- sorted_mat
    colnames(block_key) <- c('start1', 'start2', 'length', 'rel_direction_on_2', 'index1', 'index2')
    # ==== end blocks to matrix ====
    
    unique_1 <- c()
    unique_2 <- c()
    # ==== end finding unique genomic blocks ====
    
    # ==== Simplifying vector of permutations ====
    i <- 1
    indices_to_remove <- c()
    while(i < length(genome)){
      start <- j <- i
      while(i < length(genome) && (genome[i] + 1) == genome[i+1]){
        i <- i+1
        indices_to_remove <- c(indices_to_remove, i)
        block_key[j,3] <- block_key[j,3] + block_key[i,3]
        block_key[j,1] <- min(block_key[j,1], block_key[i,1])
        block_key[j,2] <- min(block_key[j,2], block_key[i,2])
      }
      if (start == i)
        i <- i+1
    }
    
    ctr <- 0
    for(ind in indices_to_remove){
      genome <- genome[-(ind-ctr)]
      ctr <- ctr + 1
    }
    
    new <- seq_len(length(genome))
    temp <- sort(abs(genome))
    
    conv <- as.list(setNames(new, temp))
    
    new_genome <- c()
    for(entry in genome){
      if(entry < 0)
        mult <- -1
      else
        mult <- 1
      
      new_genome <- c(new_genome, conv[[toString(abs(entry))]]*mult)
    }
    genome <- new_genome
    num_blocks <- length(genome)
    if (length(indices_to_remove) > 0){
      block_key <- block_key[-indices_to_remove,]
      if (!is(object = block_key,
              class2 = "array")) {
        return(list('counts'=c(0, 0, 0),
                    'scenario'='none',
                    'block_key'=as.matrix(block_key)))
      }
      # if (!('array' %in% class(block_key))){
      #   return(list('counts'=c(0, 0, 0),
      #               'scenario'='none',
      #               'block_key'=as.matrix(block_key)))
      # }
      block_key[,5] <- seq_len(nrow(block_key))
    }
    block_key <- block_key[,-6]
    # ==== end simplifying vector of permutation
    
    # ===== End synteny to vector of permutations =====
    
    # ===== Breakpoint Graph from Genome =====
    len <- 2 * length(genome)
    
    bpg <- matrix(nrow=len + 2, ncol = 3) #BreakPoint Graph
    
    #adding caps for the genomes
    #len + 1 is the end cap
    #len + 2 is the beginning cap
    
    #every iteration we create a black line between this vertex and the previous vertex in A
    prev_black <- len + 2 
    
    for (i in seq_len(length(genome))){
      val <- genome[i] #get permutation number
      
      if (val < 0){ #if negative, map permutation n to 2n and 2n-1 (in that order)
        v1 <- abs(val) * 2
        v2 <- (abs(val) * 2) - 1
        bpg[v1, 3] <- -1
        bpg[v2, 3] <- -1
        bpg <- add_connection(bpg, v2, v2 - 1, 2) #add gray line between left vertex of block and previous block's right vertex
        bpg <- add_connection(bpg, v1, v1 + 1, 2) #add gray line between right vertex of block and next block's left vertex
      }
      else { #if positive, map permutation n to 2n-1 and 2n (in that order)
        v1 <- (val * 2) - 1
        v2 <- val * 2
        bpg[v1, 3] <- 1
        bpg[v2, 3] <- 1
        bpg <- add_connection(bpg, v1, v1 - 1, 2) #add gray line between left vertex of block and previous block's right vertex
        bpg <- add_connection(bpg, v2, v2 + 1, 2) #add gray line between right vertex of block and next block's left vertex
      }
      
      bpg <- add_connection(bpg, v1, prev_black, 1) #add black line between current vertex and prev. vertex in A
      prev_black <- v2 #store right vertex to connect to the next one
    }
    #black line between last value of A and end cap of A
    bpg <- add_connection(bpg, prev_black, len+1, 1)
    
    #gray line between last value of B and end cap of B
    #and gray line between first value of B and start cap of B
    bpg <- add_connection(bpg, len, len + 1, 2)
    bpg <- add_connection(bpg, 1, len + 2, 2)
    # ===== finished creating breakpoint graph =====
    
    # ===== Rearrangement Scenarios =====
    if (NumRuns < 1)
      NumRuns <- ceiling(sqrt(length(graph_to_genome(bpg))))
    
    if(Verbose)
      progress <- txtProgressBar(max = NumRuns, char = "=", style=3)
    
    if(Mean){
      invert_count <- 0
      block_count <- 0
      seq_skip_inv_probs <- seq(0, 0.75, length.out=NumRuns)
      seq_skip_bi_probs <- seq(0.75, 0, length.out=NumRuns)
      min_total_seen <- 6 * length(graph_to_genome(bpg))
      for(i in seq_len(NumRuns)){
        rearr <- generate_rearrangements(bpg, 
                                         bi_skip_prob=seq_skip_bi_probs[i], 
                                         inv_skip_prob=seq_skip_inv_probs[i])
        invert_count <- invert_count + rearr$invert_count
        block_count <- block_count + rearr$block_count
        ctot <- rearr$invert_count + rearr$block_count
        if(ctot < min_total_seen){
          min_total_seen <- ctot
          rearrangements <- rearr$rearrangements
        }
        if (Verbose)
          setTxtProgressBar(progress, i)
      }
      invert_count <- invert_count / NumRuns
      block_count <- block_count / NumRuns
    }
    else{
      block_count <- invert_count <- 3*length(graph_to_genome(bpg))
      seq_skip_inv_probs <- seq(0, 0.75, length.out=NumRuns)
      seq_skip_bi_probs <- seq(0.75, 0, length.out=NumRuns)
      for (i in seq_len(NumRuns)){
        rearr <- generate_rearrangements(bpg, 
                                         bi_skip_prob=seq_skip_bi_probs[i], 
                                         inv_skip_prob=seq_skip_inv_probs[i])
        if((rearr$block_count + rearr$invert_count) < (block_count+invert_count)){
          block_count <- rearr$block_count
          invert_count <- rearr$invert_count
          rearrangements <- rearr$rearrangements
          graph <- rearr$rearrangements
        }
        if (Verbose)
          setTxtProgressBar(progress, i)
      }
    }
    
    # ===== end rearrangement scenarios =====
    
    # OUTPUT ---------------------------------------------------------------------
    return(list('counts'=c(invert_count, block_count, length(unique_1)/2 + length(unique_2)/2),
                'scenario'=rearrangements,
                'block_key'=block_key))
  }
  
  num_genomes <- nrow(SyntenyObject)
  rearr_mat <- matrix(data=list(), nrow=num_genomes, ncol=num_genomes)
  rearrangements <- list("Translocations"=0, "Gen1Dup"=0, "Gen2Dup"=0, 
                         "Inversions"=0, "Transpositions"=0, "InsDel"=0, 
                         "pct_hits"=0)
  row_names <- rownames(SyntenyObject)

  rownames(rearr_mat) <- row_names
  colnames(rearr_mat) <- row_names
  rates <- vector("list", num_genomes)
  for ( i in seq_len(num_genomes) ){
    for ( j in i:num_genomes ){
      if (j == i){
        rearr_mat[[i, j]] <- SyntenyObject[[i,j]]
        next
      }
        
      gen1 <- max(i, j)
      gen2 <- min(i, j)
      
      p_hits <- (sum(SyntenyObject[[gen2, gen1]][,"width"])*2) / (sum(unlist(SyntenyObject[[gen1,gen1]])) + sum(unlist(SyntenyObject[[gen2, gen2]])))
      rearrangements$pct_hits <- p_hits
      rearr_mat[[gen2, gen1]] <- paste(round(p_hits*100, 2), "% hits", sep='')
      
      gen_info <- SyntenyObject[[gen1, gen2]]
      
      num_chrom <- max(gen_info[,1],gen_info[,2]) 
      if(Verbose)
        cat("\nComputing Scenario for Genomes ", i, " and ", j, "...", sep='')
      for ( chrom in seq_len(num_chrom) ){
        rows <- gen_info[gen_info[,1] == chrom | gen_info[,2] == chrom,]
        # if only one row matches, it returns a vector
        if (is(object = rows,
               class2 = "integer")) {
          rows <- matrix(data = rows,
                         nrow = 1L)
        }
        # if ( class(rows)[1] == "integer" ) #if only one row matches, it returns a vector
        #   rows <- matrix(rows, nrow=1)

        if (nrow(rows) == 0) #if no rows match, continue
        {
          if (Verbose){
            cat("\n\tChromosome ", chrom, " of ", num_chrom, ":\n", sep='')
            progress <- txtProgressBar(max = 1, char = "=", style=3)
            setTxtProgressBar(progress, 1)
          }
          next()
        }
        
        #eliminating rows we've already seen
        rows <- rows[rows[,1] >= chrom,]
        if (is.null(nrow(rows)) || nrow(rows) == 0){ #if no rows match, continue
          if (Verbose){
            cat("\n\tChromosome ", chrom, " of ", num_chrom, ":\n", sep='')
            progress <- txtProgressBar(max = 1, char = "=", style=3)
            setTxtProgressBar(progress, 1)
          }
          next()
        }

        rows <- rows[rows[,2] >= chrom,]
        if (is.null(nrow(rows)) || nrow(rows) == 0){ #if no rows match, continue
          if ( Verbose ) { 
           cat("\n\tChromosome ", chrom, " of ", num_chrom, ":\n", sep='')
            progress <- txtProgressBar(max = 1, char = "=", style=3)
            setTxtProgressBar(progress, 1)
          }
          next()
        }

        for (k in seq_len(nrow(rows))){
          if (rows[k, 1] != rows[k, 2]){
            #determine if translocation or duplication
            gen1_match <- gen_info[gen_info[,5]==rows[k,5] & gen_info[,7]==rows[k,7],]
            gen2_match <- gen_info[gen_info[,5]==rows[k,5] & gen_info[,7]==rows[k,7],]
            found <- FALSE
            if (is(object = gen1_match,
                   class2 = "matrix")) {
              rearrangements$Gen1Dup <- rearrangements$Gen1Dup + 1
              found <- TRUE
            }
            # if(class(gen1_match)[1] == "matrix"){
            #   rearrangements$Gen1Dup <- rearrangements$Gen1Dup + 1
            #   found <- TRUE
            # }
            if (is(object = gen2_match,
                   class2 = "matrix")) {
              rearrangements$Gen2Dup <- rearrangements$Gen2Dup + 1
              found <- TRUE
            }
            # if(class(gen2_match)[1] == "matrix"){
            #   rearrangements$Gen2Dup <- rearrangements$Gen2Dup + 1
            #   found <- TRUE
            # }
            
            if (!found){
              rearrangements$Translocations <- rearrangements$Translocations + 1
            }
              
          }
        }
        matching <- rows[rows[,1] == rows[,2],]
        if (!is(object = matching,
                class2 = "matrix")) {
          matching <- matrix(data = matching,
                             nrow = 1L)
        }
        # if (class(matching)[1] != "matrix"){
        #   matching <- matrix(matching, nrow=1)
        # }
        if (nrow(matching) == 0){
          next()
        }
        if(Verbose)
          cat("\n\tChromosome ", chrom, " of ", num_chrom, ":\n", sep='')  
        chrom_results <- rearrange_chromosome(matching, SyntenyObject, Verbose=Verbose, 
                                              Mean=Mean, chrom=chrom,
                                              NumRuns=NumRuns,
                                              MinBlockLength=MinBlockLength, 
                                              first_gen = gen1, second_gen=gen2, EarlyOutput=TRUE)
        counts <- chrom_results$counts
        rearrangements$Inversions <- rearrangements$Inversions+counts[1]
        rearrangements$Transpositions <- rearrangements$Transpositions+counts[2]
        rearrangements$InsDel <- rearrangements$InsDel+counts[3]
        rearrangements$Scenario <- unlist(chrom_results$scenario)
        rearrangements$Key <- chrom_results$block_key
        
      }
        # These aren't fully fleshed out and will probably be worse to include at this point
        rearrangements$Gen1Dup <- rearrangements$Gen2Dup <- rearrangements$Translocations <- rearrangements$InsDel <- NULL
        if (is.null(rearrangements$Scenario)) rearrangements$Scenario <- 'No events'
        if (is.null(rearrangements$Key)) rearrangements$Key <- NA
      rearr_mat[[gen1, gen2]] <- rearrangements
      rearrangements <- list("Translocations"=0, "Gen1Dup"=0, "Gen2Dup"=0, "Inversions"=0, "Transpositions"=0, "InsDel"=0)
    }
  }
  if(Verbose)
    cat('\nDone.\n\nTime difference of', 
        round(difftime(Sys.time(), starttime, units = 'secs'), 2),
        'seconds.\n')
  class(rearr_mat) <- append('GenRearr', class(rearr_mat))
  return(rearr_mat)
}

print.GenRearr <- function(x, ...){
  to_print <- matrix(character(), nrow=nrow(x)+1, ncol=ncol(x)+1)
  to_print[1,] <- c('', colnames(x))
  nc <- ncol(x)
  rnames <- rownames(x)
  for (i in seq_len(nrow(x))){
    row <- x[i,]
    outvec <- character(nc)
    for ( j in seq_along(outvec) ){
      entry <- row[[j]]
      if (is(entry, 'character')) outvec[j] <- entry
      else if (is(entry, 'integer')) outvec[j] <- paste0(length(entry), ' Chromosomes')
      else {
        outvec[j] <- paste0(round(entry$Inversions, 1), 'I,', round(entry$Transpositions,1), 'T')
      }
    }
    outvec <- c(rnames[i], outvec)
    to_print[i+1,] <- outvec
  }
  lines <- format(to_print, justify='centre')
  for (i in seq_len(nrow(lines))) cat(lines[i,], '\n')
  
  return(invisible(x))
}
