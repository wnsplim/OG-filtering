function(og, threshold, n = NULL, outgroup = NULL, statistics = FALSE, quiet = FALSE) {
  ## Add "counts" column which represents a taxon coverage of each OG
  og <- og %>% mutate(counts = rowSums(. != 0)) %>%
    { 
      rownames(.) <- paste0("OG", seq_along(.$counts))
      .
    }
  
  ## If n is null, set n to min(og$counts)
  if (is.null(n)) {
    n <- min(og$counts)
  } else {
    n = n
  }
  
  ## Following loop performs weighting and filtering
  for (i in seq(max(og$counts), max(n,min(og$counts)), by = -1)) {
    ## Keep all OG with maximum coverage or if it is missing in only one taxon
    if (i == ncol(og)-1 | i == ncol(og)-2) {
      assign(paste("og", i, sep = "_"), og[og$counts == i, ] %>% select(-counts, -one_of(outgroup)))
    } else if (i == max(og$counts)){
      assign(paste("og", i, sep = "_"), og[og$counts == i, ] %>% select(-counts, -one_of(outgroup)))
    } else {
      ## Some matrices to make weighting calculation easier
      og_temp <- og[og$counts == i, ] %>% select(-counts, -one_of(outgroup)) 
      og_temp_binary <- og_temp %>% mutate(across(everything(), ~ ifelse(. != 0, 1, .)))
      
      ## For the all OGs stored in the matrix og_filtered (taxon coverage higher than i), 
      # obtain a frequency of missing data for each column (taxon) and multiply it by 100
      colsum <- colSums(
        og_filtered %>% mutate(across(everything(), ~ ifelse(. != 0, 1, .))))
      
      weights <- 100*(max(colsum)-colsum)/max(colsum)
      
      ## Above weights will be higher for the taxa with higher frequency of missing data 
      ## For each OG with a taxon coverage i, set presence as 1 and absence as -1
      ## Then, perform an elementwise multiplication between weights and presence/absence vector for each OG
      ## Thus, the final score (weightsums) of OG will be higher if the taxa that were relatively underrepresented are present
      # and penalized if they are absent
      weightsums <- as.matrix(og_temp_binary) %*% weights - as.matrix(1-og_temp_binary) %*% weights
      
      ## Filter OGs based on the predefined threshold. Threshold can be subjectively decided by the user, possibly with the aid of the threshold.R function
      assign(paste("og", i, sep = "_"), og_temp %>% cbind(weightsums) %>%
               filter(weightsums >= threshold*(ncol(og_temp)-i)) %>% select(!weightsums))
    }
    
    ## Concatenate filtered OGs to og_filtered matrix
    if (i == ncol(og) - 1) {
      og_filtered <- get(paste("og", i, sep = "_"))
    } else {
      og_filtered <- bind_rows(og_filtered, get(paste("og", i, sep = "_")))
    }
    
    ## For each OG groups with taxon coverage i, count how many OGs made it to the final matrix
    if (quiet == FALSE) {
      cat(paste("Number of OG of size", i, ":", nrow(get(paste("og", i, sep = "_"))), "\n"))
    }
  }
  
  ## Count total number of OGs per taxon and total number of OGs retained in the final matrix
  og_binary <- og_filtered %>% mutate(across(everything(), ~ ifelse(. != 0, 1, .)))
  
  taxonog <- numeric()
  taxoncov <- numeric()
  
  for (sp in colnames(og_filtered)) {
    taxonog[sp] <- sum(og_binary[[sp]])
    taxoncov[sp] <- sum(og_binary[[sp]])/nrow(og_binary)
    if (quiet == FALSE) {
      cat(paste(sp, ":", taxonog[sp], "|",round(taxoncov[sp], digits=3), "\n"))
    }
  }

  if (quiet == FALSE) {
    cat(paste("Total OGs:", nrow(og_binary), "\n"))
  }
  
  if (statistics == TRUE) {
    ## Calculate various statistics which can be used to determine a threshold value
    mean <- mean(taxoncov)
    sd <- sd(taxoncov)
    cv <- sd(taxoncov)/mean(taxoncov)
    mincov <- min(taxoncov)
    return(list(og_filtered, mean, sd, cv, mincov))
  }
   if (statistics == FALSE) {
     return(list(og_filtered)) 
   }
}
