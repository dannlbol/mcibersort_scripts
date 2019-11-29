## generate random proportions for  cell types
generate.random.prop <- function(nr = 100, ## number of rows
                                 nc = 12) ## number of columns
{
  ## create empty matrix
  prop.mat <- matrix(data = 0, nrow = nr, ncol = nc)
  dimnames(prop.mat) <- list(1:nrow(prop.mat), NULL)
  ## generate a set sequence of proportions for the population being generated
  ctoi <- seq(0, 1, 1/nr)
  ## fill in the rest with values from uniform distribution 
  for (i in 1:nrow(prop.mat)){
    tmp.set <- ctoi[i] ## set your static value
    tmp.max <- 1-tmp.set ## maximum possible value given 1-static value
    tmp.val <- rep(0, nc-1) ## populate the slots with 0
    for (j in sample(1:(nc-2))){
      tmp.val[j] <- runif(1, min = 0, max = tmp.max-sum(tmp.val)) ## generate random values
    } # end j
    tmp.val[nc-1] <- tmp.max-sum(tmp.val) ## generate the final value for the row to add up to 1
    prop.mat[i, ] <- c(sample(tmp.val, nc-1, replace = FALSE), tmp.set) ## sample without replacement then append the set value
  } # end i
  return(prop.mat)
} # end function

## generate synthetic proportions for use in benchmarking
make.synth.prop <- function(synth.colnames = NULL,      ## factor names, 1 per population
                            pop.rows = 100)             ## number of rows in matrix to populate per population
{ 
  synth.list <- list()
  ## run through the list of factors, for each factor generate mixes where 
  ## for n pop.rows the proportion is a set sequence and the rest are random
  for (i in 1:length(synth.colnames)){
    idx <- 1:(length(synth.colnames)-1)
    idx <- R.utils::insert(idx, ats = i, length(synth.colnames))
    synth.list[[synth.colnames[i]]] <- generate.random.prop(nr = pop.rows, 
                                                            nc = length(synth.colnames))[, idx]
    colnames(synth.list[[synth.colnames[i]]]) <- synth.colnames
  } # end i
  synth.prop.mat <- do.call("rbind", synth.list)
  rownames(synth.prop.mat) <- paste0("synth.", 1:nrow(synth.prop.mat))
  return(synth.prop.mat)
} # end function

make.synth.mix <- function(input.data = NULL,                                  ## data matrix to be converted to synth mixes
                           pop.factor = NULL,                                  ## factor with levels of the populations of data matrix
                           pop.rows = 100,                                     ## how many rows to generate for each population
                           output.dir = getwd(),                               ## location to save proportions table and resulting mixtures
                           output.name = gsub("-", "_", Sys.Date()),           ## name to append to the filenames
                           n.cores = 1){                                       ## do in parallel for n cores
  synth.prop.mat <- make.synth.prop(levels(pop.factor), pop.rows)
  write.table(synth.prop.mat, 
              paste0(output.dir, "/", output.name, "_synth_prop.txt"), 
              quote = FALSE, 
              row.names = FALSE)
  soi.mean <- t(apply(input.data, 1, function(x){tapply(x, pop.factor, mean)}))
  require(parallel)
  synth.list <- mclapply(1:nrow(synth.prop.mat), 
                         function(x) {apply(soi.mean[, colnames(synth.prop.mat)], 
                                            1, 
                                            function(y) weighted.mean(y, 
                                                                      w = synth.prop.mat[x, ]))},
                         mc.cores = n.cores)
  synth.df <- as.data.frame(do.call(cbind, synth.list))
  synth.df <- cbind(rownames(soi.mean), synth.df)
  colnames(synth.df) <- c("NAMES", rownames(synth.prop.mat))
  write.table(synth.df, 
              paste0(output.dir, "/", output.name, "_synth_mix.txt"), 
              quote = FALSE, 
              row.names = TRUE)
} # end function 
