## [ Edited Feature Selection Function from MethylCibersort 0.2.1 - YG, DW (2018) ] ////////////////////

## Changes made to function:
## added argument & functionality for internal conversion for M-values
## added arguments & functionality for exporting of various intermediate objects
## Fixed the way DMPs are filtered to prevent the factor order / naming order from influencing the CpG result

## Base working of the function is retained and the function takes and outputs the same objects, as before.

## FUNCTION START //////////////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////////////////////////////

feature.select.new <- function(MaxDMPs = 100,                       ## maximum differentially methylated probes to use, takes n/2 from top & n/2 from bottom
                               deltaBeta = 0.2,                     ## cutoff for the minimum difference between pairwise groups by delta-beta
                               useM = FALSE,                        ## option to use M-values not beta-values, largely not useful
                               CellLines.matrix = NULL,             ## input matrix for cell line 'cancer' data
                               export = TRUE,                       ## save a table of signature results
                               export.fit = TRUE,                   ## export the limma result
                               export.cpg = TRUE,                   ## export the CpGs selected
                               sigName = "methylCibersort",         ## name appended to start of filename
                               Stroma.matrix = NULL,                ## matrix of betas for populations
                               Phenotype.stroma = NULL,             ## pheno that corresponds to Stroma.matrix
                               FDR = 0.01,                          ## FDR cutoff
                               silent = TRUE){                      ## run function without returning output 
  
  ### import packages
  require(magrittr)
  require(matrixStats)
  require(BiocGenerics)
  require(MethylCIBERSORT)
  
  message(paste0("methylCibersort ", packageVersion("MethylCIBERSORT")))
  message(paste0("deltaBeta: ", deltaBeta, " MaxDMPs: ", MaxDMPs))
  
  ## process pheno information
  if (!is.null(ncol(CellLines.matrix))) {
    Pheno1 <- c(rep("Cancer", ncol(CellLines.matrix)))
    Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma))
  } else { 
    Pheno2 <- as.character(Phenotype.stroma)
  } # end if else 
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  if (useM){ 
    Mat3 <- minfi::logit2(Mat2) 
  } else {
    Mat3 <- Mat2
  } # end if
  
  message("Setting up for pairwise feature selection")
  ContrastMatrix <- design.pairs(levels(factor(Pheno2)))
  
  ## set up limma comparisons
  Des <- model.matrix(~0 + Pheno2)
  colnames(Des) <- rownames(ContrastMatrix)
  
  ## do limma comparisons
  Fit <- limma::lmFit(Mat3, Des) %>% 
    limma::contrasts.fit(ContrastMatrix) %>% 
    limma::eBayes()
  
  FitList <- list()
  
  ## get top results for each comparison
  for (i in 1:ncol(ContrastMatrix)) {
    tmp.name <- gsub("\\-", "\\.", colnames(ContrastMatrix)[i])
    cat(tmp.name)
    FitList[[i]] <- limma::topTable(Fit, coef = i, number = nrow(Mat2)) %>% 
      dplyr::mutate(ID = rownames(.)) %>% 
      dplyr::filter(adj.P.Val < FDR)
    cat(" ... done", "\n")
    names(FitList)[i] <- tmp.name
  } # end i
  
  if (all(c(export,export.fit))){
    message("Saving limma fit results")
  fN <- paste(sigName, deltaBeta, MaxDMPs, "Fit_list.rds", sep = "_")
  saveRDS(FitList, fN)
  } # end if
  
  message("Calculating population medians")
  ## slow
  #Split <- apply(Mat2, 1, function(x){tapply(x, as.factor(Pheno2), median)}) 
  
  Transformed <- data.frame(t(Mat2))
  
  ## get means by population
  ## faster
  Split <- split(Transformed, Pheno2) 
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split) 
  rownames(Split) <- rownames(Mat2)
  
  dbList <- list()
  message("Getting Delta Beta estimates")
  
  ## get DB estimates
  for (i in 1:ncol(ContrastMatrix)) {
    tmp.name <- gsub("\\-", "\\.", colnames(ContrastMatrix)[i])
    cat(tmp.name)
    dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
    dB <- data.frame(dB = dB, ID = rownames(Split))
    dbList[[i]] <- dB
    cat(" ... done", "\n")
    names(dbList)[i] <- tmp.name
  } # end i
  
  message("Filtering by Delta Beta")
  ## filter by DB threshold
  dbList <- lapply(dbList, function(x) dplyr::filter(x, abs(dB) > deltaBeta))
  
  message("Filtering by max DMP number")
  for (i in 1:length(FitList)) {
    A1 <- FitList[[i]]
    A1 <- dplyr::filter(A1, ID %in% dbList[[i]]$ID)
    A1 <- A1 %>% .[rev(order(.$t)), ]
    if (nrow(A1) > MaxDMPs) {
      A2 <- head(A1[, ], MaxDMPs/2) ## ///////////// Edited to correct issues with ordering of factors
      A3 <- tail(A1[, ], MaxDMPs/2) ## ///////////// affecting resulting CpG selection & erroneous 
      A1 <- rbind(A2, A3) ## /////////////////////// correlation effects, might need to be better adressed in future
    }
    FitList[[i]] <- A1
  } # end i
  
  message("Bulding unique CpG list and annotating probes")
  Nonzeros <- lapply(FitList, function(x) dplyr::select(x, ID))
  
  Nonzeros <- do.call(rbind, Nonzeros)
  #Nonzeros <- dplyr::filter(Nonzeros, !duplicated(ID))
  if (all(c(export,export.cpg))){
    message("Saving complete list of CpGs before filtering")
    fN <- paste(sigName, deltaBeta, MaxDMPs, "all_Nonzeros.txt", sep = "_")
    write.table(data.frame(ID = Nonzeros,
                           PopID = rownames(Nonzeros)), 
                file = fN, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  } # end if
  
  Nonzeros <- data.frame(ID = Nonzeros[-which(duplicated(Nonzeros$ID)), ], 
                         row.names = rownames(Nonzeros)[-which(duplicated(Nonzeros$ID))])
  
  Nonzeros <- data.frame(ID = Nonzeros$ID, PopID = gsub("\\.[0-9]+", "", rownames(Nonzeros)))
  Nonzeros.anno <- minfi::getAnnotation(minfiData::MsetEx, lociNames = Nonzeros$ID)
  Nonzeros.anno <- cbind(Nonzeros, Nonzeros.anno[as.character(Nonzeros$ID), ])
  
  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros.anno$ID, ]
  nrow(Mat3)
  Mat3 <- 100 * Mat3
  
  if (export) {
    
    message("Writing text files")
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), 
                                  stringsAsFactors = F), Collapsed)
    
    fN <- paste(sigName, deltaBeta, MaxDMPs, "Signature.txt", sep = "_")
    write.table(Collapsed, 
                file = fN, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    
    fN <- paste(sigName, deltaBeta, MaxDMPs, "CpG_Annotation.txt", sep = "_")
    write.table(Nonzeros.anno, 
                file = fN, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  } # end export
  
  if(!silent) return(list(SignatureMatrix = Mat3))
} # end function

## FUNCTION END ////////////////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////////////////////////////