#' Run Ccube analysis
#' @param sampleName sample name
#' @param dataFolder path to data folder that stores tmp result files
#' @param resultFolder path to output formats files
#' @param makeFolders flag to create data and results folders, default false
#' @param runParser run parser
#' @param variantCaller name of variant caller
#' @param cnaCaller name of CNA caller
#' @param runAnalysis run analysis
#' @param runQC run QC
#' @param runPcawgMergeQc run QC for PCAWG samples
#' @param runAnalysisSnap run a lite version of the main analysis
#' @param writeOutput write output into files
#' @param allFormats all PCAWG output formats
#' @param basicFormats basic PCAWG output formats
#' @param vcfFile path to vcf file
#' @param copyNumberFile path to copy number file
#' @param purity estimated purity of the sample
#' @param numOfClusterPool candidates of number of clusters
#' @param numOfRepeat number of repeat runs for each candidate number of clusters
#' @param epi sequencing error
#' @param tol convergence thershold
#' @param maxiter maximum iteration
#' @param multiCore use multiple cores for repeated runs
#' @param ccubeInputRDataFile path to Ccube input Rdata
#' @param ccubeResultRDataFile path to Ccube result Rdata
#' @param modelSV using double break points model
#' @param ssm Ccube input data
#' @param maxSnv maximum number of SNVs. Used in runParser mode
#' @param use option for rough estimator for ccf
#' @return returns a list including the prefered solution, res; annotated ssm, ssm; a list of all solutions, results; trace of ELBOs, lb; removed events, droppedSsm
#' @export
RunCcubePipeline <- function(sampleName = NULL, dataFolder = NULL, resultFolder = NULL, makeFolders = F,
                             runParser = F, variantCaller, cnaCaller = NULL,
                             runAnalysis = F, runQC = F, runAnalysisSnap = F, runPcawgMergeQc =F,
                             writeOutput = F, modelSV = F,
                             allFormats = F, basicFormats = T,
                             vcfFile = NULL, copyNumberFile = NULL, purity = NA,
                             numOfClusterPool = NULL, numOfRepeat = NULL,
                             epi = 1e-3, tol = 1e-8, maxiter = 1e3,
                             multiCore = F, ccubeInputRDataFile = NULL,
                             ccubeResultRDataFile = NULL,
                             ssm = NULL,
                             maxSnv = 1e7, use = "use_base"){

  # stopifnot( runParser | runAnalysis | runAnalysisSnap,
  #            runParser & !is.null(variantCaller) & !is.null(cnaCaller),
  #            ! runParser & !is.null(ccubeInputRDataFile)
  #            )


  if (is.null(sampleName) ) {
    sampleName <- "sample1"
  }

  if (is.null(dataFolder) ) {
    dataFolder <- "ccube_data/sample1/"
  }

  if (is.null(resultFolder) ) {
    resultFolder <- "ccube_res/sample1/"
  }


  if (makeFolders) {
    if (!dir.exists(dataFolder)) {
      dir.create(dataFolder, recursive = T)
    }

    if (!dir.exists(resultFolder)) {
      dir.create(resultFolder, recursive = T)
    }
  }



  if (runParser) {

    shellCommandConsensus <- paste(
      "python create_ccfclust_inputs_consensus_bb.py -v ", variantCaller,
      " --output-variants ", paste0(dataFolder, "/ssm_data.txt"),
      " ", vcfFile, sep = ""
    )

    cat(shellCommandConsensus, "\n")
    system(shellCommandConsensus, intern = TRUE)

    ssm <- read.delim(paste0(dataFolder, "/ssm_data.txt"),
                      stringsAsFactors = F)

    cna <- read.delim(batternbergFile, stringsAsFactors = F)

    if (cnaCaller == "pcawg11_std") {
      ssm = ParseSnvCnaPcawg11Format(ssm, cna)
    } else {
      ssm = ParseSnvCnaConsensus(ssm, cna)
    }


    ssm <- dplyr::filter(ssm, major_cn > 0)

    sampleSummary <- data.frame(samplename = sampleName,
                                overlap_cna = nrow(ssm),
                                overlap_cna_balance_1_1 = nrow(filter(ssm, major_cn == 1 & minor_cn ==1) ),
                                overlap_cna_balance_2_2 = nrow(filter(ssm, major_cn == 2 & minor_cn ==2) ),
                                overlap_cna_balance = nrow(filter(ssm, major_cn == minor_cn  ) ),
                                overlap_clonal_cna = nrow(filter(ssm, cn_frac == 1  ) ),
                                overlap_subclonal_cna = nrow(filter(ssm, cn_frac != 1  ) )
    )

    write.table(sampleSummary, file = paste0(dataFolder, '/', sampleName, '_pre_clustering_summary.txt'),
                sep = '\t', row.names = F, quote = F)

    if (nrow(ssm)>0) {
      if (nrow(ssm) > maxSnv) {
        ssm <- dplyr::sample_n(ssm, maxSnv)
      }
      ssm$normal_cn = 2
      ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d)
      ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts, mutation_id = gene)
      ssm$ccube_purity <- GetPurity(ssm)
      ssm$purity <- purity

      write.table(unique(ssm$ccube_purity), file = paste0(dataFolder,"/ccube_purity.txt"),
                  row.names = F, col.names=F, quote =F )
      save(ssm, file = paste0(dataFolder, "/ssm_no_chrxy.RData"))
    }
    save(ssm, file = ccubeInputRDataFile)
  }

  if (runAnalysis | runAnalysisSnap) {

    if (is.null(ssm)) {
      if (!is.null(ccubeInputRDataFile) ) {
        if (file.exists(ccubeInputRDataFile)) {
          load(ccubeInputRDataFile)
        }
      }
    }

    if (modelSV) {
      ssm <-CheckAndPrepareCcubeInupts_sv(ssm)
      ssm <- GetCcf_sv(ssm,  use=use)
    } else {
      ssm <- CheckAndPrepareCcubeInupts(ssm)
      ssm <- GetCcf(ssm,  use=use)
    }

    if (nrow(ssm) == 1) {

      message(sprintf("Only one variant in the sample! \n Exit without clustering \n Return a rough estimate. "))

      return(ssm)
    }


    # filtering
    if ( "major_cn" %in% colnames(ssm)  ) {
      droppedSsm <- dplyr::filter(ssm, major_cn <= 0)
      ssm <- dplyr::filter(ssm, major_cn > 0)
    }

    if ( "major_cn1" %in% colnames(ssm)  ) {
      droppedSsm <- dplyr::filter(ssm, major_cn1 <= 0)
      ssm <- dplyr::filter(ssm, major_cn1 > 0)
    }

    if ( "major_cn2" %in% colnames(ssm)  ) {
      droppedSsm <- dplyr::filter(ssm, major_cn2 <= 0)
      ssm <- dplyr::filter(ssm, major_cn2 > 0)
    }


    if (runAnalysisSnap) {
      iterSetting <- max(numOfClusterPool)
    } else {
      iterSetting <- sort(rep(numOfClusterPool, numOfRepeat))
    }


    if (modelSV) {
      func <- get("CcubeSVCore")
    } else {
      func <- get("CcubeCore")
    }



    if (multiCore) {
      results <- foreach::foreach(n = seq_along(iterSetting), .combine = c, .packages = "ccube") %dopar%
      {
        k <- iterSetting[n]
        list(func(mydata = ssm, epi=epi,
                  init=k, tol = tol, maxiter = maxiter,
                  fit_mult = T, fit_hyper = T, use = "use_base", verbose = F))
      }

    }else {
      results <- foreach::foreach(n = seq_along(iterSetting), .combine = c, .packages = "ccube") %do%
      {
        k <- iterSetting[n]
        list(func(mydata = ssm, epi=epi,
                  init=k, tol = tol, maxiter = maxiter,
                  fit_mult = T, fit_hyper = T, use = "use_base", verbose = F))
      }
    }

    lb <- unlist(Map( function(x) max(x$L), results))
    sortedIdx <- sort(lb, decreasing = T, index.return = T)$ix
    res <- results[[sortedIdx[1]]]

  }

  if (runQC) {

    if (!is.null(ccubeResultRDataFile) ) {
      if (file.exists(ccubeResultRDataFile)) {
        load(ccubeResultRDataFile)
        sortedIdx <- sort(lb, decreasing = T, index.return = T)$ix
      }
    }

    # Check for clonal cluster
    foundDiffRes <- F
    if (! HasClonalCluster(res) ) {
      for (ii in sortedIdx) {
        rr <- results[[ii]]
        if (HasClonalCluster(rr)) {
          foundDiffRes <- T
          break
        }
        passClonalCluterTest <- F
      }
    } else {
      passClonalCluterTest <- T
    }

    if (foundDiffRes) {
      res <- rr
      passClonalCluterTest <- T
    }


    if (modelSV) {
      func_qc1 <- get("CullEmptyClusters_sv")
      func_qc2 <- get("CullSmallClusters_sv")
      func_qc3 <- get("MergeClusters_sv")
    } else {
      func_qc1 <- get("CullEmptyClusters")
      func_qc2 <- get("CullSmallClusters")
      func_qc3 <- get("MergeClusters")
    }



    # remove empty cluster
    res <- func_qc1(res = res, ssm = ssm, epi = epi)
    passEmptyCluterTest <- T
    # remove small cluster
    res <- func_qc2(res = res, ssm = ssm, th = 1e-2, epi = epi)
    passSmallClusterTest <- T

    # Merge clusters
    res <- func_qc3(res = res, ssm = ssm, epi = epi)
    passMergeClusterTest <- T

    res$qc <- list(passClonalCluterTest = passClonalCluterTest,
                   passEmptyCluterTest = passEmptyCluterTest,
                   passSmallClusterTest = passSmallClusterTest,
                   passMergeClusterTest = passMergeClusterTest)
  }

  if (runPcawgMergeQc) {

    if (!is.null(ccubeResultRDataFile) ) {
      if (file.exists(ccubeResultRDataFile)) {
        load(ccubeResultRDataFile)
      }
    }

    if (modelSV) {
      func_qc <- get("MergeClusters_sv")
    } else {
      func_qc <- get("MergeClusters")
    }

    res = func_qc(res = res, ssm = ssm)

  }

  # Annote results
  if ( (runAnalysis | runAnalysisSnap | runQC | runPcawgMergeQc) ) {


    if (modelSV) {
      func_annote <- get("AnnotateCcubeResults_sv")
    } else {
      func_annote <- get("AnnotateCcubeResults")
    }

    ssm <- func_annote(ssm = ssm, res = res)

  }

  if ( writeOutput ) {
    # write Ccube results Rdata
    if (! is.null(ccubeResultRDataFile) ) {
      fn <-  ccubeResultRDataFile
    } else {
      fn <- paste0(dataFolder,
                   "/ccube_res.RData")
    }

    ccubeRes <- list(res = res, results = results, ssm = ssm, lb = lb)
    save(ccubeRes, file = fn)


    if (modelSV) {
      func_write <- get("WritePcawgFormats_sv")
      func_plot <- get("MakeCcubeStdPlot_sv")
    } else {
      func_write <- get("WritePcawgFormats")
      func_plot <- get("MakeCcubeStdPlot")
    }


    if (basicFormats | allFormats) {
      func_write(ssm = ssm, res = res, resultFolder = resultFolder,
                        sampleName = sampleName, allFormats = allFormats,
                        basicFormats = basicFormats)
    }

    # summary graph
    fn <- paste0(resultFolder, "/", sampleName, "_results_summary.pdf")
    func_plot(ssm = ssm, res = res, printPlot = T, fn = fn)

  }

  return(list(res = res, results = results, ssm = ssm, lb = lb, droppedSsm = droppedSsm))
}


#' A pipeline to post assign events.
#' @param snvRes A reference Ccube SNV result list
#' @param svRes A reference Ccube SV results list.
#' @param mydata A data frame of variants (SNV or SV) to be assigned. Ideally, the data has been processed by CcubeSV or Ccube model. So it should have ccube_mult1/ccube_mult2 or ccube_mult columns.
#' @return A list containing, res, the post assigned result list and, ssm, the annotated data
#' @details At least one of snvRes and svRes has to be provided. If the both lists are provided, they will be combined.
#' @export
RunPostAssignPipeline <- function(snvRes = NULL, svRes = NULL, mydata, verbose = T) {


  stopifnot( !is.null(snvRes) | !is.null(svRes)  )

  if (!is.null(snvRes) & !is.null(snvRes) ) {
    referenceRes <- CombineSNVandSVResults(snvRes, svRes)
  } else if (!is.null(snvRes) & is.null(svRes)) {
    referenceRes <- snvRes
  } else {
    referenceRes <- svRes
  }

  if ( all( c("var_counts1","ref_counts1", "var_counts2","ref_counts2") %in% names(mydata) ) &
       ! all( c("var_counts","ref_counts") %in% names(mydata) ) ) {
    modelSV = T
  } else {
    modelSV = F
  }

  if (modelSV) {
    func_assign <- get("AssignWithCcube_sv")
  } else {
    func_assign <- get("AssignWithCcube")
  }

  postAssignRes <- func_assign(referenceRes, mydata, verbose = verbose)

  if (nrow(mydata) > 1) {

    if (modelSV) {
      func_qc1 <- get("CullEmptyClusters_sv")
      func_qc2 <- get("CullSmallClusters_sv")
      func_qc3 <- get("MergeClusters_sv")
    } else {
      func_qc1 <- get("CullEmptyClusters")
      func_qc2 <- get("CullSmallClusters")
      func_qc3 <- get("MergeClusters")
    }

    postAssignRes <- func_qc1(res = postAssignRes, ssm = mydata)
    postAssignRes <- func_qc2(res = postAssignRes, ssm = mydata, th = 1e-2)
    postAssignRes <- func_qc3(res = postAssignRes, ssm = mydata)

  }

  if (modelSV) {
    func_annotate <- get("AnnotateCcubeResults_sv")
  } else {
    func_annotate <- get("AnnotateCcubeResults")
  }
  annotatedSsm <- func_annotate(mydata, postAssignRes)

  return(list(res = postAssignRes, ssm = annotatedSsm))
}
