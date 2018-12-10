rm(list = ls())
library(dplyr)
library(ccube)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

options(stringsAsFactors = F)

registerDoParallel(cores=3)
set.seed(1234)


numOfClusterPool = 1:7
numOfRepeat = 1


icgc <- "~/debug_pcawg_samples_sv/svclone_paper_20181127/"
dataFolders <- dir(icgc, full.names = T)
sampleNames <- dir(icgc)


bugSamples <- c()
bugdataFolders <- c()
for (ii in seq_along(sampleNames)) {

  sampleName <- sampleNames[ii]
  dataFolder <- dataFolders[ii]

  if (! "ccube_out/ccube_sv_results.RData" %in% dir(dataFolder, recursive = T) ) {
    bugSamples <- c(bugSamples, sampleName)
    bugdataFolders <- c(bugdataFolders, dataFolder)
  }
}


allDebugFolder <- "~/debug_pcawg_samples_sv/debuged_samples"
problemSamples <- c("0b6cd7df-6970-4d60-b7b5-85002a7d8781",
                    "f601cf2f-081f-484d-ab0e-21a8ec8d3770",
                    "fab0be4b-c84a-45cd-a76f-44d9a8bf1846")

  ii = which(bugSamples == "fab0be4b-c84a-45cd-a76f-44d9a8bf1846")
  cat(ii, "\n")

  sampleName <- bugSamples[ii]
  dataFolder <- bugdataFolders[ii]
  resultsFolder <- paste0(allDebugFolder,"/", sampleName)

  if (! dir.exists(resultsFolder) ) {
    dir.create(resultsFolder)
  }

  inputFn <- paste0(dataFolder, "/", sampleName, "_ccube_sv_input.txt")

  mydata <- try( read.delim(inputFn), T)

  if ( is.data.frame(mydata) ) {

    if (nrow(mydata) == 1) {
      singleEventSamples <- c(singleEventSamples, sampleName)
      cat("single event", file = paste0(resultsFolder, "/single_event_sample")  )
    }


    python_false_true_converter <- function(x) {
      if (x == "False") {
        x = FALSE
      }

      if (x == "True") {
        x = TRUE
      }
      x
    }


    mydata <- dplyr::mutate(rowwise(mydata),
                            subclonal_cn1 = python_false_true_converter(subclonal_cn1),
                            subclonal_cn2 = python_false_true_converter(subclonal_cn2))

    #res <- CcubeSVCore(mydata, init = numOfClusterPool, fit_mult = T, fit_hyper = T, use = "use_base", verbose = T)

    ccubeRes <-
      RunCcubePipeline(ssm = mydata, modelSV = T,
                                   numOfClusterPool = numOfClusterPool,
                                   numOfRepeat = numOfRepeat, multiCore =T,
                                   runAnalysis = T, runQC = T, returnAll = T, maxiter = 20)

    if (is.list (ccubeRes) ) {
      fn = paste0(resultsFolder, "/ccube_sv_results.RDdata")
      save(ccubeRes, file = fn)
      fn = paste0(resultsFolder, "/ccube_sv_results.pdf")
      MakeCcubeStdPlot_sv(res = ccubeRes$res, ssm = ccubeRes$ssm, printPlot = T, fn = fn)
    } else {
      cat(ccubeRes, file = paste0(resultsFolder, "/bug_info_ccube.txt"))

      if (nrow(mydata) > 1) {
        problemSamples <- c(problemSamples, sampleName)
        cat("problematic sample", file = paste0(resultsFolder, "/problematic_sample")  )
      }

    }

  } else {
    cat(mydata, file = paste0(resultsFolder, "/bug_info_mydata.txt"))

  }


