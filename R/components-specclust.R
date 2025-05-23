# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components-clust.R
NULL

#' Components based on MS/MS similarity.
#'
#' This class is derived from \code{\link{componentsClust}} and is used to store components from feature groups that
#' were clustered based on their MS/MS similarities.
#'
#' Objects from this class are generated by \code{\link{generateComponentsSpecClust}}
#'
#' @template components-altered-note
#' 
#' @seealso \code{\link{componentsClust}} for other relevant methods and \code{\link{generateComponents}}
#'
#' @templateVar class componentsSpecClust
#' @template class-hierarchy
#'
#' @export
componentsSpecClust <- setClass("componentsSpecClust", contains = "componentsClust")

#' Generate components based on MS/MS similarity
#'
#' Generates components based on MS/MS similarity between feature groups.
#'
#' @templateVar algo hierarchical clustering of MS/MS spectra
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam specclust
#' @template algo_generator
#'
#' @details The similarities are converted to a distance matrix and used as input for hierarchical clustering, and the
#'   resulting dendrogram is automatically cut with \code{\link{cutreeDynamicTree}}. The clustering is performed with
#'   \code{\link[fastcluster:hclust]{fastcluster::hclust}}.
#'
#' @param MSPeakLists The \code{\link{MSPeakLists}} object for the given feature groups that should be used for MS
#'   spectral similarity calculations.
#'
#' @templateVar noDots TRUE
#' @template compon_algo-args
#' @template compon_gen-clust
#' @template dynamictreecut
#' @template specSimParams-arg
#'
#' @inheritParams generateComponents
#'
#' @return The components are stored in objects derived from \code{\link{componentsSpecClust}}.
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow} the spectral similarities for each set are
#'   combined as is described for the \code{\link[=spectrumSimilarity,MSPeakListsSet-method]{spectrumSimilarity}} method
#'   for sets workflows.
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}> and Bas van de Velde (major contributions to spectral binning and
#'   similarity calculation).
#'
#' @templateVar what generateComponentsSpecClust
#' @template main-rd-method
#' @export
setMethod("generateComponentsSpecClust", "featureGroups", function(fGroups, MSPeakLists, method = "complete",
                                                                   specSimParams = getDefSpecSimParams(),
                                                                   maxTreeHeight = 1, deepSplit = TRUE,
                                                                   minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(componentsSpecClust(distm = NULL, method = method, gInfo = groupInfo(fGroups),
                                   properties = list(), maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                                   minModuleSize = minModuleSize, algorithm = "specclust"))

    gNames <- names(fGroups)
    
    cat("Calculating distance matrix... ")
    sims <- spectrumSimilarity(MSPeakLists, gNames, NULL, MSLevel = 2, specSimParams = specSimParams, NAToZero = TRUE,
                               drop = FALSE)
    
    # figure out fGroups with results: these must have non-zero columns (or rows), since there must be at least a 1.0
    # similarity with itself.
    grpsResults <- gNames[colSums(sims) > 0]
    sims <- sims[grpsResults, grpsResults]
    
    distm <- 1 - as.dist(sims)
    cat("Done!\n")
    
    gInfo <- groupInfo(fGroups)[grpsResults, ]

    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo,
                               properties = list(specSimParams = specSimParams),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})
