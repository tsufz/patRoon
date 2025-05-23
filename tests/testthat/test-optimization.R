# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("optimization")

initXCMS()
library(Biobase) # BUG: needed by XCMS with futures...

# disable parallelization on covr as it seems to slow things down a lot
doPar <- !requireNamespace("covr", quietly = TRUE) || !covr::in_covr()

anaInfo <- getTestAnaInfo()[1:2, ]
anaInfoOne <- getTestAnaInfo()[4, ]
epAnaInfo <- makeMZXMLs(anaInfoOne)

# HACK/UNDONE: mzR doesn't seem to be able to load mzXML files generated by
# OpenMS (as done with makeMZXMLs). Therefore, simply copy mzML
# alongside with it for now...
file.copy(patRoon:::getMzMLAnalysisPath(anaInfoOne$analysis[1], anaInfoOne$path[1]), epAnaInfo$path[1])

# BUG BUG BUG: lm seems to give slightly different results across systems, even with the same Docker image under
# different hosts!! For now just compare ref experimental data, so no results, show texts and plots... Furthermore,
# group optimizing happens with the features object from the first iteration, as these inconsistencies make the others
# are slightly different.

ffOptOpenMS <- optimizeFeatureFinding(anaInfo, "openms", list(chromFWHM = c(5, 10), mzPPM = c(5, 15),
                                                              noiseThrInt = 3E4),
                                      maxIterations = 2, parallel = doPar)
# disable 'old' xcms for now to save testing time (both interfaces are fairly similar anyway)
# ffOptXCMS <- optimizeFeatureFinding(anaInfo, "xcms", list(mzdiff = c(0.002, 0.006)))

ffOptXCMS3 <- optimizeFeatureFinding(anaInfo, "xcms3", list(mzdiff = c(0.002, 0.006), noise = 3E4), maxIterations = 2,
                                     parallel = doPar)
ffOptEnviPick <- optimizeFeatureFinding(epAnaInfo, "envipick", list(dmzgap = c(10, 20), minpeak = 25),
                                        maxIterations = 2, parallel = doPar)
ffOptKPIC2 <- suppressWarnings(optimizeFeatureFinding(anaInfo, "kpic2", list(mztol = c(0.002, 0.01), level = 2E5, kmeans = TRUE),
                                                      maxIterations = 2, parallel = doPar))

suppressWarnings(ffOptEmpty <- optimizeFeatureFinding(anaInfo, "openms", list(chromFWHM = c(5, 10), noiseThrInt = 1E9)))

ffOpenMS <- findFeatures(anaInfo, "openms", noiseThrInt = 3E4)
ffXCMS3 <- findFeatures(anaInfo, "xcms3", xcms::CentWaveParam(noise = 3E4))
ffKPIC2 <- suppressWarnings(findFeatures(anaInfo, "kpic2", level = 2E5))

fgOptOpenMS <- optimizeFeatureGrouping(ffOpenMS, "openms",
                                       list(maxGroupMZ = c(0.002, 0.007)), maxIterations = 2, parallel = doPar)
# fgOptXCMS <- optimizeFeatureGrouping(optimizedObject(ffOptXCMS), "xcms", list(groupArgs = list(bw = c(22, 28)),
#                                                                               retcorArgs = list(method = "obiwarp")))
fgOptXCMS3 <- optimizeFeatureGrouping(ffXCMS3, "xcms3",
                                      list(groupParams = list(bw = c(22, 28))), maxIterations = 2, parallel = doPar)

fgOptKPIC2 <- optimizeFeatureGrouping(ffKPIC2, "kpic2",
                                      list(groupArgs = list(mz_tolerance = c(0.002, 0.01))), maxIterations = 2,
                                      parallel = doPar)

expInfoPrepForComp <- function(...)
{
    ret <- experimentInfo(...)
    # NOTE: don't want to compare resulting object as it may be irreproducible due to file paths etc
    ret$model <- ret$max_settings <- ret$finalResult <- NULL
    return(ret)
}

test_that("verify feature optimization output", {
    expect_known_value(expInfoPrepForComp(ffOptOpenMS, 1, 1), testFile("ff-opt-oms"))
    # expect_known_show(ffOptOpenMS, testFile("ff-opt-oms-show", text = TRUE))

    # expect_known_value(expInfoPrepForComp(ffOptXCMS, 1, 1), testFile("ff-opt-xcms"))
    # expect_known_show(ffOptXCMS, testFile("ff-opt-xcms-show", text = TRUE))

    expect_length(ffOptEmpty, 1)

    skip_if(utils::packageVersion("xcms") < "3.10") # output changed a little with 3.10
    expect_known_value(expInfoPrepForComp(ffOptXCMS3, 1, 1), testFile("ff-opt-xcms3"))
    # expect_known_show(ffOptXCMS3, testFile("ff-opt-xcms3-show", text = TRUE))
    
    # disabling parallelization slightly changes the output, as this is only for coverage testing we simply skip this check
    skip_if(!doPar)
    expect_known_value(expInfoPrepForComp(ffOptKPIC2, 1, 1), testFile("ff-opt-kpic2"))    
    expect_known_value(expInfoPrepForComp(ffOptEnviPick, 1, 1), testFile("ff-opt-ep"))
    # expect_known_show(ffOptEnviPick, testFile("ff-opt-ep-show", text = TRUE))
})

test_that("verify feature group optimization output", {
    expect_known_value(expInfoPrepForComp(fgOptOpenMS, 1, 1), testFile("fg-opt-oms"))
    # expect_known_show(fgOptOpenMS, testFile("fg-opt-oms-show", text = TRUE))

    # expect_known_value(expInfoPrepForComp(fgOptXCMS, 1, 1), testFile("fg-opt-xcms"))
    # expect_known_show(fgOptXCMS, testFile("fg-opt-xcms-show", text = TRUE))
    
    skip_if(utils::packageVersion("xcms") < "4.2") # output changed a little
    expect_known_value(expInfoPrepForComp(fgOptXCMS3, 1, 1), testFile("fg-opt-xcms3"))
    # expect_known_show(fgOptXCMS3, testFile("fg-opt-xcms3-show", text = TRUE))
    
    # disabling parallelization slightly changes the output, as this is only for coverage testing we simply skip this check
    skip_if(!doPar)
    expect_known_value(expInfoPrepForComp(fgOptKPIC2, 1, 1), testFile("fg-opt-kpic2"))
})

test_that("default param generators", {
    checkmate::expect_list(generateFeatureOptPSet("xcms"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("xcms", method = "matchedFilter"),
                           min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("xcms3"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("xcms3", method = "matchedFilter"),
                           min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("openms"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("envipick"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFeatureOptPSet("kpic2"), min.len = 1, names = "unique")

    checkmate::expect_list(generateFGroupsOptPSet("xcms"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFGroupsOptPSet("xcms", groupArgs = list(method = "nearest"),
                                                  retcorArgs = list(method = "peakgroups")),
                           min.len = 1, names = "unique")
    checkmate::expect_list(generateFGroupsOptPSet("xcms3"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFGroupsOptPSet("xcms3", groupMethod = "nearest",
                                                  retAlignMethod = "peakgroups"),
                           min.len = 1, names = "unique")
    checkmate::expect_list(generateFGroupsOptPSet("openms"), min.len = 1, names = "unique")
    checkmate::expect_list(generateFGroupsOptPSet("kpic2"), min.len = 1, names = "unique")

    checkmate::expect_list(getDefFeaturesOptParamRanges("xcms"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFeaturesOptParamRanges("xcms", "matchedFilter"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFeaturesOptParamRanges("xcms3"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFeaturesOptParamRanges("xcms3", "matchedFilter"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFeaturesOptParamRanges("openms"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFeaturesOptParamRanges("envipick"))
    checkmate::expect_list(getDefFeaturesOptParamRanges("kpic2"), min.len = 1, names = "unique")

    checkmate::expect_list(getDefFGroupsOptParamRanges("xcms"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFGroupsOptParamRanges("xcms3"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFGroupsOptParamRanges("openms"), min.len = 1, names = "unique")
    checkmate::expect_list(getDefFGroupsOptParamRanges("kpic2"), min.len = 1, names = "unique")
})

test_that("plotting works", {
    # expect_doppel("fg-opt-contour", function() plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "contour"))
    expect_plot(plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "contour"))
    # expect_doppel("fg-opt-image", function() plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "image"))
    expect_plot(plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "image"))
    # expect_doppel("fg-opt-persp", function() plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "persp"))
    expect_plot(plot(ffOptOpenMS, paramSet = 1, DoEIteration = 1, type = "persp"))
    expect_error(plot(fgOptOpenMS, paramSet = 1, DoEIteration = 1)) # can't plot with single optimized variable
})

test_that("basic accessors", {
    checkmate::expect_numeric(lengths(ffOptOpenMS), min.len = 1)

    checkmate::expect_list(optimizedParameters(ffOptOpenMS, 1, 1), min.len = 1, names = "unique")
    checkmate::expect_list(optimizedParameters(ffOptOpenMS, 1), min.len = 1, names = "unique")
    checkmate::expect_list(optimizedParameters(ffOptOpenMS), min.len = 1, names = "unique")
    expect_equal(optimizedParameters(ffOptOpenMS), optimizedParameters(ffOptOpenMS)) # only 1 param set

    checkmate::expect_class(optimizedObject(ffOptOpenMS), "featuresOpenMS")
    checkmate::expect_class(optimizedObject(ffOptOpenMS, 1), "featuresOpenMS")
    expect_equal(optimizedObject(ffOptOpenMS), optimizedObject(ffOptOpenMS, 1)) # only 1 param set

    checkmate::expect_class(optimizedObject(fgOptOpenMS), "featureGroupsOpenMS")

    checkmate::expect_list(scores(ffOptOpenMS, 1, 1), min.len = 1, names = "unique")
    checkmate::expect_list(scores(ffOptOpenMS, 1), min.len = 1, names = "unique")
    checkmate::expect_list(scores(ffOptOpenMS), min.len = 1, names = "unique")
    expect_equal(scores(ffOptOpenMS), scores(ffOptOpenMS, 1)) # only 1 param set

    checkmate::expect_list(experimentInfo(ffOptOpenMS, 1, 1), min.len = 1, names = "unique")
})

verifyWithIPO <- as.logical(Sys.getenv("PATROON_VERIFY_IPO", FALSE)) && requireNamespace("IPO", quietly = TRUE)
if (verifyWithIPO)
{
    IPOParamsFeat <- IPO::getDefaultXcmsSetStartingParams("centWave")

    # only optimize one param
    IPOParamsFeat <- sapply(IPOParamsFeat, function(p) if (is.numeric(p) && length(p) == 2) mean(p) else p,
                            simplify = FALSE)
    IPOParamsFeat$min_peakwidth <- c(4, 12)
    IPOResultFeat <- IPO::optimizeXcmsSet(paste0(file.path(anaInfo$path, anaInfo$analysis), ".mzML"),
                                          IPOParamsFeat, nSlaves = 1, plot = FALSE)
    IPOOptParamsFeat <- IPOResultFeat$best_settings$parameters
    # modify output parameter format a bit so that it matches ours...
    IPOOptParamsFeat$peakwidth <- c(IPOOptParamsFeat$min_peakwidth, IPOOptParamsFeat$max_peakwidth)
    IPOOptParamsFeat$prefilter <- c(IPOOptParamsFeat$prefilter, IPOOptParamsFeat$value_of_prefilter)
    IPOOptParamsFeat[c("min_peakwidth", "max_peakwidth", "value_of_prefilter")] <- NULL

    ownResultFeat <- optimizeFeatureFinding(anaInfo, "xcms", IPOParamsFeat, maxModelDeviation = 1)
    ownOptParamsFeat <- optimizedParameters(ownResultFeat)
    ownOptParamsFeat$method <- NULL # IPO doesn't add this

    # same for groups...
    IPOParamsGroup <- IPO::getDefaultRetGroupStartingParams()
    IPOParamsGroup <- sapply(IPOParamsGroup, function(p) if (is.numeric(p) && length(p) == 2) mean(p) else p,
                             simplify = FALSE)
    IPOParamsGroup$mzwid <- c(0.01, 0.03)

    IPOResultGroup <- IPO::optimizeRetGroup(IPOResultFeat$best_settings$xset,
                                            IPOParamsGroup, nSlaves = 1, plot = FALSE)
    IPOOptParamsGroup <- IPOResultGroup$best_settings
    IPOOptParamsGroup$center <- NULL # only added by IPO. UNDONE: important?
    IPOOptParamsGroup$retcorMethod <- NULL # also only added by IPO, no need to compare

    ownParamsGroup <- list(retcorArgs = c(IPOParamsGroup[c("distFunc", "gapInit", "gapExtend", "profStep",
                                                           "plottype", "response", "factorDiag",
                                                           "factorGap", "localAlignment")],
                                          list(method = "obiwarp")),
                           groupArgs = IPOParamsGroup[c("bw", "minfrac", "mzwid", "minsamp", "max")])
    ownResultGroup <- optimizeFeatureGrouping(optimizedObject(ownResultFeat), "xcms", ownParamsGroup,
                                              maxModelDeviation = 1)
    ownOptParamsGroup <- optimizedParameters(ownResultGroup)

    # flatten and re-order so it can be compared
    ownOptParamsGroup <- c(ownOptParamsGroup$groupArgs, ownOptParamsGroup$retcorArgs)
    ownOptParamsGroup <- ownOptParamsGroup[names(IPOOptParamsGroup)]
}

test_that("IPO verification", {
    skip_if_not(verifyWithIPO)

    expect_equivalent(IPOOptParamsFeat, ownOptParamsFeat)
    expect_equivalent(unname(IPOResultFeat$best_settings$result[-1]),
                      as.vector(scores(ownResultFeat), mode = "numeric"))

    expect_equivalent(IPOOptParamsGroup, ownOptParamsGroup)
    # not present at IPO doesn't store best. Should be in first experiment for the test case.
    expect_equivalent(IPOResultGroup[[1]]$max_settings,
                      experimentInfo(ownResultGroup, 1, 1)$max_settings)
})
