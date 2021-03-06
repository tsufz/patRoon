#' @include main.R
#' @include feature_groups.R

getEICUI <- function(rtRange, mzWindow)
{
    rtRange <- round(rtRange, 1)
    showOpts <- c("Keep", "Don't keep")

    fillPage(
        tags$head(includeScript(system.file("js", "utils-EIC.js", package = "patRoon"))),
        # tags$script("setShortcuts();"),

        title = "EIC tool",

        fillCol(
            flex = c(NA, 1, NA),

            fillRow(
                height = 40,
                flex = c(NA, NA, NA, 1, NA),

                fillCol(
                    width = 45,
                    actionButton("previousGroup", "", icon("arrow-left"), onclick = "selectPrevFGroup();")
                ),
                fillCol(
                    width = 45,
                    actionButton("nextGroup", "", icon("arrow-right"),  onclick = "selectNextFGroup();")
                ),
                fillCol(
                    width = 100,
                    actionButton("toggleGroup", "Toggle group", icon("toggle-on"))
                ),

                fillCol(
                    strong(style = "font-size: 200%; text-align: center;", textOutput("groupTitle"))
                ),

                fillCol(
                    width = 150,
                    actionButton("applyClose", "Apply & Close", icon("save"))
                )
            ),

            fillRow(
                flex = c(NA, 1),

                fillCol(
                    width = 160,

                    wellPanel(
                        style = "overflow-y: auto; height: 100%;",

                        radioButtons("plotType", "Plot type", c("Interactive", "Static")),
                        radioButtons("retUnit", "Retention unit", c("Seconds", "Minutes")),
                        numericInput("mzWindow", "m/z width", mzWindow, 0.0001, 1, 0.001),

                        # UNDONE disabled for now
                        # div(
                        #     title = "Total retention time range plotted relative to group",
                        #     numericInput("rtRange", "Plot range (%)", 100, 0, 100, 10)
                        # ),
                        div(
                            title = "Retention time offset for default zoomed group",
                            # UNDONE disabled for now
                            # numericInput("rtZWindow", "Zoom window (s)", 20, 0, rtRange, step = 2)
                            numericInput("rtZWindow", "Zoom window (s)", 20, 0, step = 2)
                        ),

                        checkboxGroupInput("showWhat", "Show groups", showOpts, showOpts)
                    )
                ),

                fillCol(
                    uiOutput("plot", inline = TRUE)
                )
            ),

            fillRow(
                flex = c(1, NA, 1),
                height = 260,

                fillCol(
                    div(
                        style = "border: 1px solid black; margin: 5px;",
                        rhandsontable::rHandsontableOutput("groupHot")
                    )
                ),

                fillCol(
                    br()
                ),

                fillCol(
                    div(
                        style = "border: 1px solid black; margin: 5px;",
                        rhandsontable::rHandsontableOutput("analysesHot")
                    )
                )
            )
        )
    )
}

#' @details \code{checkChromatograms} is used to review chromatographic
#'   information for feature groups. This is especially useful to get a visual
#'   impression of the quality of detected features. In addition, this function
#'   may be used to remove unwanted (\emph{e.g.} outlier) features. Better
#'   performance is often obtained when an external browser is used to use this
#'   Shiny application. Furthermore, when a large \code{featureGroups} object is
#'   used it is recommended to limit the number of analyses/feature groups by
#'   subsetting the object.
#'
#' @param mzWindow Default \emph{m/z} window to be used for creating extracted
#'   ion chromatograms (EICs).
#' @param enabledFGroups A logical vector that states for each feature group
#'   whether it should be kept (\code{TRUE}) or not (\code{FALSE}). The order is
#'   the same as the \code{fGroups} parameter. If \code{NULL} then all feature
#'   groups are considered to be kept.
#'
#' @return \code{checkChromatograms} returns a logical vector for all feature
#'   groups that were selected to be kept (\code{TRUE}) or not (\code{FALSE}).
#'   This result can be passed to the \code{enabledFGroups} parameter for
#'   subsequent calls to \code{checkChromatograms} in order to restore the
#'   keep/not keep state from a previous call. To actually remove unwanted
#'   feature groups the object should be subset by the subsetting
#'   (\code{\link{[}}) operator to which the return value should be passed as
#'   the second parameter.
#'
#' @rdname GUI-utils
#' @aliases checkChromatograms
#' @export
setMethod("checkChromatograms", "featureGroups", function(fGroups, mzWindow, enabledFGroups)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(mzWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(enabledFGroups, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gInfo <- groupInfo(fGroups)
    gTable <- groups(fGroups)
    avgGTable <- averageGroups(fGroups)
    ftindex <- groupFeatIndex(fGroups)

    if (is.null(enabledFGroups))
        enabledFGroups <- rep(TRUE, nrow(gInfo))

    # UNDONE disabled for now
    # rtRange <- max(sapply(xrs, function(xr) xr@scantime[length(xr@scantime)]))
    rtRange <- 1

    EICPreviews <- getEICsForFGroups(fGroups, 10, mzWindow, topMost = 1, onlyPresent = FALSE)
    # format is in [[ana]][[fGroup]], since we only took top most intensive we can throw away the ana dimension
    EICPreviews <- Reduce(modifyList, EICPreviews)

    hotOpts <- list(rowHeaderWidth = 40, readOnly = TRUE, disableVisualSelection = "area",
                    columnSorting = TRUE, sortIndicator = TRUE, selectCallback = TRUE,
                    currentRowClassName = "currentRow", stretchH = "all",
                    preventOverflow = "horizontal", multiSelect = FALSE,
                    outsideClickDeselects = FALSE, manualColumnResize = TRUE,
                    rowHeaders = NULL)

    anaColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(nrow(anaInfo))
    anaColorsTrans <- adjustcolor(anaColors, alpha.f = 0.5)

    server <- function(input, output, session)
    {
        rValues <- reactiveValues(enabledFGroups = enabledFGroups,
                                  enabledHotFGroups = enabledFGroups,
                                  currentFGroup = rownames(gInfo)[1],
                                  enabledAnalyses = anaInfo$analysis,
                                  enabledHotAnalyses = rep(TRUE, nrow(anaInfo)))

        # getCurrentGroup <- reactive({
        #     if (is.null(input$groupHot_select))
        #         return(rownames(gInfo)[1])
        #     return(rownames(gInfo)[input$groupHot_select$select$r])
        # })

        getCurrentAnalysis <- reactive({
            if (is.null(input$analysesHot_select))
                return(anaInfo$analysis[1])
            return(anaInfo$analysis[input$analysesHot_select$select$rAll[1]])
        })

        plotInfo <- reactive({
            if (length(rValues$enabledAnalyses) == 0)
                return(NULL)

            g <- rValues$currentFGroup
            if (is.null(g) || !nzchar(g))
                return(NULL)

            getFTData <- function(ana)
            {
                    ftind <- ftindex[[g]][match(ana, anaInfo$analysis)]
                    if (ftind == 0)
                        return(data.table(ret=NA, retmin=NA, retmax=NA, mz=NA))
                    return(fTable[[ana]][ftind])
            }
            fts <- rbindlist(sapply(rValues$enabledAnalyses, getFTData, simplify = FALSE), fill = TRUE)

            if (all(is.na(fts$ret)))
            {
                # no features present in current set of analyses
                # try to fall back to disabled analyses
                disAna <- setdiff(anaInfo$analysis, rValues$enabledAnalyses)
                if (length(disAna) > 0)
                    fts <- rbindlist(sapply(disAna, getFTData, simplify = FALSE), fill = TRUE)

                if (length(disAna) == 0 || all(is.na(fts$ret)))
                    return(NULL) # still no luck
            }

            ret <- list(avgRt = mean(fts[, ret], na.rm = TRUE),
                        minRt = min(fts[, retmin], na.rm = TRUE),
                        maxRt = max(fts[, retmax], na.rm = TRUE),
                        maxInt = max(fts[, intensity], na.rm = TRUE))

            # UNDONE disabled for now
            # rtwin <- max(input$rtRange * rtRange / 100, input$rtZWindow)
            rtwin <- input$rtZWindow

            # subset fGroups. NOTE: don't use [ operator to avoid removing empty groups
            fg <- removeAnalyses(fGroups, match(setdiff(anaInfo$analysis, rValues$enabledAnalyses),
                                                anaInfo$analysis))
            fg <- removeGroups(fg, which(g != names(fGroups)))

            EICs <- getEICsForFGroups(fg, rtwin, input$mzWindow, topMost = NULL, onlyPresent = FALSE)
            # EICs are in [[ana]][[fgroup]] --> only have one fgroup so get rid of that dimension
            ret$data <- lapply(EICs, "[[", 1)

            ret$peaks <- fts[, c("retmin", "retmax")]
            ret$zoomRtRange <- c(ret$minRt - input$rtZWindow, ret$maxRt + input$rtZWindow)

            if (input$retUnit == "Minutes")
            {
                mod <- c("avgRt", "minRt", "maxRt", "peaks", "zoomRtRange")
                ret[mod] <- lapply(ret[mod], "/", 60)
                ret$data <- lapply(ret$data, function(d) { d$time <- d$time / 60; return(d) })
            }

            return(ret)
        })

        fGroupData <- reactive({
            # initialize data
            ret <- data.frame(group = rownames(gInfo),
                              EIC = sapply(rownames(gInfo),
                                           function(g) jsonlite::toJSON(list(values = EICPreviews[[g]]$intensity,
                                                                             xvalues = EICPreviews[[g]]$time,
                                                                             options = list(type = "line", height = 50)))),
                              keep = rValues$enabledHotFGroups, retention = gInfo$rts, mz = gInfo$mzs, stringsAsFactors = FALSE)
            ret[, unique(anaInfo$group)] <- t(avgGTable)
            ret[, "retention"] <- if (input$retUnit == "Minutes") gInfo$rts / 60 else gInfo$rts

            if (!"Keep" %in% input$showWhat)
                ret <- ret[!ret$keep, ]

            if (!"Don't keep" %in% input$showWhat)
                ret <- ret[ret$keep, ]

            return(ret)
        })

        analysesData <- reactive({
            return(data.frame(color = anaColors, analysis = anaInfo$analysis,
                              enabled = rValues$enabledHotAnalyses,
                              group = anaInfo$group, blank = anaInfo$blank, stringsAsFactors = FALSE))
        })

        observeEvent(input$toggleGroup, {
            ind <- input$groupHot_select$select$rAll[1]
            rValues$enabledHotFGroups[ind] <- !rValues$enabledHotFGroups[ind]
        })
    
        observeEvent(input$applyClose, {
            enabledFGroups <<- rValues$enabledFGroups
            stopApp()
        })

        observeEvent(input$groupHot_select$select$r, {
            rValues$currentFGroup <- rownames(gInfo)[input$groupHot_select$select$rAll[1]]
        })

        observeEvent(input$addDAEICs, {
            g <- rValues$currentFGroup
            analyses <- switch(input$addDAEICs[[1]],
                               selected = getCurrentAnalysis(),
                               all = anaInfo$analysis,
                               enabled = rValues$enabledAnalyses)

            if (!is.null(g) && nzchar(g) && length(analyses) > 0)
            {
                bgsubtr <- input$addDAEICs[[2]] == "1"
                for (f in analyses)
                    addDAEIC(f, anaInfo$path[match(f, anaInfo$analysis)], gInfo[g, "mzs"], input$mzWindow, bgsubtr = bgsubtr)
            }
        })

        observeEvent(input$enableAllGroups, {
            rValues$enabledHotFGroups = FALSE # HACK: trigger update
            rValues$enabledHotFGroups <- rep(TRUE, length(rValues$enabledHotFGroups))
        })
        observeEvent(input$disableAllGroups, {
            rValues$enabledHotFGroups = FALSE # HACK: trigger update
            rValues$enabledHotFGroups <- rep(FALSE, length(rValues$enabledHotFGroups))
        })

        observeEvent(input$enableAllAnalyses, {
            rValues$enabledHotAnalyses <- FALSE # HACK: trigger update
            rValues$enabledHotAnalyses <- rep(TRUE, nrow(anaInfo))
        })
        observeEvent(input$disableAllAnalyses, {
            rValues$enabledHotAnalyses <- FALSE # HACK: trigger update
            rValues$enabledHotAnalyses <- rep(FALSE, nrow(anaInfo))
        })

        observeEvent(input$groupHot, {
            # HACK: input$groupHot$params$maxRows: make sure we don't have empty table as hot_to_r errors otherwise
            if (input$groupHot$params$maxRows > 0)
            {
                df <- rhandsontable::hot_to_r(input$groupHot)
                rValues$enabledFGroups[match(df$group, rownames(gInfo))] <- df$keep
            }
        })

        observeEvent(input$analysesHot, {
            df <- rhandsontable::hot_to_r(input$analysesHot)
            ea <- df[df[["enabled"]] == TRUE, "analysis"]
            if (!isTRUE(all.equal(ea, rValues$enabledAnalyses)))
                rValues$enabledAnalyses <- ea
        })

        output$groupTitle <- renderText({
            rValues$currentFGroup
        })

        output$plot <- renderUI({
            if (input$plotType == "Interactive")
                tagList(plotly::plotlyOutput("plotInteractive", height = "100%"))
            else
                tagList(plotOutput("plotStatic", height = "100%"))
        })

        output$plotInteractive <- plotly::renderPlotly({
            pinfo <- plotInfo()
            p <- plotly::plot_ly(type="scatter", mode = "lines", hoverinfo = "none") %>%
                plotly::config(displaylogo = FALSE, scrollZoom = TRUE,
                               modeBarButtonsToRemove = c("hoverClosestCartesian", "hoverCompareCartesian"))

            if (!is.null(pinfo)) # NULL if no data (no active group, no enabled analyses, ...)
            {
                ea <- rValues$enabledAnalyses

                for (i in seq_along(pinfo$data))
                {
                    p <- plotly::add_trace(p, x = pinfo$data[[i]]$time, y = pinfo$data[[i]]$intensity,
                                           name = ea[i],
                                           line = list(width = if (getCurrentAnalysis() == anaInfo$analysis[i]) 2 else 1,
                                                       color = anaColors[i]))

                    if (getCurrentAnalysis() == ea[i])
                    {
                        pmin <- pinfo$peaks[["retmin"]][i]
                        pmax <- pinfo$peaks[["retmax"]][i]

                        if (!is.na(pmin))
                        {
                            sdata <- pinfo$data[[i]][numGTE(pinfo$data[[i]]$time, pmin) & numLTE(pinfo$data[[i]]$time, pmax), ]
                            p <- plotly::add_trace(p, x = sdata$time, y = sdata$intensity,
                                                   mode = "none", fill = "tozeroy",
                                                   fillcolor = anaColorsTrans[i])
                        }
                    }
                }

                p <- plotly::layout(p, showlegend = FALSE, dragmode = "pan", plot_bgcolor = "#F5FFFA",
                                    margin = list(t = 0),
                                    xaxis = list(title = "Retention time", range = pinfo$zoomRtRange),
                                    yaxis = list(title = "Intensity", exponentformat = "E", range = c(0, pinfo$maxInt * 1.1)))
            }
            return(p)
        })

        output$plotStatic <- renderPlot({
            pinfo <- plotInfo()

            if (!is.null(pinfo)) # NULL if no data (no active group, no enabled analyses, ...)
            {
                ea <- rValues$enabledAnalyses
                for (i in seq_along(pinfo$data))
                {
                    params <- list(x = pinfo$data[[i]]$time, y = pinfo$data[[i]]$intensity, type = "l",
                                   col = anaColors[i], lwd = if (getCurrentAnalysis() == ea[i]) 2 else 0.5)

                    if (i == 1)
                    {
                        do.call(plot, c(list(xlab = "Retention time", ylab = "Intensity", xlim = pinfo$zoomRtRange,
                                             ylim = c(0, pinfo$maxInt * 1.1), bty = "L"), params))
                        # bg: http://stackoverflow.com/a/7237066 (NOTE: need to plot twice!)
                        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F5FFFA", border = NA)
                    }

                    do.call(points, params)

                    if (getCurrentAnalysis() == ea[i])
                    {
                        pmin <- pinfo$peaks[["retmin"]][i]
                        pmax <- pinfo$peaks[["retmax"]][i]

                        if (!is.na(pmin))
                        {
                            sdata <- pinfo$data[[i]][numGTE(pinfo$data[[i]]$time, pmin) & numGLE(pinfo$data[[i]]$time, pmax), ]
                            polygon(c(sdata$time, rev(sdata$time)), c(sdata$intensity, rep(0, length(sdata$intensity))),
                                    col = anaColorsTrans[i], border = NA)
                        }
                    }

                    grid()
                }
            }
            else
            {
                # empty plot
                plot(1, type = "n", xlab = "Retention time", ylab = "Intensity")
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F5FFFA", border = NA)
                grid()
            }

        })

        output$groupHot <- rhandsontable::renderRHandsontable({
            gData <- fGroupData()

            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(gData, colHeaders = c("Feature group", "EIC max", "Keep", "Retention", "m/z", unique(anaInfo$group)),
                                  width = NULL, height = 250, maxRows = nrow(gData)),
                             hotOpts)) %>%
                rhandsontable::hot_cols(valign = "htMiddle", fixedColumnsLeft = 2) %>%
                rhandsontable::hot_rows(rowHeights = 50) %>%
                rhandsontable::hot_col("Keep", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_col("EIC max", renderer = htmlwidgets::JS("renderSparkline"))

            hot$x$contextMenu <- list(items = list(
                addEICsEnabled = list(
                    name = "Add DA EICs enabled",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 0, Math.random()]); }'
                    )
                ),
                addEICsEnabledBG = list(
                    name = "Add DA EICs enabled w/ bg",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 1, Math.random()]); }'
                    )
                ),
                addEICsAll = list(
                    name = "Add DA EICs all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 0, Math.random()]); }'
                    )
                ),
                addEICsAllBG = list(
                    name = "Add DA EICs all w/ bg",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 1, Math.random()]); }'
                    )
                ),
                enableAll = list(
                    name = "Enable all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("enableAllGroups", Math.random()); }'
                    )
                ),
                disableAll = list(
                    name = "Disable all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("disableAllGroups", Math.random()); }'
                    )
                )
            ))

            return(hot)
        })

        output$analysesHot <- rhandsontable::renderRHandsontable({
            hot <- do.call(rhandsontable::rhandsontable,
                           c(list(analysesData(), height = 250, maxRows = nrow(analysesData())),
                             hotOpts)) %>%
                rhandsontable::hot_col("enabled", readOnly = FALSE, halign = "htCenter") %>%
                rhandsontable::hot_col("color", renderer = "function (instance, td, row, col, prop, value, cellProperties)
                        { td.style.background = value; }")

            hot$x$contextMenu <- list(items = list(
                addEICSelected = list(
                    name = "Add DA EIC",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["selected", 0, Math.random()]); }'
                    )
                ),
                addEICSelectedBG = list(
                    name = "Add DA EIC w/ bg",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["selected", 1, Math.random()]); }'
                    )
                ),
                addEICsEnabled = list(
                    name = "Add DA EICs enabled",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 0, Math.random()]); }'
                    )
                ),
                addEICsEnabledBG = list(
                    name = "Add DA EICs enabled w/ bg",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 1, Math.random()]); }'
                    )
                ),
                addEICsAll = list(
                    name = "Add DA EICs all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 0, Math.random()]); }'
                    )
                ),
                addEICsAllBG = list(
                    name = "Add DA EICs all w/ bg",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 1, Math.random()]); }'
                    )
                ),
                enableAll = list(
                    name = "Enable all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("enableAllAnalyses", Math.random()); }'
                    )
                ),
                disableAll = list(
                    name = "Disable all",
                    callback = htmlwidgets::JS(
                        'function(key, options) { Shiny.onInputChange("disableAllAnalyses", Math.random()); }'
                    )
                )
            ))

            return(hot)
        })
    }

    runApp(shinyApp(getEICUI(rtRange, mzWindow), server))
    return(enabledFGroups)
})
