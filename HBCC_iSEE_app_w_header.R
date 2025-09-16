###################################################################
#BiocManager::install("iSEE")
library(iSEE); 
library(ggplot2)
#library("SingleCellExperiment") # dont need them as "iSEE" have these all
#library("shiny") # dont need them as "iSEE" have these all
###########################################
### Fetch the data from FigShare/ MendeleyData

# To retrieve an option
# getOption('timeout')
# To set an option
options(timeout=600)
# Download the data from FigShare:
dat <- ("https://figshare.com/ndownloader/files/39305303/sce_dlpfc_sgacc_final.RDS")
download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")
sce_small <- readRDS("sce_dlpfc_sgacc_final.RDS")
# Read tour file
tour <- read.delim("tour.txt", sep=";", stringsAsFactors = FALSE, row.names = NULL)

# ################################################
# Specify number of colors for each cell type
library(RColorBrewer)
n <- 47
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col_vector, n)
names(col_vector) <- as.vector(unique(sce_small$celltype))
#################################################
# Create list of all the panels by starting with an empty list
initial <- list()
############################################################################################
# Custom dummy plot as banner 1
############################################################################################
BANNER_FUN <- function(se, rows, columns, text = "Welcome — read the tips in the tour!") {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = text, size = 6, fontface = 2) +
    labs(title = "LAbs title ") +
    theme(plot.margin = margin(6, 6, 6, 6))
}
BannerPanel <- iSEE::createCustomPlot(
  FUN       = BANNER_FUN,
  restrict  = NULL,
  className = "BannerPanel",
  fullName  = "Demo build"
)
# Create a new class object for new panel
BannerPanel <- new("BannerPanel", text = "Demo build", PanelId = NULL, PanelWidth = 12L, PanelHeight = 400L)

############################################################################################
# Custom plot as banner 1
############################################################################################
# Set variables
title = "Human Brain Collection Core (HBCC) Data Portal"
subtitle = "Single nucleus RNA seq data from sgACC–DLPFC of Controls"
width_px = 1200
height_px = 300
bg = "white"        # slate-900
fg = "royalblue4"        # slate-200
accent = "royalblue4"    # sky-400
accent_width = 0.03  # proportion of width
logo_path =  "nimh-logo.png"      # optional: path to a PNG logo
logo_width = 0.36     # proportion of width
img <- png::readPNG(logo_path)
grob <- grid::rasterGrob(img, interpolate = TRUE)
x_max <- 1 
x_min <- 0 - 0.7
y_max <- 1 
y_min <- y_max - 0.2

# Create a function to use with custom plot
newBANNER_FUN <- function(se, rows, columns, text = "Welcome — read the tips in the tour!") {
  ggplot() +
    # background
    annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5, fill = bg, color = NA) +
    # accent bar on the left
    annotate("rect", xmin = 0, xmax = accent_width, ymin = 0, ymax = 1, fill = accent, color = NA) +
    # title
    annotate("text", x = accent_width + 0.03, y = 0.68, label = title, hjust = 0, vjust = 1, size = 20, color = fg, fontface = "bold") +
    # subtitle (optional)
    annotate("text", x = accent_width + 0.03, y = 0.38, label = subtitle, hjust = 0, vjust = 1, size = 10, color = fg) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    theme_void() + annotation_custom(grob, xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max)
}

BannerPanel2 <- iSEE::createCustomPlot(
  FUN       = newBANNER_FUN,
  restrict  = NULL,
  className = "newBannerPanel",
  fullName  = " "#"National Institutes of Health (NIH)"
)

# Subclass one panel type (e.g., ReducedDimensionPlot) and hide its Data box
setClass("BannerPanel", contains = "newBannerPanel")

setMethod(".hideInterface", "newBannerPanel",
          function(x, field) {
            if (field == "DataBoxOpen" | field == "SelectionBoxOpen") return(TRUE)  # hide the whole Data parameters box
            callNextMethod()
          }
)

initial[["newBannerPanel"]] <- new("newBannerPanel", PanelWidth = 12L, PanelHeight = 400L, PanelId = NULL)

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, 
                                          YAxis = 2L, FacetRowByColData = "Barcode", FacetColumnByColData = "Barcode", 
                                          ColorByColumnData = "ID", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_ID", 
                                          SizeByColumnData = "sum", FacetRowBy = "None", FacetColumnBy = "None", 
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                          ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---", 
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "donor4_AAACCCAAGAGTCTTC.1", 
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                          VisualChoices = c("Color", "Shape"), ContourAdd = FALSE, 
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                          CustomLabelsText = "donor4_AAACCCAAGAGTCTTC.1", FontSize = 1, 
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                          LabelCenters = FALSE, LabelCentersBy = "Barcode", LabelCentersColor = "#000000", 
                                          VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version", 
                                                                                                             "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = TRUE, 
                                          SelectionHistory = list())



################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
                                      XAxisColumnData = "broad.class", XAxisFeatureName = "SNAP25",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "SNAP25", YAxisFeatureSource = "RowDataTable1",
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "Barcode",
                                      FacetColumnByColData = "Barcode", ColorByColumnData = "broad.class",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "sample_ID", SizeByColumnData = "sum", FacetRowBy = "None",
                                      FacetColumnBy = "None", ColorBy = "Column data", ColorByDefaultColor = "#000000",
                                      ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---",
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "{{cellone}}",
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                      CustomLabels = FALSE, CustomLabelsText = "{{cellone}}",
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Barcode",
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 4L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 600L,
                                      PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                      ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                      ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                      ColumnSelectionRestrict = TRUE, SelectionHistory = list())
################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "SATB2\nGAD2\nAQP4\nMOG\nMEGF11\nPTPRC\nFLT1", 
                                        ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("neuron", 
                                                                                                           "celltype"), RowData = character(0), CustomBounds = FALSE, 
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE, 
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow", 
                                        ShowDimNames = "Rows", LegendPosition = "Right", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
                                        ShowColumnSelection = FALSE, OrderColumnSelection = TRUE, 
                                        VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", 
                                                                                                           "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 12L, 
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                        RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                        SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nFeature_RNA", 
                                    XAxisColumnData = "celltype", FacetRowByColData = "sample_ID", 
                                    FacetColumnByColData = "sample_ID", ColorByColumnData = "celltype", 
                                    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                    ShapeByColumnData = "sample_ID", SizeByColumnData = "nCount_RNA", 
                                    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                    ColorByDefaultColor = "#000000", ColorByFeatureName = "RP11-34P13.3", 
                                    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                    ColorBySampleName = "2543_sgACC_2_AAACCTGAGATAGGAG", ColorBySampleSource = "---", 
                                    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
                                    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                    CustomLabelsText = "2543_sgACC_2_AAACCTGAGATAGGAG", FontSize = 1, 
                                    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                    LabelCenters = FALSE, LabelCentersBy = "sample_ID", LabelCentersColor = "#000000", 
                                    VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", 
                                                                                                       "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L, 
                                    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "SNAP25", Search = "", SearchColumns = c("",
                                                                                                      "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                                                                                      "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0),
                                  VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
                                                                                                     "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L,
                                  PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())


######################################

sce_small <- registerAppOptions(sce_small, color.maxlevels = 47)

iSEE(
  sce_small,
  tour = tour,
  appTitle = "An official website of the United States government",
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    celltype = function(n) {
      col_vector[!grepl("drop", names(col_vector))]
    }
  ))
)

############################################################################



