###################################################################
library(iSEE); 
library(ggplot2)
library(memoise)
library(htmltools)
library(base64enc)

#library("SingleCellExperiment") # dont need them as "iSEE" have these all
#library("shiny") # dont need them as "iSEE" have these all
###########################################
### Fetch the data from FigShare/ MendeleyData
# Download the data from FigShare:
# dat <- ("https://figshare.com/ndownloader/files/39305303/sce_dlpfc_sgacc_final.RDS")
# download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")
# ##################################################################
# To retrieve an option
# getOption('timeout')
# To set an option
options(timeout=600)
# ##################################################################
# 1) Memoised loader for your SCE
get_sce <- memoise(function() {
  # OR read input from local 
  sce <- readRDS("sce_dlpfc_sgacc_final.RDS")
  sce
})
# 2) Use the cached SCE to build the app
sce_small <- get_sce()
# ##################################################################
# Read tour file
tour <- read.delim("tour.txt", sep=";", stringsAsFactors = FALSE, row.names = NULL)
# ##################################################################
# Specify number of colors for each cell type
library(RColorBrewer)
n <- 47
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col_vector, n)
names(col_vector) <- as.vector(unique(sce_small$celltype))
###################################################################

###################################################################
# Create list of all the panels by starting with an empty list
initial <- list()
############################################################################################
# Custom plot for title
############################################################################################
# Set variables
logo_path_t =  "nimh-logo.png"      # optional: path to a PNG logo
img_t <- png::readPNG("nimh-logo.png")
grob_t <- grid::rasterGrob(img_t, interpolate = TRUE)
# Create a function to use with custom plot
TITLE_FUN <- function(se, rows, columns, text = "Welcome — read the tips in the tour!") {
  ggplot() +
    # background bar
    geom_rect(aes(xmin = 0, xmax = 50, ymin = 0, ymax = 10),
              fill = "white", color = NA) +
    # Full-width bordered box for the title
    geom_rect(aes(xmin = -50, xmax = 50, ymin = 0.1, ymax = 0.5),
              fill = "white", color = "royalblue4", linewidth = 1) +
    # Left block - adjusted x position to make room for logo
    annotate("text", x = -0.03, y = 0.3, label = "Human Brain Collection Core Data Portal At NIMH",
             hjust = 0, color = "royalblue4", fontface = "bold", size = 22) +
    # # Left lower text
    # annotate("text", x = 0, y = 0.0, label = "Study: Single nucleus RNA seq data from sgACC–DLPFC of Controls",
    #          hjust = 0, color = "royalblue4", size = 12) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = TRUE) +
    # Add the logo grob at LEFT edge - adjusted coordinates
    theme_void() + 
    annotation_custom(grob_t, xmin = -0.1, xmax = 0.3, ymin = 0.6, ymax = 1)
}

TitlePanel <- iSEE::createCustomPlot(
  FUN       = TITLE_FUN,
  restrict  = NULL,
  className = "TitlePanel",
  fullName  = "iSEE App ID = HBCC"
)

# Subclass one panel type (e.g., ReducedDimensionPlot) and hide its Data box
setClass("Title", contains = "TitlePanel")

setMethod(".hideInterface", "TitlePanel",
          function(x, field) {
            if (field == "DataBoxOpen" | field == "SelectionBoxOpen") return(TRUE)  # hide the whole Data parameters box
            callNextMethod()
          }
)

# Add custom output method to include hyperlink below the plot
setMethod(".defineOutput", "TitlePanel", function(x) {
  plot_name <- .getEncodedName(x)
  tagList(
    plotOutput(plot_name, height = "300px"),  # Changed from fixed 400px to auto
    tags$div(
      style = "padding: 10px; text-align: left; padding-left: 20px; font-size: 30px; color: royalblue4",
      tags$span(
        tags$strong("Study Name:"),  # Bold PMID text
        " ",
        tags$a(
          href = "https://pubmed.ncbi.nlm.nih.gov/37037607/",
          target = "_blank",
          "Cellular Diversity in Human Subgenual Anterior Cingulate and Dorsolateral Prefrontal Cortex by Single-Nucleus RNA-Sequencing",
          style = "color: royalblue4; font-size: 28px; text-decoration: underline;"
        )
      )
    ),
    tags$div(
      style = "padding: 10px; text-align: left; padding-left: 20px; font-size: 30px; color: royalblue4",
      tags$span(
        tags$strong("PMID:"),  # Bold PMID text
        " ",
        tags$a(
          href = "https://pubmed.ncbi.nlm.nih.gov/37037607/",
          target = "_blank",
          "37037607",
          style = "color: royalblue4; font-size: 28px; text-decoration: underline;"
        )
      )
    )
  )
})

initial[["TitlePanel"]] <- new("TitlePanel", PanelWidth = 12L, PanelHeight = 400L, PanelId = 1L)
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

############################################################################################
# Custom plot for FOOTER
############################################################################################
# Create a function to use with custom plot
FOOTER_FUN <- function(se, rows, columns, text = "Welcome — read the tips in the tour!") {
  # Canvas size in arbitrary units
  # We'll treat x in [0, 100], y in [0, 10]
  ggplot() +
    # background bar
    geom_rect(aes(xmin = 0, xmax = 130, ymin = 0, ymax = 10),
              fill = "#0B3D5B", color = NA) +
    
    # Left block
    annotate("text", x = 1, y = 8.5, label = "National Institute of Mental Health",
             hjust = 0, color = "white", fontface = "bold", size = 10) +
    annotate("text", x = 1, y = 7, label = "at the National Institutes of Health",
             hjust = 0, color = "white", size = 7) +
    # Right block (top-right)
    annotate("text", x = 99, y = 8.5, label = "Contact Us",
             hjust = 1, color = "white", fontface = "bold", size = 10) +
    
    # Right block (links row)
    annotate("text", x = 99, y = 7, label = "nimhinfo@nih.gov",
             hjust = 1, color = "white", size = 7) +
    
    # Right block (bottom-right small line)
    annotate("text", x = 99, y = 5.5,
             label = "U.S. Department of Health and Human Services",
             hjust = 1, color = "white", size = 7) +
    annotate("text", x = 99, y = 4.5,
             label = "National Institutes of Health",
             hjust = 1, color = "white", size = 7) +
    annotate("text", x = 99, y = 3.5,
             label = "National Institutes of Mental Health", 
             hjust = 1, color = "white", size = 7) +
    annotate("text", x = 99, y = 2.5,
             label = "USA.gov",
             hjust = 1, color = "white", size = 7) +
    
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 10), expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
}

FooterPanel <- iSEE::createCustomPlot(
  FUN       = FOOTER_FUN,
  restrict  = NULL,
  className = "FooterPanel",
  fullName  = "iSEE App ID = HBCC"
)

# Subclass one panel type (e.g., ReducedDimensionPlot) and hide its Data box
setClass("Footer", contains = "FooterPanel")
setMethod(".hideInterface", "FooterPanel",
          function(x, field) {
            if (field == "DataBoxOpen" | field == "SelectionBoxOpen") return(TRUE)  # hide the whole Data parameters box
            callNextMethod()
          }
)

# Encode your images
socialMedia_base64 <- base64encode("socialMedia_1.png")

setMethod(".defineOutput", "FooterPanel", function(x) {
  plot_name <- .getEncodedName(x)
  tagList(
    plotOutput(plot_name, height = "200px"),
    tags$div(
      style = "padding: 1px; text-align: left; padding-left: 10px;",
      tags$div(
        # Add "Follow Us" text
        tags$span(
          "Follow Us",
          style = "color: #0B3D5B; font-weight: bold; font-size: 26px; margin-right: 0px; vertical-align: middle;"
        ),
        tags$a(
          href = "https://www.nih.gov/news-events/nih-social-media",
          target = "_blank",
          tags$img(
            src = paste0("data:image/png;base64,", socialMedia_base64), 
            height = "26px", 
            alt = "NIH social media icons",
            style = "background-color: grey;"  # Add background color here
          )
        )
      )
    )
  )
})

### Panel ends  
initial[["FooterPanel"]] <- new("FooterPanel", PanelWidth = 12L, PanelHeight = 800L, PanelId = 1L)
################################################################################

################################################################################
# # Overriding the default colors in the package
iSEEOptions$set(panel.color=c(FeatureAssayPlot="#3565AA", RowDataPlot="#3565AA", ColumnDataPlot="#3565AA", ComplexHeatmapPlot="#3565AA", RowDataTable="#3565AA", TitlePanel="#3565AA", BannerPanel1 = "white", FooterPanel = "#3565AA"))

sce_small <- registerAppOptions(sce_small, color.maxlevels = 47)
############################################################################

# Add US flag along with text title
library(base64enc)
# Encode your images
flag_base64 <- base64encode("flag.png")
logo1_base64 <- base64encode("nimh-logo.png")
link_title <- "An official website of the United States government"

# Define custom CSS for iSEE app title
custom_css <- "
  .main-header .logo {
    text-align: left;
    padding-left: 75px; /* Adjust padding as needed */;
  }
"

# Run the iSEE wrapper function to launch the app
iSEE(
  sce_small,
  tour = tour,
  appTitle = tags$img(
    height = "50px",
    tags$img(
      src = paste0("data:image/png;base64,", flag_base64), 
      height = "15px", 
      alt = "USA flag logo"), role = "img",  # Explicit role
    `aria-label` = "Application logo",  # ARIA label
    style = "float: left; margin-right: 10px",
    tags$span(
      link_title, style = "font-size: 15px; font-weight: regular; font-family: 'Open Sans', sans-serif; color: black"),
    tags$head(tags$style(HTML(custom_css)))
  ),
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    celltype = function(n) {col_vector[!grepl("drop", names(col_vector))]}))
)

############################################################################
############################################################################


