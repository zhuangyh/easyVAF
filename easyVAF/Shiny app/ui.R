library(shiny)
library(shinythemes)
library(tidyverse)
library(tibble)
library(ggplot2)
library(lme4)
library(DT)
library(ggiraph)
library(cowplot)
library(aod)
library(easyVAF)

ui <- navbarPage(
  title = "EXome Analysis",
  theme = shinytheme("flatly"),

  analysis_page <- tabPanel(
    title = "Analysis",
    titlePanel("Analysis"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1","Upload a data file", accept=".csv", width = "100%"),
        uiOutput("groups", width = "100%"),
        selectInput("method_test", label = "Select test", multiple = FALSE,
                    choices = c("Default", "Fisher's exact", "Pearson's chi-squared", "Beta-binomial"),
                    selected = "Default",
                    width = "100%"),
        selectInput("qc_method", label = "Select model for quality checking",
                    choices = c("Linear", "Linear mixed-effects"),
                    width = "100%"),

        sliderInput("digits", label ="Choose number of digits for rounding",
                    min = 1, max = 5, value = 3, step =1, ticks = FALSE),

        selectInput("table_display", label = "Display significant difference based on:", multiple = TRUE,
                    choices = c("Raw p value", "Adjusted p value", "Effect size > 0.2"),
                    selected = c("Raw p value", "Adjusted p value", "Effect size > 0.2")),

        actionButton("run_analysis", label = "Run analysis", width="100%", icon("play"),
                     style="color: #000000; background-color: #04E69E"),
        br(),
        width = 2
        ),

      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Exploratory Plots",
            br(),
            p("Displayed on this page are exploratory analysis plots to examine the distribution and variability of ", em("read depths (DP), variant counts (VC), and variant allele frequencies (VAF)"), " of loci."),
            p("Use your mouse to select an area of the plot and double click the selection to zoom in. Double-click again anywhere on the plot to return to the default view."),
            br(),
            plotOutput("exploratory_plot1", height = "700px", width = "100%", dblclick = "plot1_dblclick", brush = brushOpts(id="plot1_brush",
                                                                                                         resetOnNew = TRUE)),
            plotOutput("exploratory_plot2", height = "700px", width = "100%", dblclick = "plot2_dblclick", brush = brushOpts(id = "plot2_brush",
                                                                                                         resetOnNew = TRUE)),
            plotOutput("exploratory_plot3", height = "700px", width = "100%", dblclick = "plot3_dblclick", brush = brushOpts(id = "plot3_brush",
                                                                                                         resetOnNew = TRUE)),
            br(),
            br()
          ),
          tabPanel(
            title = "Quality of Sample",
            br(),

            p("This tab reports on the degree of biological variability among the mice, which may be a concern to the investigator if large."),
            br(),

            plotOutput("QCplot", height = 700),
            br(),

            p("Displayed below are the results of a test to assess the general variability of VAF among the experiment mice. In the side panel, you may choose to run a regular linear model or a linear mixed model, which adjusts for the non-independence of VAF within the same chromosome."),
            br(),

            h4("Summary of likelihood ratio test results:"),
            tableOutput("QC_chisquare"),
            br(),
            h4("Conclusion:"),
            htmlOutput("concl_bio_var"),
            br(),
            br()
          ),
          tabPanel(
            title = "Statistical Testing",
            br(),
            htmlOutput("stat_blurb"),
            br(),
            textOutput("stat_graphic_desc"),


            uiOutput("compare_methods"),

            br(),
            h4("Table of top selected loci"),
            DT::dataTableOutput("top_loci"),
            br(),

            downloadButton("download_toptable", "Download table of all loci as .csv"),
            br(),
            br()

          )

        )
      )
    )
  ),
  about_page <- tabPanel(
    title = "About",
    titlePanel("About this app"),
    h5(strong("Authors: "), "..."),

    hr(),

    h2("Data Requirements"),
    p("Please upload your data as a ", em("comma separated values (.csv) file"), " and make sure that it contains the following columns with these exact names:"),
    HTML("<ul><li><b>locus</b>: locus ID</li><li><b>chrom</b>: chromosome information (for linkage disequilibrium adjustment in QC test, if desired)</li><li><b>vc</b>: variant count</li><li><b>dp</b>: read depth</li><li><b>individual</b>: mouse ID (for QC test) </li><li><b>sample</b></li></ul>"),

    hr(),

    h2("Package Dependencies"),
    p("This application makes use of the following packages:"),
    HTML("<ul><li><b>tidyverse</li><li>tibble</li><li>ggplot2</li><li>lme4</li><li>DT</li><li>ggiraph</li><li>cowplot</li><li>aod</li><li>easyVAF</li></ul>"),

    hr(),
    h2("How to Navigate this App"),
    p("After uploading a data file, you may make the following selections in the sidebar panel:"),
    HTML("<ul><li>groups (treatment, generations, etc.) to compare</li><li>test to run on loci for detecting significant difference in VAF among groups</li><li>type of linear model (with random effects or without) to use for quality checking</li></ul>"),
    br(),

    br(),
    br()
  )
)
