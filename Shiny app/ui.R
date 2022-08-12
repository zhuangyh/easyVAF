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
  title = "VAF Analysis",
  theme = shinytheme("flatly"),

  analysis_page <- tabPanel(
    title = "Analysis",
    titlePanel("Analysis"),
    sidebarLayout(
      sidebarPanel(width = 3,
        fileInput("file1","Upload a data file", accept=".csv", width = "100%"),
        uiOutput("groups", width = "100%"),
        selectInput("method_test", label = "Select test", multiple = FALSE,
                    choices = c("Default", "Fisher's exact", "Pearson's chi-squared", "Beta-binomial"),
                    selected = "Default",
                    width = "100%"),
        selectInput("qc_method", label = "Select model for quality checking",
                    choices = c("Linear", "Linear mixed-effects"),
                    width = "100%"),

        numericInput("digits", label ="Choose number of digits for rounding",
                    min = 1, max = 5, value = 3, step =1,
                    width="100%"),

        selectInput("table_display", label = "Display significant difference based on:", multiple = TRUE,
                    choices = c("Raw p value", "Adjusted p value", "Effect size > 0.2")),

        actionButton("run_analysis", label = "Run analysis", width="100%", icon("play")),
        br()
        ),

      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Exploratory Plots",
            br(),
            p("Displayed on this page are exploratory analysis plots to examine the distribution and variability of ", em("read depths (DP), variant counts (VC), and variant allele frequencies (VAF)"), " of loci."),
            p("Use your mouse to select an area of the plot and double click the selection to zoom in. Double-click again anywhere on the plot to return to the default view."),
            p("Please note that the color and shape palettes for the plots will only be able to handle 6 variable levels or fewer. If your dataset has more than 6 samples or more than 6 groups, distinctions between all categories may not be reflected in the plot and legend."),
            br(),
            plotOutput("exploratory_plot1",
                       height = "700px",
                       width = "100%",
                       dblclick = "plot1_dblclick",
                       brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)),
            plotOutput("exploratory_plot2",
                       height = "700px",
                       width = "100%",
                       dblclick = "plot2_dblclick",
                       brush = brushOpts(id = "plot2_brush", resetOnNew = TRUE)),
            plotOutput("exploratory_plot3",
                       height = "700px",
                       width = "100%",
                       dblclick = "plot3_dblclick",
                       brush = brushOpts(id = "plot3_brush", resetOnNew = TRUE)),
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

            p("Displayed below are the results of a likelihood ratio test to assess the general variability of VAF among the experimental subjects. In the side panel, you may choose to run a regular linear model OR a linear mixed model, which adjusts for the non-independence of VAF within the same chromosome."),
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
            p("This tab displays relevant test results according to your selections from the sidebar panel."),

            htmlOutput("stat_blurb"),
            br(),
            textOutput("stat_graphic_desc"),

            uiOutput("compare_methods"),

            br(),
            br()

          )

        )
      )
    )
  ),
  
  about_page <- tabPanel(
    title = "About",
    titlePanel("About this software"),
    h5(strong("Authors: "), "Junxiao Hu, Vida Alami, Yonghua Zhuang, Dexiang Gao"),
    p("This app is intended to accompany the ", em("easyVAF"), " R package to help investigators analyze variant allele frequency (VAF) in their datasets. With interactive plots and tables and downloadable results, this app is a no-code, user-friendly option for applying ", em("easyVAF"), " functions to loci data."),
    p("Documentation and code for both the package and this Shiny app are available on Github at ", a(href="https://github.com/zhuangyh/easyVAF", "github.com/zhuangyh/easyVAF"), ". Please report any bugs or suggestions to the contact listed there."),
    hr(),

    h2("Data requirements"),
    p("Please upload your data as a comma separated values ", strong("(.csv)"), " file and make sure that it contains the following columns with these exact names:"),
    HTML("<ul><li><b>Locus</b>: locus ID</li><li><b>chrom</b>: chromosome information (for linkage disequilibrium adjustment in QC test, if desired)</li><li><b>vc</b>: variant count</li><li><b>dp</b>: read depth</li><li><b>sample</b></li></ul>"),

    br(),
    br()
  )
)
