server <- function(input, output) {

  # type conversion helper function

  data.type<- function(data, to_factor = NULL, to_numeric = NULL, to_character = NULL) {

    data <- data %>%
      mutate(across(all_of(to_factor), factor)) %>%
      mutate(across(all_of(to_numeric), as.numeric)) %>%
      mutate(across(all_of(to_character), as.character))

  }

  # making reactive dataset from input file (accept .csv only)

  dat <- eventReactive(input$file1, {
    dat <- req(input$file1)
    df <- read.csv(dat$datapath)
    df <- data.type(df, to_factor = c("chrom", "sample", "individual", "group"))
    df
  })

  # dynamically populate a checkbox question with the groups that are in this particular dataset

  output$groups <- renderUI({
    req(dat())
    checkboxGroupInput("groups_test", "Select groups",
                       choices = unique(dat()$group),
                       selected = unique(dat()$group),
                       inline = FALSE)
  })

  # create a subsetted version of the data with only the user-selected groups



  groups_selected <- eventReactive(input$run_analysis,{
    x <- isolate(input$groups_test)
    x
  })

  method_selected <- eventReactive(input$run_analysis,{
    x <- isolate(input$method_test)
    x
  })

  digits_round <- eventReactive(input$run_analysis,{
    x <- isolate(input$digits)
    x
  })

  quality_check <- eventReactive(input$run_analysis, {
    x <- isolate(input$qc_method)
    x
  })

  table_options <- eventReactive(input$run_analysis, {
    x <- isolate(input$table_display)
    x
  })

  dat_subset <- eventReactive(input$run_analysis, {
    df <- req(dat())
    df <- df[df$group %in% groups_selected(),]
    df
  })
########## CREATION OF RESULT DF ###############


  result <- eventReactive(input$run_analysis, {

    # adjust for difference in display test names and arguments to function
    method <- ""
    if (method_selected() == "Beta-binomial") {
      method <- "betabinom"
    } else if (method_selected() == "Fisher's exact") {
      method <- "Exact"
    } else if (method_selected() == "Pearson's chi-squared") {
      method <- "Pearson"
    } else {
      method <- NULL
    }


    # run package function to generate results dataset

    result <- VAFmain(dat_subset(), method, groups_selected(), digits_round())

    # rename columns for nicer display

    colnames(result) <- c("ID",
                          as.vector(outer(c("DP", "VC", "VAF"), groups_selected(), paste, sep=".")),
                          "P value",
                          "Test",
                          "Effect size",
                          "95% CI",
                          "Overdispersion",
                          "P value (adjusted)",
                          "Sig. diff.",
                          "Sig. diff. (FDR)",
                          "Sig. change (>0.2)",
                          "Direction")

    # convert read depth columns to numeric data type to compute total read depth for stats plot
    to_numeric <- colnames(result)[grepl("DP", colnames(result)) | grepl("VC", colnames(result)) | grepl("VAF", colnames(result)) | grepl("P value", colnames(result))]
    result <- result %>%
      mutate(across(all_of(to_numeric), as.numeric)) %>%
      mutate(`Total Read Depth` = rowSums(select(.,starts_with("DP"))))

    result
  })
  ###################################################


  # Exploratory Analysis ----
  observeEvent(input$run_analysis, {

    # INTERACTIVE SCATTER PLOT 1: READ DEPTH ----

    ranges <- reactiveValues(x = NULL, y = NULL)

    output$exploratory_plot1 <- renderPlot({
      df <- isolate(dat_subset())

      base_plot <- ggplot(data = df, aes(x=reorder(Locus, dp, mean), y=dp)) +
        geom_point(aes(shape = individual,
                       colour=group),
                   size=3)+
        xlab("Loci ordered by the mean of read depth")+
        ylab("Read depth") +
        ggtitle("Scatter plot of read depth")+
        coord_cartesian(xlim = ranges$x, ylim = ranges$y)

      if (length(groups_selected()) > 6) {
        base_plot + theme(legend.position = "none")
      } else {
        base_plot
      }
    })

    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })


    # INTERACTIVE SCATTER PLOT 2: VARIANT COUNT ----

    ranges2 <- reactiveValues(x=NULL, y=NULL)

    output$exploratory_plot2 <- renderPlot({
      df <- isolate(dat_subset())
      base_plot <- ggplot(data = df, aes(x=reorder(Locus, dp, mean), y=vc)) +
        geom_point(aes(shape = individual,
                       colour = group),
                   size=3)+
        xlab("Loci ordered by the mean of read depth")+
        ylab("Variant counts")+
        ggtitle("Scatter plot of variant count")+
        coord_cartesian(xlim = ranges2$x, ylim = ranges2$y)

      if (length(groups_selected()) > 6) {
        base_plot + theme(legend.position = "none")
      } else {
        base_plot
      }
    })

    observeEvent(input$plot2_dblclick, {
      brush <- input$plot2_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })


    # INTERACTIVE SCATTER PLOT 3: VAF ----

    ranges3 <- reactiveValues(x=NULL, y=NULL)

    output$exploratory_plot3 <- renderPlot({
      df <- isolate(dat_subset())
      base_plot <- ggplot(data = df, aes(x=reorder(Locus, dp, mean), y=vc/dp)) +
        geom_point(aes(shape = individual,
                       colour= group),
                   size=3)+
        xlab("Loci ordered by the mean of read depth")+
        ylab("VAF")+
        ggtitle("Scatter plot of VAF")+
        coord_cartesian(xlim = ranges3$x, ylim = ranges3$y)

      if (length(groups_selected()) > 6) {
        base_plot + theme(legend.position = "none")
      } else {
        base_plot
      }
    })

    observeEvent(input$plot3_dblclick, {
      brush <- input$plot3_brush
      if (!is.null(brush)) {
        ranges3$x <- c(brush$xmin, brush$xmax)
        ranges3$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges3$x <- NULL
        ranges3$y <- NULL
      }
    })
  })


  # Statistical Testing ----
observeEvent(input$run_analysis, {

  # different text to print depending on whether 2 groups have been selected (plot is displayed) or not (only table is displayed)

  output$stat_blurb <- renderText({
    paste("A total of <B>", length(unique(dat_subset()$Locus)), "</B> loci were compared independently between groups <B>", paste(groups_selected(), collapse = ", "),"</B> using <B>", method_selected(), "</B> test.")

  })

  output$stat_graphic_desc <- renderText({
    if (length(groups_selected()) == 2) {
      "The plot below shows the relationship between VAF in the selected groups for loci in this dataset. Hover over the points to see their locus IDs. Click on them or use the lasso tool in the upper right corner of the plot to select them. Selected loci will appear in the table below, which is automatically ordered by ascending p value. If no points are selected, the table displays all loci. You may download the full table at the bottom of the page."
    }
    else {
      "The table below displays loci and relevant test results. You may download the full table at the bottom of the page."
    }
  })

  # using girafe package to create interactive ggplot
    # zoom
    # hover over points and see locus ID
    # download snapshot of plot
    # lasso tool and point (de-)selection

  output$stat_plot <- renderGirafe({
    varx <- paste("VAF.", groups_selected()[1], sep = "")
    vary <- paste("VAF.", groups_selected()[2], sep = "")
    df_1 <- result() %>% select(-c("Sig. diff. (FDR)",  "Sig. change (>0.2)"))
    df_2 <- result() %>% select(-c("Sig. diff.", "Sig. change (>0.2)"))
    df_3 <- result() %>% select(-c("Sig. diff.", "Sig. diff. (FDR)"))

    stat_plot1 <- ggplot(data=df_1, aes_string(x = varx, y = vary))+
      xlab(varx)+
      ylab(vary)+
      ggtitle("Method: p value, raw")+
      geom_point_interactive(aes(tooltip = ID, data_id = ID, colour=`Sig. diff.`, size = `Total Read Depth`, alpha=0.3), show.legend=FALSE)+
      theme(text=element_text(size=30))+
      scale_size_continuous(range = c(5,20))


    stat_plot2 <- ggplot(data=df_2, aes_string(x = varx, y = vary))+
      xlab(varx)+
      theme(axis.title.y=element_blank(),
            axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            text=element_text(size=30))+
      ggtitle("Method: p value, FDR")+
      geom_point_interactive(aes(tooltip = ID, data_id = ID, colour=`Sig. diff. (FDR)`, size = `Total Read Depth`, alpha=0.3), show.legend=FALSE)+
      scale_size_continuous(range = c(5,20))

    stat_plot3 <- ggplot(data=df_3, aes_string(x = varx, y = vary))+
      xlab(varx)+
      theme(axis.title.y=element_blank(),
            axis.text.y.left = element_blank(),
            axis.ticks.y.left = element_blank(),
            legend.position = "none",
            text=element_text(size=30))+
      ggtitle("Method: change > 0.2")+
      geom_point_interactive(aes(tooltip = ID, data_id = ID, color=`Sig. change (>0.2)`, size = `Total Read Depth`, alpha=0.3))+
      scale_alpha(guide="none")+
      # scale_size(guide="none")+
      scale_size_continuous(range = c(5,20))


    # collecting all plots together in one horizontal row with functions from cowplot package

    trirow <- plot_grid(stat_plot1,stat_plot2, stat_plot3, align="h", ncol=3, rel_widths = c(1,1,1))
    legend <- get_legend(stat_plot3 + theme(legend.position = "right"))
    complete_plot <- plot_grid(trirow, legend, ncol = 2, rel_widths = c(1, 0.1))

    gx <- girafe(ggobj = complete_plot,
                 options = list(
                   opts_zoom(max=5),
                   opts_hover(css="stroke: gray;stroke-width: 7"),
                   opts_selection(css="stroke: black;stroke-width: 5"),
                   opts_resizing = FALSE
                 ),
                 width_svg = 40,
                 height_svg = 15)
    gx
  })

  # dynamically decide whether to render plot based on number of groups selected

  output$compare_methods <- renderUI ({
    req(input$run_analysis)
    if (length(groups_selected()) == 2) {
      tagList(
        br(),
        h4("Graphical comparison of selection methods"),
        h6("Point size corresponds to read depth"),
        girafeOutput("stat_plot", height = "100%", width = "100%")
      )
    }
  })


  output$top_loci <- DT::renderDataTable({

    # converting some binary variables to factor for DT column filter functionality
    fordisplay <- result() %>%
      mutate(across(all_of(c("Overdispersion" ,"Sig. diff.", "Sig. diff. (FDR)", "Sig. change (>0.2)", "Direction")), factor))

    # table display options
    if (!("Raw p value" %in% table_options())){
      fordisplay <- fordisplay %>%
        select(-c("P value", "Sig. diff."))
    }
    if (!("Adjusted p value" %in% table_options())){
      fordisplay <- fordisplay %>%
        select(-c("P value (adjusted)", "Sig. diff. (FDR)"))
    }
    if (!("Effect size > 0.2" %in% table_options()) | length(groups_selected()) > 2){
      fordisplay <- fordisplay %>%
        select(-c("Sig. change (>0.2)"))
    }

    # loci in table correspond to plot selections IF plot is displayed
    if (length(groups_selected()) ==2) {
      if (!is.null(input$stat_plot_selected)) {
        fordisplay <- fordisplay[fordisplay$ID %in% input$stat_plot_selected,]
      }
    }
    else {
      fordisplay <- fordisplay %>%
        select(-c("Effect size", "95% CI", "Direction"))
    }


    # only display test column if the test selected was "default" so user knows what was run for each locus
    if (method_selected() != "Default") {
      fordisplay <- fordisplay %>%
        select(-c("Test"))
    }

    # fordisplay <- fordisplay %>%
    #   select(-c(starts_with("VC"), starts_with("DP"), "Test", "P value (adjusted)", "Total Read Depth"))

    # have the table be ordered by ascending p value, initially
    top_dt <- datatable(fordisplay,
                        rownames = FALSE,
                        filter = "top",
                        options = list(
                                       scrollX = TRUE))
    top_dt

  })

  output$download_toptable <- downloadHandler(

    filename = paste("locitable_", Sys.Date(), ".csv", sep = ""),
    content = function(file){
      toptable <- result()[order(result()$`P value`),]
      write.csv(toptable, file, row.names = FALSE)
    }
  )
})


  # Quality of Sample ----

  # plot

    output$QCplot <- renderPlot({
      df <- req(dat())
      df$vaf <- df$vc/df$dp
      ggplot(df, aes_string(x="vaf", color="sample", fill="sample")) +
        geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
        labs(x="VAF", y = "Density", title = "Distribution of VAF")+
        facet_grid(rows = vars(groups_selected()))+
        theme_classic()
    })

  # linear or linear mixed effects model to test for significant biological variability
  test_res <- eventReactive(input$run_analysis, {
    df <- req(dat())
    test <- ""
    if(quality_check()=="Linear mixed-effects"){
      test <- "lmer"
    } else if(quality_check()=="Linear"){
      test <- "lm"
    }
    test_res <- QCchecking(df, test)
    test_res
  })

  # displaying test statistic and p value from above test
    output$QC_chisquare <- renderTable({
      if (quality_check() =="Linear") {
        table_results <- matrix(c(test_res()[2,5], test_res()[2,6]), ncol = 2)
        colnames(table_results) <- c("F statistic", "p value")
      } else {
        table_results <- matrix(c(test_res()[2,6], test_res()[2,8]), ncol = 2)
        colnames(table_results) <- c("Chi-square statistic", "p value")
      }
      table_results
    })

    # print one line about significance in summary based on the p value
    output$concl_bio_var <- renderText({
      if (quality_check() == "Linear") {
        if (test_res()[2,6] >= 0.05) {
          paste("The biological variability of the experiment mice with respect to VAF is ", strong("not significant."))
        } else {
          paste("The biological variability of the experiment mice with respect to VAF is ", strong("significant"), ". Increasing sample size or incorporating important biological covariates into analysis to increase power and accuracy is recommended.")
        }
       } else if (quality_check() == "Linear mixed-effects") {
         if (test_res()[2,8] >= 0.05) {
           paste("The biological variability of the experiment mice with respect to VAF is ", strong("not significant."))
         } else {
           paste("The biological variability of the experiment mice with respect to VAF is ", strong("significant"), ". Increasing sample size or incorporating important biological covariates into analysis to increase power and accuracy is recommended.")
         }
       }
    })

}
