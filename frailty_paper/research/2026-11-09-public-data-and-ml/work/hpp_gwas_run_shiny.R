

##' @export
##' @import fastman
shiny_manhattan <- function(x, phenoName){
  ##ylims: -log10(5*10**-n) ~ n-1, i.e for 4 put 3
  fastman::fastman(x, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8), ylim = c(min(2.5, min(-log10(x$P))), max(-log10(x$P))), main = phenoName)

}


list_features <- function(domain, options_mb){
  if (domain != "Microbiome"){
    read.csv(system.file("lists",
                         paste0(domain ,"Loader", ".csv"),
                         package = "gwasInterface"))[,2]
  }
  else{
    options_mb$species
  }
}


make_paper_figures <- function(){
  qu <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.1st_Qu..glm.linear", sep = "\t")
 # dx <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.body_head_bmd.glm.linear", sep = "\t")
  #jpeg("~/Desktop/dx.jpeg", width = 960)
  #fastman(dx, chr = "X.CHROM", bp ="POS", ylim = c(0, 11), cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex = 1, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  #dev.off()
  jpeg("~/Desktop/qu.jpeg", width = 960)
  fastman(qu, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  dev.off()
}

make_cnv_manplot <- function(){
  library(qqman)
  res <- read.csv("/net/mraid08/export/jasmine/zach/cnv/ea1c_test_formatted.csv")
  res$CHR <- as.integer(res$CHR)
  res_add_to_each_chr_start <- sapply(unique(res$CHR), function(chr){res[res$CHR == toString(as.numeric(chr) - 1), "END"][length(res[res$CHR == toString(as.numeric(chr) - 1),])]})
  res_add_to_each_chr_start[is.na(res_add_to_each_chr_start)] <- 0
  res_add_to_each_chr_start <- cumsum(res_add_to_each_chr_start)
  for (chr in names(res_add_to_each_chr_start)){
    res[res$CHR == toString(chr),"START"] <- res[res$CHR == toString(chr),"START"] + res_add_to_each_chr_start[chr]
  }
  res$BP <- (res$START + res$END)/2
  res[res$CHR == "X", "CHR"] <- 23
  res[res$CHR == "Y", "CHR"] <- 24
  res[res$CHR == "MT", "CHR"] <- 25
  res$CHR <- as.numeric(res$CHR)
  jpeg("~/Desktop/ea1c_cnv.jpeg", width = 960)
  manhattan(res[!is.na(res$P),])
  dev.off()
  }

ea1c_plot <- function(){
  library(fastman)
  ea1c <- read.csv("/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/batch0.ea1c.glm.linear", sep = "\t")
  jpeg("~/Desktop/ea1c_snps.jpeg", width = 960)
  fastman(ea1c, chr = "X.CHROM", bp ="POS", cex.axis = 1, annotatePval =  5*10**(-8)/727, snp = "ID", cex.text = 1.3, genomewideline = -log10(5*10**(-8)/727), suggestiveline = -log10(5e-8))
  dev.off()
  }

##from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
ggplot_manhattan <- function(gwas_data, theTitle = "Manhattan Plot", ymin = 4){
  colnames(gwas_data) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "AX", "TEST", "N", "BETA", "SE", "T_STAT", "P")
  data_cum <- gwas_data %>%
    group_by(CHR) %>%
    summarise(max_bp = max(BP)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR, bp_add)
  gwas_data <- gwas_data %>%
    inner_join(data_cum, by = "CHR") %>%
    mutate(bp_cum = BP + bp_add)
  axis_set <- gwas_data %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cum))

  ylim <- gwas_data %>%
    filter(P == min(P)) %>%
    mutate(ylim = abs(floor(log10(P))) + 2) %>%
    pull(ylim)

  manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P),
                                    color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0,0), limits = c(ymin, ylim)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL,
         #https://sites.google.com/view/stuck-in-the-shallow-end/home/generate-manhattan-plots-with-ggplot2-and-ggrastr
         title = theTitle) +
    ylab(expression(paste(-log[10],"(", italic(P), ")"))) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
  manhplot <- manhplot  + theme(title = element_text(size = 0), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 17), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15)) + labs(x = "CHR")
  manhplot
}

read_clumped <- function(fname){
  read_in <- do.call(rbind, lapply(strsplit(readLines(fname), "\\s+|\\t+|\\s+\\t+|\\t+\\s+"), function(x){as.data.frame(t(x))}))
  read_in <- read_in[, 1:6]
  colnames(read_in) <- read_in[1,]
  read_in <- read_in[-1,] ##Drop the column names as the first entry
  read_in$CHR <- as.numeric(read_in$CHR)
  read_in$"F" <- as.numeric(read_in$"F")
  read_in$BP <- as.numeric(read_in$BP)
  read_in$P <- as.numeric(read_in$P)
  read_in
}

format_gwas <- function(thephenoname, domain, options_mb){
  if (domain == "Microbiome"){
    thephenoname <- paste0(options_mb[options_mb$species == thephenoname, "ok"], ".glm.linear")
  }
  else{
  }
  thegwas <- readRDS(system.file(paste0("full_results_rds/", thephenoname, ".Rds"),
                                 package = "gwasInterface"))
  colnames(thegwas) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "AX", "TEST", "N", "BETA", "SE", "T_STAT", "P")
  thegwas <- thegwas[order(thegwas$P),]
  data.table(thegwas)
}

domains_list <- c("CGM", "Ultrasound", "HormonalStatus", "DEXA", "RetinaScan", "Microbiome")

##' Run the GWASInterface Shiny
##'
##' run_shiny() is an example of a single-file Shiny app in that it defines a UI object and a server method to handle that object. A benefit is that Roxygen import tags only need to be called once to become available to the entire shiny, and it's easy to run the shiny as well. After both the ui object and server function are defined, browserManager is called, which sets the viewing environment that the Shiny runs in. Anywhere with a tag$html call is either a CSS or HTML tag to change certain visual aspects of the shiny.  In addition, many ui elements are rendered in the server function and passed to ui with renderUI. RenderHTML is another wrapper for this as well. This lets us make UI elements which depend on input from the UI itself.
##'
##' @import shiny
##' @import ggplot2
##' @import shinythemes
##' @import dplyr
##' @import tibble
##' @import data.table
##' @return NULL
##' @export
run_shiny <- function(useBrowser = TRUE, usingOnline = FALSE) {
  ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
                  h1(id = "heading", "Eran Segal Human Phenotype Project Interactive GWAS Results Interface"),
                  ##Colour the top and bottom of the page appropriately.
                  tags$style(HTML("#heading {background-color: #0078a4; color: white !important;}")),
                  tags$style(HTML("#sourcelink {background-color: #0078a4; color: white !important;}")),
                  #Change the colours of text on the tab selectors to be blue.
                  tags$style(HTML("
                                  .tabbable > .nav > li > a[data-value='summarystatspanel'] {color: black}
                                  .tabbable > .nav > li > a[data-value='gencorr'] {color: black}
                                  ")),
#                  shinyWidgets::setBackgroundColor(color = "#e6ebed"),
                  sidebarLayout(
                    sidebarPanel(id = "sidebar", width = 4,
                                 ## Only show the selector to input parameters if that's selected.
                                   tabsetPanel(
                                     id = "tabs",
                                     selected = "summarystatspanel",
                                     tabPanel(
                                       uiOutput("pheno1"),
                                     ),
                                     tabPanel(
                                       title = "Feature Selection",
                                       value = "summarystatspanel",
                                       selectInput("domain",
                                                   label= "Feature Domain",
                                                   choices = domains_list),
                                       uiOutput("prs_heatmap_title_maybe"),
                                       uiOutput("maybe_heatmap"),
                                       uiOutput("featureSelect"),
                                       uiOutput("gwasDownload"),
                                       ))),
                    mainPanel(
                        h3(id = "t1", "Manhattan Plot of GWAS Results"),
                        plotOutput("plot", width  = "100%"),
                        ##Colour the top and bottom of the page appropriately.
                        htmlOutput("prs_assoc_title_maybe"),
                        dataTableOutput("prs_assoc_maybe"),
                        h3(id = "t3", "Interactive Table of GWAS Results"),
                        dataTableOutput("table")
                    )
                  ),
                  uiOutput("sourcelink")
                  )

  #Everything else.
  server <- function(input, output, session){
    prs_table_loaded <- readRDS(system.file("prs_assoc_corrected_clinical.Rds",
                                             package = "gwasInterface"))
    prs_name_dict_loaded <- readRDS(system.file("prs_name_dict.Rds",
                                            package = "gwasInterface"))
    ##do inside the server function or the shiny will break
    options_mb <- read.csv(system.file("lists",
                                      "MicrobiomeLoader.csv",package = "gwasInterface"))



    observeEvent(input[["pheno"]], {
      output$plot <- renderPlot({
        ##Force reactive loading based on these values
        if (input$domain == "Microbiome"){
          if (is.null(input$pheno)){
            input$pheno <- list_features(input$domain)[3]
          }
          actual_pheno_fname <- paste0(options_mb[options_mb$species == input$pheno, "ok"], ".glm.linear")
          shiny_manhattan(readRDS(system.file(paste0("full_results_rds/",actual_pheno_fname, ".Rds"),
                                              package = "gwasInterface")), "Selected Microbiome Species")
        }
        else{
          shiny_manhattan(readRDS(system.file(paste0("full_results_rds/", input$pheno, ".Rds"),
                                              package = "gwasInterface")), input$pheno)
        }
        })

      output$table <- renderDataTable({
        thegwas <- format_gwas(input$pheno, input$domain, options_mb)
        })
      output$prs_assoc_maybe <- renderDataTable(
        if (input$domain == "Microbiome"){
        }
        else{
          justpheno <- as.numeric(prs_table_loaded[prs_table_loaded$Phenotype == input$pheno,]) ##first two entries are NA corresponding to the index of Loader, Phenotype
          sig_ids <- justpheno < 0.05
          sig_ids[c(1,2)] <- c(FALSE, FALSE) ##set FALSE (instead of NA) to the first two columns because they are just the index
          justpheno_sig_P <- justpheno[sig_ids]
          justpheno_sig_prs_names <- colnames(prs_table_loaded)[sig_ids]
          ##use a tibble so we can have spaces in the column names
          res_table <- tibble("Phenotype" = rep(input$pheno, length(justpheno_sig_P)), "Corrected P Value" = justpheno_sig_P, "Polygenic Risk Score (PRS)" = justpheno_sig_prs_names)
          res_table["Corresponding UK Biobank Phenotype Code for PRS"] <- sapply(res_table[["Polygenic Risk Score (PRS)"]], function(thename){prs_name_dict_loaded[prs_name_dict_loaded$name == thename,]$code})
          res_table <- data.table(res_table[order(res_table["Corrected P Value"]), c("Phenotype", "Polygenic Risk Score (PRS)", "Corresponding UK Biobank Phenotype Code for PRS", "Corrected P Value")])
        }
      )
    })

    output$gwasDownload <- renderUI({
      downloadButton("downloadData", "Download GWAS Summary Statistics", class = "dbutton")
    })
    output$maybe_heatmap <- renderUI({
      if (input$domain == "Microbiome"){
      }
      else{
        renderImage({
          list(src=system.file(paste0("prs_heatmap_imgs/", input$domain, "Loader", ".png"),
                               package = "gwasInterface"),
               width = "400px", #keep square or the image will look stretched
               height = "400px")
        }, deleteFile = FALSE)
      }
    })


    output$featureSelect <- renderUI({
      selectInput("pheno",
                  label = "Phenotype:",
                  choices = list_features(input$domain, options_mb),
                  selected = list_features(input$domain, options_mb)[2])

    })
    output$prs_assoc_title_maybe <- renderUI({
      if (input$domain == "Microbiome"){

      }
      else{
        h3(id = "t2", "Top Multi-SNP Significant PRS Associations")
      }
    })

    output$prs_heatmap_title_maybe <- renderUI({
      if (input$domain == "Microbiome"){

      }
      else{
        h4(id = "t0", "PRS Associations Heatmap (Corrected P)")
      }
    })
    ##Handle downloads
    output$downloadData <- downloadHandler(
      filename = function(){
        if (input$domain == "Microbiome") {
          paste0(options_mb[options_mb$species == input$pheno, "ok"], ".csv")
        }
        else{
          paste0(input$pheno, ".csv")
        }
      },
      content = function(file) {
      write.csv(format_gwas(input$pheno, input$domain, options_mb), file, row.names = FALSE)
      }
    )
      output$sourcelink <- renderUI({tagList(
      "Source code available at ",
      a("https://github.com/ZachLevine-11/gwasInterface", href = "https://github.com/ZachLevine-11/gwasInterface"),
      "                                Based on the paper: Genome Wide Association Studies and Polygenic Risk Score Phenome Wide Association Studies across complex phenotypes in the Human Phenotype Project"
    )})
  }
    shinyApp(ui, server)
}
