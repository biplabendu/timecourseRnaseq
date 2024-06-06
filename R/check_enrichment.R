#' Title
#'
#' @param geneset A character vector of gene names of interest that you want to run functional enrichment on.
#' @param what Which functional annotation to run enrichment for? Select from "GOs", "pfams", "signalP", "TMHMM" or user-defined column name. The default is "GOs".
#' @param function.dir Specify the path/to/the/directory that contains the folder "functions". The default is "."
#' @param path_to_annot Specify the path/to/csv/file that contains the functional annotations for each gene. The file should contain at least two columns, one containing the gene names ("gene_name") and one containing the functional annotations ("GOs", "pfams", "signalP", "TMHMM").
#' @param sep Specify the separator used to put together all (GO/PFAM) annotation terms for a given gene in a single column (usually ";" or "; "). The default is ";"
#' @param gene_col What is the name of the column that contains the gene names in the user-specified annotation file? By default, the program expects your gene names to be in a column called "gene_name".
#' @param org Organisms for which annotation files are included in the package. Choose between "cflo" (camponotus floridanus), "ophio_cflo" (ophiocordyceps camponoti-floridani), "ophio_kim" (ophiocordyceps kimflemingiae), and "beau" (beauveria bassiana)
#' @param bg A character vector of gene names that make up the background universe against which enrichment is run for geneset. The default uses all genes present in the species' genome.
#' @param FDR The threshold false discovery rate to infer significance. Default is 5.
#' @param atleast Run enrichments only for annotations that are present in at least "atleast" (default is 5) of all background genes.
#' @param verbose Default is FALSE. If set to TRUE, prints summary statistics.
#' @param plot Default is TRUE. Plots the results of the enrichment analyses.
#' @param plot_file Save the generated plots to a specified path/and/filename.png. The default is "."
#' @param plot_height Height of the generated plot png in pixels. The default is 800.
#' @param plot_width Width of the generated plot png in pixels. The default is 800.
#' @param n_trunc If plot is TRUE, n_trunc specifies the number of characters in the description of the annotation that is printed on the plot.
#' @param clean If plot is TRUE, clean = F will get rid of additional labels from the plot
#' @param filter Default is TRUE. Filters to keep only significantly overrepresented terms at the specified FDR.
#' @param simple Default is TRUE. Retains only the most relevant columns for the enrichment result.
#' @param expand If set to TRUE, will return a long-formatted table of the filtered and simplified enrichment result.
#'
#' @return A table (enrichment results as a tibble) and one plot (ggplot2 object)
#' @export
#'
#' @examples
#' some.genes <- c("BBA_00100", "BBA_01000", "BBA_10000", "BBA_10001", "BBA_10002", "BBA_10003", "BBA_10004")
#' check_enrichment(some.genes, org="beau", plot=F)
check_enrichment <- function(geneset,
                             what="GOs",
                             function.dir = ".",
                             path_to_annot = "not_provided",
                             sep = ";",
                             gene_col = "gene_name",
                             org = "ophio_cflo",
                             bg = "all",
                             FDR = 5,
                             atleast = 5,
                             verbose=F,
                             plot=T,
                             plot_file=".",
                             plot_height=800,
                             plot_width=800,
                             n_trunc=40,
                             clean="no",
                             filter=T,
                             simple=T,
                             expand=F) {
  
  #-#-#-##-#-#-##-#-#-##-#-#-#
  ### STOP IF NOT, CHECKPOINT-1
  #-#-#-##-#-#-##-#-#-##-#-#-#
  
  # see here for more info: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/stopifnot
  # make sure the class of each variable is correct
  
  # save the input list of genes for enrichment test
  genes <- as.character(geneset)
  
  ## Load the required libraries
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(conflicted)
  ## set conflict preference
  conflict_prefer("filter","dplyr", quiet = T)
  conflict_prefer("select","dplyr", quiet = T)
  
  
  
  ## Define function to stop quietly
  ## From https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
  stop_quietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  
  # EDIT AlBaars 21Aug23: Changed all functions that used gene_name to use the
  # gene_col variable as to enable custom column selection in the included annotation
  # as well. The !!gene_col and .[[gene_col]] is to make sure that R/dplyr knows to
  # use the string inside the variable gene_col rather than interpreting gene_col
  # as a string (the default behaviour).
  
  if (path_to_annot=="not_provided") {
    
    ## load the selected annotation file
    if (org=="ophio_cflo"){
      
      print("Loading annotation file for Ophiocordyceps camponoti-floridani")
      # load("./data/ophio_cflo_annots.rda")
      all_genes <- ophio_cflo_annots
      
      # define the separator
      separator = "; "
      # check if annotation file is correct
      # EDIT AlBaars 21Aug23: changed gene name checking to any(x %in% y)
      if (!any(genes %in% all_genes[[gene_col]])) {
        print("Genes names do not match!")
        print("Check if provided gene names are the same as in the annotation file, or check if org is incorrect.")
        stop()
      }
      print("Done.")
      
    } else if (org=="ophio_kim"){
      
      print("Loading annotation file for Ophiocordyceps kimflemingae")
      # load("./data/ophio_kim_annots.rda")
      all_genes <- ophio_kim_annots
      
      # define the separator
      separator = ";"
      # check if annotation file is correct
      if (!any(genes %in% all_genes[[gene_col]])) {
        print("Genes names do not match!")
        print("Check if provided gene names are the same as in the annotation file, or check if org is incorrect.")
        stop()
      }
      print("Done.")
      
    } else if (org=="cflo"){
      
      print("Loading annotation file for Camponotus floridanus")
      # load("./data/cflo_annots.rda")
      all_genes <- cflo_annots
      
      # define the separator
      separator = "; "
      
      # check if annotation file is correct
      if (!any(genes %in% all_genes[[gene_col]])) {
        print("Genes names do not match!")
        print("Check if provided gene names are the same as in the annotation file, or check if org is incorrect.")
        stop()
      }
      print("Done.")
      
    } else if (org=="beau"){
      
      print("Loading annotation file for Beauveria bassiana")
      # load("./data/beau_annots.rda")
      all_genes <- beau_annots
      
      # define the separator
      separator = "; "
      # check if annotation file is correct
      if (!any(genes %in% all_genes[[gene_col]])) {
        print("Genes names do not match!")
        print("Check if provided gene names are the same as in the annotation file, or check if org is incorrect.")
        stop()
      }
      print("Done.")
      
    } else if (org=="pbar"){
      
      print("Loading annotation file for Pogonomyrmex barbatus")
      # load("./data/beau_annots.rda")
      all_genes <- pbar_annots
      
      # define the separator
      separator = "; "
      # check if annotation file is correct
      if (!any(genes %in% all_genes[[gene_col]])) {
        print("Genes names do not match!")
        print("Check if provided gene names are the same as in the annotation file, or check if org is incorrect.")
        stop()
      }
      
      print("Done.")
      
    } else {
      
      writeLines("Invalid option for argument org.")
      writeLines("Available organisms (org):")
      writeLines("Camponotus floridanus (org=cflo), Ophiocordyceps camponoti-floridani (org=ophio_cflo),")
      writeLines("Ophiocordyceps kimflemingiae (org=ophio_kim), Beauveria bassiana (org=beau)")
      writeLines("Pogonomyrmex barbatus (org=pbar)")
      writeLines("")
      writeLines("Use path_to_annot = your/path/to/annotation/csv/file to provide your own annotation file")
      writeLines("Remember to specify the separator used (sep = ?).")
      writeLines("")
      stop()
    }
    
  } else {
    print("Loading user provided annotation file...")
    all_genes <- read.csv(paste0(path_to_annot),
                          header=T, stringsAsFactors = F, na.strings = c(NA,""," ")) %>% as_tibble()
    # define the separator
    separator = sep
    
    # assign genes to all_genes variable
    all_genes <- all_genes %>%
      select(!!gene_col, everything())
      
    # check if annotation file is correct
    if (!any(genes %in% all_genes[[gene_col]])) {
      print("Genes names do not match!")
      print("Check if provided gene names are the same as in the annotation file, or if `gene_col` is incorrect.")
      stop()
    }
    print("Done.")
  }
  
  #-#-#-##-#-#-##-#-#-##-#-#-#
  ### FORMAT GENE ANNOTATIONS
  #-#-#-##-#-#-##-#-#-##-#-#-#
  
  ## make the flattened gene x annotation file,
  
  # If a gene_Id has multiple GO_terms, we want them in multiple row, instead of one
  # Select only the columns that we need
  all_genes_annots <- all_genes[,c({{gene_col}},{{what}})]
  # Let's replace the NAs in GOs and pfams with "no_annot"
  all_genes_annots[is.na(all_genes_annots)] <- "no_annot"
  
  # Let's flatten the files
  all_genes_annots <-
    all_genes_annots %>%
    # EDIT on 29Sep21: removing trailing and leading whitespace from the columns
    mutate_if(is.character, str_trim) %>%
    mutate(annot_split = str_split(all_genes_annots[[what]], separator)) %>%
    # EDIT by AlBaars on 18Aug23: unnest() syntax changed, added 'cols' argument
    unnest(cols=c(!!gene_col, annot_split)) %>%
    #dplyr::select(-GO) %>%
    separate(col=annot_split, into=c("annot","annot_desc"), sep = "([\\|])", extra = "drop") %>%
    select(!!gene_col, annot, annot_desc) %>%
    unique() %>%
    mutate(annot_desc=gsub(";$", "", annot_desc))
  
  ## all genes without an annotation (replace NAs in annot_desc to "no_annot")
  all_genes_annots[is.na(all_genes_annots)] <- "no_desc"
  
  if (org=="ophio_kim") {
    all_genes_annots %<>%
      separate(annot, c("annot","extra"), sep = "([\\.])", extra = "drop") %>%
      select(-extra)
  }
  
  #-#-#-##-#-#-##-#-#-##-#-#-#
  ### BACKGROUND GENESET
  #-#-#-##-#-#-##-#-#-##-#-#-#
  
  ## define the background geneset to run enrichment against
  if(bg == "all") {
    # save the background data frame with gene_name, GO term, and GO description in an object
    background <- all_genes_annots %>%
      arrange(!!gene_col)
  } else (
    background <- all_genes_annots %>%
      # filter and keep user specified background geneset
      filter(.[[gene_col]] %in% as.character(bg)) %>%
      arrange(!!gene_col)
  )
  
  ## Need a GO term to GO description file
  annot_to_desc <- dplyr::distinct(as.data.frame(all_genes_annots[-1]))
  
  
  #-#-#-##-#-#-##-#-#-##-#-#-#
  ### PREP FOR ENRICHMENT
  #-#-#-##-#-#-##-#-#-##-#-#-#
  
  ## Enrichment to be tested for all annotation terms that are:
  ## 1. Present in the test geneset
  ## 2. terms annotated for ≥ (at least) 5 number of genes
  annot_terms <-
    background %>%
    # Keep only the genes in my test geneset
    filter(.[[gene_col]] %in% genes) %>%
    group_by(annot) %>%
    summarize(num_genes = n()) %>%
    arrange(num_genes) %>%
    # Keep only the annotation terms that are annotated in at least 5 genes
    filter(num_genes >= atleast) %>%
    pull(annot) %>%
    .[. != ""]
  
  annot_terms_background <-
    background %>%
    group_by(annot) %>%
    summarize(num_genes=n()) %>%
    filter(num_genes >= atleast) %>%
    pull(annot) %>%
    .[. != ""]
  
  # Let's get all the genes in our geneset for each annotation term
  # Test geneset dataframe
  df.test <- background %>% filter(.[[gene_col]] %in% genes)
  annot_to_genes.test <- aggregate( .~ annot, df.test, function(x) toString(unique(x))) %>% as_tibble()
  
  ## Make an empty list that can save your results for each annotation term
  #df.list <- list()
  ## EDIT AlBaars 22Aug23: changed to rbinding dataframes as this implementation
  ## returns an empty dataframe
  df.hypertest <- data.frame(annot_term = character(),
                             annot_desc = character(),
                             sam_freq = numeric(),
                             back_freq = numeric(),
                             n_annot_test = numeric(),  ## x
                             n_test_total = numeric(),
                             n_test_annotated = numeric(), ## k
                             n_annot_bg = numeric(), ## m
                             n_not_annot_bg = numeric(), ## n
                             pVal = numeric())
  
  #-#-#-##-#-#-##-#-#-##-#-#-#
  ### SUMMARY STATS FOR DATA
  #-#-#-##-#-#-##-#-#-##-#-#-#
  
  if (verbose==T) {
    ## Print the summary stats for the enrichment test
    writeLines(paste0("Number of genes in background geneset: ", background %>% distinct(.[[gene_col]]) %>% nrow()))
    writeLines(paste0("Number of genes in the test set: ", length(genes)))
    writeLines("--------------------------------")
    writeLines(paste0("Number of ", what, " terms in background geneset: ", background %>% distinct(annot) %>% nrow()))
    writeLines(paste0("Number of ", what, " terms (at least ", atleast, " genes) in background geneset: ", background %>%
                        group_by(annot) %>% summarise(num_genes = n()) %>% filter(num_genes >= atleast) %>% nrow()))
    writeLines(paste0("Number of ", what, " terms (at least ", atleast, " genes) in test set: ",length(annot_terms)))
  }
  
  
  # EDIT by AlBaars on 18Aug023: added all() and any() to the if statements.
  # Without them, the script crashes as multiple TRUE/FALSE are submitted
  # where it expects only one.
  
  if (length(annot_terms) == 0 | all(annot_terms=="no_annot")) {
    return(paste0("There are no ", what, " terms to test enrichment for."))
    
  } else if (length(annot_terms) >= 1 & any(annot_terms!="no_annot")) {
    print("Testing for enrichment...")
    
    # Test the enrichment for each of the annot terms
    for (i in 1:length(annot_terms)) {
      
      #-#-#-##-#-#-##-#-#-##-#-#-#
      ### RUN HYPERGEOMETRIC TEST
      #-#-#-##-#-#-##-#-#-##-#-#-#
      
      ### For understanding the rationale behind the setup for the
      ### hypergeometric test, I would recommend reading the following:
      ### http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
      ### EDIT AlBaars 22Aug23: link ^ no longer works. Attaching new link:
      ### https://onlinestatbook.com/2/probability/hypergeometric.html
      
      # get the annotation term to be tested for enrichment
      annot.i <- annot_terms[i]
      
      # number of genes in the test set
      n_test <- genes %>% length()
      
      # number of genes of interest
      n_test_annotated <-
        background %>% 
        filter(.[[gene_col]] %in% genes) %>% 
        filter(annot %in% annot_terms_background) %>%
        pull(!!gene_col) %>% 
        unique() %>% 
        length()
      
      # Number of genes annotated with the annot term in the background gene set
      n_annot_background <-
        background %>% 
        filter(annot == annot.i) %>% 
        nrow()
      
      # Number of genes annotated with some annot term (at least 5) in the background gene set
      n_background_annotated <-
        background %>% 
        filter(annot %in% annot_terms_background) %>% 
        distinct(.[[gene_col]]) %>% 
        nrow()
      
      # Number of genes NOT annotated with the annot term in the background gene set (at least 5)
      n_not_annot_background <- n_background_annotated - n_annot_background
      
      # Number of genes annotated with the annot term in the test set
      n_annot_test <- 
        df.test %>% 
        filter(annot == annot.i) %>% 
        nrow()
      
      # define the number of possible genes with the given annotation
      x <- min(n_test, n_annot_background)
      
      #-#-#-##-#-#-##-#-#-##-#-#-#
      ### OBTAIN PROBABALITY
      #-#-#-##-#-#-##-#-#-##-#-#-#
      
      # calculates the total probability of obtaining at least the observed
      # number of genes annotated with the given term (GOs, pfams)
      #
      # i.e. summation of probabilities of all overlaps equal to or
      #       higher than the one observed in our test set
      #
      
      pval <-
        phyper(q=n_annot_test:x,
               m=n_annot_background,
               n=n_not_annot_background,
               k=n_test_annotated,
               lower.tail=F) %>%
        sum()
      
      #-#-#-##-#-#-##-#-#-##-#-#-#
      ### SAVE TEST RESULTS
      #-#-#-##-#-#-##-#-#-##-#-#-#
      
      df.tmp <- data.frame(annot_term = annot.i,
                                 annot_desc = annot_to_desc %>% filter(annot==annot.i) %>% pull(annot_desc),
                                 sam_freq = round(n_annot_test/n_test_annotated, 3),
                                 back_freq = round(n_annot_background/n_background_annotated, 3),
                                 n_annot_test = n_annot_test,  ## x
                                 n_test_total = n_test,
                                 n_test_annotated = n_test_annotated, ## k
                                 n_annot_bg = n_annot_background, ## m
                                 n_not_annot_bg = n_not_annot_background, ## n
                                 pVal = pval)
      df.hypertest <- rbind(df.hypertest, df.tmp)
      
    }
    
    #-#-#-##-#-#-##-#-#-##-#-#-#
    ### MAKE OUTPUT FILE
    #-#-#-##-#-#-##-#-#-##-#-#-#

    ## Make the result table:
    df.enriched <-
      df.hypertest %>%
      arrange(pVal) %>%
      mutate(adj_pVal = round(p.adjust(pVal, "BH"),5)) %>%
      # keeps only the annot terms that are found in the test set
      filter(n_annot_test != 0) %>%
      mutate(over_under = ifelse(sam_freq > back_freq, "over", "under")) %>%
      dplyr::select(annot_term, annot_desc, over_under, adj_pVal, everything()) %>%
      arrange(over_under, adj_pVal) %>%
      left_join(annot_to_genes.test[1:2], by=c("annot_term"="annot")) %>% as_tibble()
    
    ## Check if any significant enrichment has been found
    if (nrow(filter(df.enriched, adj_pVal < (FDR/100))) == 0) {
      print("No significant enrichment found.")
      stop_quietly()
    }
    
    #-#-#-##-#-#-##-#-#-##-#-#-#
    ### PLOT THE RESULTS
    #-#-#-##-#-#-##-#-#-##-#-#-#
    
    if (plot==T) {
      
      if(nrow(df.enriched)!=0) {
        
        ## build the custom theme
        theme_Publication <- function(base_size=14, base_family="Helvetica") {
          library(grid)
          library(ggthemes)
          (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.2), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(),
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.key.size= unit(0.4, "cm"),
                    legend.margin = unit(0, "cm"),
                    legend.title = element_text(face="italic"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
            ))
          
        }
        
        
        #Save the data to an object
        df <- df.enriched
        
        Title <- paste0("Enriched ", org, " ", what)
        subTitle <- paste0("FDR = ",FDR,"%")
        
        #Format the dataframe
        df <-
          df %>%
          # only keep the over-enriched terms
          filter(over_under == "over") %>%
          # remove the NA term
          filter(annot_desc != "no_annot") %>%
          # set FDR at 1%
          filter(adj_pVal < (FDR/100)) %>%
          # calculate an increasing score for pvalues (higher score == lower pvalue)
          mutate(score = -log(adj_pVal)) %>%
          # make the y-value
          mutate(percent_annot_test = round((n_annot_test/n_annot_bg*100),2))
        
        if (nrow(df)>0) {
          
          if (nrow(df)<50) {
            ## make the plot
            goplot <-
              ggplot(df) +
              # set overall appearance of the plot
              theme_Publication() +
              # Define the dependent and independent variables
              aes(x = reorder(annot_desc, score), y = percent_annot_test) +
              # From the defined variables, create a vertical bar chart
              geom_col(position = "dodge", alpha=0.4, size = 1, fill="#143740") +
              # Set main and axis titles
              labs(
                title = Title,
                subtitle = subTitle
                # caption = paste0("FDR = ",FDR,"%")
              ) +
              ylab("percent of genes") +
              # ylim(c(0,min(max(percent_annot_test)+2,100))) +
              
              # Shorten very long labels (GO descriptions)
              scale_x_discrete(label = function(x) stringr::str_trunc(x, n_trunc)) +
              
              # Flip the x and y axes
              coord_flip() +
              
              theme(legend.position = "bottom") +
              theme(strip.background = element_blank(), strip.text = element_blank(), # get rid of facet grid labels
                    plot.title = element_text(hjust = 0.5),
                    axis.line.y = element_line(colour = "transparent",
                                               size=1),
                    legend.title = element_blank(),
                    # legend.position = "None",
                    # plot.caption = element_text(hjust=0.5),
                    plot.subtitle = element_text(hjust=0.5),
                    axis.title.y = element_blank())
            # guides(fill = guide_legend(title = "Legend Title",
            #                            override.aes = aes(label = "")))
            
            if (clean == "yes") {
              goplot <- goplot
            }
            
            else if (clean != "yes") {
              goplot <-
                goplot +
                # add annotations
                ggrepel::geom_label_repel(aes(label =  paste0(" of ",n_annot_bg," ")),
                                          fill = "transparent",
                                          color = 'black',
                                          size = 5,
                                          direction = "x",
                                          # ylim=c(10,max()),
                                          point.padding = 0.25,
                                          label.padding = 0.25,
                                          segment.color = 'transparent',
                                          # get rid of the outline for the label
                                          label.size = NA)
            }
            
          } else {
            # Enriched terms word-cloud
            # (borrowed from: https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a)
            
            # load libraries
            library(tm)
            library(wordcloud)
            library(RColorBrewer)
            library(wordcloud2)
            
            # get text as a character vector
            text <- df %>% pull(annot_desc)
            # load your text data as a corpus
            docs <- Corpus(VectorSource(text)) # requires library "tm"
            # clean text (necessary?)
            docs <- docs %>%
              tm_map(removeNumbers) %>%
              tm_map(removePunctuation) %>%
              tm_map(stripWhitespace)
            docs <- tm_map(docs, content_transformer(tolower))
            docs <- tm_map(docs, removeWords, c("process", "molecular","cellular",
                                                "component", "compound", "part",
                                                "activity", "acid"
            ))
            # create document-term-matrix
            dtm <- TermDocumentMatrix(docs)
            matrix <- as.matrix(dtm)
            words <- sort(rowSums(matrix),decreasing=TRUE)
            df.wc <- data.frame(word = names(words),freq=words)
            # generate word-cloud
            goplot <- wordcloud::wordcloud(words = df.wc$word, freq = df.wc$freq, scale=c(5,0.15),
                                 min.freq = 2, max.words=200, random.order=FALSE,
                                 rot.per=0.35
                                 # colors=brewer.pal(8, "Dark2")
            # NOTE AlBaars 22Aug23: This part below breaks the whole script. I cannot
            # find col.scheme defined anywhere, which is the cause.
#                                 colors=col.scheme[[1]]
            )
            
          }
          
#         EDIT Albaars 24Aug23: added the option to save the plot to a file instead
#         rather than printing it (for non-GUI systems)
          if (plot_file == ".") {
            print(goplot)
          } else {
            png(plot_file, width=plot_width, height=plot_height, units = "px")
            print(goplot)
            dev.off()
          }
          
        }
        
      }
      
    }
    
    
    #-#-#-##-#-#-##-#-#-##-#-#-#
    ### FILTER/MODIFY OUTPUT FILE
    #-#-#-##-#-#-##-#-#-##-#-#-#
    
    if (filter == T) {
      if(nrow(df.enriched)!=0) {
        df.enriched <- df.enriched %>%
          filter(over_under == "over") %>%
          filter(adj_pVal < FDR/100)
      }
    }
    
    if (simple==T) {
      if(nrow(df.enriched)!=0) {
        df.enriched <- df.enriched %>%
          select(annot_term, annot_desc, adj_pVal, sam_freq, back_freq, n_annot_bg, !!gene_col)
      }
    }
    
    if (expand==T) {
      if(nrow(df.enriched)>=1) {
        df.enriched <-
          df.enriched %>%
          select(annot_desc,!!gene_col) %>%
          separate_rows(!!gene_col, sep=", ") %>%
          left_join(all_genes[,c(gene_col,"gene_desc")], by=gene_col) %>%
          select(!!gene_col, gene_desc, everything()) %>%
          group_by(!!gene_col,gene_desc) %>%
          summarize(annot_desc = paste(annot_desc, collapse = "; "))
      } else {
        print("No enriched terms found; can't expand.")
      }
      
    }
    
    return(df.enriched);
    
    
  }
  
  
  
  
}