### Functions for Nucleotide data processing
#' @title Import files
#'
#' @description \code{import_files} imports all files exported vie TOPPAS pipeline stored
#' at a defined input folder.
#' It also takes care of required modificiations of the input data, e.g., extraction
#' of file names and renaming column names.
#'
#' @param path_files Character.
#' @param mode Character; define either as "samples" or "cal" for calibration files
#' @param conv_list Data frame; containing conversion information, e.g., nuc_nb into nuc_id
#' @param condition_list Data frame: containing concentration level for each calibration curve, 
#' required for the mode + "calibrations". Default value - NULL.
#' @param plot_graph TRUE / FALSE for the graphical output
#'
#' @return The function returns a dataframe summarising all input information provided in the
#' single input files
#'
#' @examples import_files("input/samples/", mode = "samples")
#' @examples import_files("input/cal/", mode = "cal")
#'
#'
import_files <- function(path_files, mode, true_blank, conv_list, condition_list, plot_graph = FALSE){


  dataframe = list.files(path_files, pattern = ".unknown")

  l = list()
  temp = list()

  for (i in 1:length(dataframe)) {
    #import *.unknown file
    data = read.table(paste0(path_files, dataframe[i]), header = TRUE)
    
    #clean up header
    colnames(data)[grepl("median", colnames(data))] = "median_intensity"
    colnames(data)[grepl("nuc", colnames(data))] = "nuc_nb"

    ############## transformation for data == calibration ----------------------
    if (mode == "cal") {
      #extract: what NucleoMix
      data$calcurve = unlist(strsplit(dataframe[i], split = "_", fixed = TRUE))[4]
      
      #extract: replicate calcurve
      temp = unlist(strsplit(dataframe[i], split = "_", fixed = TRUE))[5]
      data$repl_calcurve = as.numeric(unlist(strsplit(temp, split = ".", fixed = TRUE))[1])
      
      #extract: date
      data$date = unlist(strsplit(dataframe[i], split = "_", fixed = TRUE))[1]
      
      #merge: nuc_nb matching up with nuc_id
      data_conv = merge(data, conv_list)
      
      #merge: calibration curve concentrations
      data_conv = merge(data_conv, condition_list, all.x = TRUE)
      
      #transform: transition (numeric) into factor()
      data_conv = ddply(data_conv, c("nuc_id", "transition"), transform, 
                        transition_id = paste0("t", transition))
      
      #select dataframe columns
      data_conv = data_conv[, c("nuc_id", "nuc","nuc_group" ,"transition_id", "calcurve", "level" ,"repl_calcurve", 
                                "date", "median_intensity", "cv")]
    
      
      
    }
    
    
    ############### transformation for data == samples -------------------------
    if (mode == "samples") {
      
      #extract: 
      data$file_tag = unlist(strsplit(dataframe[i], split = "_", fixed = TRUE))[3]
      data$date = unlist(strsplit(dataframe[i], split = "_", fixed = TRUE))[1]
      
      #merge: sample group annotaion
      data_conv = merge(data, conv_list, all.x = TRUE)
      
      #merge: add annotation
      data_conv = merge(data_conv, condition_list, all.x = TRUE)
      
      #transform: transition (numeric) into factor()
      data_conv = ddply(data_conv, c("nuc_id", "transition"), transform, 
                        transition_id = paste0("t", transition))
      
      #select dataframe columns
      data_conv = data_conv[, c("date","file_tag", "sample", "nuc_id", "nuc", "transition_id",  
                                "median_intensity", "cv")]
      
    }
    
    #check matching annotation file
    if (nrow(data_conv) == 0) {
      message("ERROR: Check annotation file. Empty dataframe created!")
    }
    
    #write list and transform into a dataframe
    l[[i]] = data_conv
    
  }

  data_export = do.call(rbind, l)
  
  #export into .csv file
  temp_extr = unlist(strsplit(path_files, split = "/", fixed = TRUE))[2]
  write.csv(data_export, paste0("output/", temp_extr, "_export.csv"), row.names = FALSE)
  message("Done - input data exported as .csv-file: ",temp_extr)
  
  
    
  
  if (plot_graph == TRUE & temp_extr == "calibrations") {
    
    
      
      for (var in unique(data_export$nuc_group)) {
        print(
          ggplot(subset(data_export, data_export$nuc_group == var & data_export$calcurve != "TrueBlank"), 
                 aes(x = level, y = median_intensity, colour = date)) + 
            geom_point() + 
            facet_grid(nuc~transition_id, scales = "free_y") + 
            theme_bw() + 
            theme(axis.text.x = element_text(size = 6),
                  axis.text.y = element_text(size = 8)) +
            ggtitle(paste0("Nucleotide (Calibration): ", var)) +
            theme(strip.text = element_text(size = 8), 
                  strip.background = element_rect(fill = "white")) +
            scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                               breaks = NULL) +
            scale_color_manual(values = c("dodgerblue3", "black", "red"))
        )
    }
  
    }
  }
  
  if (plot_graph == TRUE & temp_extr == "samples") {
    
    for (var in unique(data_export$nuc)) {
      print(
        ggplot(subset(data_export, data_export$nuc == var),
               aes(sample, median_intensity, fill = date)) +
          geom_boxplot() +
          facet_wrap(~transition_id) + 
          theme_bw() + 
          coord_flip() + 
          theme(axis.text.x = element_text(size = 6),
                axis.text.y = element_text(size = 8)) +
          ggtitle(paste0("Nucleotide (samples): ", var)) +
          theme(strip.text = element_text(size = 8), 
                strip.background = element_rect(fill = "white")) +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                             breaks = NULL) +
          scale_color_manual(values = c("dodgerblue3", "black", "red"))
      )
    }
  }
  
  return(data_export)
}



summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
  #' @title Create statistics
  #'
  #' @description This functions creates the statistics for each nucleotide including
  #' the mean and standard deviation.
  #' 
  #' @param data
  #' @param groupvars
  #' @param measurevar
  #' @param conf.interval
  #' @param na.rm
  #'
  #' @return Question
  #'
  #'
  #'
  
  
  require(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function(x, na.rm = FALSE){
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop = .drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm = na.rm),
                     mean = mean(xx[[col]], na.rm = na.rm),
                     sd   = sd(xx[[col]], na.rm = na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


qcurve = function(df){
  #' @title Determine calibration curve
  #'
  #' @description 
  #' 
  #' @param df dataframe
  #' 
  
  
  #here we get the r-squared for each linear regression curve
  df = ddply(df,c("nuc"), transform, adj.r.squared = summary(lm(mean ~ level))$adj.r.squared)
  
  #here we get the y-intercept for each linear regression curve
  df = ddply(df,c("nuc"), transform, intercept = coefficients(lm(mean ~ level))[1])
  
  #here we get the slope for each linear regression curve
  df = ddply(df,c("nuc"), transform, slope = coefficients(lm(mean ~ level))[2])
  
  #let's write these data into a file
  write.table(df,"Calcurve_variables.tsv", row.names = F, col.names = T, 
              sep = "\t", quote = F)
  
  #given above info, let's write it into our graphs
  ggplot(df, aes(x = level, y = mean)) + 
    geom_point(size = 3) + 
    geom_smooth(method = "lm",se = FALSE) +
    facet_wrap(~nuc, scales = 'free') + 
    theme_bw() +
    ggtitle('Calibration curves (incl. regression coefficient)')+
    scale_x_log10() +
    scale_y_log10() + 
    geom_text(aes(label = paste0("r2=", format(round(adj.r.squared,3), nsmall = 3)), x = 400, y = 10000), size = 3) + 
    theme(legend.position = "none")
  
}


islinear = function(met, test){
  #' @title Linearity check
  #' 
  #' @description This function extracts the minimum and maximum value of each calibration curve
  #' and validates the location of each measurement accordingly.
  #' 
  #' @param met Metabolite name
  #' @param test 
  #' 
  #met: metabolite name
  #test: measured intensity of a sample 
  
  if(!all(is.na(test))){
    met.ch = as.character(met)
    
    curmin = min(qt$mean[qt$nuc == met.ch], na.rm = T)  #min val of calibration
    curmax = max(qt$mean[qt$nuc == met.ch], na.rm = T)  #max val of calibration
    answer = ifelse((test >= curmin & test <= curmax), 'linear', ifelse(test < curmin, 'below','above'))
  } else {
    answer = 'na'
  }
}


absconc = function(nuc, mean){
  if(!is.na(mean)){
    intercept = qt$intercept[qt$nuc == as.character(nuc)][1]
    slope = qt$slope[qt$nuc == as.character(nuc)][1]
    y = (mean - intercept)/slope
  } else {y = NA}
  y
}
