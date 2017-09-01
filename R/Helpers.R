
#' Creates plot folder in your working directory
#'
#' @param folder
#'
#' @return
#' @export
#'
#' @examples
create_plot_dir <- function(folder){
  Data_dir <- getwd()
  # Generate the plot directory in it does not exist
  if(TRUE){
    plot_dir_path <- paste0(Data_dir, "/",folder)
    # Check the existence of the folder and if no, create it
    if (!file_test("-d", plot_dir_path)){
      dir.create(plot_dir_path)
      print(paste0("Generated ",plot_dir_path))
    }
  }

}

#' Import masslists output files from TopDownLab
#'
#' @param folder folder in which files for import are
#' @param bound whether output data_list will have files combined e.i. C_0 and N_0
#' @param mod.number number of modifications allowed in top-down lab searches and outputs
#' @param dir directory to get files if not specified working directory
#'
#' @return
#' @export
#'
#' @examples
import_masslists <- function(dir = "", folder = "/", bound = FALSE, mod.number = 0, shiny = FALSE, files = inFiles){
if(!shiny){
  if(dir == "") dir <- getwd()
  Data_dir <- dir
  TDL_out_dir <- paste0(Data_dir,folder) # Topdownlab output folder
  file_list <- list.files(TDL_out_dir, pattern = ".masslist") # Get vector with all the files

  # Leave only C_O and N_0 files
  if(bound == TRUE & mod.number != 0)
  {
    if(mod.number == 1){
      file_list_C_0 <- subset(file_list, grepl("_C_0", file_list) | grepl("_N_0", file_list))
      file_list_C_1 <- subset(file_list, grepl("_C_1", file_list) | grepl("_N_1", file_list))
      file_list <- c(file_list_C_0, file_list_C_1)
    } else if(mod.number == 2) {
      file_list_C_0 <- subset(file_list, grepl("_C_0", file_list) | grepl("_N_0", file_list))
      file_list_C_1 <- subset(file_list, grepl("_C_1", file_list) | grepl("_N_1", file_list))
      file_list_C_2 <- subset(file_list, grepl("_C_2", file_list) | grepl("_N_2", file_list))
      file_list <- c(file_list_C_0, file_list_C_1, file_list_C_2)
    }else if(mod.number == 3) {
      file_list_C_0 <- subset(file_list, grepl("_C_0", file_list) | grepl("_N_0", file_list))
      file_list_C_1 <- subset(file_list, grepl("_C_1", file_list) | grepl("_N_1", file_list))
      file_list_C_2 <- subset(file_list, grepl("_C_2", file_list) | grepl("_N_2", file_list))
      file_list_C_3 <- subset(file_list, grepl("_C_3", file_list) | grepl("_N_3", file_list))
      file_list <- c(file_list_C_0, file_list_C_1, file_list_C_2, file_list_C_3)
    }
  }
  else file_list <- subset(file_list, grepl("_C_0", file_list) | grepl("_N_0", file_list) | grepl("_C_1", file_list) | grepl("_N_1", file_list) | grepl("_C_2", file_list) | grepl("_N_2", file_list) |grepl("_C_3", file_list) | grepl("_N_3", file_list))

  # create a list with all the data
  Data_list <- lapply(file_list,function(x) {read.table(paste0(TDL_out_dir,x),
                                                         sep = "\t", stringsAsFactors = FALSE, header = TRUE) })
  names(Data_list) <- file_list # Add filenames as list entrie names


} else {
  file_list <- files$datapath[grepl(".masslist", files$datapath)]
  file_names <- files$name[grepl(".masslist", files$name)]

  if(mod.number == 1){
    file_list_C_0 <- file_list[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_list_C_1 <- file_list[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_names_C_0 <- file_names[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_names_C_1 <- file_names[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_list <- c(file_list_C_0, file_list_C_1)
    file_names <- c(file_names_C_0, file_names_C_1)
  } else if(mod.number == 2){
    file_list_C_0 <- file_list[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_list_C_1 <- file_list[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_list_C_2 <- file_list[which(grepl("C_2", file_names) | grepl("N_2", file_names))]
    file_names_C_0 <- file_names[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_names_C_1 <- file_names[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_names_C_2 <- file_names[which(grepl("C_2", file_names) | grepl("N_2", file_names))]
    file_list <- c(file_list_C_0, file_list_C_1, file_list_C_2)
    file_names <- c(file_names_C_0, file_names_C_1, file_names_C_2)
  } else if(mod.number == 3){
    file_list_C_0 <- file_list[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_list_C_1 <- file_list[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_list_C_2 <- file_list[which(grepl("C_2", file_names) | grepl("N_2", file_names))]
    file_list_C_3 <- file_list[which(grepl("C_3", file_names) | grepl("N_3", file_names))]
    file_names_C_0 <- file_names[which(grepl("C_0", file_names) | grepl("N_0", file_names))]
    file_names_C_1 <- file_names[which(grepl("C_1", file_names) | grepl("N_1", file_names))]
    file_names_C_2 <- file_names[which(grepl("C_2", file_names) | grepl("N_2", file_names))]
    file_names_C_3 <- file_names[which(grepl("C_3", file_names) | grepl("N_3", file_names))]
    file_list <- c(file_list_C_0, file_list_C_1, file_list_C_2, file_list_C_3)
    file_names <- c(file_names_C_0, file_names_C_1, file_names_C_2, file_names_C_3)
  } else {file_list <- file_list[which(grepl("_C_0", file_names) | grepl("_N_0", file_names) | grepl("_C_1", file_names) | grepl("_N_1", file_names) | grepl("_C_2", file_names) | grepl("_N_2", file_names) |grepl("_C_3", file_names) | grepl("_N_3", file_names))]
          file_names <- file_names[which(grepl("C_0|C_1|C_2|C_3|N_0|N_1|N_2|N_3", file_names))]
  }
  Data_list <- lapply(file_list, function(x) {read.table(x,
                                                       sep = "\t", stringsAsFactors = FALSE, header = TRUE) })
  names(Data_list) <- file_names
}

  Data_list

}

#' Import masslists and arrange them into the list
#'
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
combine_masslists <- function(dir = "", shiny = FALSE) {
  if(!shiny) {
  if(dir == "") dir <- getwd()
  Data_list <- import_masslists(dir)
  file_names <- list.files(dir)
  file_names <- subset(file_names, grepl(".masslist", file_names))
  } else {
          Data_list <- import_masslists(shiny = TRUE)
          #file_list <- inFiles$datapath[grepl(".masslist", inFiles$datapath)]
          file_names <- inFiles$name[grepl(".masslist", inFiles$name)]
         }
  cluster <- max(substring(regmatches(file_names,
                                      gregexpr("_C_[[:digit:]]|_N_[[:digit:]]", file_names)), 4, 5))

  cluster <- as.numeric(cluster)*2 + 2

  if((length(Data_list)/cluster) > 1)
    Data_list_new <- by(seq_along(Data_list),cut(seq_along(Data_list),
                                                 (length(Data_list)/cluster)),
                        FUN=function(x)Data_list[x]) else {
                          Data_list_new[[1]] <- Data_list
                        }

  names(Data_list_new) <- subset(file_names, grepl("_C_0", file_names))


  gList <- lapply(names(Data_list_new), function(x) {
    l  <- Data_list_new[[x]]
    #todo for various number of modifications
    if(cluster == 6){
      c <- l[[1]]
      c1 <- l[[2]]
      c2 <- l[[3]]
      n <- l[[4]]
      n1 <- l[[5]]
      n2 <- l[[6]]

      c[is.na(c$Ion.type), "Ion.type"] <- ""
      c1[is.na(c1$Ion.type), "Ion.type"] <- ""
      c2[is.na(c2$Ion.type), "Ion.type"] <- ""
      n[is.na(n$Ion.type), "Ion.type"] <- ""
      n1[is.na(n1$Ion.type), "Ion.type"] <- ""
      n2[is.na(n2$Ion.type), "Ion.type"] <- ""

      c[c$Ion.type == "",5:10] <- n[c$Ion.type == "",5:10]
      c1[c1$Ion.type == "",5:10] <- n1[c1$Ion.type == "",5:10]
      c2[c2$Ion.type == "",5:10] <- n2[c2$Ion.type == "",5:10]

      c$mod = 0
      c1$mod = 1
      c2$mod = 2

      cn <- c
      cn1 <- c1
      cn2 <- c2
      #cn <- cn[cn$Intensity > 2000,]

      df <- rbind(cn, cn1, cn2)
      df$exp <- x } else if(cluster == 4) {

        c <- l[[1]]
        c1 <- l[[2]]

        n <- l[[3]]
        n1 <- l[[4]]


        c[is.na(c$Ion.type), "Ion.type"] <- ""
        c1[is.na(c1$Ion.type), "Ion.type"] <- ""

        n[is.na(n$Ion.type), "Ion.type"] <- ""
        n1[is.na(n1$Ion.type), "Ion.type"] <- ""


        c[c$Ion.type == "",5:10] <- n[c$Ion.type == "",5:10]
        c1[c1$Ion.type == "",5:10] <- n1[c1$Ion.type == "",5:10]


        c$mod = 0
        c1$mod = 1


        cn <- c
        cn1 <- c1

        #cn <- cn[cn$Intensity > 2000,]

        df <- rbind(cn, cn1)
        df$exp <- x

      } else if(cluster == 2) {
        c <- l[[1]]
        n <- l[[2]]

        c[c$Ion.type == "", 5:10] <- n[c$Ion.type == "", 5:10]
        c$mod = 0

        cn <- c
        df <- cn
        df$exp <- x
      } else print("Incorrect input!")

    df$Intensity <- df$Intensity/max(df$Intensity)*100
    df

  })
  names(gList) <- subset(file_names, grepl("_C_0", file_names))

  gList
}


process_masslists <- function(expdef = FALSE, data = Data_list, x = "/", ppm = 3.5, sep_term = FALSE, mod.number = 0, tint = FALSE, seq = seq){

  Data_list_new <- data

  # get the sum of intensities of assigned peaks for both termini including both mod0 and mod1

  gList <- lapply(names(Data_list_new), function(x) {
    #Taking first entry of the list and working with it
    experiment <- x
    #Getting rid of unassigned masses
    ass_peak_table <- Data_list_new[[experiment]]
    #Shaping the data frame with info from the entrie's name
    if(nrow(ass_peak_table) > 0){
      ass_peak_table$exp.name <- experiment
      #source fragmentation value
      ass_peak_table$sf <- as.numeric(substr(experiment, regexpr("sf[0-9]+",experiment)+2,
                                             regexpr("sf[0-9]+",experiment) + attr(regexpr("sf[0-9]+", experiment),"match.length") -1))
      #hcd energy
      ass_peak_table$hcd <- as.numeric(substr(experiment, regexpr("hcd[0-9]+",experiment)+3,
                                              regexpr("hcd[0-9]+", experiment) + attr(regexpr("hcd[0-9]+", experiment), "match.length") -1) )
      #laser energy
      ass_peak_table$energy <- as.numeric(substr(experiment, regexpr("uvpd[0-9.]+",experiment)+4,
                                                 regexpr("uvpd[0-9.]+", experiment) + attr(regexpr("uvpd[0-9.]+", experiment), "match.length")-1))

      ass_peak_table$etd <- as.numeric(substr(experiment, regexpr("etd[0-9.]+",experiment)+3,
                                              regexpr("etd[0-9.]+", experiment) + attr(regexpr("etd[0-9.]+", experiment), "match.length")-1))

      ass_peak_table$cid <- as.numeric(substr(experiment, regexpr("cid[0-9.]+",experiment)+3,
                                              regexpr("cid[0-9.]+", experiment) + attr(regexpr("cid[0-9.]+", experiment), "match.length")-1))

      ass_peak_table$iso <- substr(experiment, regexpr("iso[0-9]+", experiment) +3,
                                   regexpr("iso[0-9]+", experiment) + attr(regexpr("iso[0-9]+", experiment), "match.length")-1)

      ass_peak_table$puls <- substr(experiment, regexpr("puls[0-9]", experiment) +4, regexpr("puls[0-9]", experiment)+4)

      ass_peak_table$enpuls <- paste0(ass_peak_table$energy, "x", ass_peak_table$puls)

      ass_peak_table$Pressure <- substr(experiment, regexpr("pres[0-9]",experiment)+4, regexpr("pres[0-9]", experiment)+4)

      ass_peak_table$nat.denat <- substr(experiment, regexpr("nat|denat", experiment),
                                         regexpr("nat|denat", experiment) + attr(regexpr("nat|denat", experiment), "match.length")-1)

      #peak intensity normalized by the sum of intensities of all peaks
      ass_peak_table$norm.intens <- ass_peak_table$Intensity/sum(ass_peak_table$Intensity)

      ass_peak_table$Log10_int <- log10(ass_peak_table$Intensity)

      ass_peak_table$log10.norm.intens <- log10(ass_peak_table$norm.intens)

      #calculate Z-scored intensities
      intens_mean <- mean(ass_peak_table$norm.intens)
      intens_sd <- sd(ass_peak_table$norm.intens)
      ass_peak_table$norm.intens.Zscored <- (ass_peak_table$norm.intens - intens_mean)/intens_sd

      ass_peak_table$seq.cov <- length(unique(ass_peak_table$Sequence.position))/(nchar(seq)-1)

      ass_peak_table$term <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1)

      ass_peak_table$mod <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3)

      ass_peak_table$chg.iso <- paste0(ass_peak_table$Charge, " ", ass_peak_table$iso)

      ass_peak_table$rep <- substr(experiment, regexpr("rep[0-9]+", experiment) +3,
                                   regexpr("rep[0-9]+", experiment) + attr(regexpr("rep[0-9]+", experiment), "match.length")-1)


      if (expdef == TRUE)
      {
        if (ass_peak_table$energy > 0 & ass_peak_table$hcd == 1) ass_peak_table$exp <- "UVPD"
        else if (ass_peak_table$hcd > 1 & !is.na(ass_peak_table$hcd))  ass_peak_table$exp <- "HCD"
        else if (ass_peak_table$energy > 0 & ass_peak_table$hcd > 1) ass_peak_table$exp <- "UVhcpD"
        else if (fusion == TRUE & ass_peak_table$etd > 0 & ass_peak_table$cid == 0 | is.na(ass_peak_table$cid)) ass_peak_table$exp <- "ETD"
        else if (fusion == TRUE & ass_peak_table$cid > 0 & ass_peak_table$etd == 0) ass_peak_table$exp <- "CID"
        else if (fusion == TRUE & ass_peak_table$cid > 0 & ass_peak_table$etd > 0) ass_peak_table$exp <- "ETciD"
        else ass_peak_table$exp <- "hz"
      }
      ass_peak_table$Mass.shift[is.na(ass_peak_table$Mass.shift)] <- ""

      #ass_peak_table$frag.type <- paste0(ass_peak_table$Ion.type," ",ass_peak_table$Mass.shift)
      #categorize to more general fragment types
      ass_peak_table <- as.data.frame(t(apply(ass_peak_table, 1, function(x) {
        if (x["Ion.type"] == "A" | x["Ion.type"] == "X") x["Ion.type.clust"] <- "a/x"
        else if (x["Ion.type"] == "B" | x["Ion.type"] == "Y") x["Ion.type.clust"] <- "b/y"
        else if (x["Ion.type"] == "C" | x["Ion.type"] == "Z") x["Ion.type.clust"] <- "c/z"
        else if (x["Ion.type"] == "D" | x["Ion.type"] == "V" | x["Ion.type"] == "W") x["Ion.type.clust"] <- "d/v/w"
        else x["Ion.type.clust"] <- ""
        #get rid of missing values in each row
        x[which(is.na(x))] <- ""
        x
      })))

      ass_peak_table$Amino.acid <- as.factor(ass_peak_table$Amino.acid)
      #ass_peak_table <- within(ass_peak_table, frag.type <- ordered(frag.type, levels = rev(sort(unique(frag.type)))))

      #get amino acid composition of the protein and make vector of frequencies with aminoacid letters as names
      aa_table <- lapply(strsplit(seq,""), table)
      aa_df <- as.data.frame(aa_table)
      aa_vector <- aa_df$Freq
      names(aa_vector) <- aa_df$Var1

      #normalize normalized intensity by aminoacid frequencies in vector aa_vector in accordance with $Amino.acid
      ass_peak_table$norm.intens.x.freq <- ass_peak_table$norm.intens/aa_vector[as.character(ass_peak_table$Amino.acid)]

      ass_peak_table$Log10.norm.intens.x.freq <- log10(ass_peak_table$norm.intens.x.freq)

      if(tint == TRUE)
      {
        fragef <- ass_peak_table %>% group_by(Sequence.position, Amino.acid) %>% summarise(fragef = sum(Intensity)/totalint[experiment])
        fragef$enpuls <- ass_peak_table$enpuls[1]
        fragef$etd <- ass_peak_table$etd[1]
        fragef$cid <- ass_peak_table$cid[1]
        fragef$hcd <- ass_peak_table$hcd[1]
        if(bound == TRUE) fragef$mod <- ass_peak_table$mod[1]
      }
    }
   ass_peak_table
  })
  names(gList) <- names(Data_list_new)
  gDF <- bind_rows(gList)
  print("Merged data frame completed")
  data_merged <- gDF
  data_merged
}


