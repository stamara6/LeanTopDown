library(scales)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
library(base)

library(tidyr)
library(compare)
library(RColorBrewer)
library(oro.nifti)

#delete all merged datasets
#rm(ass_data_merged, fragef_merged)

Data_dir <- getwd()


unique_frag_count_UVPD_HCD <- function(Plot_title_suffix = "", Plot_dir_name = "plot/", Plot_file_name = "unique frag seq cov UVPD_HCD", Plot_width = 10, Plot_height = 10) {
  #levels(ass_data_merged$hcd) <- c("20","40","60","80","100","120","140","160","180")
  #ass_data_merged <- ass_data_merged[!duplicated(ass_data_merged[6:8]),]

  ass_data_merged$fragment <- paste0(ass_data_merged$Sequence.position, " ", ass_data_merged$term, " ", ass_data_merged$Series.number)
  ass_data_dup <- ass_data_merged[duplicated(ass_data_merged[6:8]) | duplicated(ass_data_merged[6:8], fromLast = TRUE),]
  ass_data_merged <- ass_data_merged[!(ass_data_merged$fragment %in% ass_data_dup$fragment),]
  df <- count(unique(ass_data_merged[,c("exp.name","Sequence.position")]), vars = "exp.name")
  df$seq.cov <- (df$freq/(nchar(seq)-1))*100
  df["all","seq.cov"] <- (length(unique(ass_data_merged$Sequence.position))/(nchar(seq)-1))*100
  #df$cols <- c('#9e9ac8','#756bb1','#54278f')
  df$pos <- df$seq.cov + 5
  df$hcd <- as.numeric(substr(df$exp.name, regexpr("hcd[0-9]+", df$exp.name)+3,
                              regexpr("hcd[0-9]+", df$exp.name) + attr(regexpr("hcd[0-9]+", df$exp.name), "match.length") -1))
  df$energy <- as.numeric(substr(df$exp.name, regexpr("uvpd[0-9.]",df$exp.name)+4,
                                 regexpr("uvpd[0-9.]", df$exp.name) + attr(regexpr("uvpd[0-9.]", df$exp.name), "match.length")+1))
  myPalette <- colorRampPalette(hotmetal())


  plot_gg <- ggplot(ass_data_merged, aes(x = Sequence.position, y = term, fill = norm.intens.Zscored)) +
    stat_identity(geom = "tile", width  = .5) +
    #geom_hline(yintercept = 1.5) +
    #geom_histogram(aes(color = Ion.type)) +
    #geom_linerange (aes(x = Sequence.position, y = norm.intens), position = position_dodge()) +
    #stat_identity(aes(color = Ion.type)) +
    #scale_y_discrete( expand = c(0,0), limits= c(1,1)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, nchar(seq)), breaks = seq(0, nchar(seq), by = 10)) +
    theme_bw() +
    facet_grid(enpuls ~ mod) +
    #scale_fill_gradient2(limits = c(0, max(ass_data_merged$norm.intens.Zscored)), name = "Intensity", high = "black", mid = "red", low = "yellow", midpoint = median(ass_data_merged$norm.intens.Zscored), guide = "legend") +
    scale_fill_gradientn(colours = myPalette(100),name = "Normalazed\nintensity") +
    ylab("") + # x label
    xlab("Sequence position") +  # y label
    ggtitle(paste0("Sequence Coverage\n", Plot_title_suffix)) +

    #coord_flip() +
    #geom_text(aes(y=pos,label=paste0(as.integer(seq.cov), "%")), colour = "black", size=5, angle = 60, position = position_dodge(width = .55)) +
    theme(legend.position="right", title=element_text(face = "bold"), text = element_text(size=12),
          axis.text.x = element_text())

  print(plot_gg)
  Plot_dir_name <- "plot"
  plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width, Plot_height)
  save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))
}


#plot counts per fragment type cluster stacked
fragment_types_clustered_counted <- function(data = ass_data_merged, Data_list_new, Plot_dir_name = "plot", Plot_file_name = "frag_types_count_mod", Plot_title_suffix = "", Plot_width = 10, Plot_height = 10){
  #plotting of by xa and zc ions stacked horizontally and with hcd facetting (count labelling for each pair of frags)
  require(plyr)
  Data_dir <- getwd()
  seq_cov_data <- data

  #seq_cov_data$Ion.type <- factor(seq_cov_data$Ion.type)
  #seq_cov_data <- subset(seq_cov_data, seq_cov_data$rep < 3)
  cols <- c('#ef8a62','#91cf60','#67a9cf')
  seq_cov_data <- seq_cov_data[ ,c("rep", "enpuls", "Ion.type.clust", "puls", "energy", "mod")]
  seq_cov_data <- ddply(seq_cov_data, c("Ion.type.clust", "enpuls", "rep", "mod"), transform, count = count(Ion.type.clust))
  seq_cov_data <- distinct(seq_cov_data, Ion.type.clust, enpuls, rep, mod)
  seq_cov_data <- ddply(seq_cov_data, c("enpuls", "rep", "mod"), transform, pos = cumsum(count.freq) - 0.5*count.freq)


  plot_gg <- ggplot(seq_cov_data, aes(x = factor(enpuls))) +
    geom_bar(aes(y = count.freq, fill = Ion.type.clust), stat = "identity") +
    #stat_params = list(binwidth = ) +
    scale_y_continuous(expand = c(0,0)) +

    #scale_x_continuous("energy", breaks = 0.5:2.0) + #define isolated m/z
    geom_text(aes(y = pos, label = count.freq), colour = "black", size=3) +

    facet_grid(mod ~ rep, labeller = "label_both") +
    ggtitle(paste0("Fragment Ion Types\n",Plot_title_suffix)) + # plot title
    #geom_text(aes(y = ..count..), vjust = 3) +
    ylab("Count") + # x label
    xlab("Laser Energy") +  # y label
    #scale_fill_manual(values = cols) +

    coord_flip() +
    theme_bw() +
    #theme_minimal() +
    theme(legend.position="right", title=element_text(face = "bold"), text = element_text(size=14),
          axis.text.x = element_text())

  print(plot_gg)
  Plot_dir_name <- "plot"
  plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width, Plot_height)
  save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))

}

frag.type.boxplot <- function(data = ass_data_merged, Data_list_new, Plot_dir_name = "plot", Plot_file_name = "fragtypes_pin16al7_50th_annall", Plot_title_suffix = "", Plot_width = 18, Plot_height = 10){
  #plotting of by xa and zc ions stacked horizontally and with hcd facetting (count labelling for each pair of frags)
  require(plyr)
  Data_dir <- getwd()

  cols <- c('#ef8a62','#91cf60','#67a9cf')
  data$mod <- factor(data$mod,
                     levels = c(0, 1, 2),
                     labels = c("No modification", "1 phospho site", "2 phospho sites"))

  plot_gg <- ggplot(data, aes(x = frag.type, y = norm.intens.Zscored)) +
    geom_boxplot(aes(fill = Ion.type.clust), size = .5) +
    #stat_params = list(binwidth = ) +
    scale_y_continuous(expand = c(0,0)) +

    #scale_x_continuous("energy", breaks = 0.5:2.0) + #define isolated m/z
    #geom_text(aes(y = pos, label = count.freq), colour = "black", size=3) +

    facet_grid(mod ~ .) +
    ggtitle(paste0("Fragment Ion Types\n",Plot_title_suffix)) + # plot title
    #geom_text(aes(y = ..count..), vjust = 3) +
    ylab("Z-scored Intensities") + # x label
    xlab("Fragment Types") +  # y label
    #scale_fill_manual(values = cols) +
    scale_fill_discrete(name = "Ion Type\nSeries") +
    #coord_flip() +
    theme_bw() +
    #theme_minimal() +
    theme(legend.position="right", title=element_text(), text = element_text(size=14),
          axis.text.x = element_text()) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  print(plot_gg)
  Plot_dir_name <- "plot"
  plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width, Plot_height)
  save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))

}


termini_sequence_coverage <- function(data = ass_data_merged, Data_list_new, Plot_dir_name = "plot", Plot_file_name = "Pin1_ck2_term_seq_cov_mod_rep", Plot_title_suffix = "", Plot_width = 10, Plot_height = 10){

  Data_dir <- getwd()

  seq_cov_data <- data[!duplicated(ass_data_merged$exp.name),]

  seq_cov_data$energy <- as.factor(seq_cov_data$energy)

  #Determine values for text plotting
  seq_cov_data <- ddply(seq_cov_data,"seq.cov", transform, perc_cov = seq.cov*100)
  seq_cov_data$perc_cov <- as.integer(seq_cov_data$perc_cov)
  seq_cov_data$perc_cov_y <- seq_cov_data$perc_cov + 5

  seq_cov_data <- subset(seq_cov_data, seq_cov_data$rep < 3)

  cols <- c('#9e9ac8','#756bb1','#54278f')

  plot_gg <- ggplot(seq_cov_data, aes(x = enpuls, y = perc_cov)) +
    stat_identity(aes(group = interaction(enpuls,term), y = perc_cov, fill = term), geom = "bar", position = position_dodge(width = .75), width = .5) +

    scale_y_continuous(breaks=seq(10,100,by = 10), expand = c(0,0)) +

    #scale_x_discrete(limits = c("1", "100", "150")) + #define isolated m/z

    #ggtitle(paste0("Sequence Coverage with Varying Laser Energies",Plot_title_suffix)) + # plot title

    ylab("Sequence coverage") + # x label
    xlab("Energy, mj") +  # y label

    geom_text(aes(group = interaction(enpuls, term), y = perc_cov + 2, label = paste0(perc_cov, "%")), colour = "black", size=4, position = position_dodge(width = .75), angle = 70) +

    facet_grid(rep ~ mod, labeller = "label_both") +

    theme_bw() +
    #theme_minimal() +
    theme(legend.position="top", title=element_text(face = "bold"), text = element_text(size=12),
          axis.text.x = element_text())

  print(plot_gg)
  Plot_dir_name <- "plot"

  plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width, Plot_height)
  save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))

}

#Plot fragments on respective sequence positions with fill according to their z-scored intensities
#Plot Fragments' Charges
#1 choose wd
#2 create_plot_dir("plot")
#3 import_TDL_out("/")
#4 rm(ass_data_merged)
#5 set right aa sequence
#6 prepare_data_frame_cterm_nterm_sep_aa()
#7 set facetting variable, Plot_file_name, subsetting values, and point attributes in z_frags()
#8 source after changes
#9 z_frags()
#rm(ass_data_merged, fragef_merged)
fragment_charges <- function(uvhcpd = F, subset = F, Plot_title_suffix = "", Plot_dir_name = "plot/", Plot_file_name = "z_frags_pin1_z8_all_ck2", Plot_width = 10, Plot_height = 10, iso = 0, pep = FALSE) {

  Data_dir <- getwd()

  #ass_data_merged$energy <- as.factor(ass_data_merged$energy)

  ass_data_merged <- ass_data_merged[order(ass_data_merged$term, decreasing = TRUE),]

  #ass_data_merged <- subset(ass_data_merged, !ass_data_merged$iso == 2290 & !ass_data_merged$mod == 1)
  if(uvhcpd == TRUE) ass_data_merged$uvhcpd <- paste0(ass_data_merged$hcd,"Vx",ass_data_merged$energy, "mJ")

  #ass_data_merged$iso <- as.factor(ass_data_merged$iso)
  #ass_data_merged$puls <- as.numeric(ass_data_merged$puls)
  #ass_data_merged$energy <- as.factor(ass_data_merged$energy)
  if (pep == TRUE)
    {
      s <- seq(-1, 30, by = 2)
    }
  else
    {
      s <- seq(0, 500, by = 5)

    }

  ass_data_merged$exp <- regmatches(ass_data_merged$exp, gregexpr("hcd[[:digit:]]+_etd[[:digit:]]+_sa[[:digit:]]+", ass_data_merged$exp))
  ass_data_merged$exp <- unlist(ass_data_merged$exp)
  if (subset == TRUE)
  {
  ass_data_merged <- subset(ass_data_merged, ass_data_merged$Ion.type.clust == "c/z")
  #ass_data_merged$exp == "hcd26_etd75_sa20" &
  }
  sequence <- unlist(strsplit(seq, split = ""))
  names(sequence) <- c(1:nchar(seq))
  plot_gg <- ggplot(ass_data_merged[abs(ass_data_merged$calib.ppm) <= 3.5,], aes(x = Sequence.position, y = Charge)) +
    geom_point(stat = "identity", aes(shape = Ion.type.clust, color = mod, size = Intensity)) +
    scale_y_continuous(breaks = seq(0,10, by = 1), limits = c(0.5, ceiling(max(ass_data_merged$Charge)))) +
    scale_x_continuous(breaks = seq(0,length(strsplit(seq, split = "")[[1]]), by = 10)) +
    scale_shape_discrete(name = "Ion Type\nSeries") +
    scale_color_discrete(name = "Number of\nPhoshosites") +
    guides(size = FALSE) +
    geom_text(aes(y = 0.7, label = Amino.acid, x = Sequence.position), size = 1.8)  +
    #facet_grid(mod ~ .) +
    theme_bw()


  print(plot_gg)


  plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width = Plot_width, Plot_height = Plot_height)
  save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))

}



prepare_data_frame_term_sep_aa <- function(expdef = FALSE, fusion = FALSE, x = "/", ppm = 3.5, bound = FALSE, sep_term = FALSE, mod.number = 0, tint = TRUE){

  if(sep_term == TRUE)
  {
    Data_dir <- getwd()
    file_list <- list.files(paste0(Data_dir,x), pattern = ".masslist")
    #for when you have masslists with modification (e.g. _C_1, _N_1) but want only to work with those that lack modification
    # file_list <- subset(file_list, grepl("_C_0", file_list) & grepl("_N_0", file_list))
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
    Data_list_new <- Data_list

    #Bind n_o and c_o into one data_frame and assign to each list entry
    # for(ind in 1:length(Data_list_new)){
    #  experiment <- names(Data_list_new)[ind]
    # Data_list_new[[experiment]] <- rbind2(Data_list_new[[c(ind, 1)]],Data_list_new[[c(ind, 2)]]) }

    names(Data_list_new) <- file_list #names for list entries
  }
  else
  {
    Data_dir <- getwd()
    file_list <- list.files(paste0(Data_dir,x), pattern = ".masslist")
    if (bound == TRUE)
    {
      file_list_0 <- subset(file_list, grepl("_C_0", file_list))
      file_list_1 <- subset(file_list, grepl("_C_1", file_list))
      file_list <- c(file_list_0, file_list_1)
    }
    else file_list <- subset(file_list, grepl("_C_0", file_list))
    Data_list_new <- by(seq_along(Data_list),cut(seq_along(Data_list),(length(Data_list)/2)), FUN=function(x)Data_list[x]) #make list of sublists containing n_o and c_o

    for(ind in 1:length(Data_list_new)){
      experiment <- names(Data_list_new)[ind]
      Data_list_new[[experiment]] <- rbind2(Data_list_new[[c(ind, 1)]],Data_list_new[[c(ind, 2)]]) }  #bind n_o and c_o into one data_frame and assign to each list entry

    names(Data_list_new) <- file_list #name the sublists
  }
  # get the sum of intensities of assigned peaks for both termini including both mod0 and mod1
  if (tint == T)
  {
    Data_dir <- getwd()
    TDL_out_dir <- paste0(Data_dir,"/")

    fl <- list.files(TDL_out_dir, pattern = ".masslist") # Get vector with all the files

    fl <- subset(fl, grepl("_C_0", fl) | grepl("_N_0", fl) | grepl("_C_1", fl) | grepl("_N_1", fl) | grepl("_C_2", fl) | grepl("_N_2", fl) | grepl("_C_3", fl) | grepl("_N_3", fl))

    DL <- lapply(fl, function(x) {read.table(paste0(TDL_out_dir, x),
                                             sep = "\t", stringsAsFactors = FALSE, header = TRUE) })

    cluster <- max(substring(regmatches(file_list,
                                        gregexpr("_C_[[:digit:]]|_N_[[:digit:]]", file_list)), 4, 5))

    cluster <- as.numeric(cluster)*2 + 2

    if (bound == TRUE)
    {

      DL <- by(seq_along(DL), cut(seq_along(DL),(length(DL)/cluster)), function(x)DL[x])
    }
    else DL <- by(seq_along(DL), cut(seq_along(DL),(length(DL)/2)), function(x)DL[x])
    fl <- subset(fl, grepl("_C_0", fl))

    names(DL) <- fl
    totalint <- c(1:length(DL))

    for (ind in 1:length(DL))
    {
      #discard all the unidentified peaks
      DL_T <- lapply(DL[[ind]], function (x) {subset(x, !is.na(x$Sequence.position))})
      #assign summarized intensities for experiment in one named vector
      totalint[ind] <- sum(sapply(DL_T, function(x) {sum(x$Intensity)}))
    }
    if(bound == TRUE & sep_term == FALSE)
    {
      totalint <- rep(totalint, cluster)
      names(totalint) <- file_list
    }
    else if(bound == TRUE & sep_term == TRUE)
    {
      totalint <- rep(totalint, cluster)
      names(totalint) <- file_list
    }
    else  names(totalint) <- file_list
  }
  if (!exists("ass_data_merged")) {
    #Taking first entry of the list and working with it
    experiment <- names(Data_list_new)[1]
    #ass_peak_table <- data.frame(NA)
    #Getting rid of unassigned masses
    ass_peak_table <- subset(Data_list_new[[experiment]],Data_list_new[[experiment]]$Ion.type != "")

    ass_peak_table$Accuracy..PPM. <- as.numeric(ass_peak_table$Accuracy..PPM.)
    ass_peak_table$m.z.calibration <- as.numeric(ass_peak_table$m.z.calibration)
    ass_peak_table$calib.ppm <- ass_peak_table$Accuracy..PPM. - ass_peak_table$m.z.calibration
    #Setting accuracy value, discarding everything that is less accurate
    ass_peak_table <- subset(ass_peak_table, abs(ass_peak_table$calib.ppm) <= ppm) #& abs(ass_peak_table$Accuracy..PPM.) >= 5)
    #Shaping the data frame with info from the entrie's name

    if (nrow(ass_peak_table) > 0 & any(!is.na(ass_peak_table))){
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

      #for (i in 1:nrow(ass_peak_table)) {
      #  for (j in 1:nrow(ass_peak_table)) {

      #    if (ass_peak_table$Sequence.position[i])


      #  }
      #}

      #ass_peak_table$position_frag_numbers <- sapply(ass_peak_table$Sequence.position, function(x){sum(subset(ass_peak_table$norm.intens,ass_peak_table$Sequence.position==x)) })
      # ass_peak_table <- ass_peak_table[!duplicated(ass_peak_table$Sequence.position),]
      # ass_peak_table <- ass_peak_table[order(ass_peak_table$Sequence.position), ]
      # ass_peak_table$type.or.loss <- character(nrow(ass_peak_table))
      #for (ind in 1:nrow(ass_peak_table)){
      #   if(ass_peak_table$Mass.shift[ind] == "" | is.na(ass_peak_table$Mass.shift[ind])) ass_peak_table$type.or.loss[ind] <- ass_peak_table$Ion.type[ind]
      #   else if (ass_peak_table$Mass.shift[ind] == "Ammonia loss") ass_peak_table$type.or.loss[ind] <- "_NH3_Loss"
      #   else if (ass_peak_table$Mass.shift[ind] == "Water loss") ass_peak_table$type.or.loss[ind] <- "_H2O_Loss"
      #  else ass_peak_table$type.or.loss[ind] <- "WTF????"
      # }
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

      ass_peak_table$frag.type <- paste0(ass_peak_table$Ion.type," ",ass_peak_table$Mass.shift)

      for (i in 1:nrow(ass_peak_table))
      {
        if (ass_peak_table$Ion.type[i] == "A" | ass_peak_table$Ion.type[i] == "X") ass_peak_table$Ion.type.clust[i] <- "a/x"
        else if (ass_peak_table$Ion.type[i] == "B" | ass_peak_table$Ion.type[i] == "Y") ass_peak_table$Ion.type.clust[i] <- "b/y"
        else if (ass_peak_table$Ion.type[i] == "C" | ass_peak_table$Ion.type[i] == "Z") ass_peak_table$Ion.type.clust[i] <- "c/z"
        else if (ass_peak_table$Ion.type[i] == "D" | ass_peak_table$Ion.type[i] == "V" | ass_peak_table$Ion.type[i] == "W") ass_peak_table$Ion.type.clust[i] <- "d/v/w"
      }

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

      ass_peak_table$uvpd_hcd_pres_iso_nat <- paste0(ass_peak_table$energy," ",ass_peak_table$hcd," ",ass_peak_table$Pressure," ",ass_peak_table$iso, " ", ass_peak_table$nat.denat)
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
    else {ass_data_merged <<- ass_peak_table}
  }
  if(!is.null(tint) & tint == TRUE & !is.null(nrow(ass_peak_table))) fragef_merged <<- fragef
  ass_data_merged <<- ass_peak_table

  for (experiment in names(Data_list_new)[2:length(names(Data_list_new))])
  {
    #ass_peak_table <- data.frame(NA)

    ass_peak_table <- subset(Data_list_new[[experiment]], Data_list_new[[experiment]]$Ion.type != "")
    #if there are no fragments go for the next cycle iteration
    if(is.null(ass_peak_table)) next;

    if (nrow(ass_peak_table) > 0 & any(!is.na(ass_peak_table)))
    {

      ass_peak_table$exp.name <- experiment[1]

      ass_peak_table$Accuracy..PPM. <- as.numeric(ass_peak_table$Accuracy..PPM.)
      ass_peak_table$m.z.calibration <- as.numeric(ass_peak_table$m.z.calibration)
      ass_peak_table$calib.ppm <- ass_peak_table$Accuracy..PPM. - ass_peak_table$m.z.calibration
      #Setting accuracy value, discarding everything that is less accurate
      ass_peak_table <- subset(ass_peak_table, abs(ass_peak_table$calib.ppm) <= ppm) #& abs(ass_peak_table$Accuracy..PPM.) >= 5)

      if (nrow(ass_peak_table) >0){

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


        #peak intensity normalized by the sum of intensities of all peaks

        ass_peak_table$iso <- substr(experiment, regexpr("iso[0-9]+", experiment) +3,
                                     regexpr("iso[0-9]+", experiment) + attr(regexpr("iso[0-9]+", experiment), "match.length")-1)

        ass_peak_table$puls <- substr(experiment, regexpr("puls[0-9]", experiment) +4, regexpr("puls[0-9]", experiment)+4)

        ass_peak_table$enpuls <- paste0(ass_peak_table$energy, "x", ass_peak_table$puls)

        ass_peak_table$Pressure <- substr(experiment, regexpr("pres[0-9]",experiment)+4, regexpr("pres[0-9]", experiment)+4)

        ass_peak_table$nat.denat <- substr(experiment, regexpr("nat|denat", experiment),
                                           regexpr("nat|denat", experiment) + attr(regexpr("nat|denat", experiment), "match.length")-1)

        ass_peak_table$norm.intens <- ass_peak_table$Intensity/sum(ass_peak_table$Intensity)

        ass_peak_table$Log10_int <- log10(ass_peak_table$Intensity)

        ass_peak_table$log10.norm.intens <- log10(ass_peak_table$norm.intens)

        intens_mean <- mean(ass_peak_table$norm.intens)
        intens_sd <- sd(ass_peak_table$norm.intens)
        ass_peak_table$norm.intens.Zscored <- (ass_peak_table$norm.intens - intens_mean)/intens_sd

        ass_peak_table$seq.cov <- length(unique(ass_peak_table$Sequence.position))/(nchar(seq)-1)

        ass_peak_table$term <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1)

        ass_peak_table$mod <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3)

        ass_peak_table$chg.iso <- paste0(ass_peak_table$Charge, " ", ass_peak_table$iso)

        ass_peak_table$rep <- substr(experiment, regexpr("rep[0-9]+", experiment) +3,
                                     regexpr("rep[0-9]+", experiment) + attr(regexpr("rep[0-9]+", experiment), "match.length")-1)

        #ass_peak_table$position_frag_numbers <- sapply(ass_peak_table$Sequence.position, function(x){sum(subset(ass_peak_table$norm.intens,ass_peak_table$Sequence.position==x)) })
        # ass_peak_table <- ass_peak_table[!duplicated(ass_peak_table$Sequence.position),]
        # ass_peak_table <- ass_peak_table[order(ass_peak_table$Sequence.position), ]
        # ass_peak_table$type.or.loss <- character(nrow(ass_peak_table))
        #for (ind in 1:nrow(ass_peak_table)){
        #   if(ass_peak_table$Mass.shift[ind] == "" | is.na(ass_peak_table$Mass.shift[ind])) ass_peak_table$type.or.loss[ind] <- ass_peak_table$Ion.type[ind]
        #   else if (ass_peak_table$Mass.shift[ind] == "Ammonia loss") ass_peak_table$type.or.loss[ind] <- "_NH3_Loss"
        #   else if (ass_peak_table$Mass.shift[ind] == "Water loss") ass_peak_table$type.or.loss[ind] <- "_H2O_Loss"
        #  else ass_peak_table$type.or.loss[ind] <- "WTF????"
        # }
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

        ass_peak_table$frag.type <- paste0(ass_peak_table$Ion.type," ",ass_peak_table$Mass.shift)

        for (i in 1:nrow(ass_peak_table))
        {
          if (ass_peak_table$Ion.type[i] == "A" | ass_peak_table$Ion.type[i] == "X") ass_peak_table$Ion.type.clust[i] <- "a/x"
          else if (ass_peak_table$Ion.type[i] == "B" | ass_peak_table$Ion.type[i] == "Y") ass_peak_table$Ion.type.clust[i] <- "b/y"
          else if (ass_peak_table$Ion.type[i] == "C" | ass_peak_table$Ion.type[i] == "Z") ass_peak_table$Ion.type.clust[i] <- "c/z"
          else if (ass_peak_table$Ion.type[i] == "D" | ass_peak_table$Ion.type[i] == "V" | ass_peak_table$Ion.type[i] == "W") ass_peak_table$Ion.type.clust[i] <- "d/v/w"
        }

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

        ass_peak_table$uvpd_hcd_pres_iso_nat <- paste0(ass_peak_table$energy," ",ass_peak_table$hcd," ",ass_peak_table$Pressure," ",ass_peak_table$iso, " ", ass_peak_table$nat.denat)
        if(tint == TRUE)
        {
          fragef <- ass_peak_table %>% group_by(Sequence.position, Amino.acid) %>% summarise(fragef = sum(Intensity)/totalint[experiment])
          fragef$enpuls <- ass_peak_table$enpuls[1]
          fragef$etd <- ass_peak_table$etd[1]
          fragef$cid <- ass_peak_table$cid[1]
          fragef$hcd <- ass_peak_table$hcd[1]
          if(bound == TRUE) fragef$mod <- ass_peak_table$mod[1]
        }
        ass_data_merged <<- merge(ass_data_merged,ass_peak_table, all = TRUE, sort = FALSE)
        if(tint == TRUE & exists("fragef_merged") & nrow(ass_peak_table) > 0) fragef_merged <<- merge(fragef_merged, fragef, all = TRUE, sort = FALSE)
        else if (tint == TRUE & nrow(ass_peak_table) > 0) fragef_merged <<- fragef
       }
    }
  }
  print("Merged data frame completed")
}

prepare_data_frame_term_sep_aa_all <- function(expdef = FALSE, fusion = FALSE, x = "/", ppm = 3.5, bound = FALSE, sep_term = FALSE, mod.number = 0, tint = TRUE){

  if(sep_term == TRUE)
  {
    Data_dir <- getwd()
    file_list <- list.files(paste0(Data_dir,x), pattern = ".masslist")
    #for when you have masslists with modification (e.g. _C_1, _N_1) but want only to work with those that lack modification
    # file_list <- subset(file_list, grepl("_C_0", file_list) & grepl("_N_0", file_list))
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
    Data_list_new <- Data_list

    #Bind n_o and c_o into one data_frame and assign to each list entry
    # for(ind in 1:length(Data_list_new)){
    #  experiment <- names(Data_list_new)[ind]
    # Data_list_new[[experiment]] <- rbind2(Data_list_new[[c(ind, 1)]],Data_list_new[[c(ind, 2)]]) }

    names(Data_list_new) <- file_list #names for list entries
  }
  else
  {
    Data_dir <- getwd()
    file_list <- list.files(paste0(Data_dir,x), pattern = ".masslist")
    if (bound == TRUE)
    {
      file_list_0 <- subset(file_list, grepl("_C_0", file_list))
      file_list_1 <- subset(file_list, grepl("_C_1", file_list))
      file_list <- c(file_list_0, file_list_1)
    }
    else file_list <- subset(file_list, grepl("_C_0", file_list))
    Data_list_new <- by(seq_along(Data_list),cut(seq_along(Data_list),(length(Data_list)/2)), FUN=function(x)Data_list[x]) #make list of sublists containing n_o and c_o

    for(ind in 1:length(Data_list_new)){
      experiment <- names(Data_list_new)[ind]
      Data_list_new[[experiment]] <- rbind2(Data_list_new[[c(ind, 1)]],Data_list_new[[c(ind, 2)]]) }  #bind n_o and c_o into one data_frame and assign to each list entry

    names(Data_list_new) <- file_list #name the sublists
  }
  # get the sum of intensities of assigned peaks for both termini including both mod0 and mod1
  if (tint == T)
  {
    Data_dir <- getwd()
    TDL_out_dir <- paste0(Data_dir,"/")

    fl <- list.files(TDL_out_dir, pattern = ".masslist") # Get vector with all the files

    fl <- subset(fl, grepl("_C_0", fl) | grepl("_N_0", fl) | grepl("_C_1", fl) | grepl("_N_1", fl) | grepl("_C_2", fl) | grepl("_N_2", fl) | grepl("_C_3", fl) | grepl("_N_3", fl))

    DL <- lapply(fl, function(x) {read.table(paste0(TDL_out_dir, x),
                                             sep = "\t", stringsAsFactors = FALSE, header = TRUE) })

    cluster <- max(substring(regmatches(file_list,
                                        gregexpr("_C_[[:digit:]]|_N_[[:digit:]]", file_list)), 4, 5))

    cluster <- as.numeric(cluster)*2 + 2

    if (bound == TRUE)
    {

      DL <- by(seq_along(DL), cut(seq_along(DL),(length(DL)/cluster)), function(x)DL[x])
    }
    else DL <- by(seq_along(DL), cut(seq_along(DL),(length(DL)/2)), function(x)DL[x])
    fl <- subset(fl, grepl("_C_0", fl))

    names(DL) <- fl
    totalint <- c(1:length(DL))

    for (ind in 1:length(DL))
    {
      #discard all the unidentified peaks
      DL_T <- lapply(DL[[ind]], function (x) {subset(x, !is.na(x$Sequence.position))})
      #assign summarized intensities for experiment in one named vector
      totalint[ind] <- sum(sapply(DL_T, function(x) {sum(x$Intensity)}))
    }
    if(bound == TRUE & sep_term == FALSE)
    {
      totalint <- rep(totalint, cluster)
      names(totalint) <- file_list
    }
    else if(bound == TRUE & sep_term == TRUE)
    {
      totalint <- rep(totalint, cluster)
      names(totalint) <- file_list
    }
    else  names(totalint) <- file_list
  }
  if (!exists("data_merged")) {
    #Taking first entry of the list and working with it
    experiment <- names(Data_list_new)[1]
    #Getting rid of unassigned masses
    ass_peak_table <- Data_list_new[[experiment]]
    #Shaping the data frame with info from the entrie's name
    if (nrow(ass_peak_table) >0){
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

      for (i in 1:nrow(ass_peak_table))
      {
        if (ass_peak_table$Ion.type[i] == "A" | ass_peak_table$Ion.type[i] == "X") ass_peak_table$Ion.type.clust[i] <- "a/x"
        else if (ass_peak_table$Ion.type[i] == "B" | ass_peak_table$Ion.type[i] == "Y") ass_peak_table$Ion.type.clust[i] <- "b/y"
        else if (ass_peak_table$Ion.type[i] == "C" | ass_peak_table$Ion.type[i] == "Z") ass_peak_table$Ion.type.clust[i] <- "c/z"
        else if (ass_peak_table$Ion.type[i] == "D" | ass_peak_table$Ion.type[i] == "V" | ass_peak_table$Ion.type[i] == "W") ass_peak_table$Ion.type.clust[i] <- "d/v/w"
        else ass_peak_table$Ion.type.clust[i] <- ""
      }

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

      ass_peak_table$uvpd_hcd_pres_iso_nat <- paste0(ass_peak_table$energy," ",ass_peak_table$hcd," ",ass_peak_table$Pressure," ",ass_peak_table$iso, " ", ass_peak_table$nat.denat)
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
    else {data_merged <<- ass_peak_table}
  }
  if(!is.null(tint) & tint == TRUE & !is.null(nrow(ass_peak_table))) fragef_merged <<- fragef
  data_merged <<- ass_peak_table

  for (experiment in names(Data_list_new)[2:length(names(Data_list_new))])
  {
    ass_peak_table <- Data_list_new[[experiment]]

    if (nrow(ass_peak_table) > 0)
    {

      ass_peak_table$exp.name <- experiment[1]
      #Setting accuracy value, discarding everything that is less accurate

      if (nrow(ass_peak_table) >0){

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


        #peak intensity normalized by the sum of intensities of all peaks

        ass_peak_table$iso <- substr(experiment, regexpr("iso[0-9]+", experiment) +3,
                                     regexpr("iso[0-9]+", experiment) + attr(regexpr("iso[0-9]+", experiment), "match.length")-1)

        ass_peak_table$puls <- substr(experiment, regexpr("puls[0-9]", experiment) +4, regexpr("puls[0-9]", experiment)+4)

        ass_peak_table$enpuls <- paste0(ass_peak_table$energy, "x", ass_peak_table$puls)

        ass_peak_table$Pressure <- substr(experiment, regexpr("pres[0-9]",experiment)+4, regexpr("pres[0-9]", experiment)+4)

        ass_peak_table$nat.denat <- substr(experiment, regexpr("nat|denat", experiment),
                                           regexpr("nat|denat", experiment) + attr(regexpr("nat|denat", experiment), "match.length")-1)

        ass_peak_table$norm.intens <- ass_peak_table$Intensity/sum(ass_peak_table$Intensity)

        ass_peak_table$Log10_int <- log10(ass_peak_table$Intensity)

        ass_peak_table$log10.norm.intens <- log10(ass_peak_table$norm.intens)

        intens_mean <- mean(ass_peak_table$norm.intens)
        intens_sd <- sd(ass_peak_table$norm.intens)
        ass_peak_table$norm.intens.Zscored <- (ass_peak_table$norm.intens - intens_mean)/intens_sd

        ass_peak_table$seq.cov <- length(unique(ass_peak_table$Sequence.position))/(nchar(seq)-1)

        ass_peak_table$term <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+1)

        ass_peak_table$mod <- substr(experiment, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3, regexpr("_C_0|_N_0|_C_1|_N_1|_C_2|_N_2|_C_3|_N_3", experiment)+3)

        ass_peak_table$chg.iso <- paste0(ass_peak_table$Charge, " ", ass_peak_table$iso)

        ass_peak_table$rep <- substr(experiment, regexpr("rep[0-9]+", experiment) +3,
                                     regexpr("rep[0-9]+", experiment) + attr(regexpr("rep[0-9]+", experiment), "match.length")-1)

        #ass_peak_table$position_frag_numbers <- sapply(ass_peak_table$Sequence.position, function(x){sum(subset(ass_peak_table$norm.intens,ass_peak_table$Sequence.position==x)) })
        # ass_peak_table <- ass_peak_table[!duplicated(ass_peak_table$Sequence.position),]
        # ass_peak_table <- ass_peak_table[order(ass_peak_table$Sequence.position), ]
        # ass_peak_table$type.or.loss <- character(nrow(ass_peak_table))
        #for (ind in 1:nrow(ass_peak_table)){
        #   if(ass_peak_table$Mass.shift[ind] == "" | is.na(ass_peak_table$Mass.shift[ind])) ass_peak_table$type.or.loss[ind] <- ass_peak_table$Ion.type[ind]
        #   else if (ass_peak_table$Mass.shift[ind] == "Ammonia loss") ass_peak_table$type.or.loss[ind] <- "_NH3_Loss"
        #   else if (ass_peak_table$Mass.shift[ind] == "Water loss") ass_peak_table$type.or.loss[ind] <- "_H2O_Loss"
        #  else ass_peak_table$type.or.loss[ind] <- "WTF????"
        # }
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

        for (i in 1:nrow(ass_peak_table))
        {
          if (ass_peak_table$Ion.type[i] == "A" | ass_peak_table$Ion.type[i] == "X") ass_peak_table$Ion.type.clust[i] <- "a/x"
          else if (ass_peak_table$Ion.type[i] == "B" | ass_peak_table$Ion.type[i] == "Y") ass_peak_table$Ion.type.clust[i] <- "b/y"
          else if (ass_peak_table$Ion.type[i] == "C" | ass_peak_table$Ion.type[i] == "Z") ass_peak_table$Ion.type.clust[i] <- "c/z"
          else if (ass_peak_table$Ion.type[i] == "D" | ass_peak_table$Ion.type[i] == "V" | ass_peak_table$Ion.type[i] == "W") ass_peak_table$Ion.type.clust[i] <- "d/v/w"
          else ass_peak_table$Ion.type.clust[i] <- ""
        }

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

        ass_peak_table$uvpd_hcd_pres_iso_nat <- paste0(ass_peak_table$energy," ",ass_peak_table$hcd," ",ass_peak_table$Pressure," ",ass_peak_table$iso, " ", ass_peak_table$nat.denat)
        if(tint == TRUE)
        {
          fragef <- ass_peak_table %>% group_by(Sequence.position, Amino.acid) %>% summarise(fragef = sum(Intensity)/totalint[experiment])
          fragef$enpuls <- ass_peak_table$enpuls[1]
          fragef$etd <- ass_peak_table$etd[1]
          fragef$cid <- ass_peak_table$cid[1]
          fragef$hcd <- ass_peak_table$hcd[1]
          if(bound == TRUE) fragef$mod <- ass_peak_table$mod[1]
        }
        data_merged <<- merge(data_merged,ass_peak_table, all = TRUE, sort = FALSE)
        if(tint == TRUE & exists("fragef_merged") & nrow(ass_peak_table) > 0) fragef_merged <<- merge(fragef_merged, fragef, all = TRUE, sort = FALSE)
        else if (tint == TRUE & nrow(ass_peak_table) > 0) fragef_merged <<- fragef
      }
    }
  }
  print("Merged data frame completed")
}


#d <- seq_along(Data_list, 2)
stupannot <- function(subset = FALSE) {
  dir <- getwd()
  import_TDL_out()
  file_list <- list.files(dir)
  file_list <- subset(file_list, grepl(".masslist", file_list))

  cluster <- max(substring(regmatches(file_list,
                                      gregexpr("_C_[[:digit:]]|_N_[[:digit:]]", file_list)), 4, 5))

  cluster <- as.numeric(cluster)*2 + 2

  Data_list_new <- by(seq_along(Data_list),
                      cut(seq_along(Data_list),
                          (length(Data_list)/cluster)),
                      FUN=function(x)Data_list[x])


  names(Data_list_new) <- subset(file_list, grepl("_C_0", file_list))


  gList <- lapply(names(Data_list_new), function(x) {
          l  <- Data_list_new[[x]]

          c <- l[[1]]
          c1 <- l[[2]]
          c2 <- l[[3]]
          n <- l[[4]]
          n1 <- l[[5]]
          n2 <- l[[6]]


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
          df$exp <- x
          df
  })
  names(gList) <- subset(file_list, grepl("_C_0", file_list))

  gfrgsL <- lapply(names(gList), function(x) {
    cn  <- gList[[x]]

    cn_frgs <- cn[cn$Ion.type != "",]

    cn_frgs
  })

  df <- do.call("rbind", gList)
  fdf <- do.call("rbind", gfrgsL)
  fdf$mod <- as.character(fdf$mod)
  fdf[fdf$mod == "0","mod"] <- ""
  fdf[fdf$mod == "1","mod"] <- "P"
  fdf[fdf$mod == "2","mod"] <- "PP"
  df$mod <- as.character(df$mod)
  df[df$mod == "0","mod"] <- ""
  df[df$mod == "1","mod"] <- "P"
  df[df$mod == "2","mod"] <- "PP"
  df$calib.ppm <- df$Accuracy..PPM. - df$m.z.calibration
  fdf$calib.ppm <- fdf$Accuracy..PPM. - fdf$m.z.calibration
  df <- distinct(df, m.z, Intensity, exp)
  #fdf <- distinct(fdf, m.z, Intensity, exp, Sequence.position, Ion.type, Mass.shift, mod)
  #df$exp <- factor(df$exp, levels = c("sa20", "sa25"), labels = c("sa20","sa25"))
  df$exp <- regmatches(df$exp, gregexpr("hcd[[:digit:]]+_etd[[:digit:]]+_sa[[:digit:]]+", df$exp))
  fdf$exp <- regmatches(fdf$exp, gregexpr("hcd[[:digit:]]+_etd[[:digit:]]+_sa[[:digit:]]+", fdf$exp))
  fdf$exp <- unlist(fdf$exp)
  df$exp <- unlist(df$exp)

  for (i in 1:nrow(fdf))
  {
    if (fdf$Ion.type[i] == "A" | fdf$Ion.type[i] == "X") fdf$Ion.type.clust[i] <- "a/x"
    else if (fdf$Ion.type[i] == "B" | fdf$Ion.type[i] == "Y") fdf$Ion.type.clust[i] <- "b/y"
    else if (fdf$Ion.type[i] == "C" | fdf$Ion.type[i] == "Z") fdf$Ion.type.clust[i] <- "c/z"
    else if (fdf$Ion.type[i] == "D" | fdf$Ion.type[i] == "V" | fdf$Ion.type[i] == "W") fdf$Ion.type.clust[i] <- "d/v/w"
    else fdf$Ion.type.clust[i] <- ""
  }

  if(subset) {
    df <- df[df$exp == "hcd26_etd75_sa20" & df$m.z > 10000 & df$m.z < 18000,]
    fdf <- fdf[fdf$exp == "hcd26_etd75_sa20" & abs(fdf$calib.ppm) <= 3.5 & fdf$m.z > 10000 & fdf$m.z < 18000,]
    df$Intensity <- 100/max(df$Intensity)*df$Intensity
    fdf$Intensity <- 100/max(fdf$Intensity)*fdf$Intensity
  }

  ggplot(df) + theme_bw() +
    geom_bar(data = df, aes(x = m.z, y = Intensity), stat = "identity", color = "Black", size = .5) +
    geom_bar(stat = "identity", data = fdf[fdf$mod == "",],
            aes(x = m.z, y = Intensity, fill = Ion.type.clust, color = Ion.type.clust), size = 1, alpha = .5) +
    geom_bar(stat = "identity", data = fdf[fdf$mod != "",],
            aes(x = m.z, y = Intensity, fill = Ion.type.clust, color = Ion.type.clust), size = 1) +

    geom_hline(yintercept = 0) +
    #geom_label_repel(data = fdf[fdf$Intensity > 1000 & fdf$mod == "" & unlist(fdf$exp) == "hcd26_etd75_sa20" & abs(fdf$calib.ppm) < 2.5, ],
     #               min.segment.length = unit(0, "lines"),
      #              angle = 90, aes(y = Intensity, x = m.z, label = paste0(Ion.type, Sequence.position, Mass.shift), color = Ion.type.clust), size = 3) +
    #geom_label_repel(data = fdf[fdf$Intensity > 25 & fdf$mod == "P" | fdf$mod == "PP" & unlist(fdf$exp) == "hcd26_etd75_sa20" & abs(fdf$calib.ppm) < 2.5, ],
     #               min.segment.length = unit(0, "lines"),
      #              angle = 90, aes(y = Intensity, x = m.z, label = paste0(Ion.type, Sequence.position, mod), color = Ion.type.clust), size = 3) +
    geom_text_repel(data = fdf[fdf$mod != "",], segment.color = "green", angle = 90, nudge_y = 50,
              aes(y = Intensity, x = m.z, label = paste0(Ion.type, Series.number, mod, round(m.z, digits = 2),"+", Charge), color = Ion.type.clust), size = 2.5) +
    facet_grid(exp~.) +
    scale_y_continuous(limits = c(0, max(fdf$Intensity))) +
    scale_x_continuous(limits = c(14000, 18000), breaks = seq(14000, 18000, by = 500)) +
    theme(strip.text.x = element_text(size = 5))


}

plot_and_save <- function(plot_gg,Plot_dir_name,Plot_file_name, Plot_width, Plot_height){

  Data_dir <- getwd()

  print(plot_gg)
  print(paste0("Plot ", Plot_file_name, " generated"))
  Plot_dir_name <- "/plot/"
  pdf_file_name <- paste0(Data_dir, Plot_dir_name,Plot_file_name)

  pdf(file = paste0(pdf_file_name, ".pdf"), width = Plot_width, height = Plot_height)
  print(plot_gg)
  dev.off()
  print(paste0("Plot ", Plot_file_name, " saved to ", pdf_file_name))

}
