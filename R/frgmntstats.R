
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

#' Fragment Charge Plot
#'
#' @param uvhcpd
#' @param subset
#' @param Plot_title_suffix
#' @param Plot_dir_name
#' @param Plot_file_name
#' @param Plot_width
#' @param Plot_height
#' @param iso
#' @param pep
#'
#' @return
#' @export
#'
#' @examples
fragment_charges <- function(data = data_merged, seq = seq) {

  data <- data[order(data$term, decreasing = TRUE),]


  #sequence <- unlist(strsplit(seq, split = ""))
  #names(sequence) <- c(1:nchar(seq))
  plot_gg <- ggplot(data, aes(x = Sequence.position, y = Charge)) +
    geom_point(stat = "identity", aes(shape = Ion.type.clust, color = mod, size = Intensity)) +
    scale_y_continuous(breaks = seq(0,ceiling(max(as.numeric(data$Charge))), by = 1), limits = c(0.5, ceiling(max(as.numeric(data$Charge))))) +
    scale_x_continuous(breaks = seq(0,length(strsplit(seq, split = "")[[1]]), by = 10)) +
    scale_shape_discrete(name = "Ion Type\nSeries") +
    scale_color_discrete(name = "Number of\nPhoshosites") +
    guides(size = FALSE) +
    #geom_text(aes(y = 0.7, label = Amino.acid, x = Sequence.position), size = 1.8)  +
    #facet_grid(mod ~ .) +
    theme_bw()


  plot_gg
}


#' Terminal Sequence coverages
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
termini_sequence_coverage <- function(data = data_merged, ppm = 3.5){

  seq_cov_data <- data[!duplicated(data$exp),]

  #Determine values for text plotting
  seq_cov_data <- ddply(seq_cov_data,"seq.cov", transform, perc_cov = seq.cov*100)
  seq_cov_data$perc_cov <- as.integer(seq_cov_data$perc_cov)
  seq_cov_data$perc_cov_y <- seq_cov_data$perc_cov + 5

  cols <- c('#9e9ac8','#756bb1','#54278f')

  plot_gg <- ggplot(seq_cov_data[seq_cov_data$calib_ppm <= 3.5, ], aes(x = energy, y = perc_cov)) +
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

  plot_gg

}


#' Fragment Types Plotted
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fragment_types <- function(data = data_merged){
  #plotting of by xa and zc ions stacked horizontally and with hcd facetting (count labelling for each pair of frags)
  require(plyr)


  cols <- c('#ef8a62','#91cf60','#67a9cf')
  # data$mod <- factor(data$mod,
  #                    levels = c(0, 1, 2),
  #                    labels = c("No modification", "1 phospho site", "2 phospho sites"))

  plot_gg <- ggplot(data, aes(x = frag.type, y = norm.intens.Zscored)) +
    geom_boxplot(aes(fill = Ion.type.clust), size = .5) +
    #stat_params = list(binwidth = ) +
    scale_y_continuous(expand = c(0,0)) +

    #scale_x_continuous("energy", breaks = 0.5:2.0) + #define isolated m/z
    #geom_text(aes(y = pos, label = count.freq), colour = "black", size=3) +

    facet_grid(mod ~ .) +
    ggtitle(paste0("Fragment Ion Types\n")) + # plot title
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

  plot_gg


}
