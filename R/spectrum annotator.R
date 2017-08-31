#s27
#seq <- "PLAKDLLHPSPEEEKRKHKKKRLVQSPNSYFMDVKCPGCYKITTVFSHAQTVVLCVGCSTVLCQPTGGKARLTEGCSFRRKQH"
#s27_9382 <- read.table("clipboard", header = TRUE, skip = 6)

#l40
#seq <- "IIEPSLRQLAQKYNCDKMICRKCYARLHPRAVNCRKKKCGHTNNLRPKKKVK"
#writeClipboard(as.character(s27_9382$Mass))

library(ggrepel)
library(RColorBrewer)

#' Spectrum annotation
#'
#' @param dir2
#' @param dir3
#' @param subset
#' @param mergefrags
#' @param Plot_file_name
#' @param data
#' @param ranges
#'
#' @return
#' @export
#'
#' @examples
annotate_spectrum <- function(data = gList, dir2 = "", dir3 = "",
                      subset = FALSE, mergefrags = FALSE, Plot_file_name = "annotated spectrum", ranges = ranges) {

  if(mergefrags == TRUE & dir2 != "" & dir3 != "")
  {
    gList2 <- prepare_masslists(dir = "K:/mgf/topDown/l40/1ss_IIEdel")
    gList3 <- prepare_masslists(dir = "K:/mgf/topDown/l40/1ss")
  }


  gfrgsL <- lapply(names(gList), function(x) {
    cn  <- gList[[x]]


    cn_frgs <- cn[cn$Ion.type != "",]
    if(mergefrags == TRUE)
    {
      cn2 <- gList2[[x]]
      cn3 <- gList3[[x]]
      cn_frgs2 <- cn2[cn2$Ion.type != "",]
      cn_frgs3 <- cn3[cn3$Ion.type != "",]
      cn_frgs_extra <- bind_rows(anti_join(cn_frgs2, cn_frgs3, by = "m.z"),
                                 anti_join(cn_frgs3, cn_frgs2, by = "m.z"))
      cn_frgs <- bind_rows(cn_frgs, anti_join(cn_frgs_extra, cn_frgs, by = "m.z"))
    }
    cn_frgs
  })


  df <- do.call("rbind", gList)
  fdf <- do.call("rbind", gfrgsL)
  # fdf$mod <- as.character(fdf$mod)
  # fdf[fdf$mod == "0","mod"] <- ""
  # fdf[fdf$mod == "1","mod"] <- "P"
  # fdf[fdf$mod == "2","mod"] <- "PP"
  # df$mod <- as.character(df$mod)
  # df[df$mod == "0","mod"] <- ""
  # df[df$mod == "1","mod"] <- "P"
  # df[df$mod == "2","mod"] <- "PP"
  df$calib.ppm <- df$Accuracy..PPM. - df$m.z.calibration
  fdf$calib.ppm <- fdf$Accuracy..PPM. - fdf$m.z.calibration
  df <- distinct(df, m.z, Intensity, exp)
  #fdf <- distinct(fdf, m.z, Intensity, exp, Sequence.position, Ion.type, Mass.shift, mod)
  #df$exp <- factor(df$exp, levels = c("sa20", "sa25"), labels = c("sa20","sa25"))
  #df$exp <- regmatches(df$exp, gregexpr("hcd[[:digit:]]+_etd[[:digit:]]+_sa[[:digit:]]+", df$exp))
  #fdf$exp <- regmatches(fdf$exp, gregexpr("hcd[[:digit:]]+_etd[[:digit:]]+_sa[[:digit:]]+", fdf$exp))
  #fdf$exp <- unlist(fdf$exp)
  #df$exp <- unlist(df$exp)

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
  df$iso <- round(as.numeric(gsub('.*_([[:digit:]]+.[[:digit:]]+).txt.*','\\1',df$exp)),2)
  fdf$iso <- round(as.numeric(gsub('.*_([[:digit:]]+.[[:digit:]]+).txt.*','\\1',fdf$exp)),2)
  fdf$mod <- as.factor(fdf$mod)
  plot_gg <- ggplot(df) + theme_bw() +
    geom_bar(aes(x = m.z, y = Intensity), stat = "identity", color = "Black", size = .5) +
    geom_bar(stat = "identity", data = fdf[fdf$mod == "",],
             aes(x = m.z, y = Intensity, fill = mod, color = mod), size = .5, alpha = .5) +
    geom_bar(stat = "identity", data = fdf[fdf$mod != "",],
             aes(x = m.z, y = Intensity, fill = mod, color = mod), size = .5) +

    geom_hline(yintercept = 0) +
    #geom_label_repel(data = fdf[fdf$Intensity > 1000 & fdf$mod == "" & unlist(fdf$exp) == "hcd26_etd75_sa20" & abs(fdf$calib.ppm) < 2.5, ],
    #               min.segment.length = unit(0, "lines"),
    #              angle = 90, aes(y = Intensity, x = m.z, label = paste0(Ion.type, Sequence.position, Mass.shift), color = Ion.type.clust), size = 3) +
    #geom_label_repel(data = fdf[fdf$Intensity > 25 & fdf$mod == "P" | fdf$mod == "PP" & unlist(fdf$exp) == "hcd26_etd75_sa20" & abs(fdf$calib.ppm) < 2.5, ],
    #               min.segment.length = unit(0, "lines"),
    #              angle = 90, aes(y = Intensity, x = m.z, label = paste0(Ion.type, Sequence.position, mod), color = Ion.type.clust), size = 3) +
    geom_text_repel(data = fdf[fdf$Intensity > 10,], force = 2, segment.color = "green", angle = 90, nudge_y = 50,
                    aes(y = Intensity, x = m.z, label = paste0(Ion.type, Series.number), alpha = Intensity, color = mod), size = 2.5) +
    geom_text_repel(data = df[df$Intensity > 10 & df$m.z < 17000,], aes(y = Intensity, x = m.z, label = round(m.z,2), alpha = Intensity), size = 2.5, angle = 90) +
    facet_grid(exp~., scales = "free") +
    scale_y_continuous(limits = ranges$y) +
    scale_x_continuous(limits = ranges$x) +
    theme(strip.text.x = element_text(size = 5))


  Plot_dir_name <- "plot"
  #plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width = 15, Plot_height = 10)
  #save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))

  print(plot_gg)
}

#rm(ass_data_merged)
#seq <- "MKSVITTVVSAADAAGRFPSNSDLESIQGNIQRSAARLEAAEKLAGNHEAVVKEAGDACFAKYAYLKNPGEAGENQEKINKCYRDVDHYMRLVNYCLVVGGTGPLDEWGIAGAREVYRTLNLPTSAYVASIAYTRDRLCVPRDMSAQAGVEFSAYLDYLINALS"

#Data_list <- import_masslists(bound = TRUE, mod.number = 2)
#prepare_data_frame_term_sep_aa(bound = TRUE, sep_term = TRUE, mod.number = 2, ppm = 1)

#' Fragment Map Plotting function
#'
#' @param data
#' @param subset
#' @param uvhcpd
#' @param Plot_title_suffix
#' @param Plot_dir_name
#' @param Plot_file_name
#' @param Plot_width
#' @param Plot_height
#'
#' @return
#' @export
#'
#' @examples
fragment_map <- function(data = data_merged, subset = FALSE, uvhcpd = FALSE, Plot_title_suffix = "", Plot_dir_name = "plot/", Plot_file_name = "pin1_6al7_2phos_sty", Plot_width = 3, Plot_height = 10, seq = seq) {

  #Data_dir <- getwd()

  #data$energy <- as.numeric(data$energy)

  #if(subset == TRUE) data <- subset(data, abs(data$calib.ppm) <= 1.0)
  if(uvhcpd == TRUE) ass_data_merged$uvhcpd <- paste0(ass_data_merged$hcd,"Vx",ass_data_merged$energy, "mJ")
  cols <- c('#e31a1c','#33a02c','#6a3d9a','#e31a1c','#33a02c',
            '#6a3d9a')
  cols1 <- c('#000000', "#7f0000", '#d7301f','#fed976')
  #df$energy[is.na(df$energy)] <- 0

  s <- strsplit(seq, "")
  sty <- which(s[[1]] == "S" | s[[1]] == "Y" | s[[1]] == "T")
  rk <- which(s[[1]] == "R" | s[[1]] == "K" | s[[1]] == "H")
  c <- which(s[[1]] == "C")
  sty.df <- data.frame(Amino.acid = s[[1]][sty], Sequence.position = sty)
  rk.df <- data.frame(Amino.acid = s[[1]][rk], Sequence.position = rk)
  s.df <- data.frame(Amino.acid = s[[1]], Sequence.position = seq(1,length(s[[1]])))
  c.df <- data.frame(Amino.acid = s[[1]][c], Sequence.position = c)
  #myPalette <- colorRampPalette(hotmetal())
  data$term <- as.factor(data$term)
  levels(data$term) <- c("C-terminal fragments", "N-terminal fragments")
  #data$exp <- round(as.numeric(gsub('.*_([[:digit:]]+.[[:digit:]]+).txt.*','\\1',data$exp)),0)

  plot_gg <- ggplot(data) +
    stat_identity(geom = "tile", width  = 1, aes(x = Sequence.position, y = mod, fill = term, alpha = norm.intens.Zscored)) +

    geom_hline(yintercept = 0.5) +
    geom_hline(yintercept = 1.5) +
    geom_hline(yintercept = 2.5) +
    geom_hline(yintercept = 3.5) +

    geom_vline(xintercept = c.df$Sequence.position-1, color = "red", size = 1) +
    #geom_histogram(aes(color = Ion.type)) +
    #geom_linerange (aes(x = Sequence.position, y = norm.intens), position = position_dodge()) +
    #stat_identity(aes(color = Ion.type)) +
    #scale_y_discrete( expand = c(0,0), limits= c(1,1)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, nchar(seq)), breaks = seq(0, nchar(seq), by = 10)) +
    #scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    #facet_grid(exp ~ .) +
    #scale_fill_gradient2(limits = c(0, max(ass_data_merged$norm.intens.Zscored)), name = "Intensity", high = "black", mid = "red", low = "yellow", midpoint = median(ass_data_merged$norm.intens.Zscored), guide = "legend") +
    #correct gradient!!!!!
    #scale_fill_gradientn(colours = myPalette(100),name = "Z-scored\nintensity") +
    scale_fill_manual(values = c("blue", "red")) +
    ylab("Number of PEB") + # x label
    xlab("Sequence position") +  # y label
    #ggtitle(paste0("Sequence Coverage\n", Plot_title_suffix)) +
    #geom_text(aes(y = 0.5, label = Amino.acid, x = Sequence.position), size = 2)  +
    scale_y_discrete(expand = c(-0.7, .7)) +
    #coord_flip() +
    #guides(alpha = FALSE, fill = FALSE) +

    geom_text(data = c.df, aes(y = 3.6, x = Sequence.position-1, label = Amino.acid), size = 2) +
    #geom_text(data = rk.df, aes(y = 0.4, x = Sequence.position, label = Amino.acid), size = 2) +
    #geom_text(aes(y=pos,label=paste0(as.integer(seq.cov), "%")), colour = "black", size=5, angle = 60, position = position_dodge(width = .55)) +
    theme(legend.position="right", title=element_text(face = "bold"), text = element_text(size=12),
          axis.text.x = element_text()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



  #Plot_dir_name <- "plot"
  #plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width = 10, Plot_height = 4)
  #save(list = ls(), file = paste0(Data_dir,"/", Sys.Date(),"_", Plot_file_name,".RData"))
  plot_gg
}

