library(scales)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
library(base)
library(plyr)


frqfls <- function(file = input$freqfile, Plot_dir_name = "plot", Plot_file_name = "freq_flyers_l6_mod", Plot_width = 20, Plot_height = 11){

# load the design
freqflyers <- read.table(file, quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
df <- freqflyers
df$Mz.diff <- df$Mz
df$Terminus <- as.factor(df$Terminus)
df$Count[df$Terminus == 'N-term'] <- -df$Count[df$Terminus == 'N-term']
df <- ddply(df, "Mz", transform, pos = Mz + 2)
df <- ddply(df, "Count", transform, ypos = Count + 3)
#if (df$Count < 0) {df <- ddply(df, "Count", transform, ypos = Count + 2)}
#else {df <- ddply(df, "Count", transform, ypos = Count - 2)}
df <- df[df$Count > 10 | df$Count < -10,]

plot_gg <- ggplot(df) +
  theme_bw() +
  ylab("Nterm Flyers                         Cterm Flyers") +
  geom_text(aes(y = ypos, label = round(Mz, 2), x = pos, alpha = .3) )  +
  guides(alpha = FALSE) +
  geom_point(aes(x=Mz.diff, y=Count, size=Intensity, alpha = .1))


print(plot_gg)
#plot_and_save(plot_gg, Plot_dir_name, Plot_file_name, Plot_width, Plot_height)

}

#df <- df[with(df, order(abs(Mz), abs(Count))), ]

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

create_plot_dir <- function(x){
  Data_dir <- getwd()
  # Generate the plot directory in it does not exist
  if(TRUE){
    plot_dir_path <- paste0(Data_dir, "/",x)
    # Check the existence of the folder and if no, create it
    if (!file_test("-d", plot_dir_path)){
      dir.create(plot_dir_path)
      print(paste0("Generated ",plot_dir_path))
    }
  }

}
