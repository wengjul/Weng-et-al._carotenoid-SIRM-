#### Set working directory == Select data for analysis ####
#either with this code or directly from R studio
setwd(choose.dir(
  default = "D:/..."))

#### Importing Data #### 
# Extract data designation from working directory name
working_directory <- getwd()
data_designation <- tail(strsplit(working_directory, "/")[[1]], n = 1)


# importing all data file names
temp = list.files(pattern="*.txt")
n_spectra <-length(temp)

# defining data frame for raw data (i.e. spectra)
spectra <- data.frame(read.delim(temp[1], header = TRUE, dec = ".")[2])
colnames(spectra) <- "x"

# importing spectra (y values)
for(i in temp) {
  temp_y <- read.delim(i, header = TRUE, dec = ".")[3]
  name_y <- paste("y_", i, sep = "")
  spectra[,name_y] <- temp_y
}

# Designate mean spectra: manually!

colnames(spectra) <- c("x", "0% 13C", "25% 13C", "50% 13C", "75% 13C", "100% 13C")

# Shifting spectra in y-axis for visualization: first -> top
# shift_factor has to be adjusted manually

shift_factor <- 0.7
for(i in seq(3, n_spectra +1)) {
  spectra[,i] <- spectra[,i] - (i - 2) * shift_factor
}




#### Visualization #### 

# Converting data from wide to long data frame -> required for ggplot2
# melt function from "data.table" is applied
library(data.table)

# full spectra
spectra_long <- melt(data = spectra, 
                            id.vars = "x",
                            variable.name = "Spectrum",
                            value.name = "y")

# set range lower_limit1 cm-1 - upper_limit1 cm-1

lower_limit1 <- 900
upper_limit1 <- 1600

spectra_long_select1 <- melt(data = subset(spectra, x >= lower_limit1 & x <= upper_limit1), 
                            id.vars = "x",
                            variable.name = "Spectrum",
                            value.name = "y")

# set range lower_limit2 cm-1 - upper_limit2 cm-1

lower_limit2 <- 2800
upper_limit2 <- 3100

spectra_long_select2 <- melt(data = subset(spectra, x >= lower_limit2 & x <= upper_limit2), 
                            id.vars = "x",
                            variable.name = "Spectrum",
                            value.name = "y")

# Plotting using ggplot2
library(ggplot2)

# Create Output folder for plots named "Figures"

dir.create(file.path(working_directory, "Figures"), showWarnings = FALSE)

# Get colors for visualization via ColorBrewer package
#install.packages("RColorBrewer")
library(RColorBrewer)
colors_means <- brewer.pal(n_spectra, "Dark2")

# or if red / blue color code activate
colors_means <- rev(brewer.pal(n_spectra, "RdBu"))
colors_means[3] <- "#7F7F7F"
colors_means <- rev(c("#b2182b", "#ef8a62", "#7F7F7F", "#67a9cf", "#2166ac"))

# or if all spectra black activate
#colors_means <- rep("black", n_spectra)





# mean spectra in one figure, range lower_limit1 cm-1 - upper_limit1 cm-1

ggplot(spectra_long_select1, aes(x, y, color = Spectrum)) +
  geom_line(size=1.3) +
  scale_color_manual(values = colors_means) +
  xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
  ylab("Raman Intensity / a.u.") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        text = element_text(family = "sans"),
        axis.ticks.length = unit(0.15, "cm")) +
  xlim(900, 1600) +
  scale_x_continuous(breaks= seq(900,1600,by=100), 
                     labels = c("",1000,"",1200,"",1400,"",1600)) +
  scale_y_continuous(breaks = NULL) +
  geom_vline(xintercept = 1522.91, color = "black", linetype = "longdash", alpha = 0.4, size=1) +
  geom_vline(xintercept = 1160.49, color = "black", linetype = "longdash", alpha = 0.4, size=1) +
  geom_vline(xintercept = 1008.92, color = "black", linetype = "longdash", alpha = 0.4, size=1) +
  annotate(geom = "text", x = 1360, y = 0.4, label = expression(paste("0% ", ""^13, C, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1360, y = -0.3, label = expression(paste("25% ", ""^13, C, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1360, y = -0.95, label = expression(paste("50% ", ""^13, C, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1360, y = -1.7, label = expression(paste("75% ", ""^13, C, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1360, y = -2.4, label = expression(paste("100% ", ""^13, C, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  
  #adjust y (and x) values, label and number of designations manually
  # assignments of wavenumbers
  annotate(geom = "text", x = 1510, y = 0.7, 
           label = expression(paste("1522 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1150, y = 0.7, 
           label = expression(paste("1160 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1000, y = 0.4, 
           label = expression(paste("1010 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1475, y = -2.15, 
           label = expression(paste("1488 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 1120, y = -2.15, 
           label = expression(paste("1131 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans") +
  annotate(geom = "text", x = 980, y = -2.35, 
           label = expression(paste("991 ", cm^-1, sep = "")), hjust = 1, vjust = 0, size = 6, family = "sans")
  

  
ggsave(filename = paste(data_designation, "_means_spectra_range ", lower_limit1, "-", upper_limit1, ".png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)

