# Required Data Structure: spectra as .txt files in folder, 
# containing x [cm-1] and y [a.u.] column

# Names of Directories may not include special characters not suitable for
# designation of .png files

#### Select data folders / sub-directories that should be analyzed ####

# Select directory containing data folders / sub-directories 

directory <- choose.dir(
  default = "D:/...")

# list all available data folders / sub-directories in a directory
# pick a directory here (not sub directories)

sub_directories <- list.dirs(path = directory, 
                             full.names = TRUE, recursive = FALSE)

# interactively pick desired data folders / sub-directories for analysis
sub_directories <- select.list(sub_directories, 
                               graphics = TRUE, multiple = TRUE)

# Create Output folder in directory for export of mean spectra 
# named "Mean Spectra"

dir.create(file.path(directory, "Mean Spectra"), showWarnings = FALSE)

#### Calculate shift in Wavenumber axis with silicium wafer spectrum####
# Import spectrum of silicium wafer: required for shift correction
setwd(directory)
Si_waver_name <- list.files(path = "Si waver", pattern="*.txt")
Si_waver_spectrum <- data.frame(read.delim(file.path("Si waver", Si_waver_name), header = FALSE, dec = "."))
colnames(Si_waver_spectrum) <- c("x", "y")

# Cut to region of interest
Si_waver_spectrum <- subset(Si_waver_spectrum, x >= 500 & x <= 550)

# Calculate shift to standard value
# fit Gausian on Si signal around 523 cm-1 to get exact signal position
fitG <- function(x,y,mu,sig,scale){
  f = function(p){
    d = p[3]*dnorm(x,mean=p[1],sd=p[2])
    sum((d-y)^2)
  }
  optim(c(mu,sig,scale),f)
}

Si_signal_position <- fitG(Si_waver_spectrum$x, Si_waver_spectrum$y, 523, 15, 1)$par[1]

# Calculate spectral shift in cm-1
Si_signal_position_standard <- 524.5
wavenumber_shift <- Si_signal_position_standard - Si_signal_position

#### Analyze selected sub-directories each by each within loop ####

for (sub_directory in sub_directories) {
  

  #### Set working directory == Select data for analysis ####
  # either with this code or directly from R studio
  setwd(sub_directory)
  
  #### Importing Data #### 
  # Extract data designation from working directory name
  working_directory <- getwd()
  data_designation <- tail(strsplit(working_directory, "/")[[1]], n = 1)
  
  # importing all data file names
  temp = list.files(pattern="*.txt")
  
  # defining data frame for raw data (i.e. spectra)
  spectra <- data.frame(read.delim(temp[1], header = FALSE, dec = ".")[1])
  colnames(spectra) <- "x"
  
  # importing spectra (y values)
  for(i in temp) {
    temp_y <- read.delim(i, header = FALSE, dec = ".")[2]
    name_y <- paste("y_", i, sep = "")
    spectra[,name_y] <- temp_y
  }
  
  
  #### Data Pre-Processing #### 
  
  # correct for spectral shift
  
  spectra$x <- round((spectra$x + wavenumber_shift), 2)
  
  # Cropping to range selected in lower and upper limit (in cm-1)
  lower_limit <- 500
  upper_limit <- 3200
  spectra_crop <- subset(spectra, x >= lower_limit & x <= upper_limit)
  
  
  
  # Smoothing: Savitsky-Golay with p:'filter order' (same result for 2+3 and 4+5) 
  #DEACTIVATED
  # n:'filter length' (i.e. window size), m:'derivative'
  #install.packages('signal')
  library(signal)
  
  spectra_sgolay <- spectra_crop
  
  #cbind(spectra_crop[1] ,
  #                     apply(spectra_crop[-1], 2, sgolayfilt, p = 3, n = 9, m = 0))
  
  
  
  # Background correction using "baseline" package
  
  #install.packages("baseline")
  library(baseline)
  library(data.table)
  library(tibble)
  
  # Organizing data as row vectors in data frame (without x) 
  # by using "data.table" package -> required for "baseline" package
  spectra_trans <- as.matrix(transpose(spectra_sgolay[,-1]))
  rownames(spectra_trans) <- colnames(spectra_sgolay[,-1])
  colnames(spectra_trans) <- rownames(spectra_sgolay[,-1])
  
  # Rolling Ball algorithm applied
  spectra_base <- baseline(spectra_trans, 
                                wm=60, ws=5, method='rollingBall')
  
  spectra_base <- transpose(as.data.frame(spectra_base@corrected))
  spectra_base <- add_column(spectra_base, x = spectra_sgolay$x, .before = 1)
  colnames(spectra_base) <- colnames(spectra_sgolay)
  
  
  
  
  # Scaling/Normalization
  # Min-Max Scaling
  #install.packages("caret")
  library(caret)
  library(tibble)
  preProc <- preProcess(spectra_base[,-1], method = "range")
  spectra_base_norm <- predict(preProc, spectra_base[,-1])
  spectra_base_norm <- add_column(spectra_base_norm, x = spectra_crop$x, .before = 1)
  
  
  
  
  # Mean spectrum calculation
  
  # number of spectra considered
  n_mean_spectrum <- ncol(spectra_base_norm) - 1
  
  # mean calculation
  mean_base_norm <- as.data.frame(rowMeans(spectra_base_norm[,-1]))
  mean_base_norm <- add_column(mean_base_norm, x = spectra_crop$x, .before = 1)
  colnames(mean_base_norm) <- c("x", "y_mean")
  
  
  
  
  # Peak Picking based on local maxima
  
  #install.packages("IDPmisc")
  library(IDPmisc)
  
  # mean
  
  picked_peaks_mean <- peaks(x = mean_base_norm$x, y = mean_base_norm$y_mean, thr = 0.30, minPH = 0.025)
  picked_peaks_mean <- subset(picked_peaks_mean, x >= 700)
  
  
  
  
  
  
  
  #### Visualization #### 
  
  # Converting data from wide to long data frame -> required for ggplot2
  # melt function from "data.table" is applied
  library(data.table)
  
  # convert for for single spectra facet plot
  spectra_base_norm_long <- melt(data = spectra_base_norm, 
                              id.vars = "x",
                              variable.name = "Spectrum",
                              value.name = "y")
  
  # convert for single spectra facet plot 1400 cm-1 - 1600 cm-1
  
  spectra_base_norm_long_v1 <- melt(data = spectra_base_norm[340:430,], 
                                 id.vars = "x",
                                 variable.name = "Spectrum",
                                 value.name = "y")
  
  
  # convert for single spectra facet plot 1050 cm-1 - 1200 cm-1
  
  spectra_base_norm_long_v2 <- melt(data = spectra_base_norm[205:265,], 
                                    id.vars = "x",
                                    variable.name = "Spectrum",
                                    value.name = "y")
  
  # convert for single spectra facet plot 800 cm-1 - 1050 cm-1
  
  spectra_base_norm_long_v3 <- melt(data = spectra_base_norm[100:205,], 
                                    id.vars = "x",
                                    variable.name = "Spectrum",
                                    value.name = "y")
  
  # convert for single spectra facet plot 800 cm-1 - 1600 cm-1
  
  spectra_base_norm_long_v123 <- melt(data = spectra_base_norm[100:430,], 
                                    id.vars = "x",
                                    variable.name = "Spectrum",
                                    value.name = "y")
  
  # Plotting using ggplot2
  library(ggplot2)
  
  # Create Output folder for plots named "Figures"
  
  dir.create(file.path(working_directory, "Figures"), showWarnings = FALSE)
  
  
  # mean spectrum
  
  ggplot(mean_base_norm, aes(x, y_mean)) +
    geom_line() +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    geom_text(data = picked_peaks_mean, aes(x , y, label = x), 
              check_overlap = TRUE, nudge_y = 0.05) +
    annotate(geom = "text", 
             x = range(mean_base_norm$x)[2], y = range(mean_base_norm$y_mean)[2],
             label = paste("n = ", as.character(n_mean_spectrum), sep = ""), 
             hjust = 1, vjust = 0) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL)
    
  
  
  ggsave(filename = paste(data_designation, "_mean.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  
  
  # mean spectrum and shaded individual spectra
  
  ggplot(mean_base_norm, aes(x, y_mean)) +
    geom_line() +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    geom_text(data = picked_peaks_mean, aes(x , y, label = x), 
              check_overlap = FALSE, nudge_y = 0.05) +
    annotate(geom = "text", 
             x = range(mean_base_norm$x)[2], y = range(mean_base_norm$y_mean)[2],
             label = paste("n = ", as.character(n_mean_spectrum), sep = ""), 
             hjust = 1, vjust = 0) +
    geom_line(data = spectra_base_norm_long, aes(x, y), alpha = 0.1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL)
  
  
  ggsave(filename = paste(data_designation, "_mean_shaded ind.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  
  
  # single spectra facet plot
  
  ggplot(spectra_base_norm_long, aes(x, y)) +
    geom_line() +
    facet_wrap(~Spectrum, ncol = 10) +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks=c(1000, 3000))
  
  ggsave(filename = paste(data_designation, "_single spectra.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  
  
  # single spectra facet plot 1400 cm-1 - 1600 cm-1
  
  ggplot(spectra_base_norm_long_v1, aes(x, y)) +
    geom_line() +
    facet_wrap(~Spectrum, ncol = 10) +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks=c(1425, 1575)) +
  #  geom_vline(xintercept = 1485.59, color = "blue") +
    geom_vline(xintercept = 1522.1, color = "red")
  
  ggsave(filename = paste(data_designation, "_single spectra_range 1400-1600.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  # single spectra facet plot 1050 cm-1 - 1200 cm-1
  
  ggplot(spectra_base_norm_long_v2, aes(x, y)) +
    geom_line() +
    facet_wrap(~Spectrum, ncol = 10) +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks=c(1075, 1175)) +
    #  geom_vline(xintercept = 1485.59, color = "blue") +
    geom_vline(xintercept = 1157, color = "red")
  
  ggsave(filename = paste(data_designation, "_single spectra_range 1050-1200.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  
  # single spectra facet plot 800 cm-1 - 1050 cm-1
  
  ggplot(spectra_base_norm_long_v3, aes(x, y)) +
    geom_line() +
    facet_wrap(~Spectrum, ncol = 10) +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks=c(850, 1000)) +
    #  geom_vline(xintercept = 1485.59, color = "blue") +
    geom_vline(xintercept = 1005, color = "red")
  
  ggsave(filename = paste(data_designation, "_single spectra_range 800-1050.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  # single spectra facet plot 800 cm-1 - 1600 cm-1
  
  ggplot(spectra_base_norm_long_v123, aes(x, y)) +
    geom_line() +
    facet_wrap(~Spectrum, ncol = 10) +
    xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
    ylab("Raman Intensity / a.u.") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks=c(1000, 1400))
  
  ggsave(filename = paste(data_designation, "_single spectra_range 800-1600.png", sep = ""),
         path = "Figures",
         width = 20, height = 15, units = "cm",
         dpi = 600)
  
  #### Export data #### 
  # Mean spectrum into folder "Mean Spectra" in directory
  
  write.table(mean_base_norm, 
              file = file.path(directory, "Mean Spectra",
              paste(data_designation, "_mean spectrum.txt", sep = "")),
              sep = "\t",
              row.names = TRUE, col.names = NA)
}




