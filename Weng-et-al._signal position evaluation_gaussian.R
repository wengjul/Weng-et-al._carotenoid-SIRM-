# Required Data Structure: spectra as .txt files in folder, 
# containing x [cm-1] and y [a.u.] column

#### Select PCA data folders that should be analyzed (with grouped sub-directories)####
# Select directory containing data folders / sub-directories 

directory <- choose.dir(
  default = "D:/...")

# Set working directory
setwd(directory)
# Extract data designation from working directory name
working_directory <- getwd()
data_designation <- tail(strsplit(working_directory, "/")[[1]], n = 1)


sub_directories <- list.dirs(path = directory, 
                             full.names = TRUE, recursive = FALSE)

# Interactively pick desired data folders / sub-directories for analysis
sub_directories <- select.list(sub_directories, 
                               graphics = TRUE, multiple = TRUE)


# Defining data frame for raw data (i.e. spectra)
temp <- list.files(path = sub_directories[1], pattern = "*.txt")
spectra <- data.frame(read.delim(file.path(sub_directories[1], temp[1]), header = FALSE, dec = ".")[1])
colnames(spectra) <- "x"

# Initialize a vector to store number of observations in each group
n_groups <- c()

# Importing Spectra by group (single subdirectories)

for (sub_directory in sub_directories) {
  # importing all data file names of group / subdirectory
  temp <- list.files(path = sub_directory, pattern = "*.txt")
  
  # append number of observations in the group
  n_groups <- append(n_groups, length(temp))
  
  # importing spectra (y values)
  for(i in temp) {
    temp_y <- read.delim(file.path(sub_directory, i), header = FALSE, dec = ".")[2]
    name_y <- paste("y_", i, sep = "")
    spectra[,name_y] <- temp_y
  }
  
}

#### Data Pre-Processing #### 


# Cropping to range selected in lower and upper limit (in cm-1)
# Carotenoids: 900 -1600 cm-1
# Full cell: 700 - 1800 cm-1
lower_limit <- 900
upper_limit <- 1800
spectra_crop <- subset(spectra, x >= lower_limit & x <= upper_limit)



# Smoothing: Savitsky-Golay with p:'filter order' (same results for 2+3 and 4+5)
#DEACTIVATED
# n:'filter length'=window size, m:'derivative'
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

# Organizing data as row vectors in data frame (without x) by using "data.table"
# package -> required for "baseline" 
spectra_trans <- as.matrix(transpose(spectra_sgolay[,-1]))
rownames(spectra_trans) <- colnames(spectra_sgolay[,-1])
colnames(spectra_trans) <- rownames(spectra_sgolay[,-1])

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



#### Signal Position Fitting ####
# Cut to region of interest
spectra_fit_crop <- subset(spectra_base_norm, x >= 1400 & x <= 1600)

# function: fit Gaussian on signal to get exact signal position
fitgaus <- function(x,y,mu,sig,scale,y_offset){
  f = function(p){
    d = p[4]+abs(p[3])*dnorm(x,mean=p[1],sd=p[2])
    sum((d-y)^2)
  }
  optim(c(mu,sig,scale,y_offset),f)
}


# Initialize vectors to store estimated parameters of fits

fit_mu <- c()
fit_sig <- c()
fit_scale <- c()
fit_y_offset <- c()

# define function used for fitting
fitfct <- fitgaus
# execute fit function, save results for parameters in list, insert starting parameters

for(i in seq(2, length(spectra_fit_crop))) {
  fit_mu <- append(fit_mu, 
                   fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 1500, 15, 1, 0.2)$par[1])
  fit_sig <- append(fit_sig, 
                    fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 1500, 15, 1, 0.2)$par[2])
  fit_scale <- append(fit_scale, 
                      fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 1500, 15, 1, 0.2)$par[3])
  fit_y_offset <- append(fit_y_offset, 
                         fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 1500, 15, 1, 0.2)$par[4])
}



# generate vector that stores group labels of each observation

group_labels <- c(0, 25, 50, 75, 100)

group_label_vector <- c()

for (i in seq(length(n_groups))) {
  group_label_vector <- append(group_label_vector, 
                               rep(group_labels[i], n_groups[i]))
}

signal_fit_par <- as.data.frame(
  cbind("group" = group_label_vector, "mean" = fit_mu, "sd" = fit_sig, "amplitude" = fit_scale, "offset" = fit_y_offset))
signal_fit_par$group <- as.factor(signal_fit_par$group)


# calculate mean and sd of signal position grouped by isotope content

signal_positions_stat <- aggregate(signal_fit_par$mean, list(signal_fit_par$group), FUN = mean)
signal_positions_stat <- cbind(signal_positions_stat, 
                               aggregate(signal_fit_par$mean, list(signal_fit_par$group), FUN = sd)[2])
colnames(signal_positions_stat) <- c("x", "y", "sd")

signal_fit <- data.frame("x" = spectra_fit_crop$x)

for (i in seq(length(signal_fit_par$mean))) {
  signal_fit <- cbind(signal_fit, 
                       signal_fit_par[i,5]+abs(signal_fit_par[i,4])*
                         dnorm(x=spectra_fit_crop$x, mean=signal_fit_par[i,2], sd=signal_fit_par[i,3]))
}
colnames(signal_fit) <- colnames(spectra_fit_crop)

# calculate residuals of fits

residuals <- spectra_fit_crop - signal_fit
residuals$x <- spectra_fit_crop$x


# calculation of FWHM


y_max <- sapply(signal_fit[,-1], max)
x_max <- c()

for (i in seq(length(y_max))) {
  x_max <- append(x_max, 
                signal_fit$x[signal_fit[,i+1]==y_max[i]])
}




half_max <- c()

for (i in seq(length(x_max))) {
  
  half_max <- append(half_max, 
                     signal_fit[,i+1][signal_fit$x==x_max[i]]/2)
}


x_half_max_low <- c()

for (i in seq(length(half_max))) {
  
  x_half_max_low <- append(x_half_max_low, 
                           signal_fit$x[signal_fit$x < x_max[i]][which.min
                                                                 (abs(signal_fit[,i+1][signal_fit$x < x_max[i]]-half_max[i]))])
}

x_half_max_high <- c()

for (i in seq(length(half_max))) {
  
  x_half_max_high <- append(x_half_max_high, 
                           signal_fit$x[signal_fit$x > x_max[i]][which.min
                                                                 (abs(signal_fit[,i+1][signal_fit$x > x_max[i]]-half_max[i]))])
}



FWHM <- as.data.frame(
  cbind("group" = group_label_vector, "FWHM" = x_half_max_high - x_half_max_low))
FWHM$group <- as.factor(FWHM$group)






#### Visualization ####

# Create Output folder for plots named "Figures"

dir.create(file.path(working_directory, "Figures"), showWarnings = FALSE)



# single spectra facet plot

# Converting data from wide to long data frame -> required for ggplot2
# melt function from "data.table" is applied
library(data.table)

# convert for for single spectra facet plot
spectra_fit_crop_long <- melt(data = spectra_fit_crop, 
                               id.vars = "x",
                               variable.name = "Spectrum",
                               value.name = "y")
signal_fit_long <- melt(data = signal_fit, 
                              id.vars = "x",
                              variable.name = "Spectrum",
                              value.name = "y")

ggplot(spectra_fit_crop_long, aes(x, y)) +
  geom_line() +
  facet_wrap(~Spectrum, ncol = 14) +
  xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
  ylab("Raman Intensity / a.u.") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        text = element_text(family = "sans")) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks=c(1450, 1550)) +
  geom_line(data = signal_fit_long, aes(x, y), color = "darkred",size=1.15)

ggsave(filename = paste(data_designation, "_fit_facet_plot.png", sep = ""),
       path = "Figures",
       width = 40, height = 30, units = "cm",
       dpi = 600)

# single spectra facet plot: residuals

residuals_long <- melt(data = residuals, 
                        id.vars = "x",
                        variable.name = "Spectrum",
                        value.name = "y")

ggplot(residuals_long, aes(x, y)) +
  geom_line() +
  facet_wrap(~Spectrum, ncol = 14) +
  xlab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
  ylab("Raman Intensity / a.u.") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        text = element_text(family = "sans")) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks=c(1450, 1550))

ggsave(filename = paste(data_designation, "_fit_facet_plot_residuals.png", sep = ""),
       path = "Figures",
       width = 40, height = 30, units = "cm",
       dpi = 600)


# Mean signal position with error bars (sd)
ggplot(signal_positions_stat, aes(x, y, ymin = y - sd, ymax = y + sd)) +
  geom_pointrange() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Signal position as box plot
ggplot(signal_fit_par, aes(group, mean)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  xlab(expression(paste(""^13, C, "-glucose in medium / %", sep = ""))) +
  ylab(expression(paste("Wavenumber / ", cm^-1, sep = ""))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        text = element_text(family = "sans"),
        axis.ticks.length = unit(0.15, "cm")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5), linetype = "longdash", alpha = 0.4, size=1)

ggsave(filename = paste(data_designation, "_signal_positions.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)


# FWHM as box plot
ggplot(FWHM, aes(group, FWHM)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  xlab(expression(paste(""^13, C, "-glucose in medium / %", sep = ""))) +
  ylab(expression(paste("FWHM / ", cm^-1, sep = ""))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        text = element_text(family = "sans"),
        axis.ticks.length = unit(0.15, "cm")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5), linetype = "longdash", alpha = 0.4, size=1) +
  scale_y_continuous(limits = c(10,85))

ggsave(filename = paste(data_designation, "_FWHM.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)
