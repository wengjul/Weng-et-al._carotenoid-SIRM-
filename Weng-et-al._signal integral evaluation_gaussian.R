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
fit2gaus <- function(x,y,sig1,sig2,scale1,scale2,y_offset){
  f = function(p){
    d = p[5]+abs(p[3])*dnorm(x,mean=1523,sd=p[1])+abs(p[4])*dnorm(x,mean=1483.3,sd=p[2])
    sum((d-y)^2)
  }
  optim(c(sig1,sig2,scale1,scale2,y_offset),f)
}


# Initialize vectors to store estimated parameters of fits

fit_sig1 <- c()
fit_sig2 <- c()
fit_scale1 <- c()
fit_scale2 <- c()
fit_y_offset <- c()

# define function used for fitting
fitfct <- fit2gaus
# execute fit function, save results for parameters in list, insert starting parameters

for(i in seq(2, length(spectra_fit_crop))) {
  fit_sig1 <- append(fit_sig1, 
                   fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 15, 15, 1, 1, 0.2)$par[1])
  fit_sig2 <- append(fit_sig2, 
                    fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 15, 15, 1, 1, 0.2)$par[2])
  fit_scale1 <- append(fit_scale1, 
                      fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 15, 15, 1, 1, 0.2)$par[3])
  fit_scale2 <- append(fit_scale2, 
                         fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 15, 15, 1, 1, 0.2)$par[4])
  fit_y_offset <- append(fit_y_offset, 
                         fitfct(spectra_fit_crop[,1], spectra_fit_crop[,i], 15, 15, 1, 1, 0.2)$par[5])
}



# generate vector that stores group labels of each observation

group_labels <- c(0, 2, 4, 6, 8, 24)

group_label_vector <- c()

for (i in seq(length(n_groups))) {
  group_label_vector <- append(group_label_vector, 
                               rep(group_labels[i], n_groups[i]))
}

# calculate curve from fit parameters

signal_fit_par <- as.data.frame(
  cbind("group" = group_label_vector, "sig1" = fit_sig1, "sig2" = fit_sig2, 
        "scale1" = fit_scale1, "scale2" = fit_scale2, "offset" = fit_y_offset))
signal_fit_par$group <- as.factor(signal_fit_par$group)


# Both Gaussians
signal_fit <- data.frame("x" = spectra_fit_crop$x)

for (i in seq(length(signal_fit_par$group))) {
  signal_fit <- cbind(signal_fit, 
                       signal_fit_par[i,6]+
                        abs(signal_fit_par[i,4])*
                        dnorm(x=spectra_fit_crop$x, mean=1523, sd=signal_fit_par[i,2]) +
                        abs(signal_fit_par[i,5])*
                        dnorm(x=spectra_fit_crop$x, mean=1483.3, sd=signal_fit_par[i,3]))
}
colnames(signal_fit) <- colnames(spectra_fit_crop)

# 1. Gaussian

signal_fit1 <- data.frame("x" = spectra_fit_crop$x)

for (i in seq(length(signal_fit_par$group))) {
  signal_fit1 <- cbind(signal_fit1, 
                      signal_fit_par[i,6]+
                        abs(signal_fit_par[i,4])*
                        dnorm(x=spectra_fit_crop$x, mean=1523, sd=signal_fit_par[i,2]))
}
colnames(signal_fit1) <- colnames(spectra_fit_crop)
# 2. Gaussian

signal_fit2 <- data.frame("x" = spectra_fit_crop$x)

for (i in seq(length(signal_fit_par$group))) {
  signal_fit2 <- cbind(signal_fit2, 
                       signal_fit_par[i,6]+
                         abs(signal_fit_par[i,5])*
                         dnorm(x=spectra_fit_crop$x, mean=1483.3, sd=signal_fit_par[i,3]))
}
colnames(signal_fit2) <- colnames(spectra_fit_crop)


# calculate residuals of fits

residuals <- spectra_fit_crop - signal_fit
residuals$x <- spectra_fit_crop$x


# Integration of Gausians
gaussian1 <- function(x,sig1,scale1){
  abs(scale1)*dnorm(x,mean=1523,sd=sig1)
}
gaussian2 <- function(x,sig2,scale2){
  abs(scale2)*dnorm(x,mean=1483.3,sd=sig2)
}

int_gaussian1 <- c()
int_gaussian2 <- c()

for (i in seq(length(signal_fit_par$group))) {
  int_gaussian1 <- append(int_gaussian1, 
                     integrate(f=gaussian1, 
                               lower=1513, 
                               upper=1533, 
                               sig1=signal_fit_par[i,2], 
                               scale1=signal_fit_par[i,4])[1]$value)
  
  int_gaussian2 <- append(int_gaussian2, 
                     integrate(f=gaussian2, 
                               lower=1473, 
                               upper=1493, 
                               sig2=signal_fit_par[i,3], 
                               scale2=signal_fit_par[i,5])[1]$value)
}

integrals <- as.data.frame(
  cbind("group" = group_label_vector, "gauss1" = int_gaussian1, "gauss2" = int_gaussian2, 
        "ratio" = int_gaussian1/(int_gaussian1+int_gaussian2)))
integrals$group <- as.factor(integrals$group)

#### Visualization ####

# Create Output folder for plots named "Figures"

dir.create(file.path(working_directory, "Figures"), showWarnings = FALSE)

# single spectra facet plot: Both Gaussians

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
        axis.title = element_text(size = 16),
        text = element_text(family = "sans")) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks=c(1450, 1550)) +
  geom_vline(xintercept = 1524, color = "red", linetype = "longdash", alpha = 0.4) +
  geom_vline(xintercept = 1485, color = "red", linetype = "longdash", alpha = 0.4) +
  geom_line(data = signal_fit_long, aes(x, y), color = "darkred")

ggsave(filename = paste(data_designation, "_fit_facet_plot_bothgaus.png", sep = ""),
       path = "Figures",
       width = 40, height = 30, units = "cm",
       dpi = 600)

# single spectra facet plot: Single Gaussians

# Converting data from wide to long data frame -> required for ggplot2
# melt function from "data.table" is applied
library(data.table)

# convert for for single spectra facet plot
spectra_fit_crop_long <- melt(data = spectra_fit_crop, 
                              id.vars = "x",
                              variable.name = "Spectrum",
                              value.name = "y")
signal_fit1_long <- melt(data = signal_fit1, 
                        id.vars = "x",
                        variable.name = "Spectrum",
                        value.name = "y")
signal_fit2_long <- melt(data = signal_fit2, 
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
        axis.title = element_text(size = 16),
        text = element_text(family = "sans")) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks=c(1450, 1550)) +
  geom_vline(xintercept = 1524, color = "red", linetype = "longdash", alpha = 0.4) +
  geom_vline(xintercept = 1485, color = "red", linetype = "longdash", alpha = 0.4) +
  geom_line(data = signal_fit1_long, aes(x, y), color = "orange") +
  geom_line(data = signal_fit2_long, aes(x, y), color = "darkblue")

ggsave(filename = paste(data_designation, "_fit_facet_plot.png", sep = ""),
       path = "Figures",
       width = 40, height = 30, units = "cm",
       dpi = 600)


# single spectra facet plot: residuals: Both Gaussian

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
        axis.title = element_text(size = 16),
        text = element_text(family = "sans")) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks=c(1450, 1550)) +
  geom_vline(xintercept = 1524, color = "red", linetype = "longdash", alpha = 0.4) +
  geom_vline(xintercept = 1485, color = "red", linetype = "longdash", alpha = 0.4)

ggsave(filename = paste(data_designation, "_fit_facet_plot_residuals.png", sep = ""),
       path = "Figures",
       width = 40, height = 30, units = "cm",
       dpi = 600)


# Signal integral ratio as box plot
ggplot(integrals, aes(group, ratio)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  xlab("Incubation time / h") +
  ylab(expression(paste("1522 ", cm^-1, sep = "/(1522 ", cm^-1, sep = " + 1488 ", cm^-1, sep = ")"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        text = element_text(family = "sans"),
        axis.ticks.length = unit(0.15, "cm")) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5), linetype = "longdash", alpha = 0.4, size=1)


ggsave(filename = paste(data_designation, "_signal_integral.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)
