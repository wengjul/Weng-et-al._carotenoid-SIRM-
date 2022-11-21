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
upper_limit <- 1600
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






#### Principle Component Analysis (PCA) ####


#install.packages("FactoMineR") # for conducting the PCA
library(FactoMineR) # for conducting the PCA


# Bring variables to columns: transpose
spectra_base_norm_trans <- t(spectra_base_norm[-1])


# Conduct PCA with standardized data
spectra_pca <- PCA(spectra_base_norm_trans, scale.unit = TRUE, graph = FALSE)


# define data frame for coordinates of individuals
individuals_coord <- data.frame(spectra_pca$ind$coord)

# label groups and append it to the coordinates of individuals
# define labels for groups
group_labels <- c("0% 13C", "25% 13C", "50% 13C", "75% 13C", "100% 13C")

group_label_vector <- c()

for (i in seq(length(n_groups))) {
  group_label_vector <- append(group_label_vector, 
                               rep(group_labels[i], n_groups[i]))
}

group_label_factor <- factor(group_label_vector, levels = group_labels)

individuals_coord <- cbind(individuals_coord, "Samples" = group_label_factor)




# define data frame for eigenvalues -> Scree Plot
#(proportion of variances retained by the principal components )

scree_plot <- data.frame("Dimensions" = 1:length(spectra_pca$eig[, 2]),
  "Percentage" = spectra_pca$eig[, 2])
scree_plot <- scree_plot[1:10,]



#### Visualization ####

# Create Output folder for plots named "Figures"

dir.create(file.path(working_directory, "Figures"), showWarnings = FALSE)



# Visualize and safe Scree plot


ggplot(scree_plot, aes(Dimensions, Percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line() +
  geom_point() +
  xlab("Principal components") +
  ylab("Percentage of explained variance") +
  scale_x_continuous(breaks = c(1:10)) +
  geom_text(aes(label = paste(round(Percentage, 1), "%")), 
            check_overlap = TRUE, nudge_y = 0.8, nudge_x = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        text = element_text(family = "sans"))

ggsave(filename = paste(data_designation, "_scree_plot.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)



# Get colors for visualization via ColorBrewer package
#install.packages("RColorBrewer")
library(RColorBrewer)
colors_PCA <- brewer.pal(length(n_groups), "Dark2")
# or if red / blue color code activate
#colors_PCA <- rev(brewer.pal(length(n_groups), "RdBu"))
#colors_PCA[3] <- "#7F7F7F"
colors_PCA <- rev(c("#b2182b", "#ef8a62", "#7F7F7F", "#67a9cf", "#2166ac"))


# Visualize individuals, grouped by sample

# Dim1 vs. Dim2
ggplot(individuals_coord, aes(Dim.1, Dim.2)) +
  geom_point(aes(colour = Samples, shape = Samples)) +
  xlab(paste("PC1 (", round(scree_plot$Percentage[1], 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round(scree_plot$Percentage[2], 1), "%)", sep = "")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        text = element_text(family = "sans")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = colors_PCA)


ggsave(filename = paste(data_designation, "_individuals_Dim1Dim2.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)

# Dim1 vs. Dim3
ggplot(individuals_coord, aes(Dim.1, Dim.3)) +
  geom_point(aes(colour = Samples, shape = Samples)) +
  xlab(paste("PC1 (", round(scree_plot$Percentage[1], 1), "%)", sep = "")) +
  ylab(paste("PC3 (", round(scree_plot$Percentage[3], 1), "%)", sep = "")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        text = element_text(family = "sans")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = colors_PCA)


ggsave(filename = paste(data_designation, "_individuals_Dim1Dim3.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)

# Dim2 vs. Dim3
ggplot(individuals_coord, aes(Dim.2, Dim.3)) +
  geom_point(aes(colour = Samples, shape = Samples), size=2.5) +
  xlab(paste("PC 2 (", round(scree_plot$Percentage[2], 1), "%)", sep = "")) +
  ylab(paste("PC 3 (", round(scree_plot$Percentage[3], 1), "%)", sep = "")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        text = element_text(family = "sans"),
        legend.text.align = 0,
        axis.ticks.length = unit(0.15, "cm")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3, size=1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3, size=1) +
  scale_color_manual(values = colors_PCA) +
  scale_color_manual(values = colors_PCA, 
                     labels = c(expression(paste("0% ", ""^13, C, sep = "")),
                                expression(paste("25% ", ""^13, C, sep = "")),
                                expression(paste("50% ", ""^13, C, sep = "")),
                                expression(paste("75% ", ""^13, C, sep = "")),
                                expression(paste("100% ", ""^13, C, sep = "")))) +
  scale_shape_discrete(labels = c(expression(paste("0% ", ""^13, C, sep = "")),
                                  expression(paste("25% ", ""^13, C, sep = "")),
                                  expression(paste("50% ", ""^13, C, sep = "")),
                                  expression(paste("75% ", ""^13, C, sep = "")),
                                  expression(paste("100% ", ""^13, C, sep = ""))))


ggsave(filename = paste(data_designation, "_individuals_Dim2Dim3.png", sep = ""),
       path = "Figures",
       width = 20, height = 15, units = "cm",
       dpi = 600)


# 3 Dimensions : Dim1 vs. Dim2 vs. Dim3
#install.packages("rgl")
library("rgl")

get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

plot3d(x = individuals_coord$Dim.1 , y = individuals_coord$Dim.2, z = individuals_coord$Dim.3,
       xlab = "", 
       ylab = "", 
       zlab = "",
       col = get_colors(individuals_coord$Samples, colors_PCA),
       box = FALSE, size = 10)
       #type ="s", radius = 0.5) # activate for 3D spheres

mtext3d(paste("PC1 (", round(scree_plot$Percentage[1], 1), "%)", sep = ""), "x++", line = 2)
mtext3d(paste("PC2 (", round(scree_plot$Percentage[2], 1), "%)", sep = ""), "y--", line = 2)
mtext3d(paste("PC3 (", round(scree_plot$Percentage[3], 1), "%)", sep = ""), "z-+", line = 2)



rgl.snapshot(filename = paste(data_designation, "_individuals_3D_Dim1Dim2Dim3.png", sep = ""))
