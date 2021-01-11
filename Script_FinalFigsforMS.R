library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(viridis)
library(scales)
library(pals)
library(ggdendro)


##viridis colors
show_col(viridis_pal()(20))
viridis_pal()(8)
viridis_pal()(4)
show_col(viridis_pal()(3))
show_col(ocean.curl(20))

#pal colors
pal.bands(ocean.curl)
pal.bands(viridis)

#order of levels for Species will be Cp, Hz, Px, Sf
#re-ordering colors so the yellow is Cp
pal <- viridis_pal()(4)[c(4,1,2,3)]
show_col(pal)
#So Cp=yellow=pal[1], Hz=purple=pal[2], Px=blue=pal[3], Sf=green=pal[4]
pal3 <- viridis_pal()(3)


#set symbol types
symb <- c(16,16,16,16)  #all 4 are transparent dots
#symb <- values = c(8, 9, 1, 2)  #different symbol shapes

#set transparency for points
alpha_p <- 0.5

#set transparency for boxplots
alpha_b <- 0.8

#and color for any plots with just one 
dot_col <- viridis_pal()(1)


##Use function for composite figs to create legend as a separate object
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



#----------------------------------------------------------------------------------
###Figure 1A-D: A big composite fig of richness effects on performance and fungal growth
#-----------------------------

##Load workspace with standardized data
load("./Outputs/Workspaces/StandardizedData")

#For Pupal Weights
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))

p1 <- ggplot(d.temp, aes(x=Richness, y=Pupal.weight.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(Richness, 0.5), alpha = alpha_b), outlier.shape = NA)
  
p1.b <- p1 + ylab("Pupal Mass") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5)) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none") +
  labs(subtitle="(A)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1.b


#For Dev Speed
d.temp <- filter (d.rich, Treatment != "C", !is.na(Days.to.pupation.ST.inv))

p2 <- ggplot(d.temp, aes(x=Richness, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, height=0.1, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(Richness, 0.5), alpha = alpha_b), outlier.shape = NA)

p2.b <- p2 + ylab("Development Speed") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5, 10.5)) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(B)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.b

#For Survival
d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
p3 <- ggplot(d.temp, aes(x=Richness, y=PropSurv.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(Richness, 0.5), alpha = alpha_b), outlier.shape = NA)

p3.b <- p3 + ylab("Survival") + xlab("Richness") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5))+
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(C)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3.b


d.temp <- filter (d.fungi2, Treatment != "DMSO")
p4 <- ggplot(d.temp, aes(x=Richness, y=dAbs.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Fungi, col=Fungi), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(Richness, 0.5), alpha = alpha_b), outlier.shape = NA)

p4.b <- p4 + ylab("Growth Rate") + xlab("Richness") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5)) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(D)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p4.b





#To get a legend for insects
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))

p.temp <- ggplot(d.temp, aes(x=Richness, y=Pupal.weight.ST)) +
  geom_jitter(width=0.3, height=0.3, aes(shape=Species, col=Species), size=6)

p.temp.b <- p.temp + ylab("Standardized Pupal Weight") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p.temp.b)
p.temp.b



#To get a legend for fungi
d.temp <- filter (d.fungi2, Treatment != "DMSO")
labs <- c(expression(paste(italic(" B. dothidea   "))),
          expression(paste(italic(" Colletotrichum   "))),
          expression(paste(italic(" P. expansum  "))), 
          expression(paste(italic(" S. sclerotiorum  "))))
p.temp <- ggplot(d.temp, aes(x=Richness, y=dAbs.ST)) +
  geom_jitter(width=0.3, aes(shape=Fungi, col=Fungi), size=6)

p.temp.b <- p.temp + ylab("Growth Rate") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l2 <- get_legend(p.temp.b)
p.temp.b


#Make the composite fig

t1 <- textGrob("Insect Species (Plots A, B, C)", gp=gpar(fontsize=35))
t2 <- textGrob("Fungal Species (Plot D)", gp=gpar(fontsize=35))

tiff(height=1250, width=1000, filename="./Outputs/Figures/Fig1_Richness.tiff", type="cairo")
lay <- rbind(c(1,2),  #graphs A+B
             c(3,4), #graphs C+D
             c(5,5), #legend label 1
             c(6,6),  #legend 1
             c(7,7),  #legend label 2
             c(8,8)) #legend 2
grid.arrange(p1.b, p2.b, p3.b, p4.b, t1, l1, t2, l2, nrow=6, ncol=2, widths=c(500, 500),
             heights=c(500, 500, 50, 75, 50, 75), layout_matrix=lay)
dev.off()


cairo_pdf(height=18.75, width=15, family="sans",filename="./Outputs/Figures/Fig1_Richness.pdf")
grid.arrange(p1.b, p2.b, p3.b, p4.b, t1, l1, t2, l2, nrow=6, ncol=2, widths=c(500, 500),
             heights=c(500, 500, 50, 75, 50, 75), layout_matrix=lay)
dev.off()


#----------------------------------------------------------------------------------
###Figure 2A-D: A big composite fig of SD effects on performance and fungal growth
#----------------------------

##Load workspace with standardized data
load("./Outputs/Workspaces/StandardizedData")

#For Pupal Weights
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

p1 <- ggplot(d.temp, aes(x=SD, y=Pupal.weight.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(SD, 0.5), alpha = alpha_b), outlier.shape = NA)

p1.b <- p1 + ylab("Pupal Mass") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_discrete(breaks=c("L", "M", "H"), labels=c("low", "med", "high")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none") +
  labs(subtitle="(A)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1.b


#For Dev Speed
d.temp <- filter (d.rich, Treatment != "C", !is.na(Days.to.pupation.ST.inv))
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

p2 <- ggplot(d.temp, aes(x=SD, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, height=0.1, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(SD, 0.5), alpha = alpha_b), outlier.shape = NA)

p2.b <- p2 + ylab("Development Speed") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_discrete(breaks=c("L", "M", "H"), labels=c("low", "med", "high")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(B)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.b

#For Survival
d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

p3 <- ggplot(d.temp, aes(x=SD, y=PropSurv.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(SD, 0.5), alpha = alpha_b), outlier.shape = NA)

p3.b <- p3 + ylab("Survival") + xlab("Structural Diversity") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_discrete(breaks=c("L", "M", "H"), labels=c("low", "med", "high")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(C)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3.b


d.temp <- filter (d.fungi2, Treatment != "DMSO")
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

p4 <- ggplot(d.temp, aes(x=SD, y=dAbs.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=Fungi, col=Fungi), alpha=alpha_p) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(SD, 0.5), alpha = alpha_b), outlier.shape = NA)

p4.b <- p4 + ylab("Growth Rate") + xlab("Structural Diversity") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_discrete(breaks=c("L", "M", "H"), labels=c("low", "med", "high")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(D)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p4.b





#To get a legend for insects
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))

p.temp <- ggplot(d.temp, aes(x=SD, y=Pupal.weight.ST)) +
  geom_jitter(width=0.3, height=0.3, aes(shape=Species, col=Species), size=6)

p.temp.b <- p.temp + ylab("Standardized Pupal Weight") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p.temp.b)
p.temp.b



#To get a legend for fungi
d.temp <- filter (d.fungi2, Treatment != "DMSO")
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

labs <- c(expression(paste(italic(" B. dothidea   "))),
          expression(paste(italic(" Colletotrichum   "))),
          expression(paste(italic(" P. expansum  "))), 
          expression(paste(italic(" S. sclerotiorum  "))))
p.temp <- ggplot(d.temp, aes(x=SD, y=dAbs.ST)) +
  geom_jitter(width=0.3, aes(shape=Fungi, col=Fungi), size=6)

p.temp.b <- p.temp + ylab("Growth Rate") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l2 <- get_legend(p.temp.b)
p.temp.b


#Make the composite fig

t1 <- textGrob("Insect Species (Plots A, B, C)", gp=gpar(fontsize=35))
t2 <- textGrob("Fungal Species (Plot D)", gp=gpar(fontsize=35))

tiff(height=1250, width=1000, filename="./Outputs/Figures/Fig2_SD.tiff", type="cairo")
lay <- rbind(c(1,2),  #graphs A+B
             c(3,4), #graphs C+D
             c(5,5), #legend label 1
             c(6,6),  #legend 1
             c(7,7),  #legend label 2
             c(8,8)) #legend 2
grid.arrange(p1.b, p2.b, p3.b, p4.b, t1, l1, t2, l2, nrow=6, ncol=2, widths=c(500, 500),
             heights=c(500, 500, 50, 75, 50, 75), layout_matrix=lay)
dev.off()


cairo_pdf(height=18.75, width=15, family="sans",filename="./Outputs/Figures/Fig2_SD.pdf")
grid.arrange(p1.b, p2.b, p3.b, p4.b, t1, l1, t2, l2, nrow=6, ncol=2, widths=c(500, 500),
             heights=c(500, 500, 50, 75, 50, 75), layout_matrix=lay)
dev.off()

#-----------------------------------------------------
###Figure 3: effects of richness and SD on number of organisms affected 
#------------------------------

load("./Outputs/Workspaces/FinalAnalyses_Pred2a-2b_NumEffects")

#simple plot
plot(NumOrgAffected ~ jitter(Richness), data=NumEffects.rich)
abline(lm(NumOrgAffected ~ Richness, data=NumEffects.rich))

#pretty plot
d.temp <- NumEffects.rich
p1 <- ggplot(d.temp, aes(x=Richness, y=NumOrgAffected)) +
  geom_jitter(width=0.15, height=0.15, cex=4, col=pal[2]) +
  geom_smooth(method="lm", col="#6B244C", aes(group=1)) +
  ylab("Numb. Consumers Affected") + xlab("Richness") +
  labs(subtitle="(A)") +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1", "2","4","6", "8", "10")) +
  theme(text = element_text(size = 35)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1

d.temp <- NumEffects.rich
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))
p2 <- ggplot(d.temp, aes(x=SD, y=NumOrgAffected)) +
  geom_jitter(width=0.15, height=0.15, cex=4, col=pal[2]) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(SD, 0.5), alpha = alpha_b), outlier.shape = NA) +
  ylab("") + xlab("Structural Diversity") +
  labs(subtitle="(B)") +
  scale_x_discrete(breaks=c("L", "M", "H"), labels=c("low", "med", "high")) +
  scale_y_continuous(limits=c(2,7.8), breaks=c(2,3,4,5,6,7)) +
  theme(text = element_text(size = 35)) +
  annotate(geom = "text", x = 1, y = 7.6, label = "A", size=8) +
  annotate(geom = "text", x = 2, y = 7.6, label = "B", size=8) +
  annotate(geom = "text", x = 3, y = 7.6, label = "B", size=8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
p2

tiff(height=500, width=1000, filename="./Outputs/Figures/Fig3_NumbEffects.tiff", type="cairo")
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

cairo_pdf(height=7.5, width=15, family="sans",filename="./Outputs/Figures/Fig3_NumbEffects.pdf")
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

#---------------------------------------------------------------------------
## Figure 4: Effects of individual compounds, heatmaps
#------------------------------

load("./Outputs/Workspaces/FinalAnalyses_Pred2d&3a_IndivComps")

#And we also need the distance matrix for compound structural similarities
d.avg <- read.csv("./Outputs/Tables/CompoundDistances.csv")
rownames(d.avg) <- d.avg$X
d.avg <- d.avg[-1]

#Set diverging colors for heatmaps
show_col(viridis_pal()(10))
hm_high <- viridis_pal()(10)[6]
hm_low <- viridis_pal()(10)[1]
hm_mid = "gray95"

#hm_high <- ocean.curl(20)[2]
#hm_low <- ocean.curl(20)[19]
#hm_mid = ocean.curl(3)[2]


base_size <- 20
theme_update(plot.subtitle = element_text(hjust = 0.5, vjust=-5)) #originally at hjust=0 and vjust=1

#Make a function to create a character vector that indicates
#significance level with stars based on p-value
stars <- function(x){
  if(x < 0.001 & !is.na(x)){
    s <- "***"
  }else{
    if(x < 0.01 & x >= 0.001 & !is.na(x)){
      s <- "**"
    }else{
      if (x < 0.05 & x >=0.01 & !is.na(x)){
        s <- "*"
      } else{
        s <- NA
      }
    }
  }
  return (s)
}

stars <- Vectorize(stars)


####The plots


#The dendrogram

hc <- hclust(as.dist(1-d.avg), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="Avg Distance") 

dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the compound labels
comp_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, comp = as.character(label), height = 1))

# Limits for the vertical axes
comp_axis_limits <- with(
  comp_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1




# Dendrogram plot
p.dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(position="top")+   #expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = comp_pos_table$y_center, limits = comp_axis_limits, 
                     expand = c(0, 0), labels=c("Phloridzin", "Phloretin", "Rutin", 
                                                "Hyperin", "Quercetin", "Catechin", "Epicatechin",  
                                                "Gentistic Acid","Gallic Acid", "Syringic Acid",  
                                                "Chlorogenic Acid", "p-Coumaric Acid", "Caffeic Acid", 
                                                "Ferulic Acid")) + 
  labs(x = "", y = "", colour = "", size = "") +
  labs(subtitle="  \n    ", position="center") + 
  #theme_void() + 
  theme(axis.text.x = element_text(color="white"),
        panel.grid.minor = element_blank(), axis.ticks.x=element_blank(),
        panel.background= element_blank(),
        panel.grid.major.y = element_line(color="gray90"),
        panel.grid.major.x=element_blank())+
  #, panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(), axis.ticks = element_blank(),
  #panel.background = element_blank()) +
  theme(text = element_text(size = 20), plot.margin=unit(c(-0.5,0,0,0), "cm"),
        axis.text = element_text(size = 20)) 
p.dendr



#The heatmaps
d.temp <- d.sum

#re-order levels in the way we want them for the figure
d.temp$Treatment <- factor(d.temp$Treatment,levels=c("Phz", "Pht", "R", "H", "Q", "Ct",
                                                     "eCt", "GeA", "GA", "SA", "ChA",
                                                     "pCA", "CA", "FA"))

d.temp <- mutate(d.temp, stars.PW=stars(PW.p))
d.temp <- mutate(d.temp, stars.DtP=stars(DtP.p))
d.temp <- mutate(d.temp, stars.S=stars(S.pval))

p1 <- ggplot(d.temp, aes(x=Species, y=Treatment)) + 
  geom_tile(aes(fill = PropSurv.ST), colour = "white", size=0.25) + 
  #scale_fill_gradientn(colors=rev(ocean.curl(10)), guide = "colourbar", name="", limits=c(0.2, 1.8)) +
  scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high, 
                       midpoint = 1, guide = "colourbar", name="", limits=c(0.2, 1.8)) +  
  geom_text(aes(label=stars.S), col="white", size=6, vjust=0.7)
p1.b <- p1 + theme(text = element_text(size = 20), plot.margin=unit(c(-0.5,0,0,-0.85), "cm"),
                   axis.text = element_text(size = 20)) +
  labs(subtitle="Survival\n   ", position="center") + 
  labs(x="", y="") +
  #scale_y_discrete(labels=c("Phloridzin", "Phloretin", "Rutin", "Hyperin",
  #                         "Quercetin", "Catechin", "Epicatechin",  "Gentistic Acid",
  #                        "Gallic Acid", "Syringic Acid",  "Chlorogenic Acid", 
  #                       "p-Coumaric Acid", "Caffeic Acid", "Ferulic Acid")) +
  scale_x_discrete(position="top") +
  theme(legend.position="", axis.text.y=element_blank(),axis.ticks.y = element_blank())
p1.b


p2 <- ggplot(d.temp, aes(x=Species, y=Treatment))  + 
  geom_tile(aes(fill = DtP.avg), colour = "white") + 
  #scale_fill_gradientn(colors=rev(ocean.curl(10)), guide = "colourbar", name="", limits=c(0.2, 1.8)) +
  scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high, 
                       midpoint = 1, guide = "colourbar", name="", limits=c(0.2, 1.8)) +
  geom_text(aes(label=stars.DtP), col="white", size=6, vjust=0.7)
p2.b <- p2 +  theme(text = element_text(size = 20), plot.margin=unit(c(-0.5,0,0,-0.85), "cm"),
                    axis.text = element_text(size = 20),axis.text.y=element_blank(),
                    axis.ticks.y = element_blank()) +
  labs(subtitle="Development\nSpeed", x="", y="") +
  scale_x_discrete(position="top") +
  theme(legend.position="")
p2.b



p3 <- ggplot(d.temp, aes(x=Species, y=Treatment)) + 
  geom_tile(aes(fill = PW.avg), colour = "white") + 
  #scale_fill_gradientn(colors=rev(ocean.curl(10)), guide = "colourbar", name="", limits=c(0.2, 1.8)) +
  scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high, 
                       midpoint = 1, guide="colourbar", name="", limits=c(0.2, 1.8))  +
  geom_text(aes(label=stars.PW), col="white", size=6, vjust=0.7)
p3.b <- p3 + theme(text = element_text(size = 20), plot.margin=unit(c(-0.5,0,0,-0.85), "cm"),
                   axis.text = element_text(size = 20),axis.text.y=element_blank(),
                   axis.ticks.y = element_blank()) +
  labs(subtitle="Pupal\nMass", x="", y="") +
  scale_x_discrete(position="top") +
  theme(legend.position="")
p3.b


d.temp <- d.sum.F
d.temp$Treatment <- factor(d.temp$Treatment,levels=c("Phz", "Pht", "R", "H", "Q", "Ct",
                                                     "eCt", "GeA", "GA", "SA", "ChA",
                                                     "pCA", "CA", "FA"))

d.temp <- mutate(d.temp, stars.F=stars(F.p))

p4 <- ggplot(d.temp, aes(x=Fungi, y=Treatment)) + 
  geom_tile(aes(fill = F.avg), colour = "white") + 
  #scale_fill_gradientn(colors=rev(ocean.curl(10)), guide = "colourbar", name="", limits=c(0.2, 1.8)) +
  scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high, 
                       midpoint = 1, guide="colourbar", name="", limits=c(0.2, 1.8))  +
  geom_text(aes(label=stars.F), col="white", size=6, vjust=0.7)
p4.b <- p4 + theme(text = element_text(size = 20), plot.margin=unit(c(-0.5,0,0,-0.85), "cm"),
                   axis.text = element_text(size = 20),axis.text.y=element_blank(),
                   axis.ticks.y = element_blank()) +
  labs(subtitle="Growth\nRate", x="", y="") +
  scale_x_discrete(position="top", labels=c("Bd", "C", "Pe", "Ss")) +
  theme(legend.position="")
p4.b


#to get a legend
p.temp <- ggplot(d.temp, aes(x=Fungi, y=Treatment)) + 
  geom_tile(aes(fill = F.avg), colour = "white") + 
  #scale_fill_gradientn(colors=rev(ocean.curl(10)), guide = "colourbar", name="", limits=c(0.2, 1.8)) 
  scale_fill_gradient2(low = hm_low, mid = hm_mid, high = hm_high, 
                       midpoint = 1, guide="colourbar", name="", limits=c(0.2, 1.8))  
p.temp.b <- p.temp + theme(text = element_text(size = 20), plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"),
                           axis.text = element_text(size = 20),axis.text.y=element_blank(),
                           axis.ticks.y = element_blank()) +
  labs(subtitle="Growth\nRate", x="", y="") +
  scale_x_discrete(position="top", labels=c("Bd", "C", "Pe", "Ss")) +
  theme(legend.position=c(0.5,0.4), legend.direction="horizontal")+
  guides(fill=guide_colorbar(barwidth=35), ticks=TRUE)#hjust=200)
p.temp.b
l1 <- get_legend(p.temp.b)



s1 <- segmentsGrob(x0 = 0.5, y0 = 0,
                   x1 = 0.5, y1 = 1, gp = gpar(col="darkslategray"))
s2 <- segmentsGrob(x0 = 0, y0 = 0.5,
                   x1 = 1, y1 = 0.5, gp = gpar(col="darkslategray"))
t1 <- textGrob("", gp=gpar(fontsize=25))
t2 <- textGrob("Insects", gp=gpar(fontsize=25))
t3 <- textGrob("Fungi", gp=gpar(fontsize=25))


lay <- rbind(c(1,2,2,2,2,2,3,4),
             c(1,5,5,5,5,5,5,5),  
             c(6,7,8,9,10,11,12,13),   #plot, plot, line, plot, line, plot, line, plot  
             c(1,rep(14, 7)))

tiff(height=528, width=795.2, filename="./Outputs/Figures/Fig4_IndivCompHeatmaps.tiff", type="cairo") #0.8*widths
grid.arrange(t1,t2,s1,t3, s2,p.dendr, p1.b, s1,p2.b, s1, p3.b, s1, p4.b,l1, nrow=4, ncol=8, 
             widths=c(350,188,10,188,10,188,10,188),
             heights=c(50,10,600,100), layout_matrix=lay)
dev.off()



cairo_pdf(height=7.92, width=11.928, family="sans",filename="./Outputs/Figures/Fig4_IndivCompHeatmaps.pdf")
grid.arrange(t1,t2,s1,t3, s2,p.dendr, p1.b, s1,p2.b, s1, p3.b, s1, p4.b,l1, nrow=4, ncol=8, 
             widths=c(350,188,10,188,10,188,10,188),
             heights=c(50,10,600,100), layout_matrix=lay)
dev.off()


#------------------------------------------
# Fig. S1:Composite fig of all places that there are significant effects of richness or SD
#-------------------------------

##Load workspace with standardized data
load("./Outputs/Workspaces/StandardizedData")

######(A) Richness effects on Spodoptera survival

d.temp <- d.rich.PS[which(d.rich.PS$Richness != "0" & d.rich.PS$Species=="Sf"),]
p1 <- ggplot(d.temp, aes(x=Richness, y=PropSurv.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.2, color=pal[4], shape=symb[4], cex=2) +
  geom_boxplot(width=0.75, notch = FALSE, aes(group = cut_width(Richness, 0.5), alpha = 0.8), outlier.shape = NA)
p1.b <- p1 + ylab("Survival") + xlab("Richness") +
  geom_smooth(method="lm", col=pal[4], aes(group=1), fill=NA) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), labels=c("0","0.2","0.4","0.6", "0.8", "1", "1.2"), limits=c(0,1.2)) +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  labs(subtitle="(A)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p1.b




#---  (B) Sig. richness x SD interaction for Plutella
d.temp <- d.rich[which(d.rich$Richness != "0" & d.rich$Species=="Px" & d.rich$Sex=="m"),]
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

labs <- c("low", "med", "high")
p2 <- ggplot(d.temp, aes(x=Richness, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, height=0.1, aes(shape=SD, col=SD), alpha=alpha_p) +
  stat_summary(fun.y = "mean", data=d.temp, aes(color=SD), size = 4, geom = "point")

p2.b <- p2 + ylab("M Dev Speed") + xlab("Richness") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5)) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal3, labels=labs) +
  theme(legend.position="none")+
  labs(subtitle="(B)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p2.b





#---  (C) Sig. richness x SD interaction for Botrys
d.temp <- filter (d.fungi2, Fungi=="Botrys", Treatment != "DMSO")
d.temp$SD <- factor(d.temp$SD, levels=c("L", "M", "H"))

m.l <- lm(dAbs.ST ~ Richness, data=d.temp[which(d.temp$SD=="L"),])
d.new.l <- expand.grid(Richness = seq(0.5, 10.5, by = .1))
d.new.l$predict <- predict(m.l, newdata = d.new.l)

m.m <- lm(dAbs.ST ~ Richness, data=d.temp[which(d.temp$SD=="M"),])
d.new.m <- expand.grid(Richness = seq(0.5, 10.5, by = .1))
d.new.m$predict <- predict(m.m, newdata = d.new.m)

m.h <- lm(dAbs.ST ~ Richness, data=d.temp[which(d.temp$SD=="H"),])
d.new.h <- expand.grid(Richness = seq(0.5, 10.5, by = .1))
d.new.h$predict <- predict(m.h, newdata = d.new.h)

labs <- c("low", "med", "high")
p3 <- ggplot(d.temp, aes(x=Richness, y=dAbs.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.3, aes(shape=SD, col=SD), alpha=alpha_p) +
  stat_summary(fun.y = "mean", data=d.temp, aes(color=SD), size = 4, geom = "point") +
  geom_smooth(data=d.new.l, aes(y = predict), size = 1, col=pal3[1], se=FALSE) +
  geom_smooth(data=d.new.m, aes(y = predict), size = 1, col=pal3[2], se=FALSE)
  #geom_smooth(data=d.new.h, aes(y = predict), size = 1, col=pal3[3], se=FALSE)

p3.b <- p3 + ylab("Growth Rate") + xlab("Richness") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), labels=c("1","2","4", "6", "8", "10"), limits=c(0.5,10.5)) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal3, labels=labs) +
  theme(legend.position="right", legend.title=element_blank())+
  labs(subtitle="(C)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p3.b)
p3.b

#Make the composite fig

tiff(height=500, width=1600, filename="./Outputs/Figures/FigS1_IndivEffects.tiff", type="cairo")
lay <- rbind(c(1,2,3))  #graphs
grid.arrange(p1.b, p2.b, p3.b, nrow=1, ncol=3, widths=c(500, 500, 600),
             heights=c(500), layout_matrix=lay)
dev.off()


cairo_pdf(height=6.25, width=20, family="sans",filename="./Outputs/Figures/FigS1_IndivEffects.pdf")
lay <- rbind(c(1,2,3))  #graphs
grid.arrange(p1.b, p2.b, p3.b, nrow=1, ncol=3, widths=c(500, 500, 600),
             heights=c(500), layout_matrix=lay)
dev.off()

#----------------------------------------------------
### Fig. S2: Effects of richness on Zdiff (exp-obs effects of mixtures)
# For prediction 1e
#------------------------------

load("./Outputs/Workspaces/FinalAnalyses_Pred1d_synergy")

#Only case where we found any potential effects was a marginal negative effect of richness
#on Zdiff in Px pupal weights. This would indicate antagonism, or at least decreasing
#synergy with increasing richness

summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]))  #marginal antagonism

#simple plot
plot(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),])
abline(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]))


#pretty plot of just Px (only one that is even marginally signficant)

d.temp <- Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]
p <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=6, col=pal[3], shape=symb[3]) +
  #geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = 0.8), outlier.shape = NA)
  geom_smooth(method="lm", col="#6B244C", aes(group=1)) +
  #scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  ylab("Interaction Strength") +
  theme(text = element_text(size = 35))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p




#Big composite plot of significant and non-significant relationships

d.temp <- Zs.PW[which(Zs.PW$Species=="Cp" & Zs.PW$Richness !=1),]
p1 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("Interaction Strength") +
  labs(subtitle="(A)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[1]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p1

d.temp <- Zs.PW[which(Zs.PW$Species=="Hz" & Zs.PW$Richness !=1),]
p2 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(B)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[2]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p2

d.temp <- Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]
p3 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  geom_smooth(method="lm", col=pal[3], aes(group=1))+
  ylab("") +
  labs(subtitle="(C)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[3]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p3

d.temp <- Zs.PW[which(Zs.PW$Species=="Sf" & Zs.PW$Richness !=1),]
p4 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(D)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[4]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p4

d.temp <- Zs.DtP[which(Zs.PW$Species=="Cp" & Zs.PW$Richness !=1),]
p5 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("Interaction Strength") +
  labs(subtitle="(E)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[1]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p5

d.temp <- Zs.DtP[which(Zs.PW$Species=="Hz" & Zs.PW$Richness !=1),]
p6 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(F)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[2]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p6


d.temp <- Zs.DtP[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]
p7 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(G)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[3]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p7

d.temp <- Zs.DtP[which(Zs.PW$Species=="Sf" & Zs.PW$Richness !=1),]
p8 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(H)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[4]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p8

d.temp <- Zs.F[which(Zs.F$Species=="Botrys" & Zs.F$Richness !=1),]
p9 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("Interaction Strength") +
  labs(subtitle="(I)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[1]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p9

d.temp <- Zs.F[which(Zs.F$Species=="Collet" & Zs.F$Richness !=1),]
p10 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(J)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[2]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p10

d.temp <- Zs.F[which(Zs.F$Species=="Penicillium" & Zs.F$Richness !=1),]
p11 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(K)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[3]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p11

d.temp <- Zs.F[which(Zs.F$Species=="Sclerotinia" & Zs.F$Richness !=1),]
p12 <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, cex=4, aes(shape=CI.overlap, col=Species)) +
  geom_hline(yintercept=0, color="gray20", lwd=1.2, lty=3) +
  ylab("") +
  labs(subtitle="(L)") +
  scale_shape_manual(values = c(16,8)) + 
  scale_color_manual(values = pal[4]) +
  theme(text = element_text(size = 35), axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
p12





#To get a legend for insects
d.temp <- Zs.PW[which(Zs.PW$Richness !=1),]
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))
p.temp <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, aes(shape=Species, col=Species), size=8) +
  ylab("Interaction Strength") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p.temp)


#legend for fungi
d.temp <- Zs.F[which(Zs.F$Richness !=1),]
labs <- c(expression(paste(italic(" B. dothidea   "))),
          expression(paste(italic(" Colletotrichum   "))),
          expression(paste(italic(" P. expansum  "))), 
          expression(paste(italic(" S. sclerotiorum  "))))
p.temp <- ggplot(d.temp, aes(x=Richness, y=diff)) +
  geom_jitter(width=0.03, aes(shape=Species, col=Species), size=8) +
  ylab("Interaction Strength") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l2 <- get_legend(p.temp)


s1 <- segmentsGrob(x0 = 0, y0 = 0.5,
                   x1 = 1, y1 = 0.5, gp = gpar(col="darkslategray"))
t1 <- textGrob("Pupal Mass", gp=gpar(fontsize=35, fontface="bold"))
t2 <- textGrob("Development Speed", gp=gpar(fontsize=35, fontface="bold"))
t3 <- textGrob("Growth Rate", gp=gpar(fontsize=35, fontface="bold"))
t4 <- textGrob("Richness", gp=gpar(fontsize=35))
t5 <- textGrob("Insect Species (Plots A-H)", gp=gpar(fontsize=35))
t6 <- textGrob("Fungal Species (Plot I-L)", gp=gpar(fontsize=35))




lay <- rbind(c(1,1,1,1), 
             c(2,3,4,5),  
             c(6,6,6,6),  
             c(7,7,7,7),
             c(8,8,8,8),
             c(9,10,11,12),
             c(13,13,13,13),
             c(14,14,14,14),
             c(15,15,15,15),
             c(16,17,18,19),
             c(20,20,20,20),
             c(21,21,21,21),
             c(22,22,22,22),
             c(23,23,23,23),
             c(24,24,24,24)) 
tiff(height=1980, width=1800, filename="./Outputs/Figures/FigS2_RichVSIntStrength.tiff", type="cairo")
grid.arrange(t1,p1,p2,p3,p4,t4,s1,t2,p5,p6,p7,p8,t4,s1,t3,p9,p10,p11,p12,t4,t5,l1,t6,l2, nrow=15, ncol=4, 
             widths=c(300,300,300,300),
             heights=c(50, 300, 50,10,50,300,50,10,50,300,50,50,50,50,50), layout_matrix=lay)
dev.off()


cairo_pdf(height=23.76, width=21.6, family="sans",filename="./Outputs/Figures/FigS2_RichVSIntStrength.pdf")
grid.arrange(t1,p1,p2,p3,p4,t4,s1,t2,p5,p6,p7,p8,t4,s1,t3,p9,p10,p11,p12,t4,t5,l1,t6,l2, nrow=15, ncol=4, 
             widths=c(300,300,300,300),
             heights=c(50, 300, 50,10,50,300,50,10,50,300,50,50,50,50,50), layout_matrix=lay)
dev.off()


#---------------------------------------------------
### Fig. S3: Comparing mixtures to most effective singleton (Prediction 1e)
#-------------------------------

load("./Outputs/Workspaces/FinalAnalyses_Pred1e_singletonsVSmix")


#Big composite plot of significant and non-significant relationships

#  For pupal weights
d.temp <- means.PW[which(means.PW$Species=="Cp" & means.PW$Richness >1),]
p1 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[1], shape=symb[1], cex=4) +
  stat_smooth(method="glm", method.args = list(family = "binomial"), col="#6B244C", aes(group=1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(A)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p1

d.temp <- means.PW[which(means.PW$Species=="Hz" & means.PW$Richness >1),]
p2 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[2], shape=symb[2], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(B)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p2

d.temp <- means.PW[which(means.PW$Species=="Px" & means.PW$Richness >1),]
p3 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[3], shape=symb[3], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(C)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p3

d.temp <- means.PW[which(means.PW$Species=="Sf" & means.PW$Richness >1),]
p4 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[4], shape=symb[4], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(D)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p4

#For days to pupation
d.temp <- means.DtP[which(means.DtP$Species=="Cp" & means.DtP$Richness >1),]
p5 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[1], shape=symb[1], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(E)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p5

d.temp <- means.DtP[which(means.DtP$Species=="Hz" & means.DtP$Richness >1),]
p6 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[2], shape=symb[2], cex=4) +
  #stat_smooth(method="glm", method.args = list(family = "binomial"), col="#6B244C", aes(group=1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(F)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p6

d.temp <- means.DtP[which(means.DtP$Species=="Px" & means.DtP$Richness >1),]
p7 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[3], shape=symb[3], cex=4) +
  #stat_smooth(method="glm", method.args = list(family = "binomial"), col="#6B244C", aes(group=1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(G)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p7

d.temp <- means.DtP[which(means.DtP$Species=="Sf" & means.DtP$Richness >1),]
p8 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[4], shape=symb[4], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(H)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p8

#For survival

d.temp <- S[which(S$Species=="Cp" & S$Richness >1),]
p9 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[1], shape=symb[1], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(I)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p9

d.temp <- S[which(S$Species=="Hz" & S$Richness >1),]
p10 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[2], shape=symb[2], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(J)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p10

d.temp <- S[which(S$Species=="Px" & S$Richness >1),]
p11 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[3], shape=symb[3], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(K)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p11

d.temp <- S[which(S$Species=="Sf" & S$Richness >1),]
p12 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[4], shape=symb[4], cex=4) +
  stat_smooth(method="glm", method.args = list(family = "binomial"), col="#6B244C", aes(group=1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(L)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p12


#For fungal growth

d.temp <- means.F[which(means.F$Fungi=="Botrys" & means.F$Richness >1),]
p13 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[1], shape=symb[1], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(M)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p13

d.temp <- means.F[which(means.F$Fungi=="Collet" & means.F$Richness >1),]
p14 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[2], shape=symb[2], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(N)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p14

d.temp <- means.F[which(means.F$Fungi=="Penicillium" & means.F$Richness >1),]
p15 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[3], shape=symb[3], cex=4) +
  #stat_smooth(method="glm", method.args = list(family = "binomial"), col="#6B244C", aes(group=1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(O)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p15

d.temp <- means.F[which(means.F$Fungi=="Sclerotinia" & means.F$Richness >1),]
p16 <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, col=pal[4], shape=symb[4], cex=4) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0","","1"), limits = c(-0.1,1.1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10), labels=c("2","4","6", "8", "10")) +
  labs(subtitle="(P)") +
  theme(plot.margin = unit(c(0,1,0,1), "cm")) +
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p16





#To get a legend for insects
d.temp <- means.PW[which(means.PW$Richness !=1),]
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))
p.temp <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, aes(shape=Species, col=Species), size=8) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p.temp
l1 <- get_legend(p.temp)




#legend for fungi
d.temp <- means.F[which(means.F$Richness !=1),]
labs <- c(expression(paste(italic(" B. dothidea   "))),
          expression(paste(italic(" Colletotrichum   "))),
          expression(paste(italic(" P. expansum  "))), 
          expression(paste(italic(" S. sclerotiorum  "))))
p.temp <- ggplot(d.temp, aes(Richness, exceeds)) +
  geom_jitter(height=0.05, width=0.3, aes(shape=Fungi, col=Fungi), size=8) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(text = element_text(size = 35), axis.title.x=element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p.temp
l2 <- get_legend(p.temp)


s1 <- segmentsGrob(x0 = 0, y0 = 0.5,
                   x1 = 1, y1 = 0.5, gp = gpar(col="darkslategray"))
t1 <- textGrob("Insect Pupal Mass", gp=gpar(fontsize=35, fontface="bold"))
t2 <- textGrob("Insect Development Speed", gp=gpar(fontsize=35, fontface="bold"))
t3 <- textGrob("Insect Survival", gp=gpar(fontsize=35, fontface="bold"))
t4 <- textGrob("Richness", gp=gpar(fontsize=35))
t5 <- textGrob("Fungal Growth Rate", gp=gpar(fontsize=35, fontface="bold"))
t6 <- textGrob("Insect Species (Plots A-L)", gp=gpar(fontsize=35))
t7 <- textGrob("Fungal Species (Plot M-P)", gp=gpar(fontsize=35))


lay <- rbind(c(1,1,1,1), 
             c(2,3,4,5),  
             c(6,6,6,6),  
             c(7,7,7,7),
             c(8,8,8,8),
             c(9,10,11,12),
             c(13,13,13,13),
             c(14,14,14,14),
             c(15,15,15,15),
             c(16,17,18,19),
             c(20,20,20,20),
             c(21,21,21,21),
             c(22,22,22,22),
             c(23,24,25,26),
             c(27,27,27,27),
             c(28,28,28,28),
             c(29,29,29,29),
             c(30,30,30,30),
             c(31,31,31,31)) 
tiff(height=1980, width=1800, filename="./Outputs/Figures/FigS3_ProbBestSingle.tiff", type="cairo")
grid.arrange(t1,p1,p2,p3,p4,t4,s1,t2,p5,p6,p7,p8,t4,s1,t3,p9,p10,p11,p12,t4,s1,t5,
             p13,p14,p15,p16,t4,t6,l1,t7,l2, nrow=19, ncol=4, 
             widths=c(300,300,300,300),
             heights=c(50,300,50,10,50,300,50,10,50,300,50,10,50,300,50,50,50,50,50), layout_matrix=lay)
dev.off()


cairo_pdf(height=23.76, width=21.6, family="sans",filename="./Outputs/Figures/FigS3_ProbBestSingle.pdf")
grid.arrange(t1,p1,p2,p3,p4,t4,s1,t2,p5,p6,p7,p8,t4,s1,t3,p9,p10,p11,p12,t4,s1,t5,
             p13,p14,p15,p16,t4,t6,l1,t7,l2, nrow=19, ncol=4, 
             widths=c(300,300,300,300),
             heights=c(50,300,50,10,50,300,50,10,50,300,50,10,50,300,50,50,50,50,50), layout_matrix=lay)
dev.off()



#-----------------------------------------------------------------------
## Fig. S4: Effects of each treatment on each herbivore
#-------------------------------

##Load workspace with standardized data
load("./Outputs/Workspaces/StandardizedData")

#For pupal weights
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))

#re-order levels in the way we want them for the figure
d.temp$Treatment <- factor(d.temp$Treatment,levels=c("Pht", "Phz", "Q", "H", "R",  "eCt", "Ct",
                                                     "GeA",  "SA", "GA", "FA", "pCA", "CA", "ChA", 
                                                     "L2A", "L2B", "L2C", "M2A", "M2B", "M2C", "H2A", "H2B", "H2C",
                                                     "L4A", "L4B", "L4C", "M4A", "M4B", "M4C", "H4A", "H4B", "H4C",
                                                     "L6A", "L6B", "L6C", "M6A", "M6B", "M6C", "H6A", "H6B", "H6C",
                                                     "L8A", "L8B", "L8C", "M8A", "M8B", "M8C", "H8A", "H8B", "H8C",
                                                     "L10A", "L10B", "L10C", "M10A", "M10B", "M10C", "H10A", "H10B", "H10C"))
d.temp$Treatment <- factor(d.temp$Treatment,levels=rev(levels(d.temp$Treatment)))

#Pupal weights
p1 <- ggplot(d.temp, aes(x = Treatment, y = Pupal.weight.ST)) +
  geom_point(aes(col=Species), alpha=0.3, size=0.8) +
  stat_summary(fun.y = "mean", data=d.temp, aes(color=Species), size = 1.75, geom = "point") +
  #geom_boxplot(width=0.75, notch = FALSE, aes(fill=Species, alpha = alpha_b), position=position_dodge(0.9)) +
  coord_flip() 
p1b <- p1 +
  scale_color_manual(values=pal, labels=labs)+
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), panel.border = element_blank())+
  geom_hline(yintercept=1, color="gray20", lwd=0.5, lty=1)+
  labs(subtitle="(A)", y="Pupal Mass", x="Treatment") +
  theme(axis.text.y=element_text(size=8))
p1b


#Days to pupation
p2 <- ggplot(d.temp, aes(x = Treatment, y = Days.to.pupation.ST.inv)) +
  geom_point(aes(col=Species), alpha=0.3, size=0.8) +
  stat_summary(fun.y = "mean", data=d.temp, aes(color=Species), size = 1.75, geom = "point") +
  coord_flip() 
p2b <- p2 +
  scale_color_manual(values=pal, labels=labs)+
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), panel.border = element_blank())+
  geom_hline(yintercept=1, color="gray20", lwd=0.5, lty=1)+
  labs(subtitle="(B)", y="Development Speed", x="") +
  theme(axis.text.y=element_text(size=8))
p2b


#Survival
d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
#re-order levels in the way we want them for the figure
d.temp$Treatment <- factor(d.temp$Treatment,levels=c("Pht", "Phz", "Q", "H", "R",  "eCt", "Ct",
                                                     "GeA",  "SA", "GA", "FA", "pCA", "CA", "ChA", 
                                                     "L2A", "L2B", "L2C", "M2A", "M2B", "M2C", "H2A", "H2B", "H2C",
                                                     "L4A", "L4B", "L4C", "M4A", "M4B", "M4C", "H4A", "H4B", "H4C",
                                                     "L6A", "L6B", "L6C", "M6A", "M6B", "M6C", "H6A", "H6B", "H6C",
                                                     "L8A", "L8B", "L8C", "M8A", "M8B", "M8C", "H8A", "H8B", "H8C",
                                                     "L10A", "L10B", "L10C", "M10A", "M10B", "M10C", "H10A", "H10B", "H10C"))
d.temp$Treatment <- factor(d.temp$Treatment,levels=rev(levels(d.temp$Treatment)))

p3 <- ggplot(d.temp, aes(x = Treatment, y = PropSurv.ST)) +
  geom_point(aes(col=Species), size=1.75) +
  #stat_summary(fun.y = "mean", data=d.temp, aes(color=Species), size = 3, geom = "point") +
  coord_flip() 
p3b <- p3 +
  scale_color_manual(values=pal, labels=labs)+
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), panel.border = element_blank())+
  geom_hline(yintercept=1, color="gray20", lwd=0.5, lty=1)+
  labs(subtitle="(C)", y="Survival", x="") +
  theme(axis.text.y=element_text(size=8))
p3b


#Fungal Growth
d.temp <- filter (d.fungi2, Treatment != "DMSO", !is.na(dAbs.ST))

#re-order levels in the way we want them for the figure
d.temp$Treatment <- factor(d.temp$Treatment,levels=c("Pht", "Phz", "Q", "H", "R",  "eCt", "Ct",
                                                     "GeA",  "SA", "GA", "FA", "pCA", "CA", "ChA", 
                                                     "L2A", "L2B", "L2C", "M2A", "M2B", "M2C", "H2A", "H2B", "H2C",
                                                     "L4A", "L4B", "L4C", "M4A", "M4B", "M4C", "H4A", "H4B", "H4C",
                                                     "L6A", "L6B", "L6C", "M6A", "M6B", "M6C", "H6A", "H6B", "H6C",
                                                     "L8A", "L8B", "L8C", "M8A", "M8B", "M8C", "H8A", "H8B", "H8C",
                                                     "L10A", "L10B", "L10C", "M10A", "M10B", "M10C", "H10A", "H10B", "H10C"))
d.temp$Treatment <- factor(d.temp$Treatment,levels=rev(levels(d.temp$Treatment)))

p4 <- ggplot(d.temp, aes(x = Treatment, y = dAbs.ST)) +
  geom_point(aes(col=Fungi), alpha=0.3, size=0.8) +
  stat_summary(fun.y = "mean", data=d.temp, aes(color=Fungi), size = 1.75, geom = "point") +
  coord_flip() 
p4b <- p4 +
  scale_color_manual(values=pal, labels=labs)+
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), panel.border = element_blank())+
  geom_hline(yintercept=1, color="gray20", lwd=0.5, lty=1)+
  labs(subtitle="(D)", y="Growth Rate", x="") +
  theme(axis.text.y=element_text(size=8))
p4b



#to get a legend for insects
d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))
p.temp <- ggplot(d.temp, aes(x = Treatment, y = PropSurv.ST)) +
  geom_point(aes(col=Species), size=3) +
  coord_flip() 
p.temp.b <- p.temp +
  scale_color_manual(values=pal, labels=labs)+
  theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=12),
        legend.key = element_rect(colour = "transparent", fill = "white")) +
  ylab("Standardized Survival")+
  geom_hline(yintercept=1, color="gray20", lwd=1, lty=1)
l1 <- get_legend(p.temp.b)


#to get a legend for fungi
d.temp <- filter (d.fungi2, Treatment != "DMSO")
labs <- c(expression(paste(italic(" B. dothidea   "))),
          expression(paste(italic(" Colletotrichum   "))),
          expression(paste(italic(" P. expansum  "))), 
          expression(paste(italic(" S. sclerotiorum  "))))
p.temp <- ggplot(d.temp, aes(x = Treatment, y = dAbs.ST)) +
  geom_point(aes(col=Fungi), size=3) +
  coord_flip() 
p.temp.b <- p.temp +
  scale_color_manual(values=pal, labels=labs)+
  theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=12),
        legend.key = element_rect(colour = "transparent", fill = "white")) +
  ylab("Fungal Growth Rate")+
  geom_hline(yintercept=1, color="gray20", lwd=1, lty=1)
l2 <- get_legend(p.temp.b)



#Make the composite fig

t1 <- textGrob("Insect Species (Plots A, B, C)", gp=gpar(fontsize=12))
t2 <- textGrob("Fungal Species (Plot D)", gp=gpar(fontsize=12))


tiff(height=720, width=600,filename="./Outputs/Figures/FigS4_treatment_X_herbivore.tiff", type="cairo")
lay <- rbind(c(1,2,3,4),  #graphs
             c(5,5,5,5),
             c(6,6,6,6),
             c(7,7,7,7),
             c(8,8,8,8))  #legend
grid.arrange(p1b, p2b, p3b, p4b,t1,l1,t2,l2, nrow=5, ncol=4, widths=c(500, 500, 500, 500),
             heights=c(1000, 50, 50, 50, 50), layout_matrix=lay)
dev.off()




cairo_pdf(height=10.368, width=8.64, family="sans",filename="./Outputs/Figures/FigS4_treatment_X_herbivore.pdf")
lay <- rbind(c(1,2,3,4),  #graphs
             c(5,5,5,5),
             c(6,6,6,6),
             c(7,7,7,7),
             c(8,8,8,8))  #legend
grid.arrange(p1b, p2b, p3b, p4b,t1,l1,t2,l2, nrow=5, ncol=4, widths=c(500, 500, 500, 500),
             heights=c(1000, 50, 50, 50, 50), layout_matrix=lay)
dev.off()




##-----------------------------------------------
##  Fig. S5: Evenness effects on performance
##------------------------------

##Load workspace with standardized data
load("./Outputs/Workspaces/StandardizedData")

#For Pupal Weights
d.temp <- filter (d.even, Treatment != "C", !is.na(Pupal.weight.ST))

p1 <- ggplot(d.temp, aes(x=Evenness, y=Pupal.weight.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = alpha_b), outlier.shape = NA)

p1.b <- p1 + ylab("Pupal Mass") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(A)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p1.b


#For Dev Speed
d.temp <- filter (d.even, Treatment != "C", !is.na(Days.to.pupation.ST.inv))

p2 <- ggplot(d.temp, aes(x=Evenness, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, height=0.03, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = alpha_b), outlier.shape = NA)

p2.b <- p2 + ylab("Development Speed") + xlab("Evenness") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none")+
  labs(subtitle="(B)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2.b

#For Survival
d.temp <- filter (d.even.PS, Treatment != "C", !is.na(PropSurv.ST))
p3 <- ggplot(d.temp, aes(x=Evenness, y=PropSurv.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, aes(shape=Species, col=Species), alpha=alpha_p) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = alpha_b), outlier.shape = NA)

p3.b <- p3 + ylab("Survival") + xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  scale_shape_manual(values = symb) + 
  scale_color_manual(values = pal) +
  theme(legend.position="none") +
  labs(subtitle="(C)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3.b

#To get a legend for insects
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))

p.temp <- ggplot(d.temp, aes(x=Richness, y=Pupal.weight.ST)) +
  geom_jitter(width=0.3, height=0.3, aes(shape=Species, col=Species), size=6)

p.temp.b <- p.temp + ylab("Standardized Pupal Weight") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p.temp.b)
p.temp.b



#Make the composite fig

tiff(height=550, width=1500, filename="./Outputs/Figures/FigS5_Evenness.tiff", type="cairo")
lay <- rbind(c(1,2,3),  #graphs
             c(4,4,4)) #legend
grid.arrange(p1.b, p2.b, p3.b, l1, nrow=2, ncol=3, widths=c(500, 500, 500),
             heights=c(500, 50), layout_matrix=lay)
dev.off()


cairo_pdf(height=7.33, width=20, family="sans",filename="./Outputs/Figures/FigS5_Evenness.pdf")
lay <- rbind(c(1,2,3),  #graphs
             c(4,4,4)) #legend
grid.arrange(p1.b, p2.b, p3.b, l1, nrow=2, ncol=3, widths=c(500, 500, 500),
             heights=c(500, 50), layout_matrix=lay)
dev.off()




#----------------------------------------------
#  Fig. S6: Composite fig of all places there are significant effects of evenness 
#     or SD on individual herbivore species and metrics
#-------------------------------

## (A)  Evenness effects on Sf male pupal weight

d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Sf" & d.even$Sex=="m"),]
p1 <- ggplot(d.temp, aes(x=Evenness, y=Pupal.weight.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, color=pal[4], shape=symb[4], cex=2) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = 0.8), outlier.shape = NA)
p1b <- p1 + ylab("M Pupal Weight") +
  xlab("")+
  geom_smooth(method="lm", col=pal[4], aes(group=1), fill=NA) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  #scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), labels=c("0","0.2","0.4","0.6", "0.8", "1", "1.2"), limits=c(0,1.2)) +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  labs(subtitle="(A)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p1b


## (B)  Evenness effects on Px male dev speed


d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Px" & d.even$Sex=="m"),]
p2 <- ggplot(d.temp, aes(x=Evenness, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, color=pal[3], shape=symb[3], cex=2) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = 0.8), outlier.shape = NA)
p2b <- p2 + ylab("M Dev Speed") + xlab("Evenness") +
  geom_smooth(method="lm", col=pal[3], aes(group=1), fill=NA) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  #scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), labels=c("0","0.2","0.4","0.6", "0.8", "1", "1.2"), limits=c(0,1.2)) +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  labs(subtitle="(B)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p2b




## (C)  Evenness effects on Px male pupal weight


d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Px" & d.even$Sex=="m"),]
p3 <- ggplot(d.temp, aes(x=Evenness, y=Pupal.weight.ST)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.03, color=pal[3], shape=symb[3], cex=2) +
  geom_boxplot(width=0.075, notch = FALSE, aes(group = cut_width(Evenness, 0.2), alpha = 0.8), outlier.shape = NA)
p3b <- p3 + ylab("M Pupal Weight") +
  xlab("") +
  geom_smooth(method="lm", col=pal[3], aes(group=1), fill=NA) +
  scale_x_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  #scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), labels=c("0","0.2","0.4","0.6", "0.8", "1", "1.2"), limits=c(0,1.2)) +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  labs(subtitle="(C)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p3b

### (D) SD effects on Hz female development 

d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Hz" & d.even$Sex=="f"),]
d.temp$SD <- factor(d.temp$SD,levels(d.temp$SD)[c(2,3,1)]) #re-order factor levels
p4 <- ggplot(d.temp, aes(x=SD, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.2, color=pal[2], shape=symb[2], cex=2) +
  geom_boxplot(width=0.5, notch = FALSE, alpha=0.5, outlier.shape = NA)
p4b <- p4 + ylab("F Dev Speed") +
  xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Low", "Med", "High")) +
  scale_y_continuous(limits=c(0.6,1.2)) +
  annotate("text", label = "A", x = 1, y = 1.15, size=10) +
  annotate("text", label = "A", x = 2, y = 1.15, size=10) +
  annotate("text", label = "B", x = 3, y = 1.15, size=10)+
  labs(subtitle="(D)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p4b

### (E) SD effects on Sf female development 

d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Sf" & d.even$Sex=="f"),]
d.temp$SD <- factor(d.temp$SD,levels(d.temp$SD)[c(2,3,1)]) #re-order factor levels
p5 <- ggplot(d.temp, aes(x=SD, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.2, color=pal[4], shape=symb[4], cex=2) +
  geom_boxplot(width=0.5, notch = FALSE, alpha=0.5, outlier.shape = NA)
p5b <- p5 + ylab("F Dev Speed") +
  xlab("Structural Diversity") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Low", "Med", "High")) +
  scale_y_continuous(limits=c(0.25,1.25)) +
  annotate("text", label = "A", x = 1, y = 1.2, size=10) +
  annotate("text", label = "A", x = 2, y = 1.2, size=10) +
  annotate("text", label = "B", x = 3, y = 1.2, size=10)+
  labs(subtitle="(E)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p5b



### (F) SD effects on Hz male development 

d.temp <- d.even[which(d.even$Treatment != "C" & d.even$Species=="Hz" & d.even$Sex=="m"),]
d.temp$SD <- factor(d.temp$SD,levels(d.temp$SD)[c(2,3,1)]) #re-order factor levels
p6 <- ggplot(d.temp, aes(x=SD, y=Days.to.pupation.ST.inv)) +
  geom_hline(yintercept=1, color="gray20", lwd=1.2, lty=3) +
  geom_jitter(width=0.2, color=pal[2], shape=symb[2], cex=2) +
  geom_boxplot(width=0.5, notch = FALSE, alpha=0.5, outlier.shape = NA)
p6b <- p6 + ylab("M Dev Speed") +
  xlab("") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Low", "Med", "High")) +
  scale_y_continuous(limits=c(0.6,1.2)) +
  annotate("text", label = "A", x = 1, y = 1.15, size=10) +
  annotate("text", label = "A", x = 2, y = 1.15, size=10) +
  annotate("text", label = "B", x = 3, y = 1.15, size=10)+
  labs(subtitle="(F)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
p6b


#To get a legend for insects
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
labs <- c(expression(paste(italic(" C.pomonella   "))),
          expression(paste(italic(" H. zea   "))),
          expression(paste(italic(" P. xylostella  "))), 
          expression(paste(italic(" S. frugiperda  "))))

p.temp <- ggplot(d.temp, aes(x=Richness, y=Pupal.weight.ST)) +
  geom_jitter(width=0.3, height=0.3, aes(shape=Species, col=Species), size=6)

p.temp.b <- p.temp + ylab("Standardized Pupal Weight") +
  theme(text = element_text(size = 35), plot.margin=unit(c(1,1,1,1), "cm")) +
  scale_shape_manual(values = symb, labels=labs) + 
  scale_color_manual(values = pal, labels=labs) +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white"))
l1 <- get_legend(p.temp.b)
p.temp.b




tiff(height=1050, width=1500, filename="./Outputs/Figures/FigS6_EvennessIndivEffects.tiff", type="cairo")
lay <- rbind(c(1,2,3),  #graphs
             c(4,5,6),  #graphs
             c(7,7,7))   #legend
grid.arrange(p1b, p2b, p3b, p4b, p5b, p6b, l1,
             nrow=3, ncol=3, widths=c(500,500, 500), 
             heights=c(500, 500, 50),
             layout_matrix=lay)
dev.off()

cairo_pdf(height=11.9, width=17, filename="./Outputs/Figures/FigS6_EvennessIndivEffects.pdf")
lay <- rbind(c(1,2,3),  #graphs
             c(4,5,6),  #graphs
             c(7,7,7))   #legend
grid.arrange(p1b, p2b, p3b, p4b, p5b, p6b, l1,
             nrow=3, ncol=3, widths=c(500,500, 500), 
             heights=c(500, 500, 50),
             layout_matrix=lay)
dev.off()




#-----------------------------------------------------------------------------
###Figure S7: effects of evenness on number of organisms affected 
#--------------------------

load("./Outputs/Workspaces/FinalAnalyses_Evenness_NumEffects")

plot(jitter(NumHerbsAffected) ~ jitter(Evenness), data=NumEffects.even)
abline(lm(NumHerbsAffected ~ Evenness, data=NumEffects.even))

d.temp <- NumEffects.even
p1 <- ggplot(d.temp, aes(x=Evenness, y=NumHerbsAffected)) +
  geom_jitter(width=0.03, height=0.15, cex=4, col=dot_col) +
  geom_smooth(method="lm", col="#6B244C", aes(group=1)) +
  ylab("Numb. Herbivores Affected") + xlab("Evenness") +
  scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1.0), labels=c("0.2","0.4","0.6", "0.8", "1.0")) +
  theme(text = element_text(size = 35))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1




tiff(height=500, width=600, filename="./Outputs/Figures/FigS7_NumbEffects_even.tiff", type="cairo")
p1
dev.off()

cairo_pdf(height=8, width=9, family="sans",filename="./Outputs/Figures/FigS7_NumbEffects_even.pdf")
p1
dev.off()





