## Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Author: Klara Filek based on code provided by L. Kanjer

##--------------------load packages--------------------
library(maps) #v3.4.1
library(mapdata) #v2.3.1
library(tidyverse) #v2.0.0
library(ggplot2) #v3.4.2
library(ggrepel) #v0.9.3
library(ggtext) #v0.1.2

#------ Adriatic Sea maps -----
map(col="grey90", border = "grey40", fill = TRUE,
    xlim = c(10, 23), ylim = c(38, 47), mar = rep(0.1, 4))
map.scale(relwidth = 0.20, metric = TRUE, ratio = FALSE)
box()

# holobiont points
points(15.1617, 43.8895, col=1, bg = 2, pch=21, cex =1.5) #ID117 KARLO ALBANO
points(14.4720, 44.5316, col=1, bg = 2, pch=21, cex =1.5) #ID119 MARTIN
points(13.9369, 44.8217, col=1, bg = 2, pch=21, cex =1.5) #ID118 OLIVER RAUL
points(13.9225, 44.7681, col=1, bg = 2, pch=21, cex =1.5) #ID047 ZAL
points(17.6963, 42.8387, col=1, bg = 2, pch=21, cex =1.5) #ID056 SAMBA
points(15.7185, 39.9985, col=1, bg = 2, pch=21, cex =1.5) #ID057 ANGELO
points(16.5987, 41.2028, col=1, bg = 2, pch=21, cex =1.5) #ID068 KANOOH
points(16.5987, 41.2028, col=1, bg = 2, pch=21, cex =1.5) #ID069 KANFUS
points(16.5987, 41.2028, col=1, bg = 2, pch=21, cex =1.5) #ID070 FUTON
points(17.1333, 42.9625, col=1, bg = 2, pch=21, cex =1.5) #ID096 MARO
points(15.2586, 44.0931, col=1, bg = 2, pch=21, cex =1.5) #ID093 ELLA-RAVKA
points(13.8430, 44.8390, col=1, bg = 2, pch=21, cex =1.5) #ID097 FREEWINGS
points(13.8909, 44.7580, col=1, bg = 2, pch=21, cex =1.5) #ID098 MAKSIMUS
points(17.1342, 42.9584, col=1, bg = 2, pch=21, cex =1.5) #ID122 LUKA-AMADEO
points(16.1512, 41.3740, col=1, bg = 2, pch=21, cex =1.5) #ID071 COSMYN
points(16.9464, 42.9025, col=1, bg = 2, pch=21, cex =1.5) #ID010 MERRY FISHER
points(13.9930, 44.8164, col=1, bg = 2, pch=21, cex =1.5) #ID074 RYAN
points(13.9127, 44.8031, col=1, bg = 2, pch=21, cex =1.5) #ID073 MARVIN

#----- Europe map -----
map(col="grey80", border = "grey40", fill = TRUE,
    xlim = c(-25, 50), ylim = c(30, 70), mar = rep(0.1, 4))
map.scale(relwidth = 0.5, metric = TRUE, ratio = FALSE)
box()

#----- Population structure -----
# read turtle metadata
turtle_meta <- readr::read_tsv("master_metadata/collapsed_turtleid.tsv")

turtle_meta <- turtle_meta[1:18,]
turtle_meta <- turtle_meta %>% select(-c(SamplingEvent, AdmissionDate, EndoSamplingDate, EpiSamplingDate, ReleaseDate, DaysInRehabiliationBeforeSampling, HospStatus))

turtle_meta$AgeRange <- ordered(turtle_meta$AgeRange,
                                levels = c("JUVENILE", "SUB-ADULT", "ADULT"))

ccl_summary <- group_by(turtle_meta, AgeRange) %>%
  summarise(
    count = n(),
    mean = mean(CCL, na.rm = TRUE),
    sd = sd(CCL, na.rm = TRUE),
    median = median(CCL, na.rm = TRUE),
    IQR = IQR(CCL, na.rm = TRUE)
  )
write_tsv(ccl_summary, file = "r_output/turtles_CCL_summary.tsv")

weight_summary <- group_by(turtle_meta, AgeRange) %>%
  summarise(
    count = n(),
    mean = mean(Weight, na.rm = TRUE),
    sd = sd(Weight, na.rm = TRUE),
    median = median(Weight, na.rm = TRUE),
    IQR = IQR(Weight, na.rm = TRUE)
  )
write_tsv(weight_summary, file = "r_output/turtles_weight_summary.tsv")

theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

theme_pop <- theme(
  strip.text.x = element_text(size = 8, color = "black", face = "bold"),
  strip.background = element_rect(color = "black"),
  panel.grid = element_line(colour = "gray90"),
  panel.grid.major.x = element_line(colour = "gray90"),
  panel.border = element_rect(color = "black"),
  axis.ticks.x = element_line(colour = "black"), # << visibility
  axis.ticks.y = element_line(colour = "black"),
  axis.title.x = element_text(size = 9, color = "black"),
  legend.text = element_markdown(size = 8, color = "black"),
  axis.title.y = element_text(size = 9, color = "black"),
  axis.text.y = element_text(size = 8, color = "black"),
  axis.text.x = element_text(angle = 0, vjust = 0.5, 
                             hjust = 0.5, size = 8, 
                             color = "black"), 
  legend.title = element_text(size = 9, color = "black"),
  legend.spacing.y = unit(0.1, 'cm'),
  legend.key.size = unit(0.25, 'cm'),
  plot.tag = element_text(size = 10)
)
#pop_ccl_weight_summary <-
turtle_meta %>%
  ggplot(aes(x = CCL, y = Weight, fill = AgeRange, shape = TurtleSex)) + 
  geom_vline(xintercept = 60, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 70, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_point(size = 2,
             alpha = 0.8) +
  scale_shape_manual(values = c(24, 22, 21),
                     name = "Turtle sex:",
                     labels = c(FEMALE = "Female",
                                MALE = "Male",
                                ND = "Not determined")) +
  scale_fill_brewer(palette = "Set2",
                    name = "Age range:",
                    labels = c(JUVENILE = "Juvenile",
                               "SUB-ADULT" = "Subadult",
                               ADULT = "Adult")) +
  geom_text_repel(aes(label = SampleID), size = 2,force = 1, 
                  min.segment.length = unit(0, 'lines'),
                  nudge_y = 0.3) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical", override.aes=list(shape=21))) +
  labs( 
    x = "CCL (cm)",
    y = "Weight (kg)",
    tag = "b)"
  ) + 
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70), limits = c(-5, 62), expand = c(0,0)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 70, 80), limits = c(20, NA),) +
  theme_pop

ggsave(file = "r_output/ccl_weight_summary_pub.pdf",
       #plot = pop_ccl_weight_summary,
       device = cairo_pdf,
       height = 55,
       width = 85,
       units = "mm"
)

# Weight and CCL sumarry stats
metadata_turtle <- readr::read_tsv("master_metadata/2023holobiont_master_metadata.tsv")
metadata_turtle <- distinct(metadata_turtle %>% 
                              select(TurtleID, CCL, Weight, AgeRange)) %>% 
                      na.omit()

summary_stats_ccl <- group_by(metadata_turtle, AgeRange) %>%
  summarise(
    count = n(),
    mean = mean(CCL, na.rm = TRUE),
    sd = sd(CCL, na.rm = TRUE),
    median = median(CCL, na.rm = TRUE),
    IQR = IQR(CCL, na.rm = TRUE)
  )

summary_stats_weight <- group_by(metadata_turtle, AgeRange) %>%
  summarise(
    count = n(),
    mean = mean(Weight, na.rm = TRUE),
    sd = sd(Weight, na.rm = TRUE),
    median = median(Weight, na.rm = TRUE),
    IQR = IQR(Weight, na.rm = TRUE)
  )

capture.output(summary_stats_ccl, file = "r_output/turtle_pop_summary_stats_ccl.txt")
capture.output(summary_stats_weight, file = "r_output/turtle_pop_summary_stats_weight.txt")


#----- Session and citations -----
sessionInfo()
citation("maps")
citation("mapdata")
