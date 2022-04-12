
library(ggplot2)
library(tidyverse)
library(tinytex)
library(lubridate)
## plot up biomass and metal content in a readable way ----


#read in data
bio <- read.csv("2_incremental/20220404_STANDING_CROP.csv") #metal columns are duplicated and have weird numbers
metal<- read.csv("2_incremental/20220404_BB_CURRENT.csv")

#mess with data
metal_long <- pivot_longer(metal, Al:Zn_r, names_to = "Element", values_to = "concentration")
  
#only epilithic
peri <- metal %>%
  filter(SAMPLE_DESCRIPTOR == "EPIL")

peri_long <- metal_long %>%
  filter(SAMPLE_DESCRIPTOR == "EPIL")%>%
  filter(Element == "Cu"|Element =="Fe")%>%
  filter(SITE=="DL")
  
cu <- peri_long %>%
  filter(element=="Cu")

biomass <- bio %>%
  filter(SAMPLE_DESCRIPTOR=="EPIL")


peri_long<- peri_long %>%
  mutate(date = parse_date_time(peri_long$SAMPLING_DATE,dmy))

write.csv(peri_long,"2_incremental\20220404_quickplot.csv", row.names = FALSE)

peri_long_new <- read.csv("2_incremental_20404_quickplot.csv")

peri_long_new$month = factor(peri_long_new$month, levels = month.name)

#set theme
theme_set(theme_bw()+
            theme(axis.title.x = element_text( size=22),
                  panel.grid = element_blank(),
                  axis.title.y = element_text( size=22),
                  axis.text.x = element_text(size = 16),
                  axis.text.y = element_text(size = 18),
                  legend.title = element_text(size = 16),
                  legend.text=element_text(size=14),
                  strip.text = element_text(size = 16)))


####plots
#biomass

ggplot(data = biomass, mapping = aes(x = SAMPLING_DATE, y = OM.AREA.g.m2, fill=)) + 
  geom_boxplot()+
  facet_wrap(~SITE, scale = "free_y") +
  theme(axis.text.x = element_text(angle=-45, vjust = -0.5))+
  xlab("Sampling Date")+
  ylab("Periphyton Biomass (g AFDM /m^2)")





#metals

ggplot(data = peri, mapping = aes(x = SAMPLING_DATE, y = OM.AREA.g.m2)) + 
  geom_point(mapping = aes(color = SITE)) + 
  geom_smooth()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=-45, vjust = -0.5)))



ggplot(data = peri_long_new, mapping = aes(x = month, y = concentration, color = Element)) + 
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Element, scale = "free_y") +
  theme(axis.text.x = element_text(angle=-45, vjust = -0.5))+
  xlab("Month")+
  ylab("Periphyton Metal Concentration (mg/g)")
 

ggplot(data = peri_long, mapping = aes(x = SAMPLING_DATE, y = concentration, color = Element)) + 
  geom_boxplot()+
  facet_wrap(~SITE, scale = "free_y") +
  theme( axis.text.x = element_text(angle=-45, vjust = -0.5))+
  xlab("Sampling Date")+
  ylab("Metal Concentration of Periphyton (mg/g)")

ggplot(data = cu, mapping = aes(x = SAMPLING_DATE, y = concentration, color = SITE)) + 
  geom_boxplot()+
  facet_wrap(~SITE, scale = "free_y") +
  theme( axis.text.x = element_text(angle=-45, vjust = -0.5))+
  xlab("Sampling Date")+
  ylab("Copper Concentration of Periphyton (mg/g)")


#all 
ggplot(data = bio, mapping = aes(x = SAMPLING_DATE, y = OM.AREA.g.m2,color = SAMPLE_DESCRIPTOR)) + 
  geom_point() + 
  geom_smooth()


#box plots
ggplot(data = peri) + 
  geom_boxplot(mapping = aes(x = SAMPLING_DATE, y = OM.AREA.g.m2))+
  ylab("Periphyton Organic Matter (g/m$^2$)")+
  xlab("Sampling Date")


ggplot(data = bio) + 
  geom_boxplot(mapping = aes(x = SAMPLING_DATE, y = OM.AREA.g.m2))+
  facet_wrap(~SITE)
  ylab("Periphyton Organic Matter (g/m$^2$)")+
  xlab("Sampling Date")


#get rid of grid lines
  theme_set(theme_bw()+
              theme(axis.title.x = element_text( size=22),
                    panel.grid = element_blank(),
                    axis.title.y = element_text( size=22),
                    axis.text.x = element_text(size = 16),
                    axis.text.y = element_text(size = 18),
                    legend.title = element_text(size = 16),
                    legend.text=element_text(size=14),
                    strip.text = element_text(size = 16)))
  
  
  
  (size_fraction_plot_3 <- ggplot(ucfr_DTW %>% 
                                    filter(element %in% c( "Cu", "Fe")) %>%
                                    
                                    mutate(element = ordered(element, levels = c("Fe", 
                                                                                 "Cu" )))) +
      geom_boxplot(aes(x=month, y = concentration + 0.01, fill = size_fractions),  
                   outlier.shape = NA, width = 0.8) +
      scale_fill_viridis_d(name = "Size Fraction", 
                           label = c("SPM", "Colloidal", "Truly Dissolved"))+
      scale_y_log10(name = "Concentration (ppb)") +
      scale_x_discrete(name = "") +
      facet_wrap(~element, scale = "free_y") +
      theme(legend.position = "bottom", axis.text.x = element_text(angle=-45, vjust = -0.5)))
  



