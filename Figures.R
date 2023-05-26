#############################################
##                                         ##
## Reproduce figures in FishExtremeCatcher ##
##                                         ##
#############################################
library(ggplot2)
library(rgdal)
library(gridExtra)
load('datos.Rdata')
bathy = readGDAL("bathy.tiff")
load("fitstage1.Rdata")
load("fitstage2.Rdata")
load("fitstage3.Rdata")
blue = rgb(0, 18, 96, maxColorValue = 255)
# red = rgb(113,46,7, maxColorValue = 255)

ggtheme = theme(axis.title = element_text(size = 30),
                axis.text.x = element_text(size = 20, hjust = 1),
                axis.text.y = element_text(size = 20, hjust = 1),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 40, hjust = 0),
                legend.title = element_text(size = 20),
                legend.position = "right")


####################################
## Bathymetric effects (Figure 1) ##
####################################
df1 = stage1_fit$summary.random$bathy
df1$Bathymetry = mesh1d$loc

df2 = data.frame(Bathymetry   = datos$PROFUNDIDAD, 
                 mean         = datos$PROFUNDIDAD*stage2_fit$summary.fixed$mean[2],
                 `0.025quant` = datos$PROFUNDIDAD*stage2_fit$summary.fixed$`0.025quant`[2],
                 `0.975quant` = datos$PROFUNDIDAD*stage2_fit$summary.fixed$`0.975quant`[2],
                 check.names = FALSE)

df3 = stage3_fit$summary.random$bathy
df3$Bathymetry = mesh1d_ext$loc

p1 = ggplot(data = df1, aes(x = Bathymetry, y = mean)) +  geom_line(linewidth = 2, colour = blue) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), linetype = 2, alpha = 0.1, colour = blue, fill = blue) +
  ylab("Effect on linear predictor") + theme_minimal() + ggtitle("Bulk abundance") + ggtheme

p2 = ggplot(data = df2, aes(x = Bathymetry, y = mean)) +  geom_line(linewidth = 2, colour = blue) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), linetype = 2, alpha = 0.1, colour = blue, fill = blue) +
  ylab("Effect on linear predictor") + xlim(range(mesh1d$loc)) + theme_minimal() + ggtitle("Probability of ECE") + ggtheme

p3 = ggplot(data = df3, aes(x = Bathymetry, y = mean)) +  geom_line(linewidth = 2, colour = blue) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), linetype = 2, alpha = 0.1, colour = blue, fill = blue) +
  ylab("Effect on linear predictor") + xlim(range(mesh1d$loc)) + theme_minimal() + ggtitle("ECE abundance") + ggtheme

ggsave(filename = 'BathyEffect.pdf', arrangeGrob(p1, p2, p3, ncol = 3), width = 18, height = 6)

################################
## Temporal trends (Figure 2) ##
################################
df1 = stage1_fit$summary.random$year
df1$Year = unique(datos$year2)

df2 = data.frame(Year   = datos$year2, 
                 mean         = datos$year2*stage3_fit$summary.fixed$mean[2],
                 `0.025quant` = datos$year2*stage3_fit$summary.fixed$`0.025quant`[2],
                 `0.975quant` = datos$year2*stage3_fit$summary.fixed$`0.975quant`[2],
                 check.names = FALSE)

p1 = ggplot(data = df1, aes(x = Year, y = mean)) +  geom_line(linewidth = 2, colour = blue) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), linetype = 2, alpha = 0.1, colour = blue, fill = blue) +
  ylab("Effect on linear predictor") + theme_minimal() + ggtitle("Bulk abundance") + ggtheme

p2 = ggplot(data = df2, aes(x = Year, y = mean)) +  geom_line(linewidth = 2, colour = blue) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), linetype = 2, alpha = 0.1, colour = blue, fill = blue) +
  ylab("Effect on linear predictor") + theme_minimal() + ggtitle("Probability of ECE") + ggtheme


ggsave(filename = 'TemporalEffect.pdf', arrangeGrob(p1, p2, ncol = 2), width = 18, height = 6)

#########################################
## Exceedance probabilities (Figure 3) ##
#########################################

