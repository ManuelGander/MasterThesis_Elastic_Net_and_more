
# 1- swarm plot: biomarkers intensities for diffeernt normalizatio --------


# for the visualization of the Swarm and box plots for the Biomarkers
library(ggplot2)
library(ggbeeswarm)

df = read.csv('concatanated_swarm_inputs.csv')
df$type[df$type %in% "before_normalization"] = 'before'
df$type[df$type %in% "after_normalization"] = 'after'
# van der Corput sequence or Tukey texturing
ggplot(df, aes(x = factor(type,levels = c("before", "after")),
               y = intensity,
               color = patient_group)) +
  geom_quasirandom() +
  geom_boxplot(notch = TRUE) +
  facet_wrap(normalization_scenario ~ biomarker,
             ncol = length(unique(df$biomarker)),
             nrow = length(unique(df$normalization_scenario))
           )+
  theme(axis.text.x = element_text(angle = 45)) +
  xlab('Chordoma Biomarkers') + ylab('log10 (MAXLFQ_intensities)') + theme_bw()

# 2- Barplot : t_statistics --------------------------------------------------

df = read.csv('t_statistics_long.csv')
ggplot(df,aes(x=normalization_scenario,y=value,fill=normalization_scenario))+
  geom_bar(stat="identity")  +
  facet_wrap(. ~ biomarker,
               nrow = ceiling(length(unique(df$biomarker))/2),
              ncol = floor(length(unique(df$biomarker))/2)
             ) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Chordoma Biomarkers') + ylab('-log10 (p_adjusted_value)') 

