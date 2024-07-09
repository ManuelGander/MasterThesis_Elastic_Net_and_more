library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

path_to_file = '/media/kusterlab/internal_projects/active/TOPAS/Publications/Manuscript_MTB/2_figures/Figure2/Content_Construction/2J/Figure2J_input.xlsx'


input = xlsx::read.xlsx(path_to_file,sheetIndex = 2)
input = input[,c("Sample.name",
                 "num_significant_baskets_z.1.5",
             
                 "Tissue.constraints" ,
                 "num_significant_other_baskets_z.1.5"
                 
)]

input$order = 0
input$order[input$num_significant_other_baskets_z.1.5 == 'x' &
            input$num_significant_baskets_z.1.5 == 'x'       &
            input$Tissue.constraints != 'x'] = 1 #b.g


input$order[input$num_significant_other_baskets_z.1.5 == 'x' &
              input$num_significant_baskets_z.1.5 != 'x'       &
              input$Tissue.constraints != 'x'] = 2 #g


input$order[input$num_significant_other_baskets_z.1.5 == 'x' &
              input$num_significant_baskets_z.1.5 != 'x'       &
              input$Tissue.constraints == 'x'] = 3



input$order[input$num_significant_other_baskets_z.1.5 == 'x' &
              input$num_significant_baskets_z.1.5 == 'x'       &
              input$Tissue.constraints == 'x'] = 4 #r,g,b




input$order[input$num_significant_other_baskets_z.1.5 != 'x' &
              input$num_significant_baskets_z.1.5 == 'x'       &
              input$Tissue.constraints == 'x'] = 5



input$order[input$num_significant_other_baskets_z.1.5 != 'x' &
              input$num_significant_baskets_z.1.5 != 'x'       &
              input$Tissue.constraints == 'x'] = 6


input$order[input$num_significant_other_baskets_z.1.5 != 'x' &
              input$num_significant_baskets_z.1.5 != 'x'       &
              input$Tissue.constraints != 'x'] = 7   


input$order[input$num_significant_other_baskets_z.1.5 != 'x' &
              input$num_significant_baskets_z.1.5 == 'x'       &
              input$Tissue.constraints != 'x'] = 8  


input <- input[order(input$order,
                    #input$num_significant_other_baskets_z.1.5,
                    #input$num_significant_baskets_z.1.5,
                    decreasing = T
                     ), ]

#rownames(mat) = mat$Sample_name
mat = input[,! colnames(input) %in% c('Sample.name','order')]
mat$Tissue.constraints[mat$Tissue.constraints == 'x'] = 'Tissue_constrains'
mat$num_significant_baskets_z.1.5[mat$num_significant_baskets_z.1.5 == 'x'] = 'sig baskets z> 1.5'
mat$num_significant_other_baskets_z.1.5[mat$num_significant_other_baskets_z.1.5 == 'x'] = 'sig other basket z >1.5'

#mat = matrix(runif(200),4,10)
#colnames(mat) = paste0("col",1:ncol(mat))
rownames(mat) = paste0("p",1:nrow(mat))
mat  = t(mat)
df = data.frame(mat) %>% 
  rownames_to_column("row") %>% 
  pivot_longer(-row) %>%
  mutate(name=factor(name,levels=colnames(mat)),
         row=factor(row,levels=rownames(mat)))
df$value[df$value=='nan']=''
row_num = length(levels(df$row))


ggplot(df,aes(x=name,y=as.numeric(row),fill=value)) + 
  xlim(c("",colnames(mat))) + 
  ylim(c(-row_num,row_num+1))+
  geom_tile()+ ylab("")+
  #scale_fill_discrete(drop=TRUE,limits = levels(df$value)) +
  scale_fill_manual(values = c("grey",'blue','green','red')) +
  #annotate(x="",y=1:row_num,label=levels(df$row),size=2,geom="text") + 
  coord_polar(start=0, direction = 1) + 
  theme_minimal() +  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = c(0.5, 0.5),
        legend.key.size = unit(0.2, "cm"))

ggsave("~/Desktop/piechart.svg")
