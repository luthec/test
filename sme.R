library(readxl)
library(tidyverse)
library(tidyr)
library(rstatix)
library(reshape2)
library(ggplot2)
library(geepack)
library(lmerTest)
library(ggrepel)
library(ggpubr)
library(broom)
library(broom.mixed)
library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridGraphics)
library(GGally)
library(lme4)
library(compiler)
library(parallel)
library(boot)
library(lattice)



dt = read_excel("SME0510.xls",sheet = 1,col_names = T)
colnames(dt)[3]='No'

## '死亡概率','APACHEII评分' has many NA

dt_f = dt %>% select(c('No','是否诊断脓毒症','是否糖尿病','是否高血压','是否患有肿瘤','是否是ICU患者','病例分组','致病菌分类')) %>% mutate_at(vars(!matches("No")), as.factor)
colnames(dt_f)=c('No',"ifSepsis","ifDiabetes","ifHypertension","ifTumor","ifICU","Class","Pathogen")

##'presepsin' has many NA
#dt_kv = dt %>% select(c('No','PCT','IL-6','CRP','presepsin')) %>% mutate_at(vars(!matches("No")), as.double)

dt_kv = dt %>% select(c('No','PCT','IL-6','CRP')) %>% mutate_at(vars(!matches("No")), as.double)

dt_blv = dt %>% select(c('No',"WBC","NEU","NEU%","LYM","LYM%","MONO","MONO%","EOS","EOS%","BASO","BASO%","RDW-SD","RDW-CV","PLT","MPV","PDW")) %>% mutate_at(vars(!matches("No")), as.double)

dt_bcv = dt %>% select(c('No',"ALT","TB","DB","IB","TP","Alb","GLU","urea","Cr","UA","K","Na","Cl","Ca","RBC","Hb","HCT","MCV")) %>% mutate_at(vars(!matches("No")), as.double)

####basic plot

dt_test=merge(dt_kv,dt_f,by="No")

outpdf=paste("res","_pri.pdf",sep='')
pdf(outpdf, width = 16, height = 10)

ggpairs(dt_kv[, c("PCT","IL-6","CRP")])

ggplot(dt_test, aes(x = Class, y = ifSepsis)) +
  stat_sum(aes(size = ..n.., group = 1)) +
  scale_size_area(max_size=10)

tmp1 <- melt(dt_test[, c("Class","PCT","IL-6","CRP")], id.vars="Class")
ggplot(tmp1, aes(x = Class, y = value)) +
  geom_jitter(alpha = .1) +
  geom_violin(alpha = .75) +
  facet_grid(variable ~ .) +
  scale_y_sqrt()

tmp2 <- melt(dt_test[, c("ifSepsis","PCT","IL-6","CRP")],id.vars="ifSepsis")
ggplot(tmp2, aes(factor(ifSepsis), y = value, fill=factor(ifSepsis))) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free_y")

dev.off()



dt_kv2 <- gather(dt_kv, key="Marker", value='Value',-No)
p1 = ggplot(dt_kv2, aes(x=Marker,y = Value, fill=Marker ))  + 
    theme(panel.background = element_blank(),axis.line = element_line()) +
    geom_boxplot(alpha = 0.8)   + ylab(paste("Raw","value",sep="_")) + 
    geom_violin(alpha = 0.1) # + coord_cartesian(ylim = ylim1*2)


###########################data filter
z_scores<-function(data) (ifelse(!is.na(data), abs(data[!is.na(data)]-mean(data[!is.na(data)]))/sd(data[!is.na(data)]),NA ))

mad_outliers <- function(data,n){
  me = median(data, na.rm = TRUE)
  mad = mad(data, na.rm = TRUE)
  ifelse(data > me + n * mad | data < me - n * mad, NA , data)
 }

dt_kv_z <-  dt_kv %>% mutate_at(vars(!matches("No")),funs(z_score = z_scores(.)))

# dt_kv_no_outliers <- dt_kv %>% group_by(Treatment, conc) %>% identify_outliers("PCT") %>% filter(!is.extreme)

dt_kv_zf = dt_kv_z %>% select(ends_with("z_score")) %>% mutate_all(function(data)(ifelse(data>1.5,NA,data)))
dt_kv3 <- gather(dt_kv_zf, key="Marker", value='Value')
p2 = ggplot(dt_kv3, aes(x=Marker,y = Value, fill=Marker ))  + 
    theme(panel.background = element_blank(),axis.line = element_line()) +
    geom_boxplot(alpha = 0.8)   + ylab(paste("Zscore","value",sep="_")) + 
    geom_violin(alpha = 0.1) 



dt_kv_mad <- dt_kv %>% mutate_at(vars(!matches("No")),funs(mad_out = mad_outliers(.,3))) %>% select(c("No",ends_with("mad_out"))) %>% mutate_at(vars(!matches("No")),funs(z_score = z_scores(.))) 

dt_kv_madf = dt_kv_mad %>% select(c("No",ends_with("z_score"))) %>% mutate_at(vars(!matches("No")),function(data)(ifelse(data>1.5,NA,data)))

dt_kv4 <- gather(dt_kv_madf, key="Marker", value='Value',-No)
p3 = ggplot(dt_kv4, aes(x=Marker,y = Value, fill=Marker ))  + 
    theme(panel.background = element_blank(),axis.line = element_line()) +
    geom_boxplot(alpha = 0.8)   + ylab(paste("MAD_Zscore","value",sep="_")) + 
    geom_violin(alpha = 0.1) 



######################################

h0 = pheatmap(na.omit(dt_kv)[,-1])
dt_f2 = dt_f[dt_f$No%in%na.omit(dt_kv_madf)$No,-1]

row_ha = rowAnnotation(df=as.data.frame(dt_f2[,-1]),
                       col=list(Class=c('1' = 'blue','2'='red','3'='yellow','4'='blue'),
                       	        Pathogen=c('0'='white','1' = 'blue','2'='red','3'='yellow','4'='blue'),
                                ifICU=c('1'='pink','0'='darkgreen'),
                                ifTumor=c('1'='pink','0'='darkgreen'),
                                ifDiabetes=c('1'='pink','0'='darkgreen'),
                                ifHypertension=c('1'='pink','0'='darkgreen')                                
                                # Onset_admission=circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
                                )
                        )

h1 = ComplexHeatmap::Heatmap(na.omit(dt_kv_madf[,-1]),
                        column_title = paste0("Key_Value","_Heatmap_by_Sepsis"),  
                        right_annotation = row_ha, 
                        left_annotation = rowAnnotation(df=as.data.frame(dt_f2[,c('ifSepsis')]),col=list(ifSepsis=c('0' = 'blue','1'='red'))),
                        row_split = dt_f2$ifSepsis,
                        col = rev(brewer.pal(10,"RdBu"))
                        )


row_ha2 = rowAnnotation(df=as.data.frame(dt_f2[,-6]),
                       col=list(ifSepsis=c('1'='pink','0'='darkgreen'),
                       	        Pathogen=c('0'='white','1' = 'blue','2'='red','3'='yellow','4'='blue'),
                                ifICU=c('1'='pink','0'='darkgreen'),
                                ifTumor=c('1'='pink','0'='darkgreen'),
                                ifDiabetes=c('1'='pink','0'='darkgreen'),
                                ifHypertension=c('1'='pink','0'='darkgreen')                                
                                # Onset_admission=circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
                                )
                        )

h2 = ComplexHeatmap::Heatmap(na.omit(dt_kv_madf[,-1]),
                        column_title = paste0("Key_Value","_Heatmap_by_Class"),  
                        right_annotation = row_ha2, 
                        left_annotation = rowAnnotation(df=as.data.frame(dt_f2[,c('Class')]),col=list(Class=c('1' = 'blue','2'='red','3'='yellow','4'='blue'))),
                        row_split = dt_f2$Class,
                        col = rev(brewer.pal(10,"RdBu"))
                        )


###regression table

biomarker="Key_Value"

dt_select=inner_join(dt_kv_madf,dt_f,by="No")
colnames(dt_select)[2:4]=c("PCT","IL-6","CRP")
##dt_select2 <- gather(dt_select, key="Marker", value='Value',-ifSepsis,-No,-ifDiabetes,-ifHypertension,-ifTumor,-ifICU,-Class,-Pathogen)
dt_na_removed <- na.omit(dt_select)

fit.glm = glm(ifSepsis ~ ifDiabetes+ifHypertension+ifTumor+ifICU+Class+Pathogen+PCT+`IL-6`+CRP , family = binomial,data=dt_na_removed) %>%
    tidy(conf.int = TRUE) %>% 
    select(c("term","estimate","std.error","conf.low","conf.high","p.value")) %>% 
    mutate(signif = stars.pval(p.value))

main.title <- paste0(biomarker,"_GLM")
subtitle <- paste0("Generalized Linear Model") %>%
  strwrap(width = 80) %>%
  paste(collapse = "\n")
t.glm <- ggtexttable(fit.glm, theme = ttheme("light")) %>%
  tab_add_title(text = subtitle, face = "plain", size = 10) %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) %>%
  tab_add_footnote(text = "* means significant levels", size = 10, face = "italic")


fit.mem <- glmer(ifSepsis ~ Class+Pathogen+PCT+`IL-6`+CRP+ (1 | ifTumor)+ (1 | ifICU)+ (1 | ifDiabetes) + (1 | ifHypertension), data = dt_na_removed, family = binomial, control = glmerControl(optimizer = "bobyqa")) %>%
    tidy(conf.int = TRUE) %>% 
    select(c("term","estimate","std.error","conf.low","conf.high","p.value")) %>% 
    mutate(signif = stars.pval(p.value))

main.title <- paste0(biomarker,"_GMM")
subtitle <- paste0("Generalized linear mixed model") %>%
  strwrap(width = 80) %>%
  paste(collapse = "\n")
t.mem <- ggtexttable(fit.mem, theme = ttheme("light")) %>%
  tab_add_title(text = subtitle, face = "plain", size = 10) %>%
  tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) %>%
  tab_add_footnote(text = "* means significant levels", size = 10, face = "italic")

pcom = ggarrange(t.glm,t.mem, ncol = 2, nrow = 1,widths=c(1, 1))



outpdf=paste("res","_profile_new.pdf",sep='')
pdf(outpdf, width = 16, height = 10)


print(ggarrange(p1,p2, ncol = 2, nrow = 1,widths=c(1, 1)))


print(h0)


print(p3)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h1,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h2,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)

print(pcom)

dev.off()










