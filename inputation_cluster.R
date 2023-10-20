
rm(list = ls())
library(zCompositions)
library(ComplexHeatmap)
library(RColorBrewer)

lod_imput <- function(mat) ( low_dl <- c(0.1,15,0.2,0.2,0.1,2,1,9,5,35,30,3.16,2,67)
multRepl(mat,label=NA,dl=low_dl)


Progestrone_RW <- read_excel("Progestrone_RW.xlsx")


dt_s=Progestrone_RW %>% filter(性别=="女")%>% select(c("prog","E2","FSH","LH","TESTO","PRL","ALT","AST","GGT","ALB","Cr","Glu","wbc","Hb" ))
#colnames(dt_f)=c("耐受","呕吐","呕吐时期","脓毒症","鼻饲","腹泻","CRRT","死亡")

dt_f = dt_s %>% filter(!rowSums(is.na(.)) > 7) 


low_dl <- c(0.1,15,0.2,0.2,0.1,2,1,9,5,35,30,3.16,2,67)
dt_f_multRepl <- multRepl(dt_f,label=NA,dl=low_dl)

##z-score normalization

z_scores<-function(data) (ifelse(!is.na(data), abs(data[!is.na(data)]-mean(data[!is.na(data)]))/sd(data[!is.na(data)]),NA ))
dt_f_z <-  dt_f %>% mutate_all(funs(z_score = z_scores(.))) %>% select(ends_with("z_score")) %>% mutate_all(function(data)(ifelse(data>1.5,NA,data)))
dt_f_multRepl_z <-  dt_f_multRepl %>% mutate_all(funs(z_score = z_scores(.))) %>% select(ends_with("z_score")) %>% mutate_all(function(data)(ifelse(data>1.5,NA,data)))


h1 = ComplexHeatmap::Heatmap(dt_f_z,
                        column_title = paste0("Progestrone_female_all","Raw_Value"),  
                        col = rev(brewer.pal(10,"RdBu"))
                        )

h2 = ComplexHeatmap::Heatmap(dt_f_multRepl_z,
                        column_title = paste0("Progestrone_female_all","Inputation_Value"),  
                        col = rev(brewer.pal(10,"RdBu"))
                        )















outpdf=paste("res","_profile.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

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


dev.off()
