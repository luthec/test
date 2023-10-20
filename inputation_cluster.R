rm(list = ls())
library(zCompositions)
library(ComplexHeatmap)
library(RColorBrewer)

##lod & NA imputation
lod_imput <- function(mat) {  
low_dl <- c(0.1,15,0.2,0.2,0.1,2,1,9,5,35,30,3.16,2,67)
multRepl(mat,label=NA,dl=low_dl) }

##z-score normalization & heatmap
z_scores <- function(data) (ifelse(!is.na(data), abs(data[!is.na(data)]-mean(data[!is.na(data)]))/sd(data[!is.na(data)]),NA ))

z_heatmap <- function(mat,f_name,v_name){
mat_z <-  mat %>% mutate_all(z_scores) %>% mutate_all(function(data)(ifelse(data>1.5,NA,data)))
ComplexHeatmap::Heatmap(na.omit(mat_z),
                        column_title = paste0(f_name,v_name),  
                        col = rev(brewer.pal(10,"RdBu"))
                        )
}

z_heatmap_class <- function(mat,remove,f_name,v_name){
mat_z <-  mat %>% mutate(across(-remove,z_scores)) %>% mutate(across(-remove,function(data)(ifelse(data>1.5,NA,data))))
tmat=na.omit(mat_z)
ComplexHeatmap::Heatmap(tmat%>%select(-remove),
                        column_title = paste0(f_name,v_name),  
                        left_annotation = rowAnnotation(df=as.data.frame(tmat[,c('class')]),col=list(class=c('男性组' = 'blue','卵泡期'='lightblue','黄体中期'='yellow','排卵期'='red','绝经期'='grey'))),
                        row_split = tmat$class,
                        column_title_gp = gpar(fontsize = 8.2),
                        col = rev(brewer.pal(10,"RdBu"))
                        )
}


Progestrone_RW <- read_excel("Progestrone_RW.xlsx")

#  dt <- Progestrone_RW  %>% mutate(class =  ifelse(性别=="男","男性组","NA"))

dt <- Progestrone_RW  %>% mutate(class = case_when(性别=="男" ~ "男性组",
                                性别=="女" & 年龄>55 ~ "绝经期",
                                性别=="女" & !年龄>55 & LH > 19.18 ~ "排卵期",
                                性别=="女" & !年龄>55 & prog > 5 ~ "黄体中期",
                                性别=="女" & !年龄>55 & !LH > 19.18 & !prog > 5 ~ "卵泡期")) %>% 
                                select(c("性别","年龄","class","prog","E2","FSH","LH","TESTO","PRL","ALT","AST","GGT","ALB","Cr","Glu","wbc","Hb" ))


dt_f = dt %>% filter(!rowSums(is.na(.)) > 7) 

dt_f_multRepl <-cbind(dt_f$class,lod_imput(select(dt_f,-c("性别","年龄","class"))))
colnames(dt_f_multRepl)[1]="class"

h01 = z_heatmap(dt[dt$class=="男性组",-c(1:3)],"男性","_Raw_Value")
h02 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="男性组",-1],"男性","_Inputation_Value")

h1 = z_heatmap(dt[dt$class=="排卵期",-c(1:3)],"排卵期","_Raw_Value")
h2 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="排卵期",-1],"排卵期","_Inputation_Value")

h3 = z_heatmap(dt[dt$class=="黄体中期",-c(1:3)],"黄体中期","_Raw_Value")
h4 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="黄体中期",-1],"黄体中期","_Inputation_Value")

h5 = z_heatmap(dt[dt$class=="卵泡期",-c(1:3)],"卵泡期","_Raw_Value")
h6 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="卵泡期",-1],"卵泡期","_Inputation_Value")

h7 = z_heatmap(dt[dt$class=="绝经期",-c(1:3)],"绝经期","_Raw_Value")
h8 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="绝经期",-1],"绝经期","_Inputation_Value")

ho = z_heatmap_class(dt,colnames(dt)[1:3],"全图","_Raw_Value")
hi = z_heatmap_class(as_tibble(dt_f_multRepl),colnames(dt_f_multRepl)[1],"全图","_Inputation_Value")


outpdf=paste("res","_profile.pdf",sep='')
pdf(outpdf, width = 16, height = 10, family="GB1")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h01,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h02,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h7,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h8,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)

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

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h3,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h4,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(h5,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(h6,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(ho,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(hi,
        heatmap_legend_side = 'bottom',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
popViewport(0)


dev.off()
