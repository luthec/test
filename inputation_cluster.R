rm(list = ls())
library(zCompositions)
library(ComplexHeatmap)
library(RColorBrewer)

##lod
all_low_dl <- c(0.1,15,0.2,0.2,0.1,2,1,9,5,35,30,3.16,2,67)

##lod imputation
lod_imput <- function(mat) {  
tar <- colnames(mat)
low_dl <- mat %>% select(tar)%>% filter(if_any(.cols = tar, ~ ! grepl("^\\d", .))) %>% 
                 na.omit() %>% summarise(across(everything(), min)) %>%
                 slice(1) %>% unlist(., use.names=FALSE) %>% readr::parse_number()
mat1 <- mat %>% mutate_all(~ replace(., grepl("<", .), 0)) %>% mutate(across(tar, readr::parse_number))
if (any(is.na(mat1)))
  {multReplus(mat1,dl=low_dl) 
  }else{
   multRepl(mat1,label=0,dl=low_dl,closure=10^6)  
  }
}

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
                        left_annotation = rowAnnotation(df=as.data.frame(tmat[,c('class')]),col=list(class=c('男性' = 'blue','未孕女性'='yellow','绝经期'='grey'))),
                        row_split = tmat$class,
                        column_title_gp = gpar(fontsize = 8.2),
                        col = rev(brewer.pal(10,"RdBu"))
                        )
}

Progestrone_RW <- read_excel("2021P_test.xlsx")

dt <- Progestrone_RW  %>% mutate(exclusion = case_when(
                            性别=="女" & ALT>40 ~ "排除",
                            性别=="男" & ALT>50 ~ "排除",
                            性别=="女" & AST>35 ~ "排除",
                            性别=="男" & AST>40 ~ "排除",
                            性别=="女" & GGT>45 ~ "排除",
                            性别=="男" & GGT>60 ~ "排除",
                            性别=="女" & Cr>81 ~ "排除",
                            性别=="男" & Cr>111 ~ "排除",
                            HCG>5 ~ "排除",
                            Glu>6.1 ~ "排除",
                            Hb<90 ~ "排除")) %>%
                              mutate(class = case_when(性别=="男" ~ "男性",
                                性别=="女" & 年龄>55 ~ "绝经期",
                                性别=="女" & !年龄>55 ~ "未孕女性")) %>% 
                                filter(is.na(exclusion)) %>%
                                select(c("性别","年龄","class","prog","E2","FSH","LH","TESTO"))


# dt %>% mutate(lod = if_any(-c("性别","年龄","class"),\(x) str_detect(x, regex(">", ignore_case = TRUE))))
# dt %>% filter(if_any(all_of(c("prog", "E2")), ~grepl("^\\d", .)))
# dt %>% filter(if_any(.cols = c("prog", "E2"), ~ ! grepl("^\\d", .)))
# dt |> mutate(across(-a, ~ . |> str_extract("\\d+") |> as.numeric()))
# dt %>% mutate(across(c("prog", "E2"), readr::parse_number))

# tb %>% mutate(across(everything(), ~ as.numeric(na_if(.x, "<LOD"))))
# dt %>%  replace(grepl("<", .), 0)



#  dt <- Progestrone_RW  %>% mutate(class =  ifelse(性别=="男","男性组","NA"))

dt_bec <- Progestrone_RW  %>% mutate(class = case_when(性别=="男" ~ "男性组",
                                性别=="女" & 年龄>55 ~ "绝经期",
                                性别=="女" & !年龄>55 & LH > 19.18 ~ "排卵期",
                                性别=="女" & !年龄>55 & prog > 5 ~ "黄体中期",
                                性别=="女" & !年龄>55 & !LH > 19.18 & !prog > 5 ~ "卵泡期")) %>% 
                                select(c("性别","年龄","class","prog","E2","FSH","LH","TESTO","PRL","ALT","AST","GGT","ALB","Cr","Glu","wbc","Hb" ))

dt_f = dt %>% filter(!rowSums(is.na(.)) > 3) 

dt2 = dt_f %>% mutate(across(c("prog","E2","FSH","LH","TESTO"), readr::parse_number))

dt_f_multRepl <-cbind(dt_f$class,lod_imput(select(dt_f,-c("性别","年龄","class"))))
colnames(dt_f_multRepl)[1]="class"

h01 = z_heatmap(dt2[dt2$class=="男性",-c(1:3,8)],"男性","_Raw_Value")
h02 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="男性",-c(1,6)],"男性","_Inputation_Value")

h7 = z_heatmap(dt2[dt2$class=="绝经期",-c(1:3,8)],"绝经期","_Raw_Value")
h8 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="绝经期",-c(1,6)],"绝经期","_Inputation_Value")

h1 = z_heatmap(dt2[dt2$class=="未孕女性",-c(1:3,8)],"未孕女性","_Raw_Value")
h2 = z_heatmap(dt_f_multRepl[dt_f_multRepl[,1]=="未孕女性",-c(1,6)],"未孕女性","_Inputation_Value")

ho = z_heatmap_class(dt2,colnames(dt2)[1:3],"全图","_Raw_Value")
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
