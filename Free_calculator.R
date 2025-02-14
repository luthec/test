library(dplyr)
library(openxlsx)
library(readxl)
library(purrr)

cbat <- function(TT,SHBG,ALB = 43){
    # Convert TT from ng/mL to nmol/L (1 ng/mL = 3.467nmol/L)
    TT <- TT * 3.467
    Kalb <- 3.6*10^4
    Kshbg <- 10^9
    N <- 1 + Kalb*ALB/69000
    a <- N*Kshbg
    b <- N + Kshbg*(SHBG - TT)/10^9
    c <- -TT/10^9
    FT <- (-b + sqrt(b^2 - 4*a*c))/(2*a)*10^9
    # Convert results back to ng/mL (1 nmol/L = 28.84 ng/dL and 1 dL = 100 mL)
    FT <- FT * 0.2884
    cbat <- N*FT
    # return(list(free.T = FT, cbat = cbat))
    tibble::tibble(FT,cbat)
}

mat=read_excel(dir(),col_types ="text",sheet=3) %>% 
         mutate_at(2:4, as.numeric) %>% suppressWarnings() %>%
         rename_at(2,~"Alb") %>% 
         rename_at(3,~"SHBG") %>% 
         rename_at(4,~"TT")


pmap_dfr(list(mat$TT,
              mat$SHBG, 
              mat$Alb),
              cbat) %>%
    bind_cols(mat, .) %>%
    rename_at(6,~"Bioavailable Testosterone") %>%
    write.xlsx("Free & Bioavailable Testosterone_results.xlsx",  colNames = TRUE)





    


