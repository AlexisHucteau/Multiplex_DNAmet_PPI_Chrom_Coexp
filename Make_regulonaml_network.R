library(dplyr)
library(dorothea)
library(aracne.networks)
library(GeneAnswers)
library(AnnotationDbi)
library(org.Hs.eg.db)

data(regulonlaml)
viper_regulons2dorothea <- function(r) {
  res <- r %>%
    purrr::map_df(
      .f = function(i) {
        tf_target <- i$tfmode %>%
          tibble::enframe(name = "target", value = "mor") %>%
          mutate(likelihood = i$likelihood)
      },
      .id = "tf"
    )
  return(res)
}
Make_regulonaml <- function(){
  regulonaml <- viper_regulons2dorothea(regulonlaml)
  regulonaml_SYMBOL <- data.frame("source" = GeneAnswers::getSymbols(regulonaml$tf, data = "org.Hs.eg"), 
                                  "target" = GeneAnswers::getSymbols(regulonaml$target, data = "org.Hs.eg"),
                                  "mor" = regulonaml$mor,
                                  "likelihood" = regulonaml$likelihood)
    
}
