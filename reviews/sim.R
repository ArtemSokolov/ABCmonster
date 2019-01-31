## Addressing reviewers comments
## Identifies features that differ between drugs with high structural
##  similarity (Fig 4) and maps them to importance scores from the
##  predictor (Fig3C)
##
## by Artem Sokolov

library( ABCmonster )
library( tidyverse )

data(MACCSbinary)

main <- function()
{
    ## Feature importance scores
    set.seed(100)
    X1 <- MACCSbinary %>% filter( !is.na(Label) ) %>% select( -Drug, -pubchem_id )
    cv <- caret::train( Label ~ ., data=X1, method="gbm", verbose=FALSE,
                       trControl = caret::trainControl(method="cv") )
    FI <- summary( cv$finalModel, plot = FALSE ) %>% arrange( desc(rel.inf) ) %>%
        rename( Key = var ) %>% mutate_at( "Key", as.character ) %>%
        mutate_at( "Key", str_sub, 2, -2 )
    
    X <- as_tibble( MACCSbinary ) %>% mutate( Label = NA )

    ## Exemplar pairs with high structural similarity
    p1 <- c( 24360, 60700 )      # CAM    <--- MACCS(-93)
    p2 <- c( 3007984, 3792 )     # BEA    <--- MACCS(151)

    ## Compute prediction scores for everything (training and test data)
    mm <- ABCtrain( MACCSbinary )
    Y <- ABCpredict( mm, X )

    ## Identify most important features among those that differ in drug pairs
    ##   with high structural similarity
    Y %>% filter( pubchem_id %in% p1 ) %>% select( -Drug ) %>%
        gather( Key, Value, -ABCpred, -pubchem_id ) %>%
        group_by( Key ) %>% summarize( nn = length(unique(Value)) ) %>%
        filter( nn > 1 ) %>% left_join( FI ) %>% arrange( desc(rel.inf) )
}
