## Addressing reviewer comments
## Computing confusion matrix for each method evaluated in cross-validation
##
## by Artem Sokolov

library( tidyverse )
library( ABCmonster )

main <- function()
{
    set.seed(1)
    data( MACCSbinary )
    fc <- caret::trainControl( method="cv", savePred = "all" )
    XY <- select( MACCSbinary, -Drug, -pubchem_id ) %>% filter( !is.na(Label) )

    ## Method-specific caret trainer
    mytrain <- function(m) {
        cat("Training", m, "...\n")
        if (m == "gbm") 
            cv <- caret::train(Label ~ ., data = XY, method = m, 
                trControl = fc, verbose = FALSE)
        else if (m == "nnet") 
            cv <- caret::train(Label ~ ., data = XY, method = m, 
                trControl = fc, trace = FALSE, MaxNWts = 2000)
        else cv <- caret::train(Label ~ ., data = XY, method = m, 
            trControl = fc)
        cv$pred
    }

    ## Identifies the most-frequently occurring value in a vector
    mostFreq <- function( v )
    { table(v) %>% sort( decreasing=TRUE ) %>% names %>% pluck(1) }
    
    RR <- c("knn", "gbm", "glmnet", "svmLinear", "nnet") %>% 
        rlang::set_names(c("k-NN", "GBM", "Log.Reg.", "SVM", "NNet")) %>%
        purrr::map(mytrain) %>% dplyr::bind_rows(.id = "Method")

    ## AUC values from the paper
    AUC <- data_frame( Method=c("GBM", "k-NN", "Log.Reg.", "NNet", "SVM"),
                      AUC=c(0.723, 0.685, 0.669, 0.708, 0.652) )
    
    ## Aggregate predictions across parameter grid
    R1 <- RR %>% group_by( Method, rowIndex ) %>%
        summarize( pred = mostFreq(pred), obs = unique(obs) ) %>% ungroup()

    ## Compute the confusion matrix
    R2 <- R1 %>% group_by( Method ) %>%
        summarize( TP = sum(pred=="Sensitive" & obs=="Sensitive"),
                  TN = sum(pred=="Resistant" & obs=="Resistant"),
                  FP = sum(pred=="Sensitive" & obs=="Resistant"),
                  FN = sum(pred=="Resistant" & obs=="Sensitive") ) %>%
        mutate( Precision = TP / (TP+FP), Recall = TP / (TP+FN) ) %>%
        inner_join( AUC ) %>% arrange( desc(AUC) )

    write_csv( R2, "Table-Performance.csv" )
}

