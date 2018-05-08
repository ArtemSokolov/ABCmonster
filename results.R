## Analysis of results
##
## by Artem Sokolov

library( tidyverse )

## An AUC estimate that doesn't require explicit construction of an ROC curve
auc <- function( scores, lbls )
{
    stopifnot( length(scores) == length(lbls) )
    jp <- which( lbls > 0 ); np <- length( jp )
    jn <- which( lbls <= 0); nn <- length( jn )
    s0 <- sum( rank(scores)[jp] )
    (s0 - np*(np+1) / 2) / (np*nn)
}

## Given a cross-validation results data frame, computes AUC for each iteration
cv2auc <- function( X )
{ X %>% group_by( Iteration ) %>% summarize_at( vars(knn:nnet), funs(auc(.,Label)) ) %>% ungroup }

## Compares cross-validation results across methods and datasets
main <- function()
{
    ## Common set of operations applied to each file
    f <- function( fnIn, lbl )
    { read_csv( fnIn, col_types = cols() ) %>% cv2auc %>%
          mutate( Data = lbl ) %>%
          rename( `k-NN` = knn, GBM = gbm, `Log.Reg.` = glmnet, SVM = svmLinear, NNet = nnet ) %>%
          gather( Method, AUC, `k-NN`:NNet ) }
    
    ## Load individual cross-validation results files and compute AUC for each iteration/method
    XX <- bind_rows( f( "results/MACCSbinary-cv.csv.gz", "MACCS(binary)" ),
                    f( "results/MACCScount-cv.csv.gz", "MACCS(count)" ),
                    f( "results/PChem-cv.csv.gz", "PhysicoChemical" ),
                    f( "results/Morgan-cv.csv.gz", "Morgan/ECFP" ) ) %>%
        group_by( Data, Method ) %>% summarize( AUC = mean(AUC) ) %>%
        mutate( AUC = pmin( AUC, 0.73 ) )

    ## Define plot elements
    cl1 <- "#203c52"
    cl2 <- "#ccdcea"
    fpal <- colorRampPalette( c(cl1,cl2) )
    sg <- guide_legend( override.aes = list( fill = fpal(5) ) )
    etxt <- function(s, ...) element_text( face="bold", size=s, ... )
    
    ## Make a summary plot
    gg <- ggplot( XX, aes( x = Data, y = Method, size = AUC, fill = AUC ) ) +
        theme_bw() + geom_point( shape = 21, color="black" ) +
        scale_radius( range = c(2,15), breaks=seq(0.61,0.73,0.03), limits=c(0.61,0.73) ) +
        scale_fill_gradient( limits=c(0.61,0.73), low=cl1, high=cl2 ) +
        theme( axis.text.y = etxt(11), axis.title = etxt(12),
              axis.text.x = etxt(11, angle=90, hjust=1, vjust=0.5), 
              panel.grid.major = element_line( linetype = "dashed", color="gray70" ),
              legend.title = etxt(14), legend.text = etxt(12) ) +
        guides( size = sg, fill=FALSE )
    ggsave( "FigX.pdf", gg, width=5.25, height=6 )
}

