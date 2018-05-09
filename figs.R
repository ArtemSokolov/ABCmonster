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

## "Boldifies" a text element
etxt <- function(s, ...) element_text( face="bold", size=s, ... )

## Computes correlation between predictions and true labels
pred.plot <- function()
{
    ## Load predictions and true labels
    Ytrue <- read_csv( "data/test-labels.csv" )
    Ypred <- read_tsv( "results/gbm-pred.tsv" )

    ## Define plot elements
    cl1 <- "tomato"
    cl2 <- "wheat"
    fpal <- colorRampPalette( c(cl1,cl2) )
    sg <- guide_legend( override.aes = list(fill=fpal(4)), title="Correlation",
                       keyheight=c(0.2,0.3,0.4,0.6), keywidth=0.75, default.unit="inch")
    
    ## Compute Spearman correlations between true labels and predictions
    ##   for each method
    X <- inner_join( Ytrue, Ypred ) %>%
        summarize_at( vars(MACCSbin:Morgan), funs(cor(.,GMsens,method="sp")) ) %>%
        rename( `MACCS(binary)` = MACCSbin, `MACCS(count)` = MACCScnt,
               Physicochemical = PChem, `Morgan/ECFP` = Morgan ) %>%
        gather( Features, Corr. ) %>% mutate( ` ` = "GBM", Experiment = "Prospective" ) %>%
        mutate( Corr. = pmin( Corr., 0.535 ) )

    ggplot( X, aes( x = Features, y = ` `, fill = `Corr.`, size=`Corr.` ) ) +
        theme_bw() + geom_point( shape = 23, color="black" ) +
        scale_radius( range = c(2,15), breaks=seq(0.2,0.5,0.1), limits=c(0.185,0.535) ) +
        scale_fill_gradient( limits=c(0.185,0.535), low=cl1, high=cl2 ) +
        facet_wrap( ~Experiment ) +
        theme( axis.text.y = etxt(11), axis.title = etxt(12),
              axis.text.x = etxt(11, angle=90, hjust=1, vjust=0.5), 
              panel.grid.major = element_line( linetype = "dashed", color="gray70" ),
              legend.title = etxt(14), legend.text = etxt(12),
              strip.text = etxt(12), strip.background=element_rect( fill="gray95" ) ) +
        guides( fill=FALSE, size=sg )
}

## Compares cross-validation results across methods and datasets
AUC.plot <- function()
{
    ## Common set of operations applied to each file
    f <- function( fnIn, lbl )
    { read_csv( fnIn, col_types = cols() ) %>% cv2auc %>%
          mutate( Features = lbl ) %>%
          rename( `k-NN` = knn, GBM = gbm, `Log.Reg.` = glmnet, SVM = svmLinear, NNet = nnet ) %>%
          gather( Method, AUC, `k-NN`:NNet ) }
    
    ## Load individual cross-validation results files and compute AUC for each iteration/method
    XX <- bind_rows( f( "results/MACCSbinary-cv.csv.gz", "MACCS(binary)" ),
                    f( "results/MACCScount-cv.csv.gz", "MACCS(count)" ),
                    f( "results/PChem-cv.csv.gz", "PhysicoChemical" ),
                    f( "results/Morgan-cv.csv.gz", "Morgan/ECFP" ) ) %>%
        group_by( Features, Method ) %>% summarize( AUC = mean(AUC) ) %>%
        mutate( AUC = pmin( AUC, 0.73 ) ) %>%	## Clip 0.731 to 0.73 for visualization purposes
        ungroup %>% mutate( Experiment = "Cross Validation" )

    ## Define plot elements
    cl1 <- "#203c52"
    cl2 <- "#ccdcea"
    fpal <- colorRampPalette( c(cl1,cl2) )
    sg <- guide_legend( override.aes = list( fill = fpal(5) ) )
    
    ## Make a summary plot
    ggplot( XX, aes( x = Features, y = Method, size = AUC, fill = AUC ) ) +
        theme_bw() + geom_point( shape = 21, color="black" ) +
        scale_radius( range = c(2,15), breaks=seq(0.61,0.73,0.03), limits=c(0.61,0.73) ) +
        scale_fill_gradient( limits=c(0.61,0.73), low=cl1, high=cl2 ) +
        facet_wrap( ~Experiment ) +
        theme( axis.text.y = etxt(11), axis.title.y = etxt(12, hjust=0),
              axis.text.x = element_blank(), axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major = element_line( linetype = "dashed", color="gray70" ),
              legend.title = etxt(14), legend.text = etxt(12),
              strip.text = etxt(12), strip.background=element_rect( fill="gray95" ) ) +
        guides( size = sg, fill=FALSE )
##    ggsave( "FigX.pdf", gg, width=5.25, height=6 )
}

## Combines pred.plot() and AUC.plot() into the same figure
Fig.S2 <- function()
{
    ## Generate each plot and combine them into a common tableGrob
    gg1 <- AUC.plot()
    gg2 <- pred.plot()
    gt <- gridExtra::arrangeGrob( gg1, gg2, ncol=1, heights=c(1.5,1) )

    ## Adjust widths to match
    gt$grobs[[2]]$widths <- gt$grobs[[1]]$widths
    
    grid::grid.newpage()
    grid::grid.draw(gt)
    ggsave( "Fig-perf.pdf", gt, width=5.5, height=7 )
}

## Principal component analysis and multi-dimensional scaling plots
PCA.plot <- function()
{
    ## Load raw data
    X <- read_csv( "data/MACCSbinary.csv" ) %>% filter( !is.na(Label) )

    ## Compute PCA
    PCA <- X %>% select( -Label, -Drug, -pubchem_id ) %>% prcomp()
    P <- broom::augment( PCA, X ) %>% select( Label, .fittedPC1, .fittedPC2 ) %>%
        mutate( Label = factor( Label, c(1,0) ) )
    levels(P$Label) <- c( "Sensitive", "Resistant" )

    ## Compute the amount of variance explained by each component
    pcvar <- PCA$sdev / sum(PCA$sdev) * 100

    ## Plot the projection onto the first two principal components
    ggplot( P, aes( x = .fittedPC1, y = .fittedPC2, color=Label ) ) +
        geom_point() + theme_bw() + ggtitle( "Principal Component Analysis" ) +
        xlab( str_c( "PC1: ", round(pcvar[1],1), "% variance explained" ) ) +
        ylab( str_c( "PC2: ", round(pcvar[2],1), "% variance explained" ) ) +
        scale_color_manual( values=c("Sensitive"="tomato", "Resistant"="steelblue") ) +
        theme( axis.title = etxt(14), axis.text = etxt(12),
              legend.title = element_blank(), legend.text = etxt(12),
              legend.position=c(.18, .1), legend.background = element_rect( color="black" ),
              plot.title = element_text( face="bold" ) )
}

## Multi-dimensional scaling plot
MDS.plot <- function()
{
    ## Some drugs have identical binary MACCS feature profiles
    ## Remove one drug from each duplicate pair, if their labels match
    ## Remove both duplicate drugs, if their labels mismatch
    ban1 <- c( 31703, 39765, 712318, 4740, 2771 )
    ban2 <- c( c(9568614, 4594), c(33036,2725) )
    
    ## Load raw data
    X <- read_csv( "data/MACCSbinary.csv" ) %>% filter( !is.na(Label) ) %>%
        filter( !(pubchem_id %in% ban1) ) %>%
        filter( !(pubchem_id %in% ban2) ) %>%
        mutate( Label = factor( Label, c(1,0) ) )
    levels(X$Label) <- c( "Sensitive", "Resistant" )

    ## Compute MDS projections
    MDS <- X %>% select( -Label, -Drug, -pubchem_id ) %>% dist %>% MASS::isoMDS(k=2)
    M <- data_frame( Label = X$Label, MDS1 = MDS$point[,1], MDS2 = MDS$point[,2] )

    ## Plot the projections
    ggplot( M, aes( x = MDS1, y = MDS2, color = Label ) ) +
        geom_point() + theme_bw() + ggtitle( "Multi-Dimensional Scaling" ) +
        scale_color_manual( values=c("Sensitive"="tomato", "Resistant"="steelblue"), guide=FALSE ) +
        theme( axis.title = etxt(14), axis.text = etxt(12),
              plot.title = element_text( face="bold" ) )
}

## t-Stochastic Neighbor Embedding plot
tSNE.plot <- function()
{
    ## Some drugs have identical binary MACCS feature profiles
    ## Remove one drug from each duplicate pair, if their labels match
    ## Remove both duplicate drugs, if their labels mismatch
    ban1 <- c( 31703, 39765, 712318, 4740, 2771 )
    ban2 <- c( c(9568614, 4594), c(33036,2725) )
    
    ## Load raw data
    X <- read_csv( "data/MACCSbinary.csv" ) %>% filter( !is.na(Label) ) %>%
        filter( !(pubchem_id %in% ban1) ) %>%
        filter( !(pubchem_id %in% ban2) ) %>%
        mutate( Label = factor( Label, c(1,0) ) )
    levels(X$Label) <- c( "Sensitive", "Resistant" )

    ## Compute t-SNE projections
    X1 <- X %>% select( -Label, -Drug, -pubchem_id )
    TS <- Rtsne::Rtsne(X1, perplexity=30)
    P <- data_frame( Label = X$Label, `t-SNE1` = TS$Y[,1], `t-SNE2` = TS$Y[,2] )

    ## Plot the projections
    ggplot( P, aes( x = `t-SNE1`, y = `t-SNE2`, color = Label ) ) +
        geom_point() + theme_bw() + ggtitle( "t-Stochastic Neighbor Embedding" ) +
        scale_color_manual( values=c("Sensitive"="tomato", "Resistant"="steelblue"), guide=FALSE ) +
        theme( axis.title = etxt(14), axis.text = etxt(12),
              plot.title = element_text( face="bold" ) )
}

## Puts PCA, MDS and tSNE plots into the same figure
proj.plot <- function()
{
    ## Compute individual plots
    P1 <- PCA.plot()
    P2 <- MDS.plot()
    P3 <- tSNE.plot()

    ## Merge everything into a single figure
    gt <- gridExtra::arrangeGrob( P1, P2, P3, nrow=1 )
    ggsave( "Fig-proj.pdf", gt, width=14, height=5 )
}
