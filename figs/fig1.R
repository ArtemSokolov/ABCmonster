## Figure 1 for the manuscript
##
## by Artem Sokolov

library( tidyverse )
library( cba )
library( pheatmap )
library( gridExtra )

## Location of the binary MACCS features
fnMACCS <- "../data/MACCSbinary.RData"

## Given a matrix, returns hierarchical clustering with optimal leaf re-ordering
## aux.D is an auxiliary distance matrix to be added on top of the data-driven distance matrix
orderRows <- function( X, aux.D = NULL )
{
    ## Compute pairwise distances and apply standard hierarchical clustering
    d <- dist( X )
    if( !is.null(aux.D) )
        d <- as.dist( as.matrix(d) + aux.D )
    hc <- hclust(d)

    ## Apply the optimal re-ordering algorithm
    co <- order.optimal( d, hc$merge )
    hc$merge <- co$merge
    hc$order <- co$order
    hc
}

## Generates a cross-cluster distance matrix
## All pairs of samples in the sample cluster receive a distance of 0
## All pairs of samples in different clusters receive a distance of d
## Clusters are based on pam clustering
pamDistMatrix <- function( X, n, d )
{
    ## Cluster the rows
    pp <- cluster::pam( X, k = n )

    ## Compose a distance matrix
    outer( pp$clustering, pp$clustering, "!=" ) * d
}

## Splits dendrogram into n partitions
## Taken from pheatmap package
find_gaps <- function(tree, cutree_n)
{
    v <- cutree(tree, cutree_n)[tree$order]
    which((v[-1] - v[-length(v)]) != 0)
}

## Computes coordinates of the heatmap using known gaps information
find_coordinates <- function(n, gaps, m = 1:n)
{
    if(length(gaps) == 0){
        return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
    }

    if(max(gaps) > n){
        stop("Gaps do not match with matrix size")
    }
    
    size <- (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))

    gaps2 <- apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum)
    coord <- m * size + (gaps2 * unit("4", "bigpts"))

    return(list(coord = coord, size = size))
}

## Given a clustering of rows, computes coordinates for the associated cluster labels
clusCoord <- function( tr, n )
{
    gaps <- find_gaps( tr, n )
    coord <- find_coordinates( length(tr$order), gaps )$coord
    j1 <- c( 1, gaps+1 )
    j2 <- c( gaps, length(tr$order) )
    0.5*(coord[j1] + coord[j2])
}

## Defines a palette for the sensitive / resistant distinction
ABCpal <- function()
{ c("Resistant" = "tomato", "Sensitive" = "steelblue") }

## A subfigure with a single heatmap
## XY - dataset to plot
## pfx - prefix of cluster annotations
## n - total number of clusters
fig1sub <- function( XY, pfx, nClusters )
{
    ## Compose the heatmap matrix
    HH <- select( XY, -Label, -Drug) %>% as.data.frame %>%
        column_to_rownames( "pubchem_id" ) %>% as.matrix %>% t

    ## Compose the column annotations
    YY <- select( XY, pubchem_id, Label ) %>% as.data.frame %>%
        column_to_rownames( "pubchem_id" )

    ## Cluster rows and columns
    DD <- pamDistMatrix( HH, nClusters, 5 )
    hcr <- orderRows( HH, DD )
    hcc <- orderRows( t(HH) )

    ph <- pheatmap( HH, cluster_cols = hcc, cluster_rows = hcr,
                   cutree_rows = nClusters, color = c("gray90","gray30"),
                   annotation_col = YY, annotation_colors = list( Label = ABCpal() ),
                   annotation_names_col = FALSE, annotation_legend = TRUE,
                   show_rownames = FALSE, show_colnames = FALSE, silent=TRUE )
    gg <- ph$gtable

    ## Make the legend bolder
    jAnn <- which(gg$layout$name == "annotation_legend")

    ## Replicate the annotation legend to act as a binary legend for the matrix entries
    g1 <- gg$grobs[[jAnn]]
    g1$children[[1]]$label <- "Substructure"
    g1$children[[2]]$gp$fill <- c( "gray90", "gray30" )
    g1$children[[3]]$label <- c( "Absent", "Present" )
    gg$grobs[[jAnn]] <- arrangeGrob( gg$grobs[[jAnn]], g1, heights=c(1,8) )
    
    ## Replace legend with cluster annotations
    jleg <- which( gg$layout$name == "legend" )
    gg$layout$b[jleg] <- 4
    gg$grobs[[jleg]] <- textGrob( str_c( pfx, 1:nClusters ), hjust = 1,
                                  gp = gpar( fontface ="bold", fontsize = 11),
                                  y = unit(1, "npc") - clusCoord( ph$tree_row, nClusters ) )
    gg
}

## Final version of Figure 1, panel A for the manuscript
fig1A <- function()
{
    ## Load MACCS binary data and reduce to training samples
    load( fnMACCS )
    
    ## Create a heatmap over the Substructure space only
    hh <- MACCSbinary %>% filter( !is.na(Label) ) %>% fig1sub( "S", 8 )
    hh$grobs[[3]]$children[[1]]$gp$col <- NA

    ## Increase the legend width by 20%
    w <- hh$widths[6]		## Current width
    hh$widths[6] <- hh$widths[6] + 0.20 * w
    hh$widths[3] <- hh$widths[3] - 0.20 * w

##    ggsave( "Fig1A.pdf", hh, width = 6, height = 8 )
    hh
}

## "Boldifies" a text element
etxt <- function(s, ...) element_text( face="bold", size=s, ... )

## Principal component analysis plot
fig1B <- function()
{
    ## Load raw data
    load( fnMACCS )
    X <- MACCSbinary %>% filter( !is.na(Label) )

    ## Compute PCA
    PCA <- X %>% select( -Label, -Drug, -pubchem_id ) %>% prcomp()
    P <- broom::augment( PCA, X ) %>% select( Label, .fittedPC1, .fittedPC2 )

    ## Compute the amount of variance explained by each component
    pcvar <- PCA$sdev / sum(PCA$sdev) * 100

    ## Plot the projection onto the first two principal components
    gg <- ggplot( P, aes( x = .fittedPC1, y = .fittedPC2, color=Label ) ) +
        geom_point() + theme_bw() + ggtitle( "Principal Component Analysis" ) +
        xlab( str_c( "PC1: ", round(pcvar[1],1), "% variance" ) ) +
        ylab( str_c( "PC2: ", round(pcvar[2],1), "% variance" ) ) +
        scale_color_manual( values=ABCpal(), guide=FALSE ) +
        theme( axis.title = etxt(12), axis.text = etxt(10),
##              legend.title = element_blank(), legend.text = etxt(12),
##              legend.position=c(.18, .1), legend.background = element_rect( color="black" ),
              plot.title = element_text( face="bold" ) )

    ##ggsave( "Fig1B.pdf", gg, width=4, height=4 )
    gg
}

## Multi-dimensional scaling plot
fig1C <- function()
{
    ## Some drugs have identical binary MACCS feature profiles
    ## Remove one drug from each duplicate pair, if their labels match
    ## Remove both duplicate drugs, if their labels mismatch
    ban1 <- c( 31703, 39765, 712318, 4740, 2771 )
    ban2 <- c( c(9568614, 4594), c(33036,2725) )
    
    ## Load raw data
    load( fnMACCS )
    X <- MACCSbinary %>% filter( !is.na(Label) ) %>%
        filter( !(pubchem_id %in% ban1) ) %>%
        filter( !(pubchem_id %in% ban2) )

    ## Compute MDS projections
    MDS <- X %>% select( -Label, -Drug, -pubchem_id ) %>% dist %>% MASS::isoMDS(k=2)
    M <- data_frame( Label = X$Label, MDS1 = MDS$point[,1], MDS2 = MDS$point[,2] )

    ## Plot the projections
    gg <- ggplot( M, aes( x = MDS1, y = MDS2, color = Label ) ) +
        geom_point() + theme_bw() + ggtitle( "Multi-Dimensional Scaling" ) +
        scale_color_manual( values=ABCpal(), guide=FALSE ) +
        scale_y_continuous( breaks=seq(-4,4,2) ) +
        scale_x_continuous( breaks=seq(-4,4,2) ) +
        theme( axis.title = etxt(12), axis.text = etxt(10),
              plot.title = element_text( face="bold" ) )

    ##ggsave( "Fig1C.pdf", gg, width=4, height=4 )
    gg
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
    X <- read_csv( fnMACCS ) %>% filter( !is.na(Label) ) %>%
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
        scale_color_manual( values=ABCpal(), guide=FALSE ) +
        theme( axis.title = etxt(14), axis.text = etxt(12),
              plot.title = element_text( face="bold" ) )
}


## Putting all panels together
fig1 <- function()
{
    ## Compose individual panels
    f1A <- fig1A()
    f1B <- fig1B()
    f1C <- fig1C()

    ## Bring everything together into a common figure
    f1BC <- cowplot::plot_grid( f1B, f1C, labels=c("B","C"), ncol=1 )
    f1 <- cowplot::plot_grid( f1A, f1BC, labels=c("A",""), ncol=2, rel_widths=c(2,1) )
    ggsave( "Fig1.pdf", f1, width=10, height=7 )
}
