## Figure 5 for the manuscript
##
## by Artem Sokolov

library( tidyverse )

## Short-hand for boldface element_text() of desired size
etxt <- function( s, ... ) { element_text( size=s, face="bold", ... ) }

## Define the color palette
g_pal <- c( "Parental"="tomato", "Equal"="black", "ABC16"="steelblue" )
    
## Panel A
fig5a <- function( M )
{
    ## Compute summary statistics
    SS <- as_data_frame(M) %>%
        summarize( ct = list(cor.test( Response, ABCpred, method="spearman" )) ) %>%
        mutate( rho = map_dbl(ct, "estimate"), pval = map_dbl(ct, "p.value") ) %>%
        mutate( Txt1 = glue::glue( "Spearman rho = {round(rho,3)}" ),
               Txt2 = glue::glue( "p-value = {round(pval,3)}" ) ) %>%
        mutate( Txt = str_c(Txt1, "\n", Txt2, " \n") )
    
    ## Plot the correspondence between true values and estimates
    ## Include a loess fit to the data points
    ggplot( M, aes(x=ABCpred, y=Response, color=Sensitivity) ) + theme_bw() +
        geom_point() + xlab("Predicted P(ABC-16 is more sensitive)") +
        ylab( "Measured response,\n averaged across replicates" ) +
#        coord_fixed() +
        ggrepel::geom_text_repel( aes(label = Abbrev), fontface="bold", size=4 ) +
        geom_text( data = SS, aes( label = Txt ), x = Inf, y = -Inf, vjust = 0,
                  hjust = 1.1, fontface = "bold", color="black" ) +
        scale_color_manual( values=g_pal, guide=FALSE ) + scale_x_log10() +
        theme( axis.title = etxt(12), axis.text = etxt(11), aspect.ratio=1 )
}

auc <- function( pred, lbls )
{
    np <- sum( lbls == "ABC16" )
    nn <- sum( lbls == "Parental" )
    s0 <- sum( rank(pred)[lbls == "ABC16"] )
    (s0 - np*(np+1) / 2) / (np*nn)
}

## Panel B
fig5b <- function( RC )
{
    ## An extra point at (0,0)
    RC2 <- data_frame( tpr = 0, fpr = 0 ) %>% bind_rows( RC )

    ## Put together the AUC label
    Lbl <- RC %>% summarize( AUC=glue::glue("AUC = {auc(ABCpred, Sens2) %>% round(2)}\n\n") )
    
    ## Plot the ROC curve
    ggplot( RC2, aes(x = fpr, y = tpr) ) + theme_bw() +
        ## geom_smooth(color="gray70", se=FALSE) +
        geom_line( color="gray70", size=1.2 ) +
        xlim( c(0,1) ) + ylim( c(0,1) ) + coord_fixed() +
        geom_point( aes(color=Sens), size=3 ) +
        scale_color_manual( values=g_pal, guide=FALSE ) +
        geom_abline( slope=1, linetype="dashed" ) +
        xlab( "False Positive Rate" ) +
        ylab( "True Positive Rate" ) +
        geom_text( data = Lbl, aes(label = AUC), x = Inf, y = -Inf, vjust = 0,
                  hjust = 1.1, fontface = "bold", color="black" ) +
        theme( axis.text = etxt(11), axis.title = etxt(12) )
}

## Panel C
fig5c <- function( RC )
{
    ## An extra point at (0,1)
    RC2 <- data_frame( tpr = 0, prec = 1 ) %>% bind_rows( RC )
    
    ## Plot the Precision-Recall curve
    ggplot( RC2, aes(x = tpr, y = prec) ) + theme_bw() +
        xlim( c(0,1) ) + ylim( c(0,1) ) + coord_fixed() +
        geom_line( color="gray70", size=1.2 ) +
        geom_point( aes(color=Sens), size=3 ) +
        scale_color_manual( values=g_pal, breaks=c("Parental", "Equal", "ABC16"),
                           name="Sensitivity" ) +
        ylab( "Precision" ) + xlab( "Recall" ) +
        theme( axis.text = etxt(11), axis.title = etxt(12),
              legend.title = etxt(12), legend.text = etxt(11),
              legend.position = c(0.95,0.05), legend.justification = c(1,0),
              legend.background = element_rect( color="gray40", fill="white" ) )
}

## Composite figure
fig5 <- function()
{
    ## Compute predictions on the test data
    data( MACCSbinary, package="ABCmonster" )
    P <- ABCmonster::ABCtrain( MACCSbinary ) %>% ABCmonster::ABCpredict( MACCSbinary ) %>%
        select( ABCpred, Drug, pubchem_id )

    ## Load validation data and compute overall response
    data( ABCvaldata, package="ABCmonster" )
    S <- ABCvaldata %>% group_by( Drug, Abbrev, Replicate ) %>%
        summarize_at( vars(ABC16, Parental), sum ) %>%
        mutate( Value = -log2(ABC16 / Parental) ) %>% ungroup %>%
        select( -ABC16, -Parental ) %>% spread( Replicate, Value ) %>%
        mutate( Sensitivity = cut(`1`, c(-Inf, -0.1,0.1, Inf), labels=names(g_pal)) )

    ## Combine everything into a single data frame
    M <- inner_join( P, S ) %>% mutate( Response = map2_dbl(`1`, `2`, lift_vd(mean)) ) %>%
        select( -`1`, -`2`, -Drug )

    ## Create Panel A here
    f5a <- fig5a( M )

    ## Drop drugs with equal sensitivity and compute ROC and P-R curves for the rest
    RC <- M %>% arrange( desc(ABCpred) ) %>% rename( Sens = Sensitivity ) %>%
        mutate( Sens2 = recode(Sens, Equal="Parental") ) %>%
        mutate( tpr = cumsum(Sens2 == "ABC16") / sum(Sens2 == "ABC16"),
               fpr = cumsum(Sens2 == "Parental") / sum(Sens2 == "Parental"),
               prec = cumsum(Sens2 == "ABC16") / (1:length(Sens2)) )

    ## Plot the remaining panels
    f5b <- fig5b( RC )
    f5c <- fig5c( RC )

    ## Put everything together into a composite figure
    f5 <- cowplot::plot_grid( f5a, f5b, f5c, labels=LETTERS[1:3], nrow=1, label_size = 20 )
    ggsave( "Fig5.pdf", f5, width=14, height=4 )
}
