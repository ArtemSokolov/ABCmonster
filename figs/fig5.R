## Figure 5 for the manuscript
##
## by Artem Sokolov

library( tidyverse )

## Short-hand for boldface element_text() of desired size
etxt <- function( s, ... ) { element_text( size=s, face="bold", ... ) }

## Panel A
fig5a <- function()
{
    ## Define the color palette
    pal <- c( "Parental"="tomato", "Equal"="black", "ABC16"="steelblue" )
    
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
        mutate( Sensitivity = cut(`1`, c(-Inf, -0.1,0.1, Inf), labels=names(pal)) )

    ## Combine everything into a single data frame
    M <- inner_join( P, S ) %>% mutate( Response = map2_dbl(`1`, `2`, lift_vd(mean)) )

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
        ylab( "Measured response, averaged across replicates" ) +
        ggrepel::geom_text_repel( aes(label = Abbrev), fontface="bold", size=4 ) +
        geom_text( data = SS, aes( label = Txt ), x = Inf, y = -Inf, vjust = 0,
                  hjust = 1.1, fontface = "bold", color="black" ) +
        scale_color_manual( values=pal, guide=FALSE ) + scale_x_log10() +
        theme( axis.title = etxt(12), axis.text = etxt(11) )
}
