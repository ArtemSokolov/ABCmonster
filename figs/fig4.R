## Figure 4 for the manuscript
##
## by Artem Sokolov

library( tidyverse )

## Short-hand for boldface element_text() of desired size
etxt <- function( s, ... ) { element_text( size=s, face="bold", ... ) }

## Panel A
fig4a <- function()
{
    ## Load validation data
    load( "../data/ABCvaldata.RData" )

    ## Define order of labels for the "Dose" axis
    doseLbl <- c( "0", "Top/4", "Top/2", "Top" )
    ss <- c(0,0.5,0,0.5,0) %>% unit( "lines" )

    ## Prepare the data for plotting
    V <- ABCvaldata %>% gather( Strain, Growth, ABC16, Parental ) %>%
        mutate( Dose = match(Dose, doseLbl), DrugRep = str_c(Abbrev, Replicate) )
    
    gg <- ggplot( V, aes( x=Dose, y=Growth, fill=Strain ) ) + theme_bw() +
        geom_area( alpha=0.3, color="black", position="identity" ) +
        scale_fill_manual( values=c( "ABC16"="green", "Parental"="gray50" ), name="Strain:  " ) +
        scale_x_continuous( breaks = 1:4, labels = doseLbl ) +
        scale_y_continuous( breaks = c(0,0.5,1) ) +
        facet_wrap( ~DrugRep, nrow=8, ncol=6 ) +
        theme( axis.text.y = etxt(11), legend.text = etxt(11),
              axis.text.x = etxt( 11, angle=90, vjust=0.5, hjust=1 ),
              axis.title = etxt(12), legend.title = etxt(12), legend.position = "bottom",
              panel.spacing.x = ss, strip.text.x = etxt(11, margin = margin(0.1,0,0.1,0,"cm")) )

    ## Fix overlap in the x axis
    gt <- ggplotGrob(gg)
    for( j in 98:103 )
        gt$grobs[[j]]$children[[2]]$grobs[[2]]$children[[1]]$vjust <- c(1,0.5,0.5,0)

    gt
}

## Panel B
fig4b <- function()
{
    ## Load validation data
    load( "../data/ABCvaldata.RData" )
    
    ## Define the color palette and other graphical elements
    pal <- c( "Parental"="tomato", "Equal"="black", "ABC16"="steelblue" )

    ## Compute the overall response as the sum across dose-specific values
    S <- ABCvaldata %>% group_by( Drug, Abbrev, Replicate ) %>%
        summarize_at( vars(ABC16, Parental), sum ) %>%
        mutate( Value = -log2(ABC16 / Parental) ) %>% ungroup %>%
        select( -ABC16, -Parental ) %>% spread( Replicate, Value ) %>%
        mutate( Sensitivity = cut(`1`, c(-Inf, -0.1,0.1, Inf), labels=names(pal)) )

    ## Graphical elements
    rng <- select( S, `1`, `2` ) %>% range

    ## Plot the two replicates against each other
    ggplot( S, aes( x = `1`, y = `2`, color = Sensitivity ) ) +
        geom_point() + theme_bw() + scale_color_manual( values = pal ) +
        xlab( "Replicate 1 -log2(ABC-16 / parental)" ) +
        ylab( "Replicate 2 -log2(ABC-16 / parental)" ) +
        xlim( rng ) + ylim( rng ) +
        geom_abline( slope = 1, linetype = "dashed", color = "gray40" ) +
        ggrepel::geom_text_repel( aes(label = Abbrev), fontface="bold", show.legend=FALSE ) +
        theme( axis.text = etxt(11), axis.title = etxt(12),
              legend.title = etxt(12), legend.text = etxt(11),
              legend.position = c(0.95,0.05), legend.justification = c(1,0),
              legend.background = element_rect( color="black", fill="white" ) )
}
