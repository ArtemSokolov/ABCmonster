## Figure 4 for the manuscript
##
## by Artem Sokolov

library( tidyverse )
library( ggnetwork )

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

## Computes Tanimoto (Jaccard) similarity for two MACCS keys profiles
tanimoto <- function( v1, v2 )
{
    ## List-ize the inputs
    if( is.list(v1) ) return( map( v1, ~tanimoto(.x, v2) ) %>% do.call( cbind, . ) )
    if( is.list(v2) ) return( map_dbl( v2, ~tanimoto(v1, .x) ) )
    
    ## Identify MACCS keys that are present in each profile
    m1 <- v1 %>% keep( . == 1 ) %>% names
    m2 <- v2 %>% keep( . == 1 ) %>% names

    ## Compute Jaccard index
    length( intersect(m1, m2) ) / length( union(m1, m2) )
}

## Panel C
fig4c <- function()
{
    ## Load data in binary MACCS representation
    load( "../data/MACCSbinary.RData" )
    XX <- MACCSbinary %>% as_data_frame %>% nest( -Label, -Drug, -pubchem_id, .key="MACCS" ) %>%
        mutate_at( "MACCS", map, ~deframe(gather(.x)) )
    
    ## Separate into training and test
    Xtr <- filter( XX, !is.na(Label) ) %>% mutate_at( "pubchem_id", as.character )
    Xte <- filter( XX, is.na(Label) ) %>% mutate( Abbrev = str_to_upper(str_sub(Drug, 1, 3)) )

    ## Compute pair-wise similarity between the two sets
    Mtr <- with( Xtr, set_names(MACCS, pubchem_id) )
    Mte <- with( Xte, set_names(MACCS, Abbrev ) )
    S <- tanimoto( Mte, Mtr )

    ## Identify the closest drug pairs and re-annotate the match
    R <- S %>% as.data.frame %>% rownames_to_column( "pubchem_id" ) %>%
        gather( Abbrev, Tanimoto, -pubchem_id ) %>% group_by( Abbrev ) %>%
        top_n( 1, Tanimoto ) %>% ungroup %>%
        inner_join( select( Xtr, pubchem_id, Label ) ) %>%
        mutate( Sens = case_when( Abbrev %in% c("RAP", "TUN", "VAL") ~ "Parental",
                                 Abbrev %in% c("BEN", "BRO", "IMA") ~ "Equal",
                                 TRUE ~ "ABC16" ) ) %>%
        arrange( desc(Tanimoto) )

    ## Common graph elements
    sensLbls <- c( "ABC16 (Sensitive)", "Parental (Resistant)" )
    fy <- function() {rep( (8:1)/8, 3)}
    fx <- function(x) {rep(0:2,each=8)+x}

    ## Compose the bipartite graph matrix
    N1 <- R %>% select( vertex.names = Abbrev, Sensitivity = Sens ) %>%
        mutate( x=fx(0.2), xend=fx(0.2), y=fy(), yend=fy()+0.001, na.x=FALSE, na.y=NA ) %>%
        mutate_at( "Sensitivity", recode, !!!set_names(sensLbls, c("ABC16", "Parental")) )
    N2 <- R %>% select( vertex.names = pubchem_id, Sensitivity = Label ) %>%
        mutate( x=fx(0.8), xend=fx(0.8), y=fy(), yend=fy()+0.001, na.x=FALSE, na.y=NA ) %>%
        mutate_at( "Sensitivity", recode, !!!set_names(sensLbls, c("Sensitive", "Resistant")) )
    E <- R %>% select( vertex.names = Abbrev, Sensitivity = Sens, Tanimoto ) %>%
        mutate( x=fx(0.2), xend=fx(0.8), y=fy(), yend=fy()+0.001, na.x=FALSE, na.y=FALSE ) %>%
        mutate_at( "Sensitivity", recode, !!!set_names(sensLbls, c("ABC16", "Parental")) ) %>%
        mutate_at( "Tanimoto", ~as.character(round(.x,2)) )
    BG <- bind_rows( N1, N2, E )

    ## Plotting elements
    pal <- c( sensLbls, "Equal" ) %>% set_names( c("steelblue", "tomato", "black"), . )
    
    ## Plot the bipartite graph
    ggplot( BG, aes(x = x, y = y, xend = xend, yend = yend) ) + theme_blank() +
        geom_edges( color="gray40", lwd=1.25 ) +
        geom_nodes( aes(color=Sensitivity), size=8 ) +
        geom_nodetext( aes(label=vertex.names), nudge_y = -0.04, fontface="bold" ) +
        geom_edgetext( aes(label=Tanimoto), color="black", fontface="bold" ) +
        scale_color_manual( values = pal ) +
        theme( legend.position = "bottom", legend.title=etxt(12), legend.text=etxt(11) )
}

