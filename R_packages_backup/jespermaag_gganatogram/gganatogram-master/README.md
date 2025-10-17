<!-- README.md is generated from README.Rmd. Please edit that file -->

## gganatogram

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/jespermaag/gganatogram?branch=master&svg=true)](https://ci.appveyor.com/project/jespermaag/gganatogram)
[![Travis build
status](https://travis-ci.com/jespermaag/gganatogram.svg?branch=master)](https://travis-ci.com/jespermaag/gganatogram)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1434233.svg)](https://doi.org/10.5281/zenodo.1434233)

Create anatogram images for different organisms. <br/> This package uses
the tissue coordinates from the figure in Expression Atlas.
<https://www.ebi.ac.uk/gxa/home> <br/>
<https://github.com/ebi-gene-expression-group/anatomogram> <br/>

![](figure/AllSpeciesCellPlotValueTop-1.svg)

![](figure/CellPlotValueTop-1.svg)

## Citation

#### Maag JLV. gganatogram: An R package for modular visualisation of anatograms and tissues based on ggplot2 \[version 1; referees: 1 approved\]. F1000Research 2018, 7:1576 (doi: 10.12688/f1000research.16409.1)

<https://f1000research.com/articles/7-1576/v1>

``` r
citation("gganatogram")
#> To cite package 'gganatogram' in publications use:
#> 
#>   Maag J (2018). "gganatogram: An R package for modular visualisation
#>   of anatograms and tissues based on ggplot2." _f1000research_. Version
#>   1: Awaiting peer review,
#>   <https://f1000research.com/articles/7-1576/v1>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {gganatogram:  An R package for modular visualisation of anatograms and tissues based on ggplot2},
#>     author = {Jesper Maag},
#>     journal = {f1000research},
#>     year = {2018},
#>     note = {Version 1: Awaiting peer review},
#>     url = {https://f1000research.com/articles/7-1576/v1},
#>   }
```

If you use the tissue plots from gganatogram please cite Expression
Atlas as well. <br/> [Petryszak et
al. 2015](https://academic.oup.com/nar/article/44/D1/D746/2502589) <br/>
If you use the main cell figure, please cite The Protein Atlas. <br/>
[Thul PJ et
al. 2017](http://science.sciencemag.org/content/356/6340/eaal3321) <br/>

More plot examples can be found at
<https://jespermaag.github.io/blog/2018/gganatogram/>

## Install

Install from github using devtools.

``` r
## install from Github
devtools::install_github("jespermaag/gganatogram")
```

## shiny

I have now included a shiny app for gganatogram. <br/> An online version
can be found at shinapps.io. <br/>
<https://jespermaag.shinyapps.io/gganatogram/> <br/> Unfortunately,
there is a limit of 25h per month of app activity, so if you know
R/Rstudio, please run it locally. <br/> To run it locally, use the
following command.

``` r
library(shiny)
runGitHub( "gganatogram", "jespermaag",  subdir = "shiny") 
```

## Usage

This package requires `ggplot2` and `ggpolypath` which loads when
loading the package

``` r

library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)
```

``` r
hgMale <- gganatogram(data=hgMale_key, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") + theme_void()
hgFemale <- gganatogram(data=hgFemale_key, fillOutline='#a6bddb', organism='human', sex='female', fill="colour") + theme_void()
mmMale <- gganatogram(data=mmMale_key, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour") + theme_void()
mmFemale <- gganatogram(data=mmFemale_key, outline = T, fillOutline='#a6bddb', organism='mouse', sex='female', fill="colour")  +theme_void()  

grid.arrange(hgMale, hgFemale, mmMale, mmFemale, ncol=4)
```

![](figure/AllSpeciesPlot-1.svg)

``` r


hgMale <- gganatogram(data=hgMale_key, fillOutline='#440154FF', organism='human', sex='male', fill="value") + theme_void() +  scale_fill_viridis()
hgFemale <- gganatogram(data=hgFemale_key, fillOutline='#440154FF', organism='human', sex='female', fill="value") + theme_void() +  scale_fill_viridis()
mmMale <- gganatogram(data=mmMale_key, fillOutline='#440154FF', organism='mouse', sex='male', fill="value") + theme_void() +  scale_fill_viridis()
mmFemale <- gganatogram(data=mmFemale_key, outline = T, fillOutline='#440154FF', organism='mouse', sex='female', fill="value")  +theme_void()   +  scale_fill_viridis()

grid.arrange(hgMale, hgFemale, mmMale, mmFemale, ncol=2)
```

![](figure/AllSpeciesPlotValue-1.svg)

In order to use the function gganatogram, you need to have a data frame
with organ, colour, and value if you want to.

``` r
organPlot <- data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"), 
 type = c("circulation", "circulation",  "nervous system", "nervous system", "digestion", "digestion", "digestion"), 
 colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"), 
 value = c(10, 5, 1, 8, 2, 5, 5), 
 stringsAsFactors=F)

 head(organPlot)
#>       organ           type colour value
#> 1     heart    circulation    red    10
#> 2 leukocyte    circulation    red     5
#> 3     nerve nervous system purple     1
#> 4     brain nervous system purple     8
#> 5     liver      digestion orange     2
#> 6   stomach      digestion orange     5
```

Using the function gganatogram with the filling the organs based on
colour.

``` r
gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="colour")
```

![](figure/organPlot-1.svg)

Of course, we can use the ggplot themes and functions to adjust the
plots

``` r
gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") + 
theme_void()
```

![](figure/organPlotvoid-1.svg)

We can also plot all tissues available using hgMale_key

``` r
hgMale_key$organ
#>  [1] "thyroid_gland"             "bone_marrow"              
#>  [3] "frontal_cortex"            "prefrontal_cortex"        
#>  [5] "pituitary_gland"           "aorta"                    
#>  [7] "gastroesophageal_junction" "left_ventricle"           
#>  [9] "caecum"                    "ileum"                    
#> [11] "rectum"                    "nose"                     
#> [13] "breast"                    "tongue"                   
#> [15] "left_atrium"               "pulmonary_valve"          
#> [17] "mitral_valve"              "penis"                    
#> [19] "nasal_pharynx"             "spinal_cord"              
#> [21] "throat"                    "tricuspid_valve"          
#> [23] "diaphragm"                 "liver"                    
#> [25] "stomach"                   "spleen"                   
#> [27] "duodenum"                  "gall_bladder"             
#> [29] "pancreas"                  "colon"                    
#> [31] "small_intestine"           "appendix"                 
#> [33] "smooth_muscle"             "urinary_bladder"          
#> [35] "bone"                      "cartilage"                
#> [37] "esophagus"                 "salivary_gland"           
#> [39] "parotid_gland"             "submandibular_gland"      
#> [41] "skin"                      "pleura"                   
#> [43] "brain"                     "heart"                    
#> [45] "adrenal_gland"             "lymph_node"               
#> [47] "adipose_tissue"            "skeletal_muscle"          
#> [49] "leukocyte"                 "temporal_lobe"            
#> [51] "atrial_appendage"          "coronary_artery"          
#> [53] "hippocampus"               "vas_deferens"             
#> [55] "seminal_vesicle"           "epididymis"               
#> [57] "tonsil"                    "lung"                     
#> [59] "amygdala"                  "trachea"                  
#> [61] "bronchus"                  "nerve"                    
#> [63] "cerebellum"                "cerebellar_hemisphere"    
#> [65] "kidney"                    "renal_cortex"             
#> [67] "testis"                    "prostate"
gganatogram(data=hgMale_key, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") +theme_void()
```

![](figure/organPlotAll-1.svg)

We can also skip plotting the outline of the graph

``` r
organPlot %>%
    dplyr::filter(type %in% c('circulation', 'nervous system')) %>%
gganatogram(outline=F, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") + 
theme_void()
```

![](figure/organPlotSubset-1.svg)

We can fill the tissues based on the values given to each organ

``` r
gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="value") + 
theme_void() +
scale_fill_gradient(low = "white", high = "red")
```

![](figure/organPlotValue-1.svg)

We can also use facet_wrap to compare groups.  
First create add two data frames together with different values and the
conditions in the type column

``` r
compareGroups <- rbind(data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"), 
  colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"), 
 value = c(10, 5, 1, 8, 2, 5, 5), 
 type = rep('Normal', 7), 
 stringsAsFactors=F),
 data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"), 
  colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"), 
 value = c(5, 5, 10, 8, 2, 5, 5), 
 type = rep('Cancer', 7), 
 stringsAsFactors=F))
```

``` r
gganatogram(data=compareGroups, fillOutline='#a6bddb', organism='human', sex='male', fill="value") + 
theme_void() +
facet_wrap(~type) +
scale_fill_gradient(low = "white", high = "red") 
```

![](figure/Condition-1.svg)

You can also split the tissues into types while retaining the outline

``` r
gganatogram(data=hgMale_key, outline = T, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") +
facet_wrap(~type, ncol=4) +
theme_void()
```

![](figure/maleOrgans-1.svg)

## Added female option

All female tissues

``` r
hgFemale_key$organ
#>  [1] "atrial_appendage"          "ectocervix"               
#>  [3] "hippocampus"               "pleura"                   
#>  [5] "bronchus"                  "trachea"                  
#>  [7] "lung"                      "tonsil"                   
#>  [9] "submandibular_gland"       "breast"                   
#> [11] "spinal_cord"               "pancreas"                 
#> [13] "liver"                     "colon"                    
#> [15] "bone_marrow"               "urinary_bladder"          
#> [17] "stomach"                   "duodenum"                 
#> [19] "esophagus"                 "gall_bladder"             
#> [21] "spleen"                    "small_intestine"          
#> [23] "placenta"                  "endometrium"              
#> [25] "vagina"                    "aorta"                    
#> [27] "pituitary_gland"           "gastroesophageal_junction"
#> [29] "caecum"                    "appendix"                 
#> [31] "ileum"                     "left_atrium"              
#> [33] "left_ventricle"            "pulmonary_valve"          
#> [35] "mitral_valve"              "diaphragm"                
#> [37] "bone"                      "cartilage"                
#> [39] "throat"                    "rectum"                   
#> [41] "nasal_septum"              "nasal_pharynx"            
#> [43] "cerebellum"                "cerebellar_hemisphere"    
#> [45] "prefrontal_cortex"         "frontal_cortex"           
#> [47] "nose"                      "temporal_lobe"            
#> [49] "cerebral_cortex"           "kidney"                   
#> [51] "renal_cortex"              "coronary_artery"          
#> [53] "tricuspid_valve"           "thyroid_gland"            
#> [55] "skin"                      "parotid_gland"            
#> [57] "adipose_tissue"            "heart"                    
#> [59] "smooth_muscle"             "brain"                    
#> [61] "adrenal_gland"             "lymph_node"               
#> [63] "skeletal_muscle"           "ovary"                    
#> [65] "leukocyte"                 "salivary_gland"           
#> [67] "fallopian_tube"            "uterus"                   
#> [69] "uterine_cervix"            "nerve"
gganatogram(data=hgFemale_key, outline = T, fillOutline='#a6bddb', organism='human', sex='female', fill="colour")  +theme_void()
```

![](figure/female-1.svg)

You can also split the tissues into types while retaining the outline

``` r
gganatogram(data=hgFemale_key, outline = T, fillOutline='#a6bddb', organism='human', sex='female', fill="colour") +
facet_wrap(~type, ncol=4) +
theme_void()
```

![](figure/femaleOrgans-1.svg)

To display the female reproductive system with outline.

``` r
hgFemale_key %>%
    dplyr::filter(type=='reproductive') %>%
    gganatogram( outline = T, fillOutline='#a6bddb', organism='human', sex='female', fill="colour")  +
    theme_void()  +
    coord_cartesian(xlim = c(30, 75), ylim = c(-110, -80))
```

![](figure/femaleRep-1.svg)

# Added mouse

## Male

``` r
mmMale_key$organ
#>  [1] "aorta"                     "brown_adipose_tissue"     
#>  [3] "stomach"                   "duodenum"                 
#>  [5] "pancreas"                  "spleen"                   
#>  [7] "adrenal_gland"             "kidney"                   
#>  [9] "colon"                     "small_intestine"          
#> [11] "caecum"                    "jejunum"                  
#> [13] "ileum"                     "esophagus"                
#> [15] "gall_bladder"              "lymph_node"               
#> [17] "seminal_vesicle"           "penis"                    
#> [19] "femur"                     "bone_marrow"              
#> [21] "cartilage"                 "quadriceps_femoris"       
#> [23] "spinal_cord"               "lung"                     
#> [25] "diaphragm"                 "trachea"                  
#> [27] "hindlimb"                  "trigeminal_nerve"         
#> [29] "sciatic_nerve"             "intestinal_mucosa"        
#> [31] "liver"                     "heart"                    
#> [33] "brain"                     "skeletal_muscle"          
#> [35] "circulatory_system"        "blood_vessel"             
#> [37] "skin"                      "prostate_gland"           
#> [39] "vas_deferens"              "epididymis"               
#> [41] "testis"                    "urinary_bladder"          
#> [43] "thymus"                    "peripheral_nervous_system"
#> [45] "eye"
gganatogram(data=mmMale_key, outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +theme_void()  +coord_fixed()
```

![](figure/maleMouse-1.svg)

``` r

gganatogram(data=mmMale_key, outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +theme_void()+facet_wrap(~type, ncol=4)
```

![](figure/maleMouseOrgan-1.svg)

## Female

``` r
mmFemale_key$organ
#>  [1] "aorta"                     "circulatory_system"       
#>  [3] "brown_adipose_tissue"      "stomach"                  
#>  [5] "duodenum"                  "pancreas"                 
#>  [7] "spleen"                    "adrenal_gland"            
#>  [9] "kidney"                    "colon"                    
#> [11] "small_intestine"           "caecum"                   
#> [13] "jejunum"                   "ileum"                    
#> [15] "esophagus"                 "gall_bladder"             
#> [17] "vagina"                    "uterus"                   
#> [19] "urinary_bladder"           "tongue"                   
#> [21] "Peyer's_patch"             "femur"                    
#> [23] "bone_marrow"               "cartilage"                
#> [25] "quadriceps_femoris"        "skeletal_muscle"          
#> [27] "spinal_cord"               "diaphragm"                
#> [29] "hindlimb"                  "trigeminal_nerve"         
#> [31] "eye"                       "intestinal_mucosa"        
#> [33] "brain"                     "heart"                    
#> [35] "liver"                     "sciatic_nerve"            
#> [37] "blood_vessel"              "skin"                     
#> [39] "mammary_gland"             "title8178"                
#> [41] "reproductive_system"       "lymph_node"               
#> [43] "thymus"                    "thyroid_gland"            
#> [45] "lung"                      "peripheral_nervous_system"
#> [47] "trachea"
gganatogram(data=mmFemale_key, outline = T, fillOutline='#a6bddb', organism='mouse', sex='female', fill="colour")  +theme_void()    +coord_fixed()
```

![](figure/femaleMouse-1.svg)

``` r

gganatogram(data=mmFemale_key, outline = T, fillOutline='#a6bddb', organism='mouse', sex='female', fill="colour")  +theme_void()+facet_wrap(~type, ncol=4)
```

![](figure/femaleMouseOrgan-1.svg)

## Cellular structures

I have now included cellular substructures, using the cell.svg from the
Protein Atlas. If you use the main cell figure (hopefully more will be
added), please cite [Thul PJ et
al. 2017](http://science.sciencemag.org/content/356/6340/eaal3321) <br/>

The cellular data can be access using cell_key

``` r
length(cell_key)
#> [1] 1
cell_key
#> $cell
#>                            organ  type    colour       value
#> 1                        cytosol other steelblue  2.07159434
#> 4         intermediate_filaments other   #984EA3 14.89497057
#> 6                actin_filaments other   #FFFF33  5.87440944
#> 8           focal_adhesion_sites other   #F781BF  8.12483660
#> 10 microtubule_organizing_center other   #66C2A5  8.67564889
#> 12                    centrosome other   #8DA0CB  1.02852838
#> 13                  microtubules other   #E78AC3  9.48882657
#> 16              microtubule_ends other   #E5C494  4.80457195
#> 18             secreted_proteins other   #8DD3C7  9.20191105
#> 20                lipid_droplets other   #BEBADA  3.48903574
#> 22                     lysosomes other   #80B1D3  3.73790434
#> 24                   peroxisomes other   #B3DE69  6.79465458
#> 26                     endosomes other   #D9D9D9 13.48636296
#> 28         endoplasmic_reticulum other   #CCEBC5 11.36654344
#> 30               golgi_apparatus other   #7FC97F 11.29225961
#> 32                   nucleoplasm other   #FDC086  2.07964782
#> 34              nuclear_membrane other   #386CB0  7.98595837
#> 36                nuclear_bodies other   #BF5B17  0.05868359
#> 38              nuclear_speckles other   #1B9E77  0.61672243
#> 40                      nucleoli other   #7570B3 14.96900579
#> 42     nucleoli_fibrillar_center other   #66A61E  8.72324527
#> 44                rods_and_rings other   #A6761D  9.53194209
#> 46                  mitochondria other   #A6CEE3  1.29396698
#> 48               plasma_membrane other   #B2DF8A 13.45657571
```

To plot the whole cell with colours or values, use the following
command. If you want to specify a background colour, you either have to
remove the cytosol or change the colour of cytosol to the desired
colour.

``` r
gganatogram(data=cell_key[['cell']], outline = T, fillOutline='steelblue', organism="cell", fill="colour")  +theme_void()   + coord_fixed()
```

![](figure/unnamed-chunk-10-1.svg)

``` r


gganatogram(data=cell_key[['cell']], outline = T, fillOutline='lightgray', organism="cell", fill="value")  +theme_void() +  coord_fixed() +  scale_fill_viridis()
```

![](figure/unnamed-chunk-10-2.svg)

To see all the subsstructures individually, you can plot the data one at
a time

``` r
figureList <- list()
for (i in 1:nrow(cell_key[['cell']])) {
    figureList[[i]] <- gganatogram(data=cell_key[['cell']][i,], outline = T, fillOutline='steelblue', organism="cell", fill="colour")  +theme_void() +ggtitle(cell_key[['cell']][i,]$organ) + theme(plot.title = element_text(hjust=0.5, size=16)) + coord_fixed()
}

do.call(grid.arrange,  c(figureList[1:4], ncol=2))
```

![](figure/unnamed-chunk-11-1.svg)

``` r

do.call(grid.arrange,  c(figureList[5:8], ncol=2))
```

![](figure/unnamed-chunk-11-2.svg)

``` r

do.call(grid.arrange,  c(figureList[9:12], ncol=2))
```

![](figure/unnamed-chunk-11-3.svg)

``` r

do.call(grid.arrange,  c(figureList[13:16], ncol=2))
```

![](figure/unnamed-chunk-11-4.svg)

``` r

do.call(grid.arrange,  c(figureList[17:20], ncol=2))
```

![](figure/unnamed-chunk-11-5.svg)

``` r

do.call(grid.arrange,  c(figureList[21:24], ncol=2))
```

![](figure/unnamed-chunk-11-6.svg)

## Other organisms i.e. tier 2 organisms

Expression atlas contains other organisms than human and mice, however,
these are not as well anotated. <br/> All the expression atlas
anatograms can be found here
<https://ebi-gene-expression-group.github.io/anatomogram/> <br/>
Unfortunately, I won’t be able to add other organs to these since I’m
neither an anatomist nor artist. <br/> If anyone would like to add more
organs, I would love for you to contribute. <br/> <br/> To create these
plots, I have added two other objects other_key and other_list. <br/>
These are lists within lists, and to plot all the organs from an
organisms use other_key\[\[“organism”\]\] as data, and “organism” as
organism. <br/> Also, the organ names are so far a mix of UBERON and
plant ids.

``` r
length(other_key)
#> [1] 24
names(other_key)
#>  [1] "anolis_carolinensis"                 
#>  [2] "arabidopsis_thaliana"                
#>  [3] "bos_taurus"                          
#>  [4] "brachypodium_distachyon.flower_parts"
#>  [5] "brachypodium_distachyon.whole_plant" 
#>  [6] "gallus_gallus"                       
#>  [7] "hordeum_vulgare.flower_parts"        
#>  [8] "hordeum_vulgare.whole_plant"         
#>  [9] "macaca_mulatta"                      
#> [10] "monodelphis_domestica"               
#> [11] "oryza_sativa.flower_parts"           
#> [12] "oryza_sativa.whole_plant"            
#> [13] "papio_anubis"                        
#> [14] "rattus_norvegicus"                   
#> [15] "solanum_lycopersicum.flower_parts"   
#> [16] "solanum_lycopersicum.whole_plant"    
#> [17] "sorghum_bicolor.flower_parts"        
#> [18] "sorghum_bicolor.whole_plant"         
#> [19] "tetraodon_nigroviridis"              
#> [20] "triticum_aestivum.flower_parts"      
#> [21] "triticum_aestivum.whole_plant"       
#> [22] "xenopus_tropicalis"                  
#> [23] "zea_mays.flower_parts"               
#> [24] "zea_mays.whole_plant"
```

To plot bos_taurus use the following command. Unfortunately, I have not
managed to add the correct names yet.

``` r
other_key[["bos_taurus"]]
#>             organ  type  colour     value
#> 2        duodenum other #E41A1C 11.381132
#> 3           brain other #377EB8  2.264810
#> 4          kidney other #4DAF4A  4.131599
#> 5            lung other #984EA3  3.182946
#> 6           colon other #FF7F00  3.114481
#> 7           heart other #FFFF33 13.141334
#> 8           liver other #A65628 17.251310
#> 9  pulmonary vein other #F781BF 13.414659
#> 19 UBERON_0001013 other #999999 12.126515
#> 20 UBERON_0001013 other #66C2A5  1.898023
#> 21 UBERON_0001013 other #FC8D62 19.290389
#> 22 UBERON_0014892 other #8DA0CB 10.994221
#> 23 UBERON_0014892 other #E78AC3 16.761115
#> 24 UBERON_0014892 other #A6D854  2.468627
#> 25 UBERON_0014892 other #FFD92F  1.556285
#> 26 UBERON_0014892 other #E5C494  3.461740
#> 27 UBERON_0014892 other #B3B3B3 18.595027

gganatogram(data=other_key[["bos_taurus"]], outline = T, fillOutline='white', organism="bos_taurus", sex='female', fill="colour")  +
        theme_void() +
        ggtitle("bos_taurus") + 
        theme(plot.title = element_text(hjust=0.5)) + 
        coord_fixed()
```

![](figure/bosTaurus-1.svg)

Here is a way to loop through all the other organisms and plot their
organs.

``` r
library(gridExtra)
plotList <- list()
for (organism in names(other_key)) {
    plotList[[organism]] <- gganatogram(data=other_key[[organism]], outline = T, fillOutline='white', organism=organism, sex='female', fill="colour")  +
                theme_void() +
                ggtitle(organism) + 
                theme(plot.title = element_text(hjust=0.5, size=9)) + 
                coord_fixed()
}

do.call(grid.arrange,  c(plotList[1:4], ncol=2))
```

![](figure/othersFirst12-1.svg)

``` r

do.call(grid.arrange,  c(plotList[5:8], ncol=2))
```

![](figure/othersFirst12-2.svg)

``` r

do.call(grid.arrange,  c(plotList[9:12], ncol=2))
```

![](figure/othersFirst12-3.svg)

``` r

do.call(grid.arrange,  c(plotList[13:16], ncol=2))
```

![](figure/othersFirst12-4.svg)

``` r

do.call(grid.arrange,  c(plotList[17:20], ncol=2))
```

![](figure/othersFirst12-5.svg)

``` r

do.call(grid.arrange,  c(plotList[21:24], ncol=2))
```

![](figure/othersFirst12-6.svg)

# Added mouse brain

``` r
mouse_brain_key <- data.frame(organ = names(gganatogram::other_list[['mus_musculus.brain']])[1:12], 
          type = 'brain', 
          colour = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
                     "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
          value = rnorm(n =12, mean = 100, sd  =50))


gganatogram(data=mouse_brain_key, outline = T, fillOutline='white', organism="mus_musculus.brain", sex='female', fill="colour") +
  theme_void() 
```

![](figure/unnamed-chunk-13-1.svg)

``` r

mouse_brain_key$type <- mouse_brain_key$organ



gganatogram(data=mouse_brain_key, outline = T, fillOutline='white', organism="mus_musculus.brain", sex='female', fill="colour") +
  theme_void() +
  facet_wrap(~type)
```

![](figure/unnamed-chunk-13-2.svg)
