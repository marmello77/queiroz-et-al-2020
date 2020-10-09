# queiroz_et_al_2020

Supplement to the paper Queiroz *et al*. (2020, Biotropica).

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Joel A. Queiroz, Ugo M. Diniz, Diego P. Vázquez, Zelma M. Quirino, Francisco A.R. Santos, Marco A.R. Mello, Isabel C. Machado.

E-mail: marmello@usp.br. 

Published on October 2nd, 2020 (English version).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4064012.svg)](https://doi.org/10.5281/zenodo.4064012)


Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again".

Disclaimer: You may freely use the software and data provided here for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this material, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the material without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this material helps you produce any academic work (paper, book, chapter, monograph, dissertation, report or similar), please acknowledge the authors and cite the source.


## Functionality and origin

The data, scripts, and functions provided here aim at making or study fully reproducible. You will find code to reproduce both the analysis and the figures, as well as the main supplementary material. We have also added some bonus material related to data visualization.


## List of folders and files

1. data (folder)

    a. morph_pla.xlsx -> morphometric data of plants.
    
    b. morph_pol.xlsx - > morphometric data of pollinators.
    
    c. network.txt -> the nocturnal pollination network composed of plants, bats, and hawkmoths.
    
    d. plants.xlsx -> guild data of plants.
    
    e. sampbat.xlsx -> capture data of bats.
    
    f. samphawk.xlsx -> capture data of hawkmoths.
    
    g. sampling_bats.xlsx -> sampling data of bats at the individual level.
    
    h. sampling_hawkmoths.xlsx -> sampling data of hawkmoths at the individual level.
    
    i. vertices.csv and edges.csv -> additional data to reproduce Figure 2 as it is presented in the paper. See "tutorial.pdf".
    

2. figures (folder)

    a. centrality.tiff -> boxplots of centrality by pollinatiom syndrome.
    
    b. compound.png -> matrix of the nocturnal pollination network.

    c. morph.tiff -> boxplot comparing morphometric data between network modules.

    d. network.tif -> graph of the nocturnal pollination network.

    e. sampling.tiff -> rarefaction curves representing interaction sampling completeness.


3. results (folder)

    a. dplants.csv -> scores of specialization (d') for the plants.

    b. BCplants.csv -> scores of betweenness centrality for the plants.

    c. NDplants.csv -> scores of normalized degree for the plants.

    d. results_centrality.txt -> results of the GLMs comparing centraliy by pollination syndrome.
    
    e. results_topology.txt -> results of the topological analyses.


4. analysis.R -> main script to run the analyses and plot the graphs. Best choice if you want not only to reproduce our analyses and figures, but also to aplly our solutions to your own data.

5. MyDiamond.R -> additional function to create the diamond shape used in "network.tif". We found this script on the internet some years ago, and could not identify its authors. We thank the person who wrote it, whomever it may be.

6. PosteriorProb.R -> additional function to calculate restricted probabilities to be used in the compound topology analysis. See the original [repo](https://github.com/gabrielmfelix/Restricted-Null-Model).

7. RestNullModel.R -> additional function to run the compound topology analysis. See the original [repo](https://github.com/gabrielmfelix/Restricted-Null-Model).

8. tutorial.pdf -> knit version of the tutorial provided in "tutorial.Rmd". Best choice if are not familiar with R, and would like to just take a look at how our paper was cooked.

9. tutorial.Rmd -> a tutorial written in RMarkDown to help you reproduce the analyses and figures. Best choice if you just want to reproduce our results, but does not want to apply them to your own data.


## Instructions

1. Run the script "analysis.R" and follow the steps provided in it;

2. You can also run specific sections of the script to reproduce specific results, but beware of dependences between some sections.

3. The analyses based on null models can be quite resource-consuming. Check your computer's power and memory before setting the number of permutations. If you want to analyze your own data for real, consider that you will need at least 1,000 permutations. If you want only to explore your data, 10 permutations will do it.

4. If you have any questions, corrections, or suggestions, please feel free to open an [issue](https://github.com/marmello77/queiroz-et-al-2020/issues) or make a [pull request](https://github.com/marmello77/queiroz-et-al-2020/pulls).

5. Have fun!


## Original sources

The codes used in this repo have borrowed solutions from other repos developed by or lab:

1. [network-significance](https://github.com/marmello77/network-significance)

2. [Restricted-Null-Model](https://github.com/gabrielmfelix/Restricted-Null-Model)

3. [ihsmodel](https://github.com/pinheirorbp/ihsmodel)


## Acknowledgments

We thank our colleagues, who helped us at different stages of this project. Special thanks go to [Renata Muylaert](https://github.com/renatamuy), [Pavel Dodonov](https://github.com/pdodonov), [Gabriel Felix](https://github.com/gabrielmfelix), and [Rafael Pinheiro](https://github.com/pinheirorbp), who, together with us, developed the additional functions and scripts used here. We also thank our sponsors, especially the Alexander von Humboldt-Stiftung, CNPq, CAPES, and FAPESP, who gave us grants, fellowships, and scholarships. Last, but not least, we thank the [Stack Overflow Community](https://stackoverflow.com), where we solve most of our coding dilemmas.


## Reference

Queiroz JA, Diniz UM, Vázquez DP, Quirino ZM, Santos FAR, Mello MAR, Machado IC. 2020. Bats and hawkmoths form mixed modules with flowering plants in a nocturnal interaction network. Biotropica, *accepted*.
