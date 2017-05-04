# observer-free-morphotype-characterization
These R scripts accompany the manuscript "Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape"

A snapshot of the code as used for the manuscript is reposited at 10.6084/m9.figshare.4968170. Analyses are centered around the script `sensitivity_base.R`, sourcing necessary internal and external functions from `sensitivity_dependencies.R` . 

As part of the code for data analysis, the repository contains functions to:
* calculate per-principal component repeatibility of an auto3dgm alignment based on duplicate specimens (find_repeatablePCs.R)
* determine whether the principal components of one analyses should be inverted in order to better correspond to a second set of principal components
* perform the Anderson test on the standard deviations resulting from a principal components analysis
* create heat maps visualizing the differences between principal component extremes or between two shapes

If you use part of this code in published materials, please cite:

Vitek, N.S., Manz, C.L., Gao, T. Bloch, J.I., Strait, S.G., Boyer, D.M. In Press. Semi-supervised determination of pseudocryptic morphotypes using observer-free characterizations of anatomical alignment and shape. Ecology and Evolution. 
