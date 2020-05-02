# hubViz
<!--
hubViz: A Novel Tool for Hub-centric Visualization
-->

hubViz (a novel tool for **hub**-centric **Vis**ualization) is a novel tool for hub-centric visualization which is based on a latent space joint model (LSJM) for visualization.
'hubviz' package provides computationally efficient and user friendly interface to fit the hubViz models. 
The following help pages provide a good start point for the genetic analysis using the 'hubviz' package, including the overview of 'hubviz' package and the example command lines:

```
library(hubviz)
package?hubviz
class?hubviz
```

Installation
============ 

The stable versions of 'hubviz' package can be obtained from the following URLs:

Package source: [https://github.com/dongjunchung/chunglab_binary_packages/blob/master/hubviz_0.1.tar.gz](https://github.com/dongjunchung/chunglab_binary_packages/blob/master/hubviz_0.1.tar.gz?raw=true)

Windows binary: [https://github.com/dongjunchung/chunglab_binary_packages/blob/master/hubviz_0.1.zip](https://github.com/dongjunchung/chunglab_binary_packages/blob/master/hubviz_0.1.zip?raw=true)

Mac OS/X binary: [comming soon](https://)

To install the developmental versions of 'hubviz' package, it's easiest to use the 'devtools' package. Note that the ‘hubviz’ package depends on the ‘Rcpp’ package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/hubviz")
```

References
==========
Nam JH, Yun J, Jin IH and Chung D (2020) ''hubViz: A novel tool for hub-centric visualization''.
