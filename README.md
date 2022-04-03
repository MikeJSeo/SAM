Web-enabled and Cross-platform SAM via Shiny
============================================

First make sure you have a very recent version of R or RStudio.

Next install required packages. Cut and paste what’s below in an R
session. You only need to do this once.

``` r
install.packages(c("matrixStats", "GSA", "shiny", "openxlsx", "Rcpp"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")

#install.packages("devtools") if devtools not installed
#library(devtools)
#install_github("cran/samr")
install.packages("samr")
```

You may need to download R tools:

<a href="https://cran.r-project.org/bin/windows/Rtools/" class="uri">https://cran.r-project.org/bin/windows/Rtools/</a>

Then, you may run SAM any time in an R session as follows.

``` r
library(shiny)
library(impute)
runGitHub("SAM", "MikeJSeo")
```

That will bring up a browser window with a user interface. More details
are provided in the sam-manual.pdf in this directory. Data examples and
sample output are in a seperate folder to view.

GSA.correlate.revised.R, GSA.listsets.revised.R, and GSA.plot.revised.R
are revised version of the function in GSA. server.r and ui.r contains
the main code that generates the GUI.

Please post to the group regarding any issues. This will help us ensure
we have all the kinks ironed out before merging the code into the next
version of the `samr` package.

Note that users should use Firefox or Chrome as the \_default\_browser
on Windows: IE will not work.

Thank you,  
Michael Seo
