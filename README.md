## TransTEx (Transcript-level Tissue Expression)
### Novel tissue-specificity scoring method for grouping human transcriptome into different expression groups

`TransTEx` is a comprehensive R package developed by the Davuluri Lab for analyzing and exploring tissue-specific transcript/gene expression patterns. The package facilitates the identification of tissue-specific and tissue-enhanced transcripts across various conditions, leveraging bulk-transcriptomic data to uncover insights into gene expression dynamics.

### Important links:
* [Paper](https://doi.org/10.1093/bioinformatics/btae475)
* [Database](https://bmi.cewit.stonybrook.edu/transtexdb/)
  ##### Note: The database works on LTE or Wi-Fi (not both sometimes). Please try other networks in case one does not work.
* [Talk](https://youtu.be/p4S01JPJHsU?si=P8Jq7A6VVY1on-sQ)

### Installation
Before installing `TransTEx`, ensure you have R installed on your system. `TransTEx` can be installed directly from GitHub using the `remotes` package. If you do not have `remotes` installed, you can install it using the following command:

``` r
install.packages("remotes")
```

With `remotes` installed, proceed to install `TransTEx`:

``` r
BiocManager::install("edgeR")
remotes::install_github("pallavisurana1/TransTEx")
```

### Quick Start
Load `TransTEx` into your R session:

``` r
library(TransTEx)
```

Explore the functions and datasets provided by `TransTEx`:

``` r
help(package = "TransTEx")
```

### Usage Examples
For detailed examples of using `TransTEx` to analyze tissue-specific gene expression, please refer to the vignettes:

``` r
browseVignettes("TransTEx")
```

### Reporting Issues
If you encounter any bugs or issues while using `TransTEx`, please report them on the [Issues](https://github.com/pallavisurana1/TransTEx/issues) page of the GitHub repository. We appreciate your contributions to improving `TransTEx`.

### Citation

If you use TransTEx in your research, please cite:
```bibtex
@article{Surana2024,
  author    = {Pallavi Surana and Pratik Dutta and Ramana V Davuluri},
  title     = {TransTEx: novel tissue-specificity scoring method for grouping human transcriptome into different expression groups},
  journal   = {Bioinformatics},
  volume    = {40},
  number    = {8},
  year      = {2024},
  article   = {btae475},
  doi       = {10.1093/bioinformatics/btae475},
  url       = {https://doi.org/10.1093/bioinformatics/btae475}
}
```

### Contact
For additional information or inquiries about `TransTEx`, please contact Pallavi Surana (pallavi.surana\@stonybrook.edu)

### License
`TransTEx` is made available under the (specify license, e.g., MIT License). See the LICENSE file for more details.
