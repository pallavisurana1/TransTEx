# TransTEx (Transcript-level Tissue Expression)

**TransTEx** is a comprehensive R package developed by the Davuluri Lab for analyzing and exploring tissue-specific transcript/gene expression patterns. This tool introduces a novel tissue-specificity scoring method to group the human transcriptome into distinct expression categories. TransTEx enables the discovery of tissue-specific and tissue-enhanced transcripts across conditions using bulk-transcriptomic data.

<img width="903" height="509" alt="image" src="https://github.com/user-attachments/assets/f80a3461-ec3b-42d2-a712-8281b9f67ec9" />

---

## 🔗 Important Links

* 📄 [Published Paper](https://doi.org/10.1093/bioinformatics/btae475)
* 🌐 [Online Database](https://bmi.cewit.stonybrook.edu/transtexdb/)
  *Note: The database may work best on either LTE or Wi-Fi — try switching networks if it doesn’t load.*
* 🎥 [Introductory Talk](https://youtu.be/p4S01JPJHsU?si=P8Jq7A6VVY1on-sQ)

---

## 📦 Installation

Before installing `TransTEx`, make sure R is installed on your system. Then install required dependencies and the package itself:

```r
install.packages("remotes")
BiocManager::install("edgeR")
remotes::install_github("pallavisurana1/TransTEx")
```

---

## 🚀 Getting Started

Load the package in R:

```r
library(TransTEx)
```

Explore available functions and documentation:

```r
help(package = "TransTEx")
```

---

## 📘 Usage Examples

Access detailed usage examples and workflows via the package vignettes:

```r
browseVignettes("TransTEx")
```

---

## 🛠 Reporting Issues

Found a bug or issue? Report it on the official GitHub repository:

👉 [GitHub Issues](https://github.com/pallavisurana1/TransTEx/issues)

We welcome community feedback and contributions!

---

## 📖 Citation

If you use **TransTEx** in your work, please cite the following:

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

---

## 📬 Contact

For inquiries or more information about TransTEx, please contact:
**Pallavi Surana**
📧 [pallavi.surana@stonybrook.edu](mailto:pallavi.surana@stonybrook.edu)

---

## 📄 License

TransTEx is made available under the **[MIT License](https://opensource.org/licenses/MIT)**. Please see the LICENSE file for details.
