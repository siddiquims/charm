# CHarm
[![Build Status](https://travis-ci.org/Athemis/charm.png?branch=master)](https://travis-ci.org/Athemis/charm)

## About

CHarm is a tool harmonizing the codon usage in a DNA sequence prior to heterologous expression.
This is done by matching the codon usage of the target host organism to the codon usage of the
native origin organism.

The underlying algorithm has been described by Angov et al. (2008)<sup>[1](#angov2009)</sup>
The current implementation differs from the described algorithm, though:

 1. It implements by default a **lower threshold** for the codon usage in the target host. This threshold is only
 fallen below if the codon usage in the origin organism is also below this threshold. The reference behavior can be
 restored by explicitly choosing a threshold of zero.
 2. In its current state, CHarm does **not search for putative link/end segments** that require slow translational
 progress as described by Thanaraj & Argos (1996)<sup>[2](#thanaraj1996)</sup>. If no structural data for the
 translated protein is available, this might not have a big impact, as those predicted link/end segments might not
 represent the real protein structure.

### Features

 1. **Easy to use:** Self-explaning command line frontend: Obtains codon usage table directly from
 http://www.kazusa.or.jp/codon. You only have to enter the species id (e.g. 83333 for *E. coli K12*) and the path to
 your input sequence.
 2. **Cross platform:** CHarm is written in [Python][1] and uses wide-spread modules like [matplotlib][2], [NumPy][3]
 and [Biopython][4] for data processing and visualization. It will run on any platform that is supported by Python.
 3. **Open source:** CHarm is licensed under the [MIT License][5]. You can freely alter its codebase as it fits
 your needs.

## Get started
### Installation

Make sure you have all the necessary dependencies in place.

**GNU/Linux users:** Use your distribution's package manager to download and install Python and the dependency modules.

**Mac OS X users:** Download and install installer packages for Python and the dependencies from the respective
download sites (see below).

**Windows users:** Installers for Python and the necessary modules are available at the respective download sites
(see below). Alternatively, if you are member of a degree granting educational institution (e.g. university), you can
apply for a free license of [Enthought Canopy][6] which will provide you with a full-flavored Python environment
including all necessary additional packages.

**Dependencies**

 - [Python 3][7]      (tested with Python 3.3)
 - [NumPy][8]         (tested with NumPy 1.8.0)
 - [Biopython][9]     (tested with Biopython 1.63)
 - [matplotlib][10]   (tested with matplotlib 1.3.1)

### Usage of charm-cli.py

 1. Open a terminal. You can check whether the python executable is available by executing

 ```bash
 python --version
 ```

 This should result in e.g. 'Python 3.4.0'
 2. Open http://www.kazusa.or.jp/codon in a web browser of your choice and navigate to the codon usage tables both of
 the organism of origin and the expression host. Note down the id numbers as they appear in the navigation bar
 ("[..]?species=XXXX"). For *E. coli (K12)* the id is '83333'.
 3. Save the sequence to be harmonized (DNA/RNA) in FASTA format.
 4. Run the harmonization with default options:

 ```bash
 python ./charm-cli.py <id origin> <id host> <path to sequence file>
 ```

 To see all available options, run

 ```bash
 python ./charm-cli.py --help
 ```
  

----------
## References
1. <a name="angov2009">**Angov, E., Hillier, C. J., Kincaid, R. L., & Lyon, J. a. (2008)**</a>. Heterologous
protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the
expression host. PloS one, 3(5), e2189. doi:10.1371/journal.pone.0002189
2. <a name="thanaraj1996">**Thanaraj, T. a, & Argos, P. (1996)**</a>. Protein secondary structural types are
differentially coded on messenger RNA. Protein science : a publication of the Protein Society, 5(10), 1973–83.
doi:10.1002/pro.5560051003


  [1]: http://www.python.org "Python"
  [2]: http://www.matplotlib.org "Matplotlib"
  [3]: http://www.numpy.org "NumPy"
  [4]: http://www.biopython.org "Biopython"
  [5]: http://opensource.org/licenses/MIT "MIT License"
  [6]: https://www.enthought.com/products/canopy/academic/ "Enthought Canopy"
  [7]: http://www.python.org "Python"
  [8]: http://www.numpy.org "NumPy"
  [9]: http://www.biopython.org "Biopython"
  [10]: http://www.matplotlib.org "Matplotlib"
