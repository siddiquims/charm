# CHarm

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
 translated protein is available, this might not have a big impact, as those predicted link/end segments might
 not represent the real protein structure.

### Features

 1. **Easy to use:** Self-explaning command line frontend: Obtains codon usage table directly from
 http://www.kazusa.or.jp/codon. You only have to enter the species id (e.g. 83333 for *E. coli K12*) and the path to
 your input sequence.
 2. **Cross platform:** CHarm is written in [Python][1] and uses wide-spread modules like [matplotlib][2],
 [NumPy][3] and [Biopython][4] for data processing and visualization. It will run on any platform that is supported
 by Python.
 3. **Open source:** CHarm is licensed under the [MIT License][5].
 You can freely alter its codebase as it fits your needs.


----------
## References
1. <a name="angov2009">**Angov, E., Hillier, C. J., Kincaid, R. L., & Lyon, J. a. (2008)**</a>.
Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those
of the expression host. PloS one, 3(5), e2189. doi:10.1371/journal.pone.0002189
2. <a name="thanaraj1996">**Thanaraj, T. a, & Argos, P. (1996)**</a>. Protein secondary structural types are
differentially coded on messenger RNA. Protein science : a publication of the Protein Society, 5(10),
1973–83. doi:10.1002/pro.5560051003


  [1]: http://www.python.org "Python"
  [2]: http://www.matplotlib.org "Matplotlib"
  [3]: http://www.numpy.org "NumPy"
  [4]: http://www.biopython.org "Biopython"
  [5]: http://opensource.org/licenses/MIT "MIT License"
