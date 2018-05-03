[![Build Status](https://travis-ci.org/BackofenLab/CAM.svg?branch=master)](https://travis-ci.org/BackofenLab/CAM)

# CAM

CAM is a constraint programming based approach to compute atom-atom maps for chemical reactions based on the identification of chemically feasible imaginary transition states (ITS). 

Chemical reactions are rearrangements of chemical bonds. Each atom in an educt molecule thus appears again in a specific position of one of the reaction products. This bijection between educt and product atoms is not reported by chemical reaction databases, however, so that the “Atom Mapping Problem” of finding this bijection is left as an important computational task for many practical applications in computational chemistry and systems biology. Elementary chemical reactions feature a cyclic imaginary transition state (ITS) that imposes additional restrictions on the bijection between educt and product atoms that are not taken into account by previous approaches. We demonstrate that Constraint Programming is well-suited to solving the Atom Mapping Problem in this setting. The performance of our approach is evaluated for a manually curated subset of chemical reactions from the [KEGG database](http://www.genome.jp/kegg/) featuring various ITS cycle layouts and reaction mechanisms.

## Dependencies

To compile CAM you need

- [boost](http://www.boost.org/) library >= 1.55
- [Graph Grammar Library](https://github.com/BackofenLab/GGL) (ggl) >= 4.1.1
- [Gecode](http://www.gecode.org/) library >= 4.0 < 5.0 *NOTE: CAM requires version 4!*


## Contribution

Feel free to contribute to this project by raising [Issues](https://github.com/BackofenLab/CAM/issues) with feature requests or bug reports.

## Cite
If you use CAM, please cite our articles:

- [Atom mapping with constraint programming](http://almob.biomedcentral.com/articles/10.1186/s13015-014-0023-3),
  Martin Mann, Feras Nahar, Norah Schnorr, Rolf Backofen, Peter F Stadler and Christoph Flamm,
  Algorithms for Molecular Biology, 2014, 9:23, [DOI:10.1186/s13015-014-0023-3](http://almob.biomedcentral.com/articles/10.1186/s13015-014-0023-3)
  = `cam` binary
- [Kekule structure enumeration yields unique SMILES](http://www.bioinf.uni-freiburg.de/Publications/Mann_kekule_13.pdf),
   Martin Mann and Bernhard Thiel,
   In Proceedings of the Workshop on Constraint Based Methods for Bioinformatics (WCB 2013), pages 1-9, 2013
   = `molKekule` binary
