# CAM

CAM is a constraint programming based approach to compute atom-atom maps for chemical reactions based on the identification of chemically feasible imaginary transition states (ITS). 

Chemical reactions are rearrangements of chemical bonds. Each atom in an educt molecule thus appears again in a specific position of one of the reaction products. This bijection between educt and product atoms is not reported by chemical reaction databases, however, so that the “Atom Mapping Problem” of finding this bijection is left as an important computational task for many practical applications in computational chemistry and systems biology. Elementary chemical reactions feature a cyclic imaginary transition state (ITS) that imposes additional restrictions on the bijection between educt and product atoms that are not taken into account by previous approaches. We demonstrate that Constraint Programming is well-suited to solving the Atom Mapping Problem in this setting. The performance of our approach is evaluated for a manually curated subset of chemical reactions from the KEGG database featuring various ITS cycle layouts and reaction mechanisms.


## Contribution

Feel free to contribute to this project by wirting [Issues](https://github.com/BackofenLab/CAM/issues) with feature requests or bug reports.

## Cite
If you use IntaRNA, please cite our [article](http://almob.biomedcentral.com/articles/10.1186/s13015-014-0023-3):
```
DOI: 10.1186/s13015-014-0023-3
```
