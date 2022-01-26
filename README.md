# contacTrees case study on Celtic, Germanic and Romance languages
The BEAST 2 package _contacTrees_ (https://github.com/NicoNeureiter/contacTrees) offers a model for horizontal transfer (borrowing) in linguistic phylogenies. This repository contains data and python scripts to set up a case study on Indo-European languages (Celtic, Germanic, Romance), as described in the Neureiter et al. (2022).

The data is provided in the `ressources` directory and all scripts are provided in the `contacTreesIE` directory. You can reproduce the case study using the following steps:
1. Run `python -m contacTreesIE.compile_beast_xmls` to generate the XML files for the BEAST analysis.
2. Run BEAST 2 using an XML file in one of the generated sub-directories of `runs`.
3. The results will contain a `.log` and a `.trees` file. The trees can be summarised using the ContactreesAnnotator app:
```applauncher ContactreesAnnotator -threshold 50 -burnin 10 CT_full_1.trees CT_full_1.summary.tree```
4. The python scripts in `contacTreesIE.plotting` provides various ways of visualising the results. E.g. `plot_contact_edges.py` will generate a plot for each summarised contact edge including the data for the corresponding loanwords reconstructed at this edge.

The results of the analysis in Neureiter et al. (2022) are avaialbe in the `results` directory and visualised in the `visualisations` directory.


> Neureiter N, Ranacher P, Efrat-Kowalsky N, Kaiping G, Weibel R, Widmer P, and Bouckaert R R. "Detecting contact in language trees: a Bayesian phylogenetic model with horizontal transfer", 2022 (in review).
