This repository contains the model and data used in the paper "The Role of Power-to-Gas and Carbon Capture Technologies in Cross-Sector Decarbonisation Strategies" ([DOI](https://doi.org/10.1016/j.epsr.2019.106039)) published in Electric Power Systems Research in 2020.

The model was originally implemented in Pyomo but has now been replaced by a GBOML implementation. The GBOML implementation (belgian_energy_system_model.txt) is much  simpler to use and analyse than the original Pyomo model. GBOML is pip-installable and the model can be run by calling the GBOML executable on the input file directly in the command line (e.g., by typing 'gboml belgian_energy_system_model.txt --gurobi --transposed_csv' to run the model, solve it with Gurobi and print the results to file). More information about GBOML can be found in the [online documentation](https://gboml.readthedocs.io/en/latest/) or in its [host repository](https://gitlab.uliege.be/smart_grids/public/gboml).

The repository also includes the full data required to run the GBOML model. Specifically, time series are stored in the data/ folder, while scalar parameters are directly encoded into the GBOML input file.

Based on the input data provided in this repository, the different scenarios should lead to the following objectives (system cost):

1. Scenario 1 - 3.3501e+05
2. Scenario 2 - 2.5315e+05
3. Scenario 3 - 2.0508e+05
4. Scenario 4 - 6.1521e+04
5. Scenario 5 - 4.0455e+04

Before running the model with GBOML, make sure that the data paths specified in the input file are consistent - currently, the input file must be one level above the data folder in the file hierarchy.

Any enquiries or comments can be emailed to mathias.berger@alumni.duke.edu. Alternatively, an issue can be opened in this repository directly.