This repository initially contained a Pyomo implementation of the model developed in the paper "Centralised Planning of National Integrated Energy System with Power-to-Gas and Gas Storages" ([DOI](https://doi.org/10.1049/cp.2018.1912)) published at MEDPOWER 2018.

The original model was extended and the updated model was discussed in a follow-up paper titled "The Role of Power-to-Gas and Carbon Capture Technologies in Cross-Sector Decarbonisation Strategies" ([DOI](https://doi.org/10.1016/j.epsr.2019.106039)) that was published in Electric Power Systems Research in 2020. A Pyomo implementation of the updated model then replaced the original one in this repository.

The Pyomo implementation of the updated model has now been replaced by a GBOML implementation. The GBOML implementation (model/belgian_energy_system_model.txt) is much  simpler to use and analyse than the original Pyomo model. GBOML is pip-installable and the model can be run by calling the GBOML executable on the input file directly in the command line (e.g., by typing 'gboml model/belgian_energy_system_model.txt --cplex --opt solver_params/cplex.opt --col_csv --detailed solver_params/query_attrs_cplex.txt' to run the model, solve it with CPLEX - barrier algorithm and no cross-over - and print the results to file, including the dual variables associated with the hyperedges). More information about GBOML can be found in the [online documentation](https://gboml.readthedocs.io/en/latest/) or in its [host repository](https://gitlab.uliege.be/smart_grids/public/gboml).

The repository also includes the full data required to run the GBOML model as well as basic plotting and postprocessing scripts. Specifically, time series are stored in the data/ folder, while scalar parameters are directly encoded into the GBOML input file. Postprocessing scripts, on the other hand, can be found in the postprocessing_scripts/ folder.

Based on the input data provided in this repository, the different scenarios should lead to the following objectives (system cost over five years):

1. Scenario 1 - 3.3501e+05
2. Scenario 2 - 2.5315e+05
3. Scenario 3 - 2.0508e+05
4. Scenario 4 - 6.1521e+04
5. Scenario 5 - 4.0455e+04

Before running the model with GBOML, make sure that the data paths specified in the input file are consistent.

Any enquiries or comments can be emailed to mathias.berger@alumni.duke.edu. Alternatively, an issue can be opened in this repository directly.