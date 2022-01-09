#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa

from optlang.exceptions import SolverError

#from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_matlab_model, load_json_model


from pytfa.io import import_matlab_model, load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo
from pytfa.analysis.variability import variability_analysis
import pickle

thermo_database = '/projectnb2/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/data/thermo_data.thermodb'
out_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/tests/singleCoenzymeModel.01092021.fva'

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = GUROBI

case = 'full' # 'reduced' or full'

# Load reaction DB
print("Loading thermo data...")
thermo_data = load_thermoDB(thermo_database)

print("Done !")
#biomass_rxn = 'BIOMASS_Ec_iJO1366_WT_53p95M'
biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'

# We import pre-compiled data as it is faster for bigger models
model_path = '/projectnb2/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/models'
cobra_model = load_json_model(model_path + '/iJO1366_NAD_semi-unconstrained.json')
lexicon = read_lexicon(model_path + '/iJO1366/lexicon.csv')
compartment_data = read_compartment_data(model_path + '/iJO1366/compartment_data.json')

# Initialize the cobra_model
mytfa = pytfa.ThermoModel(thermo_data, cobra_model)

# Annotate the cobra_model
annotate_from_lexicon(mytfa, lexicon)
apply_compartment_data(mytfa, compartment_data)


mytfa.name = 'iJO1366[NAD]'
mytfa.solver = solver
mytfa.objective = biomass_rxn

# Solver settings

def apply_solver_settings(model, solver = solver):
    model.solver = solver
    # model.solver.configuration.verbosity = 1
    model.solver.configuration.tolerances.feasibility = 1e-9
    if solver == 'optlang_gurobi':
        model.solver.problem.Params.NumericFocus = 3
    model.solver.configuration.presolve = True

apply_solver_settings(mytfa)


## FBA
fba_solution = cobra_model.optimize()
fba_value = fba_solution.objective_value
print('FBA Solution found : {0:.5g}'.format(fba_value))


## TFA conversion
mytfa.prepare()
mytfa.convert()#add_displacement = True)

## Info on the cobra_model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()
tfa_value = tfa_solution.objective_value
print('TFA Solution found : {0:.5g}'.format(tfa_value))


# It might happen that the model is infeasible. In this case, we can relax 
# thermodynamics constraints:

data_for_later = {}

if tfa_value < 0.1:
    from pytfa.optim.relaxation import relax_dgo

    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5*fba_value
    relaxed_model, slack_model, relax_table = relax_dgo(mytfa)

    original_model, mytfa = mytfa, relaxed_model

    print('Relaxation: ')
    print(relax_table)
    
    tfa_solution = mytfa.optimize()
    tfa_value = tfa_solution.objective_value
    
    data_for_later['original_model'] = original_model
    data_for_later['slack_model'] = slack_model
    data_for_later['relaxed_model'] = mytfa
    data_for_later['relax_table'] = relax_table
    data_for_later['tfa_solution'] = tfa_solution


    print('relaxed TFA Solution found : {0:.5g}'.format(tfa_value))

# Report
print('FBA Solution found : {0:.5g}'.format(fba_value))
print('TFA Solution found : {0:.5g}'.format(tfa_value))


# perform variability analysis on reactions
mytfa.reactions.get_by_id(biomass_rxn).lower_bound = tfa_value - 1e-6
varib_results_ = variability_analysis(mytfa,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.sc.results.Optimal.31Aug2021.csv')


# perform variability analysis on reactions
mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.1
varib_results_ = variability_analysis(mytfa,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.sc.results.gtzp1.31Aug2021.csv')

