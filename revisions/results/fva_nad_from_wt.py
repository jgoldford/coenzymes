#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa
import cobra
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
out_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/tests/FVA/12Oct2021'

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
cobra_model = load_json_model(model_path + '/iJO1366_WT_semi-unconstrained.11Oct2021.json')
lexicon = read_lexicon(model_path + '/iJO1366/lexicon.csv')
compartment_data = read_compartment_data(model_path + '/iJO1366/compartment_data.json')

# Initialize the cobra_model
mytfa = pytfa.ThermoModel(thermo_data, cobra_model)

# Annotate the cobra_model
annotate_from_lexicon(mytfa, lexicon)
apply_compartment_data(mytfa, compartment_data)


mytfa.name = 'iJO1366[WT]'
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

    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5
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

double_coenzyme_models = {}
z = mytfa.copy()
double_coenzyme_models['tmodel'] = z
double_coenzyme_models['tmodel-maxGrowth'] = tfa_value
double_coenzyme_models['tmodel-solution'] = tfa_solution



def single_coenzyme_transform(model,reaction_id):
    
    met = {}
    met['nad'] = np.where([x.id == 'nad_c' for x in model.metabolites])[0][0]
    met['nadp'] = np.where([x.id == 'nadp_c' for x in model.metabolites])[0][0]
    met['nadh'] = np.where([x.id == 'nadh_c' for x in model.metabolites])[0][0]
    met['nadph'] = np.where([x.id == 'nadph_c' for x in model.metabolites])[0][0]

    met_objs = {}
    met_objs['nad'] = [x for x in model.metabolites if x.id == 'nad_c'][0]
    met_objs['nadh'] = [x for x in model.metabolites if x.id == 'nadh_c'][0]
    met_objs['nadp'] = [x for x in model.metabolites if x.id == 'nadp_c'][0]
    met_objs['nadph'] = [x for x in model.metabolites if x.id == 'nadph_c'][0]
    
    rxn = [x for x in model.reactions if x.id == reaction_id][0].copy()
    # make a new dictionary with coenzyme swapped
    v = {x:y for x,y in rxn.metabolites.items() if x.id in [x + '_c' for x in list(met)]}
    
    nad_stoich = 0;
    nadh_stoich  = 0;
    set([x.id for x in m.reactions.get_by_id(rxnid).metabolites])
    if len(v) > 1:
        v2 = {}
        for x,y in v.items():
            if x.id == 'nadph_c':
                nadh_stoich = nadh_stoich + y 
                #v2[met_objs['nadh']] = y
            elif x.id == 'nadp_c':
                nad_stoich = nad_stoich + y 
                #v2[met_objs['nad']] = y
            elif x.id == 'nadh_c':
                nadh_stoich = nadh_stoich + y 
                #v2[met_objs['nadph']] = y
            elif x.id == 'nad_c':
                nad_stoich = nad_stoich + y 
                #v2[met_objs['nadp']] = y

        v2[met_objs['nad']] = nad_stoich
        v2[met_objs['nadh']] = nadh_stoich

        rxn.subtract_metabolites(v)
        rxn.add_metabolites(v2)
        #rxn.id = rxn.id + '[condensed]'
        # keep old reaction ID
        #print('adding new reaction '+ rxn.id)
        model.remove_reactions([x for x in model.reactions if x.id == reaction_id][0])
        model.add_reaction(rxn)
        coverted = True
        return model,coverted
    else:
        coverted = False
        return model,coverted



m  = z.copy()

rxns_to_remove = ['NADTRHD','NADPPPS','NADK','THD2pp']
rxn_ids = [x.id for x in m.reactions]

#m.remove_reactions([x for x in m.reactions if x.id == 'NADTRHD'][0])
#m.remove_reactions([x for x in m.reactions if x.id == 'NADPPPS'][0])
#m.remove_reactions([x for x in m.reactions if x.id == 'NADK'][0])
#m.remove_reactions([x for x in m.reactions if x.id in rxns_to_remove][0])

m.remove_reactions(rxns_to_remove)

rxn_ids = [x for x in rxn_ids if x not in rxns_to_remove]
#rxn_ids = [x for x in rxn_ids if x not in ['Ec_biomass_iJO1366_WT_53p95M','Ec_biomass_iJO1366_core_53p95M']]

converted_rxns = []
for rxnid in rxn_ids:
    m,conv = single_coenzyme_transform(m,rxnid)
    if conv:
        converted_rxns.append(rxnid)

        
#m.remove_reactions([x for x in m.reactions if x.id == 'NADTRHD[condensed]'][0])
#m.remove_reactions([x for x in m.reactions if x.id == 'NADPPPS[condensed]'][0])
#m.remove_reactions([x for x in m.reactions if x.id == 'NADK[condensed]'][0])
#m.remove_reactions([x for x in m.reactions if x.id == 'NADDP'][0])

m.id = 'iJO1366[NAD]'
#m.objective = 'BIOMASS_Ec_iJO1366_WT_53p95M'
m.objective = 'Ec_biomass_iJO1366_WT_53p95M'
m.name = 'RelaxedModel iJO1366[NAD]'
m.description = 'RelaxedModel iJO1366[NAD]'

met_objs = {}
met_objs['nad'] = [x for x in m.metabolites if x.id == 'nad_c'][0]
met_objs['nadh'] = [x for x in m.metabolites if x.id == 'nadh_c'][0]
met_objs['nadp'] = [x for x in m.metabolites if x.id == 'nadp_c'][0]
met_objs['nadph'] = [x for x in m.metabolites if x.id == 'nadph_c'][0]
#m.remove_metabolites([met_objs['nadp'],met_objs['nadph']])
m,unusedmets = cobra.manipulation.delete.prune_unused_metabolites(m)

## re-prepare the TFA model
m.prepare()
m.convert()#add_displacement = True)
## Info on the cobra_model
m.print_info()

tfa_solution_sc = m.optimize()
tfa_value_sc = tfa_solution_sc.objective_value
# Report
print('TFA Solution found for NAD-only metabolic model: {0:.5g}'.format(tfa_value_sc))



# perform variability analysis on reactions
mytfa.reactions.get_by_id(biomass_rxn).lower_bound = tfa_value - 1e-6
varib_results = variability_analysis(mytfa,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.wt.results.optimal.csv')
print("finished variability analysis on optimal solution for WT model")

# perform variability analysis on reactions
mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.1
varib_results = variability_analysis(mytfa,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.wt.results.gtzp1.csv')
print("finished variability analysis with growth rate = 0.1 for WT model")

# perform variability analysis on reactions
mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5
varib_results = variability_analysis(mytfa,kind='reactions')
print("finished variability analysis with growth rate = 0.5 for WT model")
varib_results.to_csv(out_dir + '/FVA.wt.results.gtzp5.csv')


# perform variability analysis on reactions
m.reactions.get_by_id(biomass_rxn).lower_bound = tfa_value_sc - 1e-6
varib_results = variability_analysis(m,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.sc.results.optimal.csv')
print("finished variability analysis on optimal solution for SC model")

# perform variability analysis on reactions
m.reactions.get_by_id(biomass_rxn).lower_bound = 0.1
varib_results = variability_analysis(m,kind='reactions')
varib_results.to_csv(out_dir + '/FVA.sc.results.gtzp1.csv')
print("finished variability analysis with growth rate = 0.1 for SC model")

# perform variability analysis on reactions
m.reactions.get_by_id(biomass_rxn).lower_bound = 0.5
varib_results = variability_analysis(m,kind='reactions')
print("finished variability analysis with growth rate = 0.5 for SC model")
varib_results.to_csv(out_dir + '/FVA.sc.results.gtzp5.csv')


