#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import cobra
import pytfa
import numpy as np

from optlang.exceptions import SolverError

#from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_matlab_model, load_json_model


from pytfa.io import import_matlab_model, load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

import pickle
from pytfa.optim.variables import LogConcentration

thermo_database = '/projectnb2/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/data/thermo_data.thermodb'
out_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/tfa_sc_results/iJO1366/glucose_oxygen'
#root_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/tests/singleCoenzymeModel.08272021.v2'
#thermo_database = '../assets/thermo_data.thermodb'
model_path = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/models'

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = GUROBI

# Load reaction DB
print("Loading thermo data...")
thermo_data = load_thermoDB(thermo_database)

print("Done !")
#biomass_rxn = 'BIOMASS_Ec_iJO1366_WT_53p95M'
biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'


model_id = 'iJO1366_WT_semi-unconstrained.json'
#model_id = 'iJO1366.json'

# We import pre-compiled data as it is faster for bigger models
cobra_model = load_json_model(model_path + '/' + model_id)


from copy import deepcopy
def removeDuplicateRxn(model):
	model2 = deepcopy(model)
	toRemove = []
	for eachReaction in model2.reactions:

		if eachReaction.id in toRemove:
			continue

		ids = []
		stechiometry = []
		for eachMet in eachReaction.metabolites:
			ids.append(eachMet.id)
			stechiometry.append(eachReaction.metabolites[eachMet])

		removedEachReaction = 0
		for eachReaction2 in model2.reactions:
			if eachReaction2 == eachReaction or eachReaction2.id in toRemove or len(eachReaction2.metabolites) != len(eachReaction.metabolites):
				continue
			
			duplicate = 1
			for eachMet in eachReaction2.metabolites:
			
				if eachMet.id not in ids:
					duplicate = 0
					break
					
				if abs(eachReaction2.metabolites[eachMet]) != abs(stechiometry[ids.index(eachMet.id)]):
					duplicate = 0
					break
			
			if duplicate == 1:
				
				if len(eachReaction2.genes) == 0 and len(eachReaction.genes) == 0:
					if eachReaction2.lower_bound < eachReaction.lower_bound:
						toRemove.append(eachReaction.id)
						removedEachReaction = 1
					else:
						toRemove.append(eachReaction2.id)
				
				elif len(eachReaction2.genes) == 0 and len(eachReaction.genes) != 0:
					toRemove.append(eachReaction.id)
					removedEachReaction = 1
					
				elif len(eachReaction2.genes) != 0 and len(eachReaction.genes) == 0:
					toRemove.append(eachReaction2.id)
				
				elif len(eachReaction2.genes) != 0 and len(eachReaction.genes) != 0:
					
					totalGenes = []
					for eachGene in eachReaction2.genes:
						totalGenes.append(eachGene.id)
					for eachGene in eachReaction.genes:
						if eachGene.id not in totalGenes:
							totalGenes.append(eachGene.id)
				
					if eachReaction2.lower_bound < eachReaction.lower_bound:
						toRemove.append(eachReaction.id)
						if len(eachReaction2.genes) != len(totalGenes):
							string = ''
							for i in range(len(totalGenes)):
								string = string + totalGenes[i]
								if i != len(totalGenes) -1:
									string = string + " or "
							eachReaction2.gene_reaction_rule = '( ' + string +' )'
						removedEachReaction = 1
					else:
						toRemove.append(eachReaction2.id)
						if len(eachReaction.genes) != len(totalGenes):
							string = ''
							for i in range(len(totalGenes)):
								string = string + totalGenes[i]
								if i != len(totalGenes) -1:
									string = string + " or "
							eachReaction.gene_reaction_rule = '( ' + string +' )'
							
				
			if removedEachReaction == 1:
				break
		

	model2.remove_reactions(toRemove)
	
	return [model2, toRemove]

cobra_model = load_json_model(model_path + '/' + model_id)
exchange_rxns = [x.id for x in  cobra_model.reactions if x.id[0:3] == 'DM_']
lb = []
ub = []
for e in exchange_rxns:
    #cobra_model.reactions.get_by_id(e).lower_bound = 0
    lb.append(cobra_model.reactions.get_by_id(e).lower_bound )
    ub.append(cobra_model.reactions.get_by_id(e).upper_bound )
    
bounds = pd.DataFrame({'exchange': exchange_rxns, 'lb': lb, 'ub':ub})


# change bounds to glucose consuming, oxygen growth

carbon_source = 'DM_glc_e'
oxygen = True
carbon_flux = 8.16
cobra_model.reactions.get_by_id('DM_glc_e').lower_bound = 0

cobra_model.reactions.get_by_id(carbon_source).lower_bound = -carbon_flux
if oxygen:
    cobra_model.reactions.get_by_id('DM_o2_e').lower_bound = -1000
else:
    cobra_model.reactions.get_by_id('DM_o2_e').lower_bound = 0



# find solution
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

# fva = flux_variability_analysis(mytfa)

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

if tfa_value < 0.1:
    
    from pytfa.optim.relaxation import relax_dgo
    #mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5*fba_value
    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5
    relaxed_model, slack_model, relax_table = relax_dgo(mytfa)
    original_model, mytfa = mytfa, relaxed_model
    print('Relaxation: ')
    print(relax_table)
    
    tfa_solution = mytfa.optimize()
    tfa_value = tfa_solution.objective_value

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

# copy model to mutant model
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
print('TFA Solution found : {0:.5g}'.format(tfa_value_sc))

tfa_solution_wt = z.optimize()
tfa_value_wt = tfa_solution_wt.objective_value
# Report
print('TFA Solution found : {0:.5g}'.format(tfa_value_wt))

tfa_res = pd.DataFrame({'model': ['WT','SC'], 'max growth rate' : [tfa_value_wt,tfa_value_sc]})

tfa_res.to_csv(out_dir + '/tfa.wt.sc.growth-rates.csv')
tfa_solution_wt.fluxes.to_csv(out_dir + '/tfa.wt.fluxes.csv')
tfa_solution_sc.fluxes.to_csv(out_dir + '/tfa.sc.fluxes.csv')


