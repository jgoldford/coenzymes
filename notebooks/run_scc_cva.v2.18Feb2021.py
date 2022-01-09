from cobra.io import load_matlab_model

from pytfa.io import load_thermoDB,                    \
                        read_lexicon, annotate_from_lexicon,            \
                        read_compartment_data, apply_compartment_data
import pytfa
import cobra
import pandas as pd
import numpy as np
from cobra.manipulation.delete import prune_unused_metabolites

solver = 'optlang-gurobi'
root_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/assets/pytfba/'
#root_dir = '../../pytfa/'
fileOut = root_dir + 'results.RVA_ANG.v2.18Feb2021.csv'
fileOut_temp = root_dir + 'temp.results.RVA_ANG.v2.18Feb2021.csv'
# Load the model
cobra_model = cobra.io.load_matlab_model(root_dir + 'models/small_ecoli.mat')

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
    
    if len(v) > 0:
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
        rxn.id = rxn.id + '[condensed]'
        model.remove_reactions([x for x in model.reactions if x.id == reaction_id][0])
        model.add_reaction(rxn)
        
    return model

m  = cobra_model.copy()
rxn_ids = [x.id for x in cobra_model.reactions]
m  = cobra_model.copy()
rxn_ids = [x.id for x in cobra_model.reactions]
for rxnid in rxn_ids:
    m = single_coenzyme_transform(m,rxnid)
    
m.remove_reactions([x for x in m.reactions if x.id == 'NADTRHD[condensed]'][0])
m.objective = 'Ec_biomass_iJO1366_WT_53p95M[condensed]'


# Load reaction DB
thermo_data = load_thermoDB(root_dir + 'data/thermo_data.thermodb')
lexicon = read_lexicon(root_dir + 'models/small_ecoli/lexicon.csv')
compartment_data = read_compartment_data(root_dir + 'models/small_ecoli/compartment_data.json')

def convert2thermo(model,name):
    # Initialize the model
    tmodel = pytfa.ThermoModel(thermo_data,model)
    tmodel.name = name

    # Annotate the model
    annotate_from_lexicon(tmodel, lexicon)
    apply_compartment_data(tmodel, compartment_data)

    ## TFA conversion
    tmodel.prepare()
    tmodel.convert()
    tmodel.solver = solver
    return tmodel

tmodel_wt = convert2thermo(cobra_model,'wt')
solution_wt = tmodel_wt.optimize()

tmodel_mut = convert2thermo(m,'mut')
solution_mut = tmodel_mut.optimize()


# do this for 250 metabolites in model

results = {'model': [], 'min/max': [], 'met': [], 'value': []}
mets = [x.id for x in tmodel_wt.LC_vars]


# perform metabolite variability analysis

m_i = cobra_model.copy()
lb_growth_rate = solution_wt.objective_value
m_i.reactions.get_by_id('Ec_biomass_iJO1366_WT_53p95M').lower_bound = lb_growth_rate
tmodel_obj = convert2thermo(m_i,'wt_model_lb_gr')



glutamate = [x for x in tmodel_obj.LC_vars if 'glu-L_c' in x.id][0]
nh4 = [x for x in tmodel_obj.LC_vars if 'nh4_c' in x.id][0]
akg = [x for x in tmodel_obj.LC_vars if 'akg_c' in x.id][0]
exp = tmodel_obj.LC_vars[akg] + tmodel_obj.LC_vars[nh4] - tmodel_obj.LC_vars[glutamate]
tmodel_obj.objective = exp
tmodel_obj.objective.direction = 'max'
max_met = tmodel_obj.slim_optimize() 
results['model'].append('two-coenzyme')
results['min/max'].append('max')
results['met'].append('akg+nh3-glt')
results['value'].append(max_met)

tmodel_obj.objective.direction = 'min'
min_met = tmodel_obj.slim_optimize()
results['model'].append('two-coenzyme')
results['min/max'].append('min')
results['met'].append('akg+nh3-glt')
results['value'].append(min_met)
pd.DataFrame(results).to_csv(fileOut_temp)

m_i = m.copy()
lb_growth_rate = solution_mut.objective_value
m_i.reactions.get_by_id('Ec_biomass_iJO1366_WT_53p95M[condensed]').lower_bound = lb_growth_rate
tmodel_obj = convert2thermo(m_i,'mut_model_lb_gr')
glutamate = [x for x in tmodel_obj.LC_vars if 'glu-L_c' in x.id][0]
nh4 = [x for x in tmodel_obj.LC_vars if 'nh4_c' in x.id][0]
akg = [x for x in tmodel_obj.LC_vars if 'akg_c' in x.id][0]
exp = tmodel_obj.LC_vars[akg] + tmodel_obj.LC_vars[nh4] - tmodel_obj.LC_vars[glutamate]
tmodel_obj.objective = exp
tmodel_obj.objective.direction = 'max'
max_met = tmodel_obj.slim_optimize() 
results['model'].append('one-coenzyme')
results['min/max'].append('max')
results['met'].append('akg+nh3-glt')
results['value'].append(max_met)

tmodel_obj.objective.direction = 'min'
min_met = tmodel_obj.slim_optimize()
results['model'].append('one-coenzyme')
results['min/max'].append('min')
results['met'].append('akg+nh3-glt')
results['value'].append(min_met)
pd.DataFrame(results).to_csv(fileOut_temp)
                            
results = pd.DataFrame(results)
results.to_csv(fileOut)