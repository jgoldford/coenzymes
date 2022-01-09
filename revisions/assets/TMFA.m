function [results] = TMFA(irr_model,params)

mets_lb = params.met.lb;
mets_ub = params.met.ub;
err_bound = params.deltag_err_bound;
T = params.T;
global RT K faraday dPsi dPH
RT = 0.008309424 * (273.15+T);
K = 1e8;
% farady in kJ/mol mV
faraday = 0.09648304;
dPsi = -130;
dPH = 0.4;

% load model with KEGG IDs for compounds
%load('ecoli_core_keggids.mat', 'model');

% make irreversible model
S = num2cell(irr_model.S,1);

% collapse all data in G into a simple table for convience
% the rxn index corresponds to which irreversible reaction we are modeling
% in the model - corresponds to column in S
g = irr_model.deltaG;

free_energy_struct.rxn_index = find(~isnan(g));
% delta G for each reaction
free_energy_struct.deltaG = g(~isnan(g));
free_energy_struct.deltaGerr = irr_model.deltaGerr(~isnan(g));
free_energy_struct.rxn = irr_model.rxns(free_energy_struct.rxn_index);
[unique_rxns,ri,rj]= unique(cellfun(@(x) x(1:end-2),free_energy_struct.rxn,'uni',0));
free_energy_struct.urxns = unique_rxns;
free_energy_struct.urxns_ri = ri;
free_energy_struct.urxns_rj = rj;
direction = cellfun(@(x) x(end), free_energy_struct.rxn,'uni',0);
ds = [];
for i = 1:length(direction)
    if direction{i} == 'f'
        ds(i) = 1;
    else
        ds(i) = -1;
    end
end
free_energy_struct.rxn_dir = ds';
    

% error in delta G

% direction (+1 is forward and -1 is backwards)
% index for which error varaiable should be appended to the TMFA
% formulation - i.e. error in  free energy estimates for forward and
% backward reactions are coupled.

% specificty variable dimensions
dims.v = size(S,2);
dims.m = length(irr_model.b);
dims.z = length(free_energy_struct.rxn_index);
%dims.e = length(free_energy_struct.rxn_index);
dims.e = max(free_energy_struct.urxns_rj);
dims.x = dims.v + dims.z + dims.e + dims.m;

f = @(LETTER, DIM) arrayfun(@(x) LETTER,1:DIM);
varIdx = @(x) find(strcmp(num2cell([f('v',dims.v),f('z',dims.z),f('m',dims.m),f('e',dims.e)]),x));

% set variable type
vartype{1} = arrayfun(@(x) 'C', 1:dims.v,'uni',1);
vartype{2} = arrayfun(@(x) 'B', 1:dims.z,'uni',1);
vartype{3} = arrayfun(@(x) 'C', 1:dims.m,'uni',1);
vartype{4} = arrayfun(@(x) 'C', 1:dims.e,'uni',1);


%% construct constraint 1: mass balance constraints
C1 = [cell2mat(S), zeros(dims.m,dims.x-dims.v)];
B1 = zeros(dims.m,1);
csense1 = arrayfun(@(x) 'E', 1:dims.m,'uni',1);

%% construct constraint 2: reaction feasibility constraints
D = -diag(irr_model.ub(free_energy_struct.rxn_index));
U = eye(dims.v);
U = U(free_energy_struct.rxn_index,:);

C2 = [U,D,zeros(dims.z,dims.m+dims.e)];
%C2 = C2(free_energy_struct.rxn_index,:);
B2 = zeros(dims.z,1);
csense2 = arrayfun(@(x) 'L', 1:dims.z,'uni',1);


%% construct constraint 3: thermodynamic constraints feasibility constraints

% set up thermodynamic constraints
%I = eye(dims.v);
%I = I(free_energy_struct.rxn_index,:);
s = cell2mat(S);


%E = diag(free_energy_struct.deltaGerr);
% construct error matrix +1 for forward reaction -1 for reverse
E = zeros(dims.z,dims.e);
for zi =1:dims.z
        rxn = cellfun(@(x) x{1} ,splitString(free_energy_struct.rxn(zi),'_'),'uni',0);
        ei = find(strcmp(free_energy_struct.urxns,rxn));
        if free_energy_struct.rxn_dir(zi) > 0
            E(zi,ei) = free_energy_struct.deltaGerr(zi);
        elseif free_energy_struct.rxn_dir(zi) < 0
            E(zi,ei) = -free_energy_struct.deltaGerr(zi);
        end
end

C3 = [zeros(dims.z,dims.v),K.*eye(dims.z), RT.*s(:,free_energy_struct.rxn_index)',E];
B3 = [K - free_energy_struct.deltaG];
csense3 = arrayfun(@(x) 'L', 1:dims.z,'uni',1);


%intracelluarl metabolite range is from 0.100 M to 1nM
mets.ub = log(mets_ub.*ones(length(irr_model.mets),1));
%mets.lb = log(10^(-9).*ones(length(irr_model.mets),1));
mets.lb = log(mets_lb.*ones(length(irr_model.mets),1));

% set bounds of intracellular and extracellular metabolites
h2o = find(strcmp(irr_model.mets,'C00001'));
mets.lb(h2o) = log(1);
mets.ub(h2o) = log(1);

hplus = find(strcmp(irr_model.mets,'C00080'));
mets.lb(hplus) = log(1);
mets.ub(hplus) = log(1);

% construct upper and lower bounds for all variables (fluxes, reaction
% feasibility, metabolites and error function

ub = [irr_model.ub;ones(dims.z,1);mets.ub;err_bound.*ones(dims.e,1)];
lb = [irr_model.lb;zeros(dims.z,1);mets.lb;-err_bound.*ones(dims.e,1)];     


MILP.A = [C1;C2;C3];
MILP.b = [B1;B2;B3];
MILP.csense = [csense1,csense2,csense3];
MILP.vartype = [vartype{:}];
MILP.lb = lb;
MILP.ub = ub;
MILP.c = zeros(dims.x, 1); MILP.c(~~(irr_model.c)) = 1;
MILP.osense = -1;

x0 = zeros(dims.x,1);
x0(varIdx('v')) = zeros(dims.v,1);
x0(varIdx('m')) = mets.ub;
x0(varIdx('z')) = 1;
MILP.x0 = x0; 

x = solveCobraLP(MILP);

MILP.x0 = x.full; 

out = solveCobraMILP(MILP);

if strcmp(out.stat,'INFEASIBLE')
    results=[];
else
    
       % solve the MILP that minimizes the norm of fluxes while
       % constraining the objective value
       epsilon = 1e-3.*out.obj;
       c = find(MILP.c);
       MILP.lb(c) = out.obj - epsilon;
       MILP.ub(c) = out.obj + epsilon;
       MILP.x0 = out.full;
       MILP.osense = +1;
       MILP.c(varIdx('v'))=1;
       out2 = solveCobraMILP(MILP);
    
    
results = out;
results.blocked_reactions = irr_model.rxns(free_energy_struct.rxn_index(out.full(varIdx('z')) ==0));
results.fluxes = out2.full(varIdx('v'));
results.concentrations = out.full(varIdx('m'));

end

end

    



