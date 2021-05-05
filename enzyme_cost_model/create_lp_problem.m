function [LPproblem] = create_lp_problem(N,C,v_r,s,K,kappa,kappa_e,coenzyme_ratios)
% function to generate the linear program for the enzyme cost model
% returns a LPproblem stucture that can be used directly by the Cobra Toolbox

% create stoichiometric matrix for transhydrognase reactions
thd_idx = nchoosek(1:C,2);
Ne = size(thd_idx,1);

s_thd = zeros(Ne,C);

for i =1:Ne
    alpha = min(thd_idx(i,:));
    beta = max(thd_idx(i,:));
    s_thd(i,alpha) = 1;
    s_thd(i,beta) = -1;
end



A0 = repmat(diag(ones(N,1)),[1,C]);

% construct the constraint matrix, each row corresponds to a new
% stoichiometric constraint for each coenzyme in C.  The additional stoichiometric
% matrix appendeded to the end ecodes the transhydrogenase reactions variables.

A1 = zeros(C,Ne+(N*C));
for c =1:C
    A1(c,(c-1) * N +  1 : c * N) = s';
end
A1(:,(C*N+1):(Ne + N*C)) = s_thd';

% now construct with 
LPproblem.A = [A0,zeros(N,Ne);[A1]];
LPproblem.b = [v_r;zeros(C,1)];

% define weight funcrtion
% compute wights for coenzyme ratios

w_e = cellfun(@(s,ke) 1./ke * (1+prod(coenzyme_ratios.^s))./(1-prod(coenzyme_ratios.^s)), num2cell(s_thd,2),num2cell(kappa_e),'uni',1);
weight_r = @(x) (1./ kappa) .* (K + x.^s) ./ (K - x.^s);

weights = [];

for c = 1:C
    weights = [weights;weight_r(coenzyme_ratios(c))];
end

w = [weights;w_e];

LPproblem.c = w;

direct = 0;
for j =1:length(w)
    if w(i) == Inf;
       direct(j) = 0;
    else
       direct(j) = sign(w(j));
    end
end


LPproblem.lb = zeros(C*N+Ne,1);
LPproblem.lb(direct < 0 ) = -10000;
LPproblem.ub = zeros(C*N+Ne,1);
LPproblem.ub(direct > 0 ) = 10000;
            
LPproblem.csense = [arrayfun(@(x) 'E',LPproblem.b)];
LPproblem.osense = 1;
end