
% initialize cobra toolbox for LP solver
initCobraToolbox
changeCobraSolver('glpk')

clear
% This script models Scenario 1, where reaction directionality is not dictated by coenzyme choise.
% firt sample effective equilibrium constants, and catalytic
% rates

% randomly sample net reaction directionality and coenzyme stoichiometry
% metabolism sets the demands of these sources and sink fluxes, and are
% typically constratined by the environment and the definition of biomass
% equation

% for the stoichiometry and K values, total flux has to sum to v_r.  
% define number of reactions and coenzymes

% define number of simulations
numReals = 1;
mu_k = 6;
var_k = 1;
min_ratio = nan(numReals,1);
max_ratio = nan(numReals,1);
min_ratio = nan(numReals,1);
minRelativeCosts = nan(numReals,1);
params = {};
fluxes = {};
for m =1:numReals
    
    N = 100;
    C = 2;

    % randomly sample the stoichiometry:
   
    s = 2*binornd(1,0.5,[N,1]) - 1;

    % sample the flux distribution
    TotalFlux = 100;
    N_b = sum(s<0);
    N_f = sum(s>0);
    amp = 1;
    v_f = drchrnd(amp*ones(N_f,1)',1) * TotalFlux;
    v_b = drchrnd(amp*ones(N_b,1)',1) * TotalFlux;
    v_r = zeros(N,1);
    v_r(s<0) = v_b;
    v_r(s>0) = v_f;

    kappa = lognrnd(1,5,[N,1]);
    kappa_e = lognrnd(1,5,1);

    % sample the effective thermodynamic drive with
    
    mu_kp = mu_k; var_kp = var_k;
    mu_kc = mu_k; var_kc = var_k;

    N_p = sum(s>0);
    N_c = sum(s<0);

    cond = true;
    while cond
        
        Kp = lognrnd(mu_kp,var_kp,[N_p,1]);
        Kc = lognrnd(mu_kc,var_kc,[N_c,1]);
           
        K = zeros(N,1);

        K(s >0) = Kp;
        K(s <0) = Kc;
        %$K = lognrnd(mu_k,var_k,[N,1]);
        %K(s<0) = 1./K(s<0);

        minGamma = max(1./Kc);
        %minGamma = max(K(s<0));
        %minGamma = max(1./K(s<0));
        maxGamma = min(Kp);
        %maxGamma = min(K(s>0));
        
        if minGamma > maxGamma
            cond = true;
        else
            
            weight = @(x) (1./ kappa) .* (K + x.^s) ./ (K - x.^s);
            gammas = linspace(minGamma,maxGamma,1000);
            gammas = gammas(2:end-1);
            costs = arrayfun(@(x) weight(x)' * v_r,gammas);
            minCost = min(costs);
            if minCost < 0
                cond = true;
            else
                cond = false;
            end
        
        end
    end
    
    weight = @(x) (1./ kappa) .* (K + x.^s) ./ (K - x.^s);


    %gammas = logspace(minGamma,maxGamma,100);
    gammas = 10.^linspace(log10(minGamma),log10(maxGamma),100);
    Gas = [];
    Gbs = [];
    Feasible = [];
    F = nan(length(gammas));
    T = zeros(length(gammas));
    Fr = cell(length(gammas));
    Entropy = zeros(length(gammas));

    a = 1;
    
    for i = 1:length(gammas)
        for j = (i+1):length(gammas)

            G_a = gammas(i);
            G_b = gammas(j);
                clear LP
                gradients = [G_a,G_b];
                LP  = create_lp_problem(N,C,v_r,s,K,kappa,kappa_e,gradients);
                

                y = solveCobraLP(LP);

                if strcmp('INFEASIBLE',y.origStat) || isempty(y.obj)
                    F(i,j) = nan;
                    Fr{i,j} = [nan,nan];
                    Entropy(i,j) = nan;
                    fluxes{i,j}= [];
                else
                    F(i,j) = y.obj;
                    t = reshape(y.full(1:end-1),[N,C]);
                    fluxes{i,j} = reshape(y.full(1:end-1),[N,C]);
                    fr = sum(t);
                    Fr{i,j} = fr;
                    psuedo = 1e-9;
                    T(i,j) = y.full(end);
                    t = cellfun(@(x) (x + psuedo) ./ sum((x + psuedo)), num2cell(t,2),'uni',0);
                    Entropy(i,j) = mean(cellfun(@(x) -sum(x.*log(x)),t));
                end

            

        end
    end
    
    
    
    clear param
    param.mu_k = mu_k;
    param.var_k = var_k;
    param.K = K;
    param.v_r = v_r;
    param.s = s;
    param.N = N;
    param.C = C;
    param.kappa_r = kappa;
    param.kappa_e = kappa_e;
    param.singleCoenzymeCost = minCost;
    param.minGamma_single = minGamma;
    param.maxGamma_single = maxGamma;
    %parmas{m} = param;
    
    g = min(min(F));    
    param.doubleCoenzymeCost = g;
    param.relativeCostReduction = 1 - g / minCost;
    [ii,jj] = find(F == g);
    ga_min = gammas(ii);
    ga_max = gammas(jj);
    ss = sort([ga_min,ga_max]);
    
    param.gamma_alpha_double = ss(1);
    param.gamma_beta_double = ss(2);
    param.fluxes_opt = fluxes{ii,jj};
    param.cost_gradients = F;
    param.gradients = gammas;
    param.entropy = Entropy;
    %minRelativeCosts(m) = costRelative;
    params{m} = param;
    
    disp(['finished simulation: ' , int2str(m)])
end

