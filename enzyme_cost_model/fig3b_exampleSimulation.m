clear
load('sims_gridSearch_August10_muk_6_vark_1.mat');
param = params{3};


%gammas = logspace(minGamma,maxGamma,100);
gammas = logspace(-5,5,100);
Gas = [];
Gbs = [];
Feasible = [];
F = nan(length(gammas));
T = zeros(length(gammas));
Fr = cell(length(gammas));
Entropy = zeros(length(gammas));

a = 1;
N = param.N;
C = param.C;
F2 = nan(length(gammas));
F2Single = nan(length(gammas));

for i = 1:length(gammas)
    for j = (i+1):length(gammas)

        G_a = gammas(i);
        G_b = gammas(j);
            clear LP
            gradients = [G_a,G_b];
            LP  = create_lp_problem(param.N,param.C,param.v_r,param.s,param.K,param.kappa_r,param.kappa_e,gradients);
            y = solveCobraLP(LP);

            if strcmp('INFEASIBLE',y.origStat) || isempty(y.obj)
                F(i,j) = nan;
                F2(i,j) = nan;
                Fr{i,j} = [nan,nan];
                Entropy(i,j) = nan;
                fluxes{i,j}= [];
            else
                F(i,j) = y.obj/param.singleCoenzymeCost;
                F2(i,j) = 1;
                t = reshape(y.full(1:end-1),[N,C]);
                fluxes{i,j} = reshape(y.full(1:end-1),[N,C]);
                fr = sum(t);
                Fr{i,j} = fr;
                psuedo = 1e-9;
                T(i,j) = y.full(end);
                t = cellfun(@(x) (x + psuedo) ./ sum((x + psuedo)), num2cell(t,2),'uni',0);
                Entropy(i,j) = mean(cellfun(@(x) -sum(x.*log(x)),t));
            end

            if ((G_a > param.minGamma_single) && (G_a < param.maxGamma_single)) || ((G_b > param.minGamma_single) && (G_b < param.maxGamma_single))
                F2Single(i,j) = 1;
            else
                F2Single(i,j) = nan;
            
            end
    end
end

figure()
hold on
[nr,nc] = size(F);
%plot figure
r = [log10(gammas),nan];
pcolor(r,r,[F2' nan(nr,1); nan(1,nc+1)]);
%pcolor(r,r,[z nan(nr,1); nan(1,nc+1)]);
%[X,Y] = meshgrid(output(3).gradients,output(3).gradients);
%v = [1,1];
%Z = output(k).cost_gradients/output(k).singleCoenzymeCost;
%contour(X,Y,Z)

shading flat;
%set(gca, 'ydir', 'reverse');
%cmap = cbrewer('div','Spectral',100);
%cmap = flipud(cmap);
%colormap(cmap);
%caxis([0,1]);
xlabel('log (\Gamma_\alpha)')
ylabel('log (\Gamma_\beta)')
title('Enzyme cost')
xlim([min(log10(gammas)) max(log10(gammas))])
ylim([min(log10(gammas)) max(log10(gammas))])
%colorbar()

pcolor(r,r,[F2Single' nan(nr,1); nan(1,nc+1)]);

shading flat;

g = min(min(F));    
minCost = param.singleCoenzymeCost;
[ii,jj] = find(F == g);
x = log10(gammas);
hold on; 
scatter(x(ii),x(jj),'s')

