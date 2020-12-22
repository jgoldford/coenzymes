
mets.nad = find(strcmp(model.mets,'nad[c]'));
mets.nadh = find(strcmp(model.mets,'nadh[c]'));
mets.nadp = find(strcmp(model.mets,'nadp[c]'));
mets.nadph = find(strcmp(model.mets,'nadph[c]'));
nad_rxns = find((sum(model.S([mets.nad,mets.nadh],:)) == 0) & any(model.S([mets.nad,mets.nadh],:)));
nadp_rxns = find((sum(model.S([mets.nadp,mets.nadph],:)) == 0) & any(model.S([mets.nadp,mets.nadph],:)));

ntemp = setdiff(nad_rxns,nadp_rxns);
ptemp = setdiff(nadp_rxns,nad_rxns);

nad_rxns = ntemp;
nadp_rxns = ptemp;


l = {};

for i = 1:length(model.genes)

    gene = model.genes{i};
    if ~isempty(find(strcmp(model.genes,gene)))
        z = deleteModelGenes(model,gene);
        dels = setdiff(find((z.lb == 0) & (z.ub == 0)),find((model.lb == 0) & (model.ub == 0)));
        dels = intersect(dels,[nad_rxns,nadp_rxns]);
        l{end+1} = strjoin(model.rxns(dels),';');
    else
        l{end+1} = '';
    end
end
