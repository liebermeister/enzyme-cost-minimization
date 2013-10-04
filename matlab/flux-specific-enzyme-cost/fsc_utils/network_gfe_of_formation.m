function [G0,GO_sbtab] = network_gfe_of_formation(network)

% [G0,GO_sbtab] = network_gfe_of_formation(network)
% 
% network must contain the field metabolite_KEGGID with valid KEGG IDs

if ~exist('parseKeggModel','file'),
  error('Please install matlab code for component contribution method');
end

kegg_formulae = network_print_formulae(network,network.actions,network.metabolite_KEGGID,1);


% -----------------------------------------------
% compute GFE of formation (call component contribution method)

[S, cids] = parseKeggModel(kegg_formulae);

% Use the standard parameters

pH = 7.5;
I  = 0.2;
T  = 300; 

[DfG0_prime, covf] = getGibbsForKeggModel(S, cids, pH, I, T);


% -----------------------------------------------
% map GFE of formation back onto model

KEGG_ids = {};

for it = 1:length(cids),
  KEGG_ids{it,1} = ['C' repmat('0',1,5-length(num2str(cids(it)))), num2str(cids(it))];
end

ll = label_names(network.metabolite_KEGGID,KEGG_ids);
G0 = DfG0_prime(ll);


% -----------------------------------------------
% output in SBtab format (THIS CAN LEAD TO PROBLEMS)

if nargout>1,

  nm = length(network.metabolites);
  clear s
  s.QuantityType = repmat({'standard chemical potential'},nm,1);
  s.Unit         = repmat({'kJ/mol'},nm,1);
  s.Compound_MiriamID__urn_miriam_kegg_compound = network.metabolite_KEGGID;
  s.SBMLSpeciesID = network.metabolite_KEGGID;
  if isfield(network,'metabolite_names'),
    s.CompoundName = network.metabolite_names;
  else,
    s.CompoundName = network.metabolites;
  end
  s.Value        = G0;
  s.Reference    = repmat({'Component contribution method'},nm,1);
  G0_sbtab = sbtab_table_construct_from_struct(struct,s,{'QuantityType','Unit','Compound MiriamID::urn.miriam.kegg:compound','SBMLSpeciesID','CompoundName','Value','Reference'}');
  
end
