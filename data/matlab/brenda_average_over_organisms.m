M = sbtab_table_load('/home/wolfram/projekte/flux_specific_enzyme_cost/flux-specific-enzyme-cost/data/data-kinetic/kinetic_constants_brenda.tsv');

QuantityType  = sbtab_table_get_column(M,'QuantityType');
Value         = sbtab_table_get_column(M,'Value',1);
Unit          = sbtab_table_get_column(M,'Unit');
KEGG_Reaction = sbtab_table_get_column(M,'Reaction MiriamID::urn:miriam:kegg.reaction');
KEGG_Compound = sbtab_table_get_column(M,'Compound MiriamID::urn:miriam:kegg.compound');
KEGG_Combined = cellstr([char(KEGG_Reaction) char(KEGG_Compound)]);

parameter_types = {'substrate catalytic rate constant','product catalytic rate constant','Michaelis constant'};

QuantityTypeList  = {};
KEGG_CompoundList = {};
KEGG_ReactionList = {};
Value_List        = [];
Unit_List         = {};

for it = 1:length(parameter_types),

  my_QuantityTypeList  = {};
  my_KEGG_CompoundList = {};
  my_KEGG_ReactionList = {};
  my_Value_List        = [];
  my_Unit_List         = {};
  
  my_quantity_type = parameter_types{it};
  ind = find(strcmp(parameter_types{it},QuantityType));
  combinations = unique(KEGG_Combined(ind));
  for itt = 1:length(combinations),
    my_combination = combinations{itt}; 
    my_ind = ind(find(strcmp(my_combination, KEGG_Combined(ind))));
    my_QuantityTypeList  = [my_QuantityTypeList; {my_quantity_type}];
    my_KEGG_CompoundList = [my_KEGG_CompoundList; KEGG_Compound(my_ind(1))];
    my_KEGG_ReactionList = [my_KEGG_ReactionList; KEGG_Reaction(my_ind(1))];
    my_Value_List        = [my_Value_List; exp(mean(log(Value(my_ind))))];
    my_Unit_List         = [my_Unit_List; Unit(my_ind(1))];
  end

  QuantityTypeList  = [QuantityTypeList ; my_QuantityTypeList  ];
  KEGG_CompoundList = [KEGG_CompoundList; my_KEGG_CompoundList ];
  KEGG_ReactionList = [KEGG_ReactionList; my_KEGG_ReactionList ];
  Value_List        = [Value_List       ; my_Value_List        ];
  Unit_List         = [Unit_List        ; my_Unit_List         ];

end

clear ts
ts.QuantityType = QuantityTypeList;
ts.Reaction_MiriamID__urn_miriam_kegg_reaction =  KEGG_ReactionList;
ts.Compound_MiriamID__urn_miriam_kegg_compound =  KEGG_CompoundList;
ts.Value = Value_List;
ts.Unit  = Unit_List;

my_sbtab_table = sbtab_table_construct_from_struct({},ts, {'QuantityType','Reaction MiriamID::urn:miriam:kegg.reaction','Compound MiriamID::urn:miriam:kegg.compound','Value','Unit'});

sbtab_table_save(my_sbtab_table,'/home/wolfram/projekte/flux_specific_enzyme_cost/flux-specific-enzyme-cost/data/data-kinetic/average_kinetic_constants_brenda.tsv');