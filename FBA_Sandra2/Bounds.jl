function Bounds(DF,TXTL)

  FB = DF["default_flux_bounds_array"]

  FB[258,2] = DF["Oxygen"]; #[]-> O2
  FB[270,2] = DF["GlcUptake"]; #[]-> GLC
  FB[274,2] = 0; #[]-> PYR

  #uncomment the following line if running "Ensemble.jl"
  #FB[251,2] = DF["HeroinUptake"]

  org = collect(277:1:281)
  for i=1:length(org)
    FB[org[i],2] = 0; #succ, Mal, fum, etoh, mglx
  end

  abc = collect(287:2:325)
  for i=1:length(abc)
      FB[abc[i],2] = DF["AAUptake"]
  end

  turt = collect(288:2:326)
  for i=1:length(turt)
      FB[turt[i],2] = DF["AASecretion"]
  end

  FB[252,2] = 0; #no heroin secretion

  #FB[287:2:325,2] = DF["AAUptake"]; #AA UPTAKE
  #FB[288:2:326,2] = DF["AASecretion"]; #AA ->[]
  #FB[264,2] = 0; #ATP ->[]
  #FB[265,2] = 0; #[]-> ADP
  #AASyn = TXTL["AA_synthesis_rxn"]
  #AADeg = TXTL["AA_degradation_rxn"]
  #FB[AASyn,2] = DF["AASyn"];
  #FB[AADeg,2] = 0;

#==============================================TXTL=====================================================#
  RNAP_concentration_nM = TXTL["RNAP_concentration_nM"];
  RNAP_elongation_rate = TXTL["RNAP_elongation_rate"];
  RIBOSOME_concentration = TXTL["RIBOSOME_concentration"];
  RIBOSOME_elongation_rate = TXTL["RIBOSOME_elongation_rate"];
  kd = TXTL["mRNA_degradation_rate"];
  mRNA_length_1 = TXTL["mRNA_length_1"];
  protein_length_1 = TXTL["protein_length_1"];
  mRNA_length_2 = TXTL["mRNA_length_2"];
  protein_length_2 = TXTL["protein_length_2"];
  mRNA_length_3 = TXTL["mRNA_length_3"];
  protein_length_3 = TXTL["protein_length_3"];
  gene_copies = TXTL["gene_copies"];
  volume = TXTL["volume"];
  polysome_amplification = TXTL["polysome_gain"];
  plasmid_saturation_coefficient = TXTL["plasmid_saturation_coefficient"];
  mRNA_saturation_coefficient = TXTL["mRNA_saturation_coefficient"];
  Promoter = TXTL["Promoter"]
  inducer = TXTL["inducer"]
  protein_degradation_rate = TXTL["protein_degradation_rate"]
#====================================Transcription===================================================#
  #Compute the promoter strength P -
  n = Promoter[1]
  KD = Promoter[2]
  K1 = Promoter[3]
  K2 = Promoter[4]
  K1_T7 = Promoter[5]
  f = inducer^n/(KD^n+inducer^n)
  if Promoter_model == 1
    P = (K1_T7)/(1+K1_T7)
  elseif Promoter_model == 2
    P = (K1+K2*f)/(1+K1+K2*f);
  end

  gene_concentration = gene_copies*(1e9/6.02e23)*(1/volume);
  saturation_term = (gene_concentration)/(plasmid_saturation_coefficient+gene_concentration);
  RNAP_concentration = RNAP_concentration_nM/1e6; #nM to mM
  TX_1 = (RNAP_elongation_rate*(1/mRNA_length_1)*(RNAP_concentration)*(saturation_term)*3600)*P;
  TX_2 = (RNAP_elongation_rate*(1/mRNA_length_2)*(RNAP_concentration)*(saturation_term)*3600)*P;
  TX_3 = (RNAP_elongation_rate*(1/mRNA_length_3)*(RNAP_concentration)*(saturation_term)*3600)*P;

#====================================Translation===================================================#
  mRNA_steady_state_1 = (TX_1/kd);
  translation_rate_constant_1 = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length_1)*3600;
  TL_1 = translation_rate_constant_1*RIBOSOME_concentration*mRNA_steady_state_1/(mRNA_saturation_coefficient+mRNA_steady_state_1);

  mRNA_steady_state_2 = (TX_2/kd);
  translation_rate_constant_2 = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length_2)*3600;
  TL_2 = translation_rate_constant_2*RIBOSOME_concentration*mRNA_steady_state_2/(mRNA_saturation_coefficient+mRNA_steady_state_2);

  mRNA_steady_state_3 = (TX_3/kd);
  translation_rate_constant_3 = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length_3)*3600;
  TL_3 = translation_rate_constant_3*RIBOSOME_concentration*mRNA_steady_state_3/(mRNA_saturation_coefficient+mRNA_steady_state_3);

#==============================Steady State Protein Concentrations==================================#
Prot_1 = TL_1/protein_degradation_rate;
Prot_2 = TL_2/protein_degradation_rate;
Prot_3 = TL_3/protein_degradation_rate;

#==========================Vmax for Biosensor Reactions===========================================#
her_spec_act = 8.9      #specific activity - units of umol/min/mg
her_MW = 34240              #g/mol

morA_spec_act = 230
morA_MW = 32124

kcat_her = her_spec_act*her_MW*(1/60)*10^-3         #converting to kcat units = s^-1
kcat_morA = morA_spec_act*morA_MW*(1/60)*10^-3
kcat_fmn = 1.5
kcat_lux = 0.1

her_Vmax = kcat_her*Prot_1
morA_Vmax = kcat_morA*Prot_2
fmn_Vmax = kcat_fmn*1               #enzyme concentration unknown, assumed to be 1mM
lux_Vmax = kcat_lux*Prot_3

#update bounds of biosensor reactions
FB[242,2] = her_Vmax;
FB[243,2] = morA_Vmax;
FB[244,2] = fmn_Vmax
FB[245,2] = lux_Vmax;

#===================================================================================================#

  FB[167,1] = TX_1 #transcriptional initiation
  FB[167,2] = TX_1 #transcriptional initiation
  #FB[169,1] = TX_1 #transcriptional initiation
  #FB[169,2] = TX_1 #mRNA_degradation
  FB[170,1] = TX_2
  FB[170,2] = TX_2 #translations initiation
  #FB[172,1] = TX_2
  #FB[172,2] = TX_2
  FB[173,1] = TX_3
  FB[173,2] = TX_3
  #FB[175,1] = TX_3
  #FB[175,2] = TX_3
  FB[176,1] = TL_1
  FB[198,1] = TL_2
  FB[220,1] = TL_3
  FB[176,2] = TL_1
  FB[198,2] = TL_2
  FB[220,2] = TL_3

  DF["default_flux_bounds_array"] = FB
  return DF

end
