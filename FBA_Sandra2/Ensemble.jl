include("Include.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")
include("DataDictionary.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])

#Set objective reaction
data_dictionary["objective_coefficient_array"][242] = -1;
data_dictionary["objective_coefficient_array"][243] = -1;
data_dictionary["objective_coefficient_array"][244] = -1;
data_dictionary["objective_coefficient_array"][245] = -1;

#=============================Cases=========================================#
#Define case number
# 1 = Amino Acid Uptake & Synthesis
# 2 = Amino Acid Uptake w/o Synthesis
# 3 = Amino Acid Synthesis w/o Uptake
Case = 1
if Case == 1
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if Case == 2
  data_dictionary["AASyn"] = 0;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if Case == 3
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 0
  data_dictionary["AASecretion"] = 100;
end

#Set Promoter
#1 = T7 Promoter model
#2 = P70a Promoter model
Promoter_model = 1
#===========================================================================#

#Set Plasmid Dose (nM)
plasmid_concentration = 5;
volume = TXTL_dictionary["volume"]
gene_copy_number = (volume/1e9)*(6.02e23)*plasmid_concentration;
TXTL_dictionary["gene_copies"] = gene_copy_number

#Set Glucose and Oxygen (mM/h)
data_dictionary["GlcUptake"] = 30
data_dictionary["Oxygen"] = 100;


runs = 100
index = 1
flux_ensemble = zeros(number_of_fluxes,runs)
#uptake_ensemble = zeros(number_of_species,runs)
while index <= runs
  data_dictionary["HeroinUptake"] = 0 + 10*rand(1)[1];
  # solve the lp problem -
  turkey = Bounds(data_dictionary,TXTL_dictionary);
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(turkey)
  if exit_flag == 5
    flux_ensemble[:,index] = flux_array
    #uptake_ensemble[:,index] = uptake_array
    global index += 1;
  end
end

x = flux_ensemble[251,:]
y = flux_ensemble[242,:]
z = flux_ensemble[243,:]
a = flux_ensemble[244,:]
b = flux_ensemble[245,:]

her_up_max = maximum(x)
her_deg_max = maximum(y)
mor_deg_max = maximum(z)
fmn_max = maximum(a)
lux_max = maximum(b)

using PyPlot
plot(x,y,label="Heroin Degradation",linewidth=3,linestyle="-")
plot(x,z,label="Morphine Degradation",linewidth=1,linestyle=":")
plot(x,a,label="FMN:NADPH Redox",linewidth=3,linestyle="-")
plot(x,b,label="Bioluminescence",linewidth=1,linestyle=":")

xlabel("Heroin Uptake (mM/hr)")
ylabel("Reaction Flux (mM/hr)")
legend()
savefig("heroin_dep1.png")

return (her_up_max,her_deg_max,mor_deg_max,fmn_max,lux_max)
