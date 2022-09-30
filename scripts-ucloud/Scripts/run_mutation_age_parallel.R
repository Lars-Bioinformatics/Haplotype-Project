setwd("/work/sduvarcall/haplotype-project")
source("Scripts/Mutation-Age.R")

# run_compareAgeMethod_generateFiles_parallel(g_start=10, g_stop=510, g_step=10, sim=20, s_type = "famHaplotypesNoHet", 
#                                             ancestral_method="simulatedFounder", cs_correction = "none")

# run_compareAgeMethod_generateFiles_parallel(g_start=10, g_stop=1210, g_step=10, sim=50, s_type = "knownBreaks", 
#                                             ancestral_method="simulatedFounder", cs_correction = "none")

for (sim_type in c("haplotypes", "famHaplotypes", "famHaplotypesNoHet")){
# for (sim_type in c("famHaplotypes","haplotypes")){
# for (sim_type in c("famHaplotypes")){
# for (sim_type in c("haplotypes")){
    # for (ancestral in c("simulatedFounder", "branchBoundIndep", "branchBound", "mostFreqBase")){
    # for (ancestral in c("branchBoundIndep","mostFreqBase")){
    for (ancestral in c("branchBoundIndep")){
        run_compareAgeMethod_generateFiles_parallel(g_start=10, g_stop=510, g_step=10, sim=50, minSamplesBB='auto3',
                                                    s_type = sim_type, ancestral_method=ancestral, samples = 100,
                                                    cs_correction = "none", eps = 0.01)
                                                    # cs_correction = "gandolfo", eps = 0.01)
    }
}

### COMPARE SAMPLE SIZES
# run_compareSampleSize_generateFiles_parallel(sample_sizes = c(seq(10,100,10), seq(150,1000,50)), gen = 100, sim = 50,
#                                              minSamplesBB = 0.25, s_type = "knownBreaks", ancestral_method = "simulatedFounder",
#                                              cs_correction = "none", eps = 0.01)

# for (s_type in c("famHaplotypes")){
# # for (s_type in c("haplotypes")){
# # for (s_type in c("famHaplotypes", "haplotypes")){
#     # for (gen in c(50,100)){
#     for (gen in c(100)){
#         for (ancestral in c("simulatedFounder", "branchBoundIndep", "mostFreqBase")){
#             run_compareSampleSize_generateFiles_parallel(sample_sizes = c(seq(10,100,10), seq(150,1000,50)), gen = gen, sim = 50,
#                                                          minSamplesBB = 0.25, s_type = s_type, ancestral_method = ancestral,
#                                                          cs_correction = "none", eps = 0.01)
#         }
#     }
# }
