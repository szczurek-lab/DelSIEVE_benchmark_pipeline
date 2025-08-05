configfile: "config.yaml"
include: "scripts/constants.py"
include: "scripts/utils.py"


#########################################################
#                Generate simulated data                #
#########################################################

generate_simulated_data(
    SIMUDATAROOT,
    config["tools"]["indelsieve_simulator"],
    config["configFiles"]["indelsieve_simulator"],
    SIMUDATAROOT + "simulation.log",
    config["benchmark"]["simulation"]["dataDir"],
    config["fineTune"]["indelsieve_simulator"]["overwrite"]
)


########################################################
#          Include predefined snakemake files          #
########################################################

include: "simulated_data.snake"
include: "monovar.snake"
include: "sciphi.snake"
include: "sciphin.snake"
include: "delsieve_2stage.snake"

# Monovar must be included for SiFit.
# Three options for SiFit is available. Choose one, comment out the other two.
# 1. Run SiFit on the local machine.
include: "sifit_local.snake"
# 2. Results from SiFit are available. Only collect results without running SiFit locally. This mostly happens when SiFit has successfully been run.
# include: "sifit_fixed.snake"
# 3. Run SiFit on a remote server.
# include: "sifit_remote.snake"


#########################################################
#                       All rules                       #
#########################################################

rule all:
    input:
        INDELSIEVEDIR + "CANDIDATE_COLLECTED"

        # Inferred trees
        ANALYSISTREECOMPDIR + "trees_info_updated.tsv",

        # Parameter estimates
        ANALYSISPARAMSDIR + "params_info.tsv",

        # Variant calling results
        ANALYSISVARCALLDIR + "variants_info.tsv",

        # ADO calling results
        ANALYSISADOCALLDIR + "ado_info.tsv",

        # Allelic information
        ANALYSISALLELICINFO + "allelic_info.rds",

        # Sites information
        ANALYSISSITESINFODIR + "sites_info.tsv",

        # Summary for parallel mutations
        SIMUPARALLELMUTDIR + "summary.tsv"


#######################################################
#                                                     #
#                Gather inferred trees                #
#                                                     #
#######################################################

# Rule: gather all inferred trees from different tools
rule gatherInferredTrees:
    input:
        ANALYSISTREECOMPDIR + "from_sciphi.tsv",
        ANALYSISTREECOMPDIR + "from_sciphin.tsv",
        ANALYSISTREECOMPDIR + "from_indelsieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISTREECOMPDIR + "from_indelsieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISTREECOMPDIR + "from_sifit_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISTREECOMPDIR + "trees_info.tsv")
    script:
        "scripts/gather_inferred_trees.py"

# Rule: get tree comparison results
rule getTreeComparisonResults:
    input:
        ANALYSISTREECOMPDIR + "trees_info.tsv"
    output:
        treeInfo=protected(ANALYSISTREECOMPDIR + "trees_info_updated.tsv")
    params:
        log=ANALYSISTREECOMPDIR + "out"
    shell:
        "Rscript --vanilla scripts/compare_trees.R "
        "-i {input} "
        "-o {output.treeInfo} &>{params.log}"


#######################################################
#                                                     #
#               Get parameter estimates               #
#                                                     #
#######################################################

# Rule: get all parameter estimates from different tools
rule getParams:
    input:
        ANALYSISPARAMSDIR + "from_sciphi.tsv",
        ANALYSISPARAMSDIR + "from_sciphin.tsv",
        ANALYSISPARAMSDIR + "from_indelsieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISPARAMSDIR + "from_indelsieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISPARAMSDIR + "from_sifit_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISPARAMSDIR + "params_info.tsv")
    script:
        "scripts/get_param_estimates.py"


#######################################################
#                                                     #
#             Get variant calling results             #
#                                                     #
#######################################################

# Rule: gather variant calling results from different tools
rule gatherVarResults:
    input:
        ANALYSISVARCALLDIR + "from_monovar_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv",
        ANALYSISVARCALLDIR + "from_sciphi.tsv",
        ANALYSISVARCALLDIR + "from_sciphin.tsv",
        ANALYSISVARCALLDIR + "from_indelsieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISVARCALLDIR + "from_indelsieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv"
    output:
        protected(ANALYSISVARCALLDIR + "variants_info.tsv")
    script:
        "scripts/gather_variant_calling_results.py"


#######################################################
#                                                     #
#               Get ADO calling results               #
#                                                     #
#######################################################

# Rule: gather ADO calling results from different tools
rule gatherAdoResults:
    input:
        ANALYSISADOCALLDIR + "from_indelsieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISADOCALLDIR + "from_indelsieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv"
    output:
        protected(ANALYSISADOCALLDIR + "ado_info.tsv")
    script:
        "scripts/gather_ado_calling_results.py"


#######################################################
#                                                     #
#                Get sites information                #
#                                                     #
#######################################################

# Rule: gather sites information from different tools
rule gatherSitesInfo:
    input:
        ANALYSISSITESINFODIR + "from_sciphi.tsv",
        ANALYSISSITESINFODIR + "from_sciphin.tsv",
        ANALYSISSITESINFODIR + "from_indelsieve.tsv",
        ANALYSISSITESINFODIR + "from_monovar_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISSITESINFODIR + "sites_info.tsv")
    script:
        "scripts/gather_sites_info.py"

