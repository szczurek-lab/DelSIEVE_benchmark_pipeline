# This file will be loaded by snakemake files.


configfile: "config.yaml"


##############################################
#         Set up constant variables          #
##############################################

# Constants for simulated data; 'SIMU' for 'SIMULATION'
SIMUDATAROOT = config["benchmark"]["baseDir"] + config["benchmark"]["simulation"]["dir"]
SIMUDATADIR = SIMUDATAROOT + config["benchmark"]["simulation"]["dataDir"]
SIMUTUMORCELLNAMESDIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["tumorCellNamesDir"]
)
SIMUTUMORCELLNAMES = (
    SIMUDATADIR + config["benchmark"]["simulation"]["tumorCellNamesPre"]
)
SIMUSNVSITESDIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVSitesNamesDir"]
)
SIMUTRUESNVDIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVDir"]
SIMUPARALLELMUTDIR = SIMUDATADIR + config["benchmark"]["simulation"]["parallelMutDir"]
# SIMUTRUESNVADODIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVAdoDir"]
SIMUALLELICSEQINFODIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["trueAllelicSeqInfoSNVsDir"]
)
SIMUTRUEADOSTATESDIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["trueAdoStatesDir"]
)
SIMUTRUESIZEFACTORSDIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["trueSizeFactorsDir"]
)

# Simulated true trees
SIMUTRUETREEDIR = SIMUDATADIR + config["benchmark"]["simulation"]["treesDir"]
SIMUCONVERTEDTRUETREEDIR = (
    SIMUDATADIR + config["benchmark"]["simulation"]["convertedTreesDir"]
)

# Simulated raw data
SIMURAWDATADIR = SIMUDATADIR + config["benchmark"]["simulation"]["rawDataDir"]
SIMUDUMMYBAMDATADIR = SIMURAWDATADIR + config["benchmark"]["simulation"]["dummyBamDir"]

# Simulated aligned reads (MPILEUP)
SIMUMPILEUPDATADIR = SIMUDATADIR + config["benchmark"]["simulation"]["MPileUpDataDir"]
SIMUMPILEUPMERGED = (
    config["benchmark"]["simulation"]["MPileUpDataMergedPre"]
    + config["benchmark"]["simulation"]["MPileUpDataSuf"]
)
SIMUMPILEUPUNMERGED = (
    config["benchmark"]["simulation"]["MPileUpDataUnmergedPre"]
    + config["benchmark"]["simulation"]["MPileUpDataSuf"]
)
SIMUMPILEUPNORMALMERGED = (
    config["benchmark"]["simulation"]["MPileUpDataNormalMergedPre"]
    + config["benchmark"]["simulation"]["MPileUpDataSuf"]
)
SIMUMPILEUPNORMALUNMERGED = (
    config["benchmark"]["simulation"]["MPileUpDataNormalUnMergedPre"]
    + config["benchmark"]["simulation"]["MPileUpDataSuf"]
)

# Add CNVs to simulated data
SIMUDATATYPESDIR = SIMUDATADIR + config["benchmark"]["simulation"]["dataTypes"]["dir"]

# Results of INDELSIEVE
INDELSIEVEDIR = (
    config["benchmark"]["baseDir"] + config["benchmark"]["indelsieve"]["dir"]
)
INDELSIEVEFROMCANDIDATEDIR = (
    INDELSIEVEDIR + config["benchmark"]["phyloInf"]["fromCandidate"]
)
INDELSIEVEFROMCANDIDATEPHYINFDIR = (
    INDELSIEVEFROMCANDIDATEDIR + config["benchmark"]["indelsieve"]["phyloInfDir"]
)
INDELSIEVEFROMCANDIDATEPHYINFSTAGE2DIR = (
    INDELSIEVEFROMCANDIDATEPHYINFDIR + config["benchmark"]["indelsieve"]["stage2Dir"]
)
INDELSIEVEFROMCANDIDATEVARCALDIR = (
    INDELSIEVEFROMCANDIDATEDIR + config["benchmark"]["indelsieve"]["variantCallingDir"]
)
INDELSIEVEFROMCANDIDATEVARCALSTAGE2DIR = (
    INDELSIEVEFROMCANDIDATEVARCALDIR + config["benchmark"]["indelsieve"]["stage2Dir"]
)

# Results of SIEVE
# SIEVEDIR = config["benchmark"]["baseDir"] + \
#     config["benchmark"]["sieve"]["dir"]
# SIEVEFROMCANDIDATEDIR = SIEVEDIR + config["benchmark"]["phyloInf"]["fromCandidate"]
# SIEVEFROMCANDIDATEPHYINFDIR = SIEVEFROMCANDIDATEDIR + config["benchmark"]["sieve"]["phyloInfDir"]
# SIEVEFROMCANDIDATEPHYINFSTAGE2DIR = SIEVEFROMCANDIDATEPHYINFDIR + config["benchmark"]["sieve"]["stage2Dir"]
# SIEVEFROMCANDIDATEVARCALDIR = SIEVEFROMCANDIDATEDIR + config["benchmark"]["sieve"]["variantCallingDir"]
# SIEVEFROMCANDIDATEVARCALSTAGE2DIR = SIEVEFROMCANDIDATEVARCALDIR + config["benchmark"]["sieve"]["stage2Dir"]

# # Results of CELLPHY
# CELLPHYDIR = config["benchmark"]["baseDir"] + \
#     config["benchmark"]["cellphy"]["dir"]
# CELLPHYTRUEMONOVARSNVSDIR = CELLPHYDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]

# Git manager
GITSKELETONLOCALDIR = (
    config["benchmark"]["baseDir"] + config["benchmark"]["git"]["skeletonDir"]
)
GITREMOTEDIR = (
    config["servers"]["remoteServer"]["rootPath"]
    + config["benchmark"]["git"]["skeletonDir"]
)

# Results of SIFIT
# Skeleton
SIFITGITSKELETONLOCALDIR = GITSKELETONLOCALDIR + config["benchmark"]["sifit"]["dir"]
SIFITGITSKELETONTRUEMONOVARLOCALDIR = (
    SIFITGITSKELETONLOCALDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]
)

SIFITGITREMOTEDIR = GITREMOTEDIR + config["benchmark"]["sifit"]["dir"]
SIFITGITSKELETONTRUEMONOVARREMOTEDIR = (
    SIFITGITREMOTEDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]
)

# Local
SIFITLOCALDIR = config["benchmark"]["baseDir"] + config["benchmark"]["sifit"]["dir"]

SIFITTRUEMONOVARSNVSLOCALDIR = (
    SIFITLOCALDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]
)
SIFITLOCALFLAGFILETRUEMONOVARINITSKELETONIN = (
    SIFITTRUEMONOVARSNVSLOCALDIR
    + config["benchmark"]["sifit"]["flag"]["dir"]
    + config["benchmark"]["sifit"]["flag"]["initSkeletonIn"]
)
SIFITLOCALFLAGFILETRUEMONOVARINITSKELETONOUT = (
    SIFITTRUEMONOVARSNVSLOCALDIR
    + config["benchmark"]["sifit"]["flag"]["dir"]
    + config["benchmark"]["sifit"]["flag"]["initSkeletonOut"]
)
SIFITLOCALFLAGFILEPUSHEDTRUEMONOVARDATA = (
    SIFITTRUEMONOVARSNVSLOCALDIR
    + config["benchmark"]["sifit"]["flag"]["dir"]
    + config["benchmark"]["sifit"]["flag"]["pushedData"]
)
SIFITTRUEMONOVARSNVSTRUEPARAMSLOCALDIR = (
    SIFITTRUEMONOVARSNVSLOCALDIR
    + config["benchmark"]["paramTypes"]["fromTrueParamsDir"]
)

# Remote, if applicable
SIFITREMOTEDIR = (
    config["servers"]["remoteServer"]["rootPath"] + config["benchmark"]["sifit"]["dir"]
)
SIFITTRUEMONOVARSNVSREMOTEDIR = (
    SIFITREMOTEDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]
)

# Results of MONOVAR
MONOVARDIR = config["benchmark"]["baseDir"] + config["benchmark"]["monovar"]["dir"]
MONOVARTRUEPARAMSDIR = (
    MONOVARDIR + config["benchmark"]["paramTypes"]["fromTrueParamsDir"]
)

# Results of SCIPHI
SCIPHIDIR = config["benchmark"]["baseDir"] + config["benchmark"]["sciphi"]["dir"]

# Results of SCIPHIN
SCIPHINDIR = config["benchmark"]["baseDir"] + config["benchmark"]["sciphin"]["dir"]

# # Results of GATK
# GATKDIR = config["benchmark"]["baseDir"] + \
#     config["benchmark"]["gatk"]["dir"]

# # Results of BCFTOOLS
# BCFTOOLSDIR = config["benchmark"]["baseDir"] + \
#     config["benchmark"]["bcftools"]["dir"]

# # Results of Sccaller
# SCCALLERDIR = config["benchmark"]["baseDir"] + \
#     config["benchmark"]["sccaller"]["dir"]

# Analysis results
ANALYSISDIR = config["benchmark"]["baseDir"] + config["benchmark"]["analysis"]["dir"]

# Analysis results of tree comparison
ANALYSISTREECOMPDIR = ANALYSISDIR + config["benchmark"]["analysis"]["treeComparisonDir"]

# Analysis results of parameter estimates
ANALYSISPARAMSDIR = ANALYSISDIR + config["benchmark"]["analysis"]["paramsEstimatesDir"]

# Analysis results of variant calling
ANALYSISVARCALLDIR = ANALYSISDIR + config["benchmark"]["analysis"]["varCallResultsDir"]

# Analysis results of ado calling
ANALYSISADOCALLDIR = ANALYSISDIR + config["benchmark"]["analysis"]["adoCallResultsDir"]

# Analysis results of sites info
ANALYSISSITESINFODIR = ANALYSISDIR + config["benchmark"]["analysis"]["sitesInfoDir"]

# Analysis results of allelic information
ANALYSISALLELICINFO = (
    ANALYSISDIR + config["benchmark"]["analysis"]["allelicInfoPerCellDir"]
)

# Efficiency benchmarking results
EFFICIENCYDIR = (
    config["benchmark"]["baseDir"] + config["benchmark"]["efficiency"]["dir"]
)

EFFICIENCYMONOVARDIR = EFFICIENCYDIR + config["benchmark"]["monovar"]["dir"]

EFFICIENCYSCIPHIDIR = EFFICIENCYDIR + config["benchmark"]["sciphi"]["dir"]

EFFICIENCYSCIPHINDIR = EFFICIENCYDIR + config["benchmark"]["sciphin"]["dir"]

# EFFICIENCYSIEVEDIR = EFFICIENCYDIR + config["benchmark"]["sieve"]["dir"]

EFFICIENCYINDELSIEVEDIR = EFFICIENCYDIR + config["benchmark"]["indelsieve"]["dir"]

EFFICIENCYCELLPHYDIR = EFFICIENCYDIR + config["benchmark"]["cellphy"]["dir"]

EFFICIENCYSIFITDIR = EFFICIENCYDIR + config["benchmark"]["sifit"]["dir"]
