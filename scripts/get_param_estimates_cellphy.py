import sys
import pandas as pd
import numpy as np
import os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from scripts.utils import (
    NON_WORDS_PATTERN,
    which_patterns,
    vec_pattern,
    get_matched_files,
    get_cnv_rate,
)


def parse_params(params_file: str) -> tuple[float, float, float]:
    if os.path.isfile(params_file):
        eff_seq_err_rate = ado = gamma = np.nan
        with open(params_file, "r") as fh:
            for line in fh:
                if "ERR_P17" in line:
                    index = line.index("ERR_P17")
                    start_index = end_index = index
                    for i in range(index, len(line)):
                        if line[i] == "{":
                            start_index = i + 1
                        elif line[i] == "}":
                            end_index = i
                            break
                    if start_index >= end_index:
                        raise ValueError(
                            "Error! Invalid format of error parameters in this file: "
                            + params_file
                        )
                    else:
                        params = line[start_index:end_index].strip().split("/")
                        eff_seq_err_rate = float(params[0].strip())
                        ado = float(params[1].strip())

                if "G4m" in line:
                    index = line.index("G4m")
                    start_index = end_index = index
                    for i in range(index, len(line)):
                        if line[i] == "{":
                            start_index = i + 1
                        elif line[i] == "}":
                            end_index = i
                            break
                    if start_index >= end_index:
                        raise ValueError(
                            "Error! Invalid format of error parameters in this file: "
                            + params_file
                        )
                    else:
                        gamma = float(line[start_index:end_index].strip())

                return eff_seq_err_rate, ado, gamma


search_for_tool_setup_flag = False
if snakemake.params["toolSetup"] is None or type(snakemake.params["toolSetup"]) == str:
    tool_setups = [snakemake.params["toolSetup"]]
elif type(snakemake.params["toolSetup"]) == list:
    search_for_tool_setup_flag = True
    tool_setups = sorted(snakemake.params["toolSetup"])
else:
    raise ValueError(
        f'toolSetup must be a string or a list of strings ({snakemake.params["toolSetup"]} was given).'
    )

fine_tune_type_flag = False
if snakemake.params["fineTuneType"] is None:
    fine_tune_types = np.repeat(np.nan, len(snakemake.input))
elif len(snakemake.params["fineTuneType"]) == 1:
    fine_tune_types = np.repeat(snakemake.params["fineTuneType"], len(snakemake.input))
else:
    fine_tune_types = []
    fine_tune_type_flag = True

data_type_flag = False
if snakemake.params["dataType"] is None:
    data_types = np.repeat(np.nan, len(snakemake.input))
elif len(snakemake.params["dataType"]) == 1:
    data_types = np.repeat(snakemake.params["dataType"], len(snakemake.input))
else:
    data_types = []
    data_type_flag = True


tool_setup = []
dataset = []
eff_seq_err_rates = []
ado_rates = []
gamma_shape = []

for __tool_setup in tool_setups:
    for dataset_name in sorted(snakemake.params["datasetNames"]):
        if search_for_tool_setup_flag:
            __pat = "".join(
                [
                    f"(?=.*{i})"
                    for i in vec_pattern(
                        [__tool_setup, dataset_name], NON_WORDS_PATTERN
                    )
                ]
            )
        else:
            __pat = vec_pattern(dataset_name, NON_WORDS_PATTERN)

        inferred_idx, __inferred_params = get_matched_files(snakemake.input, __pat)
        cnt = np.sum(inferred_idx)

        if cnt == 0:
            print(f"Warning: no matches found for {__pat} in {snakemake.input}.")

        tool_setup.extend(np.repeat(__tool_setup, cnt))
        dataset.extend(np.repeat(dataset_name, cnt))

        for __params in __inferred_params:
            eff_seq_err_rate, ado, gamma = parse_params(__params)
            eff_seq_err_rates.append(eff_seq_err_rate)
            ado_rates.append(ado)
            gamma_shape.append(gamma)

        if fine_tune_type_flag:
            fine_tune_type_idx = which_patterns(
                __inferred_params,
                vec_pattern(snakemake.params["fineTuneType"], NON_WORDS_PATTERN),
            )
            fine_tune_types.extend(
                [
                    snakemake.params["fineTuneType"][j]
                    for i in fine_tune_type_idx
                    for j in i
                    if j >= 0
                ]
            )

        if data_type_flag:
            data_type_idx = which_patterns(
                __inferred_params,
                vec_pattern(snakemake.params["dataType"], NON_WORDS_PATTERN),
            )
            data_types.extend(
                [
                    snakemake.params["dataType"][j]
                    for i in data_type_idx
                    for j in i
                    if j >= 0
                ]
            )


params_collection = pd.DataFrame(
    {
        "cell_num": snakemake.params["cellNum"],
        "coverage_mean": snakemake.params["covMean"],
        "coverage_variance": snakemake.params["covVariance"],
        "mutation_rate": snakemake.params["mutationRate"],
        "doublet_rate": snakemake.params["doubletRate"],
        "cnv_rate": get_cnv_rate(
            data_types,
            snakemake.params["cnvRate"],
        ),
        "dataset": pd.Series(dataset),
        "tool": snakemake.params["tool"],
        "snv_type": snakemake.params["snvType"],
        "tool_setup": snakemake.params["toolSetup"],
        "fine_tune_type": pd.Series(fine_tune_types),
        "data_type": pd.Series(data_types),
        "eff_seq_err_rate": pd.Series(eff_seq_err_rates),
        "ado_rate": pd.Series(ado_rates),
        "gamma_shape": pd.Series(gamma_shape),
    }
)

params_collection.to_csv(snakemake.output[0], sep="\t", na_rep="NA", index=False)
