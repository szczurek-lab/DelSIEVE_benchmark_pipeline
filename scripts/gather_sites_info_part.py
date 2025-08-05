import sys
from typing import Union
import pandas as pd
import numpy as np
import os
import re

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from scripts.utils import (
    NON_WORDS_PATTERN,
    which_patterns,
    vec_pattern,
    get_matched_files,
    get_cnv_rate,
)


def get_sites_from_sieve(file: str) -> tuple[int, int, set[tuple[str, int]]]:
    num_mu_sites: int = 0
    num_background_sites: int = 0
    candidate_sites: set[tuple[str, int]] = set()

    with open(file, "r") as fh:
        flag_1: int = 0
        flag_2: int = 0
        flag_3: int = 0
        flag_4: int = 0
        for line in fh:
            if line:
                _line = line.strip()
                if flag_1 == 0 and _line == "=numCandidateMutatedSites=":
                    flag_1 = 1
                    continue
                elif flag_2 == 0 and _line == "=numBackgroundSites=":
                    flag_2 = 1
                    continue
                elif flag_3 == 0 and _line == "=mutations=":
                    flag_3 = 1
                    continue
                elif flag_4 == 0 and _line == "=background=":
                    flag_3 = 2
                    flag_4 = 1
                    continue
                elif flag_1 == 1:
                    flag_1 = 2
                    num_mu_sites = int(_line)
                elif flag_2 == 1:
                    flag_2 = 2
                    num_background_sites = int(_line)
                elif flag_3 == 1:
                    _comp: list[str] = _line.split()
                    candidate_sites.add((_comp[0], int(_comp[1])))
                elif flag_4 == 1:
                    flag_4 = 2
                    break
    return num_mu_sites, num_background_sites, candidate_sites


def get_sites_from_monovar(
    file: str,
) -> tuple[int, Union[int, None], set[tuple[str, int]]]:
    num_mu_sites: int = 0
    candidate_sites: set[tuple[str, int]] = set()

    with open(file, "r") as fh:
        for line in fh:
            if line and not line.strip().startswith("#"):
                _comp: list[str] = line.strip().split()
                candidate_sites.add((_comp[0], int(_comp[1])))
    return num_mu_sites, None, candidate_sites


def get_sites_from_sciphi(file: str) -> tuple[int, int, set[tuple[str, int]]]:
    num_samples: int = 0
    num_background_sites: int = 0
    candidate_sites: set[tuple[str, int]] = set()

    with open(file, "r") as fh:
        flag_1: int = 0
        flag_2: int = 0
        flag_3: int = 0
        for line in fh:
            if line:
                _line = line.strip()
                if flag_1 == 0 and _line == "=numSamples=":
                    flag_1 = 1
                    continue
                elif flag_2 == 0 and _line == "=mutations=":
                    flag_2 = 1
                    continue
                elif flag_3 == 0 and _line == "=background=":
                    flag_2 = 2
                    flag_3 = 1
                    continue
                elif flag_1 == 1:
                    flag_1 = 2
                    num_samples = int(_line)
                elif flag_2 == 1:
                    _comp: list[str] = _line.split()
                    candidate_sites.add((_comp[0], int(_comp[1])))
                elif flag_3 == 1:
                    _comp: list[str] = _line.split()
                    idx: int = 0
                    while idx < len(_comp):
                        num_background_sites += int(_comp[idx + 1])
                        idx += 2
                    break

    assert num_samples > 0
    if num_background_sites % num_samples != 0:
        raise ValueError(f"Error! the number of background sites is not integer.")
    else:
        num_background_sites = int(num_background_sites / num_samples)

    return len(candidate_sites), num_background_sites, candidate_sites


def get_sites_info(
    row: pd.Series, extra_eval: Union[None, list[str]] = ["m&d", "m&i", "i", "d", "m|d"]
) -> pd.Series:
    true_sites: set[tuple[str, int]] = set()

    extra_true_sites: Union[None, list[set[tuple[str, int]]]] = (
        None if extra_eval is None else list()
    )
    extra_tp: Union[None, list[set[tuple[str, int]]]] = (
        None if extra_eval is None else list()
    )
    extra_fn: Union[None, list[set[tuple[str, int]]]] = (
        None if extra_eval is None else list()
    )

    _extra_eval: Union[None, list[str]] = (
        None if extra_eval is None else [i.upper().replace(" ", "") for i in extra_eval]
    )
    extra_tp_name: Union[None, list[str]] = (
        None if _extra_eval is None else ["_".join(["tp", i]) for i in _extra_eval]
    )
    extra_fn_name: Union[None, list[str]] = (
        None if _extra_eval is None else ["_".join(["fn", i]) for i in _extra_eval]
    )
    extra_recall_name: Union[None, list[str]] = (
        None if _extra_eval is None else ["_".join(["recall", i]) for i in _extra_eval]
    )

    ret_name: list[str] = ["tp", "fp", "fn", "recall", "precision", "f1_score"]

    with open(row.true_sites_file, "r") as fh:
        _map: dict[str, str] = {
            "M": "(int(_line[2]) > 0)",
            "D": "(int(_line[3]) > 0)",
            "I": "(int(_line[4]) > 0)",
        }
        __extra_eval: Union[None, list[str]] = None if _extra_eval is None else list()
        if _extra_eval is not None:
            for i in _extra_eval:
                ret: str = i
                for k, v in _map.items():
                    ret = ret.replace(k, v)
                __extra_eval.append(ret)
                extra_true_sites.append(set())
                extra_tp.append(set())
                extra_fn.append(set())

        for line in fh:
            if line and not line.strip().startswith("#"):
                _line = line.strip().split()

                true_sites.add((_line[0], int(_line[1])))

                if _extra_eval is not None:
                    for id, _eval in enumerate(__extra_eval):
                        if eval(_eval):
                            extra_true_sites[id].add((_line[0], int(_line[1])))

    num_mu_sites = None
    num_background_sites = None
    candidate_sites = set()
    if re.search(r"sieve", row.tool) is not None:
        num_mu_sites, num_background_sites, candidate_sites = get_sites_from_sieve(
            row.candidate_sites_file
        )
    elif re.search(r"monovar", row.tool) is not None:
        num_mu_sites, num_background_sites, candidate_sites = get_sites_from_monovar(
            row.candidate_sites_file
        )
    elif re.search(r"sciphi", row.tool) is not None:
        num_mu_sites, num_background_sites, candidate_sites = get_sites_from_sciphi(
            row.candidate_sites_file
        )
    else:
        raise ValueError(f"Error! Unknown tool {row.tool}.")

    tp: set[tuple[str, int]] = true_sites.intersection(candidate_sites)
    fp: set[tuple[str, int]] = candidate_sites.difference(true_sites)
    fn: set[tuple[str, int]] = true_sites.difference(candidate_sites)
    ret_val: list[float] = [
        len(tp),
        len(fp),
        len(fn),
        len(tp) / (len(tp) + len(fn)),
        len(tp) / (len(tp) + len(fp)),
        len(tp) / (len(tp) + (len(fp) + len(fn)) / 2),
    ]

    if extra_eval is not None:
        for id, _true_sites in enumerate(extra_true_sites):
            extra_tp[id] = _true_sites.intersection(candidate_sites)
            extra_fn[id] = _true_sites.difference(candidate_sites)

            ret_name.append(extra_tp_name[id])
            ret_val.append(len(extra_tp[id]))
            ret_name.append(extra_fn_name[id])
            ret_val.append(len(extra_fn[id]))
            ret_name.append(extra_recall_name[id])

            if (len(extra_tp[id]) + len(extra_fn[id])) == 0:
                ret_val.append(np.nan)
            else:
                ret_val.append(
                    len(extra_tp[id]) / (len(extra_tp[id]) + len(extra_fn[id]))
                )

    return pd.Series(
        [
            row.candidate_sites_file,
            row.true_sites_file,
            num_mu_sites,
            num_background_sites,
            *ret_val,
        ],
        index=[
            "candidate_sites_file",
            "true_sites_file",
            "num_candidate_mutated_sites",
            "num_background_sites",
            *ret_name,
        ],
    )


search_for_tool_setup_flag = False
if snakemake.params["toolSetup"] is None or type(snakemake.params["toolSetup"]) == str:
    tool_setups = [snakemake.params["toolSetup"]]
elif type(snakemake.params["toolSetup"]) == list:
    search_for_tool_setup_flag = True
    tool_setups = sorted(set(snakemake.params["toolSetup"]))
else:
    raise ValueError(
        f'toolSetup must be a string or a list of strings ({snakemake.params["toolSetup"]} was given).'
    )

fine_tune_type_flag = False
if snakemake.params["fineTuneType"] is None:
    fine_tune_types = np.repeat(np.nan, len(snakemake.input["candidateSites"]))
elif len(snakemake.params["fineTuneType"]) == 1:
    fine_tune_types = np.repeat(
        snakemake.params["fineTuneType"], len(snakemake.input["candidateSites"])
    )
else:
    fine_tune_types = []
    fine_tune_type_flag = True

tool = []
tool_setup = []
data_type = []
dataset = []
candidate_sites_file_names = []
true_sites_file_names = []


for __tool_setup in tool_setups:
    for __data_type in sorted(snakemake.params["dataType"]):
        for dataset_name in sorted(snakemake.params["datasetNames"]):
            if search_for_tool_setup_flag:
                __pat = "".join(
                    [
                        f"(?=.*{i})"
                        for i in vec_pattern(
                            [__tool_setup, __data_type, dataset_name], NON_WORDS_PATTERN
                        )
                    ]
                )
            else:
                __pat = "".join(
                    [
                        f"(?=.*{i})"
                        for i in vec_pattern(
                            [__data_type, dataset_name], NON_WORDS_PATTERN
                        )
                    ]
                )

            inferred_idx, __candidate_sites = get_matched_files(
                snakemake.input["candidateSites"], __pat
            )
            cnt = np.sum(inferred_idx)

            if cnt == 0:
                print(
                    f"Warning: no matches found for {__pat} in {snakemake.input['candidateSites']}."
                )

            if type(snakemake.params["tool"]) == dict:
                tool.extend(
                    np.repeat([snakemake.params["tool"][__tool_setup]["tool"]], cnt)
                )
            else:
                tool.extend(np.repeat([snakemake.params["tool"]], cnt))
            tool_setup.extend(np.repeat(__tool_setup, cnt))
            data_type.extend(np.repeat(__data_type, cnt))
            dataset.extend(np.repeat(dataset_name, cnt))
            candidate_sites_file_names.extend(__candidate_sites)

            __true_sites_pat = "".join(
                [
                    f"(?=.*{i})"
                    for i in vec_pattern([__data_type, dataset_name], NON_WORDS_PATTERN)
                ]
            )
            true_idx, __true_sites = get_matched_files(
                snakemake.input["trueSites"], __true_sites_pat
            )
            if np.sum(true_idx) == 0 or np.sum(true_idx) > 1:
                raise ValueError(
                    f"{__true_sites_pat} should have only a single file for true sites, but {np.sum(true_idx)} found."
                )
            true_sites_file_names.extend(np.repeat(__true_sites, cnt))

            if fine_tune_type_flag:
                fine_tune_type_idx = which_patterns(
                    __candidate_sites,
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


sites_collection = pd.DataFrame(
    {
        "cell_num": snakemake.params["cellNum"],
        "coverage_mean": snakemake.params["covMean"],
        "coverage_variance": snakemake.params["covVariance"],
        "eff_seq_err_rate": snakemake.params["effSeqErrRate"],
        "ado_rate": snakemake.params["adoRate"],
        "deletion_rate": snakemake.params["deletionRate"],
        "insertion_rate": snakemake.params["insertionRate"],
        "doublet_rate": snakemake.params["doubletRate"],
        "cnv_rate": get_cnv_rate(
            data_type,
            snakemake.params["cnvRate"],
        ),
        "dataset": pd.Series(dataset),
        "tool": pd.Series(tool),
        "snv_type": snakemake.params["snvType"],
        "tool_setup": pd.Series(tool_setup),
        "fine_tune_type": pd.Series(fine_tune_types),
        "data_type": pd.Series(data_type),
        "candidate_sites_file": pd.Series(candidate_sites_file_names),
        "true_sites_file": pd.Series(true_sites_file_names),
    }
).drop_duplicates()

sites_collection = sites_collection.merge(
    sites_collection.apply(
        lambda x: get_sites_info(x, snakemake.params["typeOfQuery"]), axis=1
    ),
    on=["candidate_sites_file", "true_sites_file"],
).to_csv(snakemake.output[0], sep="\t", na_rep="NA", index=False)
