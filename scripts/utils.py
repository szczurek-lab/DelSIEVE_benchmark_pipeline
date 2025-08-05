import os
import numpy as np
import re

from snakemake import shell


NON_WORDS_PATTERN = "[^a-zA-Z0-9_]+"
NONE_PATTERN = r"(?!.*)"


def which_patterns(
    s: list[str], patterns: list[str], min_pat_num: int = 0, max_pat_num: int = 1
) -> list[list[int]]:
    if max_pat_num <= 0:
        raise ValueError("max_pat_num must be greater than 0.")
    ret: list[list[int]] = []
    for __s in s:
        __ret: list[int] = []
        for id, pat in enumerate(patterns):
            if re.search(re.compile(pat), __s):
                __ret.append(id)
        if len(__ret) > max_pat_num:
            raise ValueError(
                f"Error! More than {max_pat_num} pattern(s) ({__ret}) is found in {__s}."
            )
        elif len(__ret) < min_pat_num:
            raise ValueError(
                f"Error! Less than {min_pat_num} pattern(s) ({__ret}) is found in {__s}."
            )

        if len(__ret) > 0:
            ret.append(__ret)
        else:
            ret.append([-1])

    return ret


vec_search = np.vectorize(lambda s, pat: re.search(pat, s), otypes=[bool])
vec_abspath = np.vectorize(lambda i: os.path.abspath(i), otypes=[str])
vec_pattern = np.vectorize(lambda i, pat: f"{pat}{i}{pat}", otypes=[str])
vec_pattern_2 = np.vectorize(lambda i, pat: f"{pat}{i}", otypes=[str])


def get_matched_files(files: list[str], pattern: str) -> tuple[list[bool], list[str]]:
    idx = vec_search(files, pattern)
    matched = vec_abspath(sorted(np.array(files)[idx]))
    return idx, matched


def get_fine_tune_types(tool_fine_tune: dict, fine_tune_types: dict) -> list[str]:
    ret: list[str] = []
    for i in get_dict_keys(tool_fine_tune):
        if i in fine_tune_types.keys():
            ret.append(fine_tune_types[i])
        else:
            raise ValueError(
                f'Error! {i} is not in {fine_tune_types}. Make sure you have listed all fine tune types under "benchmark - fineTuneTypes" in "config.yaml".'
            )
    return ret


def get_data_types(
    used_data_types: list[str],
    available_data_types: dict,
    require_cnv: bool | None = None,
) -> list[str]:
    ret: list[str] = []
    for i in used_data_types:
        if i in available_data_types.keys():
            if (
                require_cnv is None
                or available_data_types[i]["requireCNV"] == require_cnv
            ):
                ret.append(available_data_types[i]["dir"].rstrip("/"))
        else:
            raise ValueError(
                f'Error! {i} is not in {available_data_types.keys()}. Do not change any values under "benchmark - simulatedDataTypes"!'
            )
    if require_cnv is not None and not require_cnv and len(ret) > 1:
        raise ValueError(
            f'Error! At most one CNV type is allowed to have "requireCNV" being false.'
        )
    return ret


def get_cnv_keep_type(key: str, options: dict) -> int:
    ret: list[int] = []
    for i in options.keys():
        if key == options[i]["dir"].rstrip("/"):
            if options[i]["requireCNV"]:
                if options[i]["includeCNV"]:
                    ret.append(1)
                else:
                    ret.append(2)
    if len(ret) == 0:
        raise ValueError(f"Error! {key} is not in {options.keys()}.")
    elif len(ret) > 1:
        raise ValueError(
            f"Error! More than one key in {{options.keys()}} is found for {key}."
        )
    return ret[0]


def get_subset_option(
    targetedTrueVarSitesNum: int | None, targetedBg2varRatio: float | None
) -> str:
    cmd1 = cmd2 = ""
    if targetedTrueVarSitesNum is None and targetedBg2varRatio is None:
        raise ValueError(
            f'At least one of "targetedTrueVarSitesNum" (={targetedTrueVarSitesNum}) and "targetedBg2varRatio" (={targetedBg2varRatio}) should be defined.'
        )
    if targetedTrueVarSitesNum is not None:
        cmd1 = f"--targetedTrueVarSitesNum {targetedTrueVarSitesNum}"
    if targetedBg2varRatio is not None:
        cmd2 = f"--targetedBg2varRatio {targetedBg2varRatio}"
    return " ".join([cmd1, cmd2])


def get_dict_keys(dict: dict) -> list:
    return list(dict.keys())


def get_dict_values(dict: dict, keys: list[str]):
    tmp = dict
    for key in keys:
        tmp = tmp[key]
    return tmp


def create_sifit_flag_file(file_name: str):
    if os.path.isfile(file_name):
        shell(
            'echo "Existing flag file detected: %s; skipping creation..." >&2'
            % file_name
        )
        return
    else:
        shell("mkdir -p %s" % os.path.abspath(os.path.dirname(file_name)))
        shell("touch %s" % file_name)
        shell('echo "Created: %s" >&2' % file_name)


def generate_simulated_data(
    base_path, tool_path, config_file_path, log_path, results_path, overwrite_flag
) -> None:
    """Generate simulated data.

    Args:
        base_path (string): Everything should be moved here. If exists, choose to overwrite or not.
        tool_path (string): Path to simulator.
        config_file_path (string): Path to simulation configuration file.
        log_path (string): Path to running log.
        results_path (string): Path to generated data defined in the configuration file.
        overwrite_flag (bool): TRUE or FALSE.
    """
    _overwrite_flag = False
    if os.path.isdir(results_path):
        # This is a sign of a failure previous run.
        shell('echo "Existing directory detected: %s" >&2' % results_path)
        shell('echo "Removing..." >&2')
        shell("rm -r %s" % results_path)
        _overwrite_flag = True
    if os.path.isdir(base_path):
        shell('echo "Existing directory detected: %s" >&2' % base_path)
        if _overwrite_flag | overwrite_flag:
            shell('echo "Overwriting..." >&2')
            shell("rm -r %s" % base_path)
            _generate_simulated_data(
                base_path, tool_path, config_file_path, log_path, results_path
            )
        else:
            shell('echo "Skipping simulated data generation..." >&2')
            return
    else:
        shell('echo "Generating simulated data..." >&2')
        _generate_simulated_data(
            base_path, tool_path, config_file_path, log_path, results_path
        )


def _generate_simulated_data(
    base_path, tool_path, config_file_path, log_path, results_path
):
    """Generate simulated data. Should only be called by generate_simulated_data internally.

    Args:
        base_path (string): Everything should be moved here.
        tool_path (string): Path to simulator.
        config_file_path (string): Path to simulation configuration file.
        log_path (string): Path to running log.
        results_path (string): Path to generated data defined in the configuration file.
    """
    shell("mkdir -p %s" % base_path)
    shell("%s -F%s &>%s" % (tool_path, config_file_path, log_path))
    shell("cp %s %s" % (config_file_path, base_path))
    shell("mv %s %s" % (results_path, base_path))


def get_dataset_names(path: str, num: int):
    file_full_names = [
        f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
    ]
    dataset_names = []
    for f in sorted(file_full_names):
        comp = f.split(".")
        dataset_names.append(comp[-1])

    if 0 < num < len(dataset_names):
        shell(
            'echo "Only using the first %d simulated datasets among the %d simulated ones in total." >&2'
            % (num, len(dataset_names))
        )
        return dataset_names[:num]
    else:
        return dataset_names


def get_abs_working_dir():
    return os.getcwd() + "/"


def get_matched_file_names(root: str, pattern: str, group: int = 0) -> list[str]:
    if not os.path.isdir(root):
        raise ValueError(f"Error! {root} is not a directory.")

    vec_match = np.vectorize(
        lambda pat, s, group=0: (
            "" if re.match(pat, s) is None else re.match(pat, s)[group]
        ),
        otypes=[str],
    )

    return [i for i in vec_match(pattern, os.listdir(root), group).tolist() if i != ""]


def get_tumor_cell_names(
    file_name: str,
    additional_names: list[str] | None = None,
    suffix: str | None = None,
) -> list[str] | None:
    if suffix is not None and suffix == "":
        suffix = None
    if os.path.isfile(file_name):
        ret: list[str] = []
        with open(
            file_name,
            "r",
            encoding="utf-8",
        ) as fh:
            for line in fh:
                ret.append(line.strip() if suffix is None else line.strip() + suffix)
        if additional_names is not None:
            ret.extend([i if suffix is None else i + suffix for i in additional_names])
        return ret
    return None


def get_parameter_estimate_type(input):
    if input == "all":
        return ["median", "mean", "mode_gaussian"]
    elif input == "median":
        return ["median"]
    elif input == "mean":
        return ["mean"]
    elif input == "mode":
        return ["mode_gaussian"]
    else:
        shell(
            "echo \"Unrecognized parameter estimater type. Should be one of 'all', 'median', 'mean', and 'mode'.\" >&2"
        )


def get_line_num(file: str) -> int:
    if not os.path.isfile(os.path.abspath(file)):
        return -1

    line_num: int = 0
    with open(file, "r") as fh:
        for line in fh:
            if not line.strip().startswith(" "):
                line_num += 1

    return line_num


def get_item_num(file: str) -> int:
    if not os.path.isfile(os.path.abspath(file)):
        return -1

    item_num: int = 0
    with open(file, "r") as fh:
        for line in fh:
            if not line.strip().startswith(" "):
                item_num += len(line.strip().split())

    return item_num


def get_dir_name(file: str) -> str:
    return os.path.dirname(file) + "/"


def get_content_files(file: str) -> list[str]:
    if os.path.isfile(file):
        results = []
        root = os.path.abspath(os.path.dirname(file))
        with open(file, "r") as fh:
            for line in fh:
                content = os.path.join(root, line.strip())
                if os.path.isfile(content):
                    results.append(content)
        return results


def get_thread_num_for_sieve(rawDataFile: str, nrOfSitesPerThread: int) -> None | int:
    if not os.path.isfile(os.path.abspath(rawDataFile)):
        return -1

    target_line: bool = False
    with open(rawDataFile, "r") as fh:
        for line in fh:
            if target_line:
                sites_num: int = int(line.strip())
                break
            if line.strip() == "=numCandidateMutatedSites=":
                target_line = True

    rawNrOfThreads: float = sites_num / nrOfSitesPerThread
    if int(rawNrOfThreads) == 0:
        # shell("echo \"%f\"" % rawNrOfThreads)
        return 1
    elif 0 <= rawNrOfThreads - int(rawNrOfThreads) < 0.5:
        # shell("echo \"%d\"" % int(rawNrOfThreads))
        return int(rawNrOfThreads)
    elif 0.5 <= rawNrOfThreads - int(rawNrOfThreads) <= 1:
        # shell("echo \"%d\"" % (int(rawNrOfThreads) + 1))
        return int(rawNrOfThreads) + 1


def get_thread_num_for_cellphy(dataFile: str, nrOfSitesPerThread: int) -> int:
    if not os.path.isfile(os.path.abspath(dataFile)):
        return -1

    target_line: bool = False
    sites_num: int = 0
    with open(dataFile, "r") as fh:
        for line in fh:
            if line.strip():
                if target_line:
                    sites_num += 1
                elif line.strip().upper().startswith("#CHROM"):
                    target_line = True

    rawNrOfThreads: float = sites_num / nrOfSitesPerThread
    if int(rawNrOfThreads) == 0:
        return 1
    elif 0 <= rawNrOfThreads - int(rawNrOfThreads) < 0.5:
        return int(rawNrOfThreads)
    elif 0.5 <= rawNrOfThreads - int(rawNrOfThreads) <= 1:
        return int(rawNrOfThreads) + 1


def get_dict_value(
    dict_str: dict,
    key: str,
    prefix: None | str = None,
    suffix: None | str = None,
    default: None | str | int | float = "",
) -> None | str | int | float:
    if key in dict_str.keys():
        ret: str = dict_str[key]

        if prefix is not None:
            ret = f"{prefix}{dict_str[key]}"
        if suffix is not None:
            ret = f"{ret}{suffix}"

        return ret
    else:
        return default


def sieve_config_files_sanity_check(
    setup: dict[str, dict], two_stage_sep: str = "-"
) -> tuple[list[str], list[str], list[str], dict[str, dict[str, str | float | int]]]:
    stage1_all: list[str] = []
    stage1: list[str] = []
    stage2: list[str] = []
    config2key: dict[str, dict[str, str | float | int]] = {}

    for tool_name, cases in setup.items():
        stage2_names: list[str] = []

        for case_name, case in cases.items():
            if "stage1" not in case.keys():
                raise ValueError(
                    f"Error! Missing `stage1` for {tool_name}:{case_name}."
                )

            if case["stage1"] not in stage1_all:
                stage1_all.append(case["stage1"])
                config2key[case["stage1"]] = {"file": case["stage1"], "tool": tool_name}

            if "stage2" in case.keys():
                if (
                    "name" not in case["stage2"].keys()
                    or "file" not in case["stage2"].keys()
                ):
                    raise ValueError(
                        f"Error! Both `name` and `file` should be defined for {tool_name}:{case_name}:stage2."
                    )

                if case["stage2"]["name"] in stage2_names:
                    raise ValueError(
                        f'Error! Duplicate name of `stage2` ({case["stage2"]["name"]}) found for {tool_name}:{case_name}.'
                    )
                stage2_names.append(case["stage2"]["name"])

                __name: str = two_stage_sep.join(
                    [case["stage1"], case["stage2"]["name"]]
                )

                if __name in config2key.keys():
                    raise ValueError(
                        f'Error! Duplicate names for stage1 ({case["stage1"]}) and stage2 ({case["stage2"]["name"]}) of {tool_name} found.'
                    )

                stage1.append(case["stage1"])
                stage2.append(case["stage2"]["name"])
                config2key[__name] = {"file": case["stage2"]["file"], "tool": tool_name}

                if (
                    "bgNum" in case["stage2"].keys()
                    and "bg2varRatio" in case["stage2"].keys()
                ):
                    raise ValueError(
                        f"Error! `bgNum` and `bg2varRatio` for {tool_name}:{case_name}:stage2 are mutually exclusive."
                    )

                if "bgNum" in case["stage2"].keys():
                    if case["stage2"]["bgNum"] <= 0:
                        raise ValueError(
                            f'Error! The number of background sites ({case["stage2"]["bgNum"]}) of {tool_name}:{case_name}:stage2:bgNum should be larger than 0.'
                        )
                    config2key[__name]["bgNum"] = case["stage2"]["bgNum"]

                if "bg2varRatio" in case["stage2"].keys():
                    if case["stage2"]["bg2varRatio"] <= 0:
                        raise ValueError(
                            f'Error! The ratio between the number of background sites and the number of variant sites ({case["stage2"]["bg2varRatio"]}) of {tool_name}:{case_name}:stage2:bg2varRatio should be larger than 0.'
                        )
                    config2key[__name]["bg2varRatio"] = case["stage2"]["bg2varRatio"]

    return stage1_all, stage1, stage2, config2key


def get_cnv_rate(data_types: list[str], cnv_rate: str | float | None) -> list[str]:
    if cnv_rate is None:
        cnv_rate = 0
    return ["0" if i == "no_cnv" else str(cnv_rate) for i in data_types]
