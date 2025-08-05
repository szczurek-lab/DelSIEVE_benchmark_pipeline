import argparse
import os
import subprocess
import numpy as np


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(
        description="Subset variant sites in pileup data"
    )
    sys_args_parser.add_argument(
        "--targetedTrueVarSitesNum",
        type=int,
        default=None,
        help="the targeted number of true variant sites after subsetting",
    )
    sys_args_parser.add_argument(
        "--targetedBg2varRatio",
        type=float,
        default=None,
        help="the targeted ratio between the number of background sites and the number of variant sites after subsetting",
    )
    sys_args_parser.add_argument(
        "--prefix",
        type=str,
        required=True,
        help="prefix to output files",
    )
    sys_args_parser.add_argument(
        "--inPileupDir",
        type=str,
        required=True,
        help="input directory to pileup files",
    )
    sys_args_parser.add_argument(
        "--trueVarSites",
        type=str,
        required=True,
        help="input true variant sites",
    )
    sys_args_parser.add_argument(
        "--addTrueInfo",
        type=str,
        help="additional information of true sites (comma-separated)",
    )
    sys_args_parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="random seed",
    )

    sys_args = sys_args_parser.parse_args()
    if not os.path.isdir(sys_args.inPileupDir):
        raise ValueError(
            f"Error! Input directory {sys_args.inPileupDir} does not exist."
        )
    if (
        sys_args.targetedTrueVarSitesNum is None
        and sys_args.targetedBg2varRatio is None
    ):
        raise ValueError(
            f'At least one of "targetedTrueVarSitesNum" (={sys_args.targetedTrueVarSitesNum}) and "targetedBg2varRatio" (={sys_args.targetedBg2varRatio}) should be defined.'
        )
    if not os.path.isdir(os.path.dirname(sys_args.prefix)):
        os.makedirs(
            os.path.dirname(sys_args.prefix),
            exist_ok=True,
        )
    return sys_args


def parse_pileup(file: str) -> dict[tuple[str, int], str]:
    if not os.path.isfile(file):
        raise ValueError(f"Error! File {file} does not exist.")

    ret: dict[tuple[str, int], str] = dict()
    with open(
        file,
        "r",
        encoding="utf-8",
    ) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("#"):
                comp: list[str] = line.split("\t")
                key: tuple[str, int] = (comp[0], int(comp[1]))
                if key in ret:
                    raise ValueError(f"Error! Duplicate position found: {key}.")
                ret[key] = line
    return ret


def parse_true_var_sites(
    file: str,
) -> None | tuple[list[str], dict[int, tuple[tuple[str, int], bool, str]]]:
    if not os.path.isfile(file):
        raise ValueError(f"Error! File {file} does not exist.")
    # for comments
    ret1: list[str] = list()
    # relative index, chromosome label, position, has SNVs, original string
    ret2: dict[int, tuple[tuple[str, int], bool, str]] = dict()
    with open(
        file,
        "r",
        encoding="utf-8",
    ) as fh:
        idx: int = 0
        for line in fh:
            line = line.strip()
            if line.startswith("#"):
                ret1.append(line)
            else:
                comp: list[str] = line.split("\t")
                ret2[idx] = ((comp[0], int(comp[1])), int(comp[2]) == 1, line)
                idx += 1
    return (ret1, ret2)


def get_targeted_var_sites_num(
    total_sites_num: int,
    true_var_sites_num: int,
    targeted_true_var_sites_num: None | int,
    targeted_bg2var_ratio: None | float,
) -> None | int:
    if targeted_true_var_sites_num is None and targeted_bg2var_ratio is None:
        return None
    if targeted_true_var_sites_num is None:
        targeted_true_var_sites_num = true_var_sites_num
    if targeted_bg2var_ratio is not None:
        tmp: int = int((total_sites_num - true_var_sites_num) / targeted_bg2var_ratio)
        targeted_true_var_sites_num = (
            targeted_true_var_sites_num if tmp >= targeted_true_var_sites_num else tmp
        )
    if targeted_true_var_sites_num >= true_var_sites_num:
        return None
    else:
        return targeted_true_var_sites_num


def get_file_name(prefix: str, file: str) -> str:
    return "".join([prefix, file.split(os.path.sep)[-1]])


def main():
    sys_args: argparse.Namespace = parse_system_args()

    reads: list[tuple[str, dict[tuple[str, int], str]]] = []
    total_num_sites: int = 0

    for i in os.listdir(sys_args.inPileupDir):
        if i.endswith(".pileup"):
            reads.append(
                (
                    i,
                    parse_pileup(
                        os.path.join(
                            sys_args.inPileupDir,
                            i,
                        ),
                    ),
                )
            )

            if total_num_sites == 0:
                total_num_sites = len(reads[-1][-1])
            else:
                assert total_num_sites == len(reads[-1][-1])

    true_var_sites: (
        None | tuple[list[str], dict[int, tuple[tuple[str, int], bool, str]]]
    ) = parse_true_var_sites(sys_args.trueVarSites)
    assert true_var_sites is not None
    subsetable_idx: list[int] = [k for k, v in true_var_sites[1].items() if v[1]]
    unsubsetable_idx: list[int] = [k for k, v in true_var_sites[1].items() if not v[1]]

    target_var_num: None | int = get_targeted_var_sites_num(
        total_num_sites,
        len(true_var_sites[1]),
        sys_args.targetedTrueVarSitesNum,
        sys_args.targetedBg2varRatio,
    )
    if target_var_num is None:
        print("No need to subset. Copy files to destination...")
        for file in [
            sys_args.inPileupMerged,
            sys_args.inPileupUnmerged,
            sys_args.trueVarSites,
            *sys_args.addTrueInfo.split(","),
        ]:
            if file is not None:
                subprocess.call(
                    [
                        "cp",
                        file,
                        get_file_name(
                            sys_args.prefix,
                            file,
                        ),
                    ]
                )
        return

    assert target_var_num is not None
    unsubsetable_sites_num_to_keep: int = int(
        len(unsubsetable_idx) * target_var_num / len(subsetable_idx)
    )

    rng = np.random.default_rng(sys_args.seed)
    idx_to_be_removed: list[int] = sorted(
        [
            *rng.choice(
                a=subsetable_idx,
                size=len(subsetable_idx) - target_var_num,
                replace=False,
            ),
            *rng.choice(
                a=unsubsetable_idx,
                size=len(unsubsetable_idx) - unsubsetable_sites_num_to_keep,
                replace=False,
            ),
        ]
    )
    print(f"{len(idx_to_be_removed)} variant sites will be removed.")

    keys_to_be_removed: list[tuple[str, int]] = [
        true_var_sites[1][i][0] for i in idx_to_be_removed
    ]

    for i in reads:
        with open(
            os.path.join(
                sys_args.prefix,
                i[0],
            ),
            "w",
            encoding="utf-8",
        ) as fh:
            for key in sorted(i[1].keys()):
                if key not in keys_to_be_removed:
                    fh.write(i[1][key])
                    fh.write("\n")

    with open(
        get_file_name(
            sys_args.prefix,
            sys_args.trueVarSites,
        ),
        "w",
        encoding="utf-8",
    ) as fh:
        for i in true_var_sites[0]:
            fh.write(i)
            fh.write("\n")
        for key in sorted(true_var_sites[1].keys()):
            if key not in idx_to_be_removed:
                fh.write(true_var_sites[1][key][2])
                fh.write("\n")

    for file in sys_args.addTrueInfo.split(","):
        if not os.path.isfile(file):
            raise ValueError(f"Error! File {file} does not exist.")
        __comments: list[str] = []
        __contents: list[str] = []

        with open(
            file,
            "r",
            encoding="utf-8",
        ) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#"):
                    __comments.append(line)
                else:
                    __contents.append(line)

        with open(
            get_file_name(
                sys_args.prefix,
                file,
            ),
            "w",
            encoding="utf-8",
        ) as fh:
            if len(__comments) > 0:
                for i in __comments:
                    fh.write(i)
                    fh.write("\n")
            for i in range(len(__contents)):
                if i not in idx_to_be_removed:
                    fh.write(__contents[i])
                    fh.write("\n")


if __name__ == "__main__":
    main()
