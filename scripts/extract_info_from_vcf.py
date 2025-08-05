import argparse
import os
import re
from enum import Enum, IntEnum, unique
from typing import Iterable, Union

TERNARY_COMMENT = """#-3: .
#-2: 1/.
#-1: 0/.
#0: 0/0
#1: 0/1
#2: 1/1
#3: 1/2
#4: 000
#5: 001
#6: 011
#7: 012
#8: 111
#9: 112
#10: 123
#11: deletions in general
#12: amplifications in general
"""


@unique
class GenotypeProbComment(Enum):
    UNIV = "#Genotypes: 0/0, 0/1, 1/1\n"
    SCIPHI_NO_PROB = "#Genotypes: 0/0, 0/1\n"
    SCIPHIN = "#Genotypes: 0/0, 0/1, 1/1, 1/-, 0/-, parallel\n"


@unique
class ExtendedTernaryGenotype(IntEnum):
    MISSING = -3
    LOH_REF = -2
    LOH_ALT = -1
    HOMO_REF = 0
    HET = 1
    HOMO_ALT = 2
    HET_ALT = 3
    TRIPLE_000 = 4
    TRIPLE_001 = 5
    TRIPLE_011 = 6
    TRIPLE_012 = 7
    TRIPLE_111 = 8
    TRIPLE_112 = 9
    TRIPLE_123 = 10
    DEL = 11
    INS = 12


def parse_system_args() -> tuple[argparse.Namespace, dict[str, str]]:
    sys_args_parser = argparse.ArgumentParser(
        description="Extract information from vcf file."
    )
    sys_args_parser.add_argument(
        "input",
        type=str,
        help="a vcf file or a file list containing a vcf file (for sieve2)",
    )
    sys_args_parser.add_argument(
        "--mode",
        type=str,
        choices=["monovar", "sciphi", "sciphin", "sieve1", "sieve2"],
        help="which type of the input vcf file is? For monovar and sieve, a standard vcf file is the input; for sciphi and sciphin, the input vcf file has no PL filed. For sieve1, --input would be vcf; for sieve2, --input should be a file list containing a vcf file name.",
    )
    sys_args_parser.add_argument(
        "--prefix",
        type=str,
        help="prefix of output files; if not provided, prefix of input file will be used",
    )
    sys_args_parser.add_argument(
        "--probs",
        type=str,
        default=None,
        help="a .probs file should be specified when --mode is sciphi or sciphin",
    )
    sys_args_parser.add_argument(
        "--sifit",
        type=str,
        choices=["n", "b", "t"],
        default="n",
        help="save input data to SIFIT: binary (b), ternary (t), or not saving (n)?",
    )
    sys_args = sys_args_parser.parse_args()

    if (
        sys_args.mode == "monovar" or re.match(r"sieve\d", sys_args.mode)
    ) and sys_args.probs is not None:
        raise ValueError(
            "Error! When in mode 'monovar', 'sieve1' or 'sieve2', no --probs file should be specified."
        )
    if re.match(r"sieve\d", sys_args.mode) and sys_args.sifit == "n":
        raise ValueError(
            "Error! When in mode 'sieve1' or 'sieve2', --sifit should be set to 'b' or 't'."
        )
    if not os.path.isfile(sys_args.input):
        raise ValueError(f"Error! File {sys_args.input} does not exist.")

    if sys_args.prefix is not None:
        prefix = sys_args.prefix
    else:
        prefix = os.path.splitext(sys_args.input)[0]

    output_files = {
        "cell_names": prefix + ".cell_names",
        "loci_info": prefix + ".loci_info",
        "genotypes": prefix + ".genotypes",
        "ternary": prefix + ".ternary",
        "genotype_probs": prefix + ".genotype_probs",
        "sifit_data": prefix + ".sifit_data",
        "sifit_cell_names": prefix + ".sifit_cell_names",
    }
    return sys_args, output_files


def is_all_at_least_one(values: list[str]) -> bool:
    if values is None or len(values) == 0:
        return False

    for i in values:
        try:
            value = int(i)
        except ValueError:
            return False

        if value < 1:
            return False

    return True


def is_zero_and_more(values: list[str]) -> bool:
    if values is None or len(values) == 0:
        return False

    results: list[int] = []
    for i in values:
        try:
            value = int(i)
        except ValueError:
            return False

        results.append(value)

    results.sort()

    if results[0] == 0 and results[-1] > 0:
        return True
    else:
        return False


"""
-3: .
-2: 1/.
-1: 0/.
0: 0/0
1: 0/1
2: 1/1
3: 1/2
4: 000
5: 001
6: 011
7: 012
8: 111
9: 112
10: 123
11: deletions in general
12: amplifications in general
"""


def genotype2ternary(
    genotype: str, allele2nuc_map: Union[dict[int, str], None]
) -> Union[int, None]:
    if re.match(r"([.-])([/|]\1)?", genotype):
        return ExtendedTernaryGenotype.MISSING
    elif re.match(r"^([.-]|[1-9]\d*)[/|](?!\1$)([.-]|[1-9]\d*)$", genotype):
        return ExtendedTernaryGenotype.LOH_REF
    elif re.match(r"^([.-]|0)[/|](?!\1$)([.-]|0)$", genotype):
        return ExtendedTernaryGenotype.LOH_ALT
    elif re.match(r"^0[/|]0$", genotype):
        return ExtendedTernaryGenotype.HOMO_REF
    elif genotype == "000":
        return ExtendedTernaryGenotype.TRIPLE_000
    elif genotype == "001":
        return ExtendedTernaryGenotype.TRIPLE_001
    elif genotype == "011":
        return ExtendedTernaryGenotype.TRIPLE_011
    elif genotype == "012":
        return ExtendedTernaryGenotype.TRIPLE_012
    elif genotype == "111":
        return ExtendedTernaryGenotype.TRIPLE_111
    elif genotype == "112":
        return ExtendedTernaryGenotype.TRIPLE_112
    elif genotype == "123":
        return ExtendedTernaryGenotype.TRIPLE_123

    if allele2nuc_map is None or all(
        [len(val) == 1 for _, val in allele2nuc_map.items()]
    ):
        if re.match(r"(?=^0|0$)(\d+)[|/](?!\1)\d+$", genotype):
            return ExtendedTernaryGenotype.HET
        elif re.match(r"^([1-9]\d*)[/|]\1$", genotype):
            return ExtendedTernaryGenotype.HOMO_ALT
        elif re.match(r"^([1-9]\d*)[/|](?!\1$)[1-9]\d*$", genotype):
            return ExtendedTernaryGenotype.HET_ALT
    else:
        raise ValueError(
            f"Error! Unhandled genotype {genotype} and nucleotides {allele2nuc_map}."
        )


def ternary_loose2strict(ter: ExtendedTernaryGenotype) -> int:
    if ter == ExtendedTernaryGenotype.MISSING:
        return 3
    elif (
        ter == ExtendedTernaryGenotype.LOH_ALT
        or ter == ExtendedTernaryGenotype.HOMO_REF
        or ter == ExtendedTernaryGenotype.TRIPLE_000
    ):
        return 0
    elif (
        ter == ExtendedTernaryGenotype.LOH_REF
        or ter == ExtendedTernaryGenotype.HOMO_ALT
        or ter == ExtendedTernaryGenotype.HET_ALT
        or ter == ExtendedTernaryGenotype.TRIPLE_001
        or ter == ExtendedTernaryGenotype.TRIPLE_011
        or ter == ExtendedTernaryGenotype.TRIPLE_012
    ):
        return 2
    else:
        return 1


def get_genotypes_probs(pls: Iterable) -> tuple[float]:
    gls: list[float] = [pow(10, -1 * i / 10) for i in pls]
    __sum: float = sum(gls)
    return tuple(
        gls[i] / __sum if i < len(gls) - 1 else 1 - sum(gls[:-1]) / __sum
        for i in range(len(gls))
    )


def get_sciphi_gt(index: int) -> Union[str, None]:
    if index == 0:
        return "0/0"
    elif index == 1:
        return "0/1"
    elif index == 2:
        return "1/1"
    else:
        raise ValueError(f"Error! Unknown index {index}.")


def get_sciphin_gt(index: int) -> Union[str, None]:
    if index == 0:
        return "0/0"
    elif index == 1:
        return "0/1"
    elif index == 2:
        return "1/1"
    elif index == 3:
        return "1/."
    elif index == 4:
        return "0/."
    elif index == 5:
        return "0/1"
    else:
        raise ValueError(f"Error! Unknown genotype index {index}.")


def get_vcf_file(input_file: str) -> str:
    with open(input_file, "r") as fh:
        for line in fh:
            if line.strip().endswith(".vcf"):
                return os.path.join(os.path.dirname(input_file), line.strip())
    raise ValueError("Error! File " + input_file + " does not contain a vcf file.")


def parse_vcf_file(input_file: str, mode: str, probs_file: Union[str, None]) -> tuple[
    list[str],
    list[tuple[str, int]],
    dict[tuple[str, int], list[str]],
    dict[tuple[str, int], list[ExtendedTernaryGenotype]],
    dict[tuple[str, int], list[tuple[float]]],
]:
    tmp_input_file: str = input_file

    if mode == "sieve2":
        tmp_input_file = get_vcf_file(input_file)

    if os.path.isfile(tmp_input_file):
        cell_names = []
        loci_info = []
        genotypes = {}
        ternary = {}
        genotype_probs = {}
        allele2nuc_map = {}

        chrom_index = pos_index = ref_index = alt_index = gt_index = pl_index = -1
        cell_info_index = 0

        with open(tmp_input_file, "r") as fh:
            for line in fh:
                if line.strip().startswith("#"):
                    if line.strip().startswith("#CHROM"):
                        comp = line.strip().split("\t")
                        reach_cells = False
                        for item in comp:
                            if reach_cells is True:
                                cell_name = item.strip().replace(".bam", "")
                                cell_names.append(cell_name)
                            else:
                                _item = item.strip().upper()
                                if "CHROM" in _item:
                                    chrom_index = cell_info_index
                                elif _item == "POS":
                                    pos_index = cell_info_index
                                elif _item == "REF":
                                    ref_index = cell_info_index
                                elif _item == "ALT":
                                    alt_index = cell_info_index
                                elif _item == "FORMAT":
                                    reach_cells = True
                                cell_info_index += 1
                else:
                    duplicate_pos = False
                    comp = line.strip().split("\t")

                    __loci_info = (
                        comp[chrom_index].strip(),
                        int(comp[pos_index].strip()),
                    )

                    if __loci_info not in loci_info:
                        loci_info.append(__loci_info)
                        genotypes[__loci_info] = ["0/0" for _ in cell_names]
                        ternary[__loci_info] = [0 for _ in cell_names]

                        if mode.lower() == "monovar":
                            genotype_probs[__loci_info] = [
                                (1.0, 0.0, 0.0) for _ in cell_names
                            ]
                        elif (
                            mode.lower() in ["sciphi", "sciphin"] and probs_file is None
                        ):
                            genotype_probs[__loci_info] = [
                                (1.0, 0.0) for _ in cell_names
                            ]
                    else:
                        duplicate_pos = True

                    allele2nuc_map[__loci_info] = {}
                    allele2nuc_map[__loci_info].update(
                        (
                            (idx, i)
                            for idx, i in enumerate(
                                [
                                    comp[ref_index].strip(),
                                    *comp[alt_index].strip().split(","),
                                ]
                            )
                        )
                    )
                    if len(allele2nuc_map[__loci_info][0]) == 0:
                        raise ValueError(
                            f"Error! Reference allele is empty for chromosome {loci_info[-1][0]} at position {loci_info[-1][1]}."
                        )

                    if mode.lower() == "monovar" and pl_index == -1:
                        pl_index = 0
                        _comp = comp[cell_info_index - 1].split(":")
                        for item in _comp:
                            if not item.strip().upper() == "PL":
                                pl_index += 1
                            else:
                                break

                    if gt_index == -1:
                        gt_index = 0
                        _comp = comp[cell_info_index - 1].split(":")
                        for item in _comp:
                            if not item.strip().upper() == "GT":
                                gt_index += 1
                            else:
                                break

                    for cell_index in range(0, len(cell_names)):
                        cell_info = (
                            comp[cell_info_index + cell_index].strip().split(":")
                        )

                        if not duplicate_pos or genotypes[__loci_info][cell_index] in [
                            "0/0",
                            "0|0",
                        ]:
                            genotypes[__loci_info][cell_index] = cell_info[
                                gt_index
                            ].strip()

                            _ternary = genotype2ternary(
                                cell_info[gt_index].strip(), allele2nuc_map[__loci_info]
                            )
                            ternary[__loci_info][cell_index] = _ternary

                            if mode.lower() == "monovar":
                                if _ternary == ExtendedTernaryGenotype.MISSING:
                                    genotype_probs[__loci_info][cell_index] = (
                                        0.0,
                                        0.0,
                                        0.0,
                                    )
                                else:
                                    genotype_probs[__loci_info][cell_index] = (
                                        get_genotypes_probs(
                                            [
                                                int(i.strip())
                                                for i in cell_info[pl_index]
                                                .strip()
                                                .split(",")
                                            ]
                                        )
                                    )
                            elif (
                                mode.lower() in ["sciphi", "sciphin"]
                                and probs_file is None
                                and genotypes[__loci_info][cell_index]
                                not in [
                                    "0/0",
                                    "0|0",
                                ]
                            ):
                                genotype_probs[__loci_info][cell_index] = (0.0, 1.0)

        if mode.lower() in ["sciphi", "sciphin"] and probs_file is not None:
            assert os.path.isfile(
                probs_file
            ), f"Error! File {probs_file} does not exist."

            unprocessed_loci_info = loci_info.copy()

            cell_index_in_vcf = []
            with open(probs_file, "r") as fh:
                for line in fh:
                    if line.strip().startswith("chrom"):
                        comp = line.strip().split("\t")

                        if len(comp) - 2 != len(cell_names):
                            raise ValueError(
                                "Error! The number of cells in vcf file is different from that in probs file."
                            )

                        for index in range(2, len(comp)):
                            __cell_name = comp[index].strip().replace(".bam", "")
                            cell_index_in_vcf.append(cell_names.index(__cell_name))
                    else:
                        comp = line.strip().split("\t")
                        __loci_info = (comp[0].strip(), int(comp[1].strip()))

                        if __loci_info in unprocessed_loci_info:
                            unprocessed_loci_info.remove(__loci_info)
                        else:
                            continue

                        for index in range(2, len(comp)):
                            _comp = comp[index].strip().split("|")

                            if mode.lower() == "sciphi":
                                assert (
                                    len(_comp) == 3
                                ), f"Error! The number of probabilities for cell {cell_names[cell_index_in_vcf[index - 2]]} at {__loci_info[0]}, {__loci_info[1]} is not 3."

                                # 0/0, 0/1, 1/1
                                genotype_probs[__loci_info][
                                    cell_index_in_vcf[index - 2]
                                ] = (
                                    1 - float(_comp[2].strip()),
                                    float(_comp[0].strip()),
                                    float(_comp[1].strip()),
                                )
                            elif mode.lower() == "sciphin":
                                assert (
                                    len(_comp) == 6
                                ), f"Error! The number of probabilities for cell {cell_names[cell_index_in_vcf[index - 2]]} at {__loci_info[0]}, {__loci_info[1]} is not 6."

                                # 0/0, 0/1, 1/1, 1/-, 0/-, parallel
                                genotype_probs[__loci_info][
                                    cell_index_in_vcf[index - 2]
                                ] = (
                                    1 - float(_comp[5].strip()),
                                    float(_comp[0].strip()),
                                    float(_comp[1].strip()),
                                    float(_comp[2].strip()),
                                    float(_comp[3].strip()),
                                    float(_comp[4].strip()),
                                )
                            else:
                                raise ValueError(f"Error! Unknown mode {mode}.")

                            if (
                                genotype_probs[__loci_info][
                                    cell_index_in_vcf[index - 2]
                                ]
                                is not None
                            ):
                                if (
                                    ternary[__loci_info][cell_index_in_vcf[index - 2]]
                                    != 0
                                ):
                                    if mode.lower() == "sciphi":
                                        genotypes[__loci_info][
                                            cell_index_in_vcf[index - 2]
                                        ] = get_sciphi_gt(
                                            genotype_probs[__loci_info][
                                                cell_index_in_vcf[index - 2]
                                            ].index(
                                                max(
                                                    genotype_probs[__loci_info][
                                                        cell_index_in_vcf[index - 2]
                                                    ]
                                                )
                                            )
                                        )
                                        ternary[__loci_info][
                                            cell_index_in_vcf[index - 2]
                                        ] = genotype2ternary(
                                            genotypes[__loci_info][
                                                cell_index_in_vcf[index - 2]
                                            ],
                                            None,
                                        )
                                    elif mode.lower() == "sciphin":
                                        genotypes[__loci_info][
                                            cell_index_in_vcf[index - 2]
                                        ] = get_sciphin_gt(
                                            genotype_probs[__loci_info][
                                                cell_index_in_vcf[index - 2]
                                            ].index(
                                                max(
                                                    genotype_probs[__loci_info][
                                                        cell_index_in_vcf[index - 2]
                                                    ]
                                                )
                                            )
                                        )
                                        ternary[__loci_info][
                                            cell_index_in_vcf[index - 2]
                                        ] = genotype2ternary(
                                            genotypes[__loci_info][
                                                cell_index_in_vcf[index - 2]
                                            ],
                                            None,
                                        )

        return cell_names, loci_info, genotypes, ternary, genotype_probs


def print_cell_names(out_file: str, cell_names: list[str]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, "w") as fh:
        for index in range(0, len(cell_names)):
            fh.write(cell_names[index])
            if index < len(cell_names) - 1:
                fh.write("\n")


def print_loci_info(out_file: str, loci_info: list[tuple[str, int]]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, "w") as fh:
        for index in range(0, len(loci_info)):
            fh.write(loci_info[index][0] + "\t" + str(loci_info[index][1]))
            if index < len(loci_info) - 1:
                fh.write("\n")


def print_genotypes(
    out_file: str,
    loci_info: list[tuple[str, int]],
    cell_names: list[str],
    genotypes: dict[tuple[str, int], list[str]],
):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, "w") as fh:
        for loci_index in range(0, len(loci_info)):
            for cell_index in range(0, len(cell_names)):
                fh.write(genotypes[loci_info[loci_index]][cell_index])
                if cell_index < len(cell_names) - 1:
                    fh.write("\t")
            if loci_index < len(loci_info) - 1:
                fh.write("\n")


def print_ternary(
    out_file: str,
    loci_info: list[tuple[str, int]],
    cell_names: list[str],
    ternary: dict[tuple[str, int], list[ExtendedTernaryGenotype]],
):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, "w") as fh:
        fh.write(TERNARY_COMMENT)
        for loci_index in range(0, len(loci_info)):
            for cell_index in range(0, len(cell_names)):
                fh.write(str(ternary[loci_info[loci_index]][cell_index].value))
                if cell_index < len(cell_names) - 1:
                    fh.write("\t")
            if loci_index < len(loci_info) - 1:
                fh.write("\n")


def print_genotype_probs(
    out_file: str,
    loci_info: list[tuple[str, int]],
    cell_names: list[str],
    genotype_probs: dict[tuple[str, int], list[tuple[float]]],
    mode: Union[str, None] = None,
    provide_genotype_probs: bool = False,
):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, "w") as fh:
        if mode is not None:
            if mode == "sciphin" and provide_genotype_probs:
                fh.write(GenotypeProbComment.SCIPHIN.value)
            if mode == "sciphin" and not provide_genotype_probs:
                fh.write(GenotypeProbComment.SCIPHI_NO_PROB.value)
            else:
                fh.write(GenotypeProbComment.UNIV.value)

        for loci_index in range(0, len(loci_info)):
            for cell_index in range(0, len(cell_names)):
                probs = genotype_probs[loci_info[loci_index]][cell_index]
                for genotype_index in range(0, len(probs)):
                    fh.write(str(probs[genotype_index]))
                    if genotype_index < len(probs) - 1:
                        fh.write(",")
                if cell_index < len(cell_names) - 1:
                    fh.write("\t")
            if loci_index < len(loci_info) - 1:
                fh.write("\n")


def print_for_sifit(
    out_data_file: str,
    out_cell_names_file: str,
    mode: str,
    loci_info: list[tuple[str, int]],
    cell_names: list[str],
    ternary: dict[tuple[str, int], list[ExtendedTernaryGenotype]],
):
    if not os.path.exists(os.path.dirname(out_data_file)):
        os.makedirs(os.path.dirname(out_data_file))
    with open(out_data_file, "w") as fh:
        for loci_index in range(0, len(loci_info)):
            fh.write(str(loci_info[loci_index][1]))
            for cell_index in range(0, len(cell_names)):
                fh.write(" ")
                code = ternary_loose2strict(ternary[loci_info[loci_index]][cell_index])
                if mode == "b" and code == 2:
                    code = 1
                fh.write(str(code))
            if loci_index < len(loci_info) - 1:
                fh.write("\n")

    if not os.path.exists(os.path.dirname(out_cell_names_file)):
        os.makedirs(os.path.dirname(out_cell_names_file))
    with open(out_cell_names_file, "w") as fh:
        for index in range(0, len(cell_names)):
            fh.write(cell_names[index])
            if index < len(cell_names) - 1:
                fh.write(" ")


def main():
    sys_args, output_files = parse_system_args()
    cell_names, loci_info, genotypes, ternary, genotype_probs = parse_vcf_file(
        sys_args.input, sys_args.mode, sys_args.probs
    )
    if sys_args.mode.lower() in ["monovar", "sciphi", "sciphin"]:
        print_cell_names(output_files["cell_names"], cell_names)
        print_loci_info(output_files["loci_info"], loci_info)
        print_genotypes(output_files["genotypes"], loci_info, cell_names, genotypes)
        print_ternary(output_files["ternary"], loci_info, cell_names, ternary)
        print_genotype_probs(
            output_files["genotype_probs"],
            loci_info,
            cell_names,
            genotype_probs,
            sys_args.mode.lower(),
            sys_args.probs,
        )

    if sys_args.sifit != "n":
        print_for_sifit(
            output_files["sifit_data"],
            output_files["sifit_cell_names"],
            sys_args.sifit,
            loci_info,
            cell_names,
            ternary,
        )


if __name__ == "__main__":
    main()
