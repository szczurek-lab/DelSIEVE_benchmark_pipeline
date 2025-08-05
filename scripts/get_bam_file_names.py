import argparse
import os
import re

import pandas as pd


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(
        description="Collect dummy bam file names for MonoVar and SCIPhI."
    )
    sys_args_parser.add_argument(
        "--ibam",
        required=True,
        type=str,
        help="a directory to simulated dummy bam files",
    )
    sys_args_parser.add_argument(
        "--icell",
        required=True,
        type=str,
        help="cell names (the order is relevant)",
    )
    sys_args_parser.add_argument(
        "--type",
        required=True,
        type=int,
        choices=[1, 2],
        help="for MonoVar (1) or SCIPhI (2)",
    )
    sys_args_parser.add_argument(
        "--prefix",
        type=str,
        default="tumcell",
        help="prefix of bam files for tumor cells",
    )
    sys_args_parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="output file for bam file names",
    )
    args = sys_args_parser.parse_args()

    if not os.path.isdir(args.ibam):
        raise ValueError("Directory for dummy bam files does not exist.")

    if not os.path.isfile(args.icell):
        raise ValueError("Input file for cell names does not exist.")

    os.makedirs(os.path.dirname(args.o), exist_ok=True)

    if os.path.isfile(args.o):
        print(f"Output file {args.o} exists; will be overwritten.")

    return args


def get_cell_names(file_name: str) -> list[tuple[str, int]]:
    ret: list[tuple[str, int]] = []
    cell_names: pd.DataFrame = pd.read_csv(file_name)
    for i in range(cell_names.shape[0]):
        ret.append((cell_names.iloc[i, 0], i))
    return ret


def get_dummy_bam_file_names(dir: str, prefix: str) -> list[str]:
    ret: list[str] = list()

    for dir_path, _, file_names in os.walk(dir):
        for file_name in file_names:
            file_name = file_name.strip()
            if file_name.startswith(prefix):
                ret.append(os.path.join(dir_path, file_name))
    return ret


def match_file_names(
    ref: list[tuple[str, int]],
    target: list[str],
) -> list[tuple[str, str, int]]:
    ret: list[tuple[str, str, int]] = list()
    for ref_name, id in ref:
        _ret: list[str] = list(
            filter(re.compile(f"{ref_name}[^a-zA-Z0-9]").search, target)
        )
        if len(_ret) == 1:
            ret.append((*_ret, ref_name, id))
        elif len(_ret) > 1:
            raise ValueError(f"Error! More than one cell {ref} found: {_ret}")
    return ret


def write_to_output(
    file_names: list[tuple[str, str, int]],
    output: str,
    type: int,
) -> None:
    with open(
        output,
        "w",
        encoding="utf-8",
    ) as fh:
        for file_name, _, _ in file_names:
            fh.write(file_name)

            if type == 2:
                fh.write("\tCT")
            fh.write("\n")


def main() -> None:
    args: argparse.Namespace = parse_system_args()
    cell_names: list[tuple[str, int]] = get_cell_names(args.icell)
    bam_file_names: list[str] = get_dummy_bam_file_names(
        args.ibam,
        args.prefix,
    )
    matched_bam_files: list[tuple[str, str, int]] = match_file_names(
        cell_names,
        bam_file_names,
    )
    write_to_output(
        sorted(matched_bam_files, key=lambda x: x[2]),
        args.o,
        args.type,
    )


if __name__ == "__main__":
    main()
