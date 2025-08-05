import sys
import pandas as pd
import numpy as np
import os
from typing import Union

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from scripts.utils import (
    NON_WORDS_PATTERN,
    which_patterns,
    vec_pattern,
    vec_pattern_2,
    get_matched_files,
    get_cnv_rate,
)


def _get_content_files(file: str) -> Union[list[str], None]:
    if os.path.isfile(file):
        results = []
        root = os.path.abspath(os.path.dirname(file))
        with open(file, "r") as fh:
            for line in fh:
                content = os.path.join(root, line.strip())
                if content != os.path.abspath(file) and os.path.isfile(content):
                    results.append(content)
        return results
    return None


vec_get_content_files = np.vectorize(_get_content_files, otypes=[list])


def _get_cell_names(input: str, is_true: bool) -> Union[pd.DataFrame, None]:
    if os.path.isfile(input):
        tmp: pd.DataFrame = pd.read_csv(
            input,
            header=0 if is_true else None,
        )

        if is_true:
            return (
                tmp.assign(
                    is_doublet=lambda x: ~x.doublet_with.isna(),
                )
                .filter(
                    ["cell", "is_doublet"],
                )
                .reset_index(
                    names="index_in_true",
                )
            )
        else:
            return tmp.rename(
                columns={0: "cell"},
            ).reset_index(
                names="index_in_pred",
            )
    return None


def _get_loci_info(
    input: str, is_truth: bool = False
) -> Union[tuple[list[tuple[str, int]], Union[dict[str, list[bool]], None]], None]:
    if os.path.isfile(input):
        ret = []
        ret_truth = {"mu": [], "del": [], "amp": []}
        with open(input, "r") as fh:
            for line in fh:
                if line.strip().startswith("#"):
                    continue
                comp = line.strip().split()
                ret.append((comp[0].strip(), int(comp[1].strip())))
                if is_truth:
                    ret_truth["mu"].append(bool(comp[2].strip()))
                    ret_truth["del"].append(bool(comp[3].strip()))
                    ret_truth["amp"].append(bool(comp[4].strip()))
        return ret, ret_truth if is_truth else None
    return None


def _get_input_data(input: str) -> list[tuple[str, int]] | None:
    if os.path.isfile(input):
        ret = []
        with open(input, "r") as fh:
            is_start: bool = False
            for line in fh:
                if line:
                    if line.strip() == "=mutations=":
                        is_start = True
                    elif line.strip() == "=background=":
                        is_start = False
                        break
                    elif is_start:
                        comp: list[str] = line.strip().split()
                        ret.append((comp[0].strip(), int(comp[1].strip())))
        return ret
    return None


def _get_ternary(row_num: int, col_num: int, input: str) -> Union[np.ndarray, None]:
    if os.path.isfile(input):
        results = np.empty((row_num, col_num), dtype=np.int32)
        with open(input, "r") as fh:
            row_index = 0
            for line in fh:
                if line.strip().startswith("#"):
                    continue
                if row_num <= row_index:
                    raise ValueError(
                        "Error! Unmatched size. Expecting "
                        + str(row_num)
                        + " rows, detecting "
                        + str(row_index + 1)
                        + " rows."
                    )
                comp = line.strip().split()
                if col_num != len(comp):
                    raise ValueError(
                        "Error! Unmatched size. Expecting "
                        + str(col_num)
                        + " columns, detecting "
                        + str(len(comp))
                        + " columns."
                    )
                for col_index in range(0, len(comp)):
                    results[row_index, col_index] = int(comp[col_index])
                row_index += 1
        return results
    return None


"""
mode: 0 - mu, 1 - del, 2 - amp

mode = 0:
-3: -
-2: 1/-
-1: 0/-
0: 0/0
1: 0/1
2: 1/1
3: 1/2 (namely 1/1')
4: 000
5: 001
6: 011
7: 012
8: 111
9: 112
10: 123
"""


def _get_strict_ternary(ternary: int, mode: int = 0) -> Union[None, int]:
    if mode == -1:
        if ternary == 0:
            return 0
        else:
            return 1
    elif mode == 0:
        if ternary < 0:
            return ternary
        elif ternary in [0, 4]:
            return 0
        elif ternary in [1, 5, 6, 7]:
            return 1
        elif ternary in [2, 3, 8, 9, 10]:
            return 2
        else:
            raise ValueError(f"Error! Unmatched ternary value: {ternary}")
    elif mode == 1:
        pass
    elif mode == 2:
        if ternary < 4:
            return 0
        else:
            return 1
    else:
        raise ValueError(f"Error! Unknown mode: {mode}")
    return None


def _get_properties_all(
    loci_info: list[tuple[str, int]],
    true_loci_info: list[tuple[str, int]],
    cell_names: pd.DataFrame,
    ternary_value: np.ndarray,
    true_ternary_value: np.ndarray,
):
    detected_true_loci_info = []
    tp = fp = tn = fn = 0

    for _loci_info_index in range(0, len(loci_info)):
        if loci_info[_loci_info_index] in true_loci_info:
            detected_true_loci_info.append(loci_info[_loci_info_index])
            _true_loci_info_index = true_loci_info.index(loci_info[_loci_info_index])

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[_true_loci_info_index, row.index_in_true],
                    -1,
                )
                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred], -1
                )

                if _true_ternary_value == 0:
                    if _ternary_value == 0:
                        tn += 1
                    elif _ternary_value == 1:
                        fp += 1
                elif _true_ternary_value == 1:
                    if _ternary_value == 0:
                        fn += 1
                    elif _ternary_value == 1:
                        tp += 1
        else:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred], -1
                )
                if _ternary_value == 0:
                    tn += 1
                elif _ternary_value == 1:
                    fp += 1

    for _true_loci_info_index in range(0, len(true_loci_info)):
        if true_loci_info[_true_loci_info_index] not in detected_true_loci_info:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[_true_loci_info_index, row.index_in_true], -1
                )
                if _true_ternary_value == 0:
                    tn += 1
                elif _true_ternary_value == 1:
                    fn += 1

    return (
        tp,
        fp,
        tn,
        fn,
    )


def _get_properties_mu(
    loci_info: list[tuple[str, int]],
    true_loci_info: list[tuple[str, int]],
    cell_names: pd.DataFrame,
    ternary_value: np.ndarray,
    true_ternary_value: np.ndarray,
):
    detected_true_loci_info = []
    tp = fp = tn = fn = tp_hetero_mu = fp_hetero_mu = tn_hetero_mu = fn_hetero_mu = (
        tp_homo_mu
    ) = fp_homo_mu = tn_homo_mu = fn_homo_mu = t_homo_ref_as_hetero_mu = (
        t_homo_mu_as_homo_ref
    ) = t_homo_mu_as_hetero_mu = 0

    for _loci_info_index in range(0, len(loci_info)):
        if loci_info[_loci_info_index] in true_loci_info:
            detected_true_loci_info.append(loci_info[_loci_info_index])
            _true_loci_info_index = true_loci_info.index(loci_info[_loci_info_index])

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[_true_loci_info_index, row.index_in_true]
                )
                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred]
                )

                if _true_ternary_value in [0, -1, -3]:
                    if _ternary_value in [0, -1, -3]:
                        tn += 1
                        tn_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 1:
                        fp += 1
                        fp_hetero_mu += 1
                        tn_homo_mu += 1
                        if _true_ternary_value == 0:
                            t_homo_ref_as_hetero_mu += 1
                    elif _ternary_value == 2:
                        fp += 1
                        tn_hetero_mu += 1
                        fp_homo_mu += 1
                    elif _ternary_value == -2:
                        fp += 1
                        tn_hetero_mu += 1
                        tn_homo_mu += 1
                elif _true_ternary_value == 1:
                    if _ternary_value in [0, -1, -3]:
                        fn += 1
                        fn_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 1:
                        tp += 1
                        tp_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 2:
                        tp += 1
                        fn_hetero_mu += 1
                        fp_homo_mu += 1
                    elif _ternary_value == -2:
                        tp += 1
                        fn_hetero_mu += 1
                        tn_homo_mu += 1
                elif _true_ternary_value == 2:
                    if _ternary_value in [0, -1, -3]:
                        fn += 1
                        tn_hetero_mu += 1
                        fn_homo_mu += 1
                        if _ternary_value == 0:
                            t_homo_mu_as_homo_ref += 1
                    elif _ternary_value == 1:
                        tp += 1
                        fp_hetero_mu += 1
                        fn_homo_mu += 1
                        t_homo_mu_as_hetero_mu += 1
                    elif _ternary_value == 2:
                        tp += 1
                        tn_hetero_mu += 1
                        tp_homo_mu += 1
                    elif _ternary_value == -2:
                        tp += 1
                        fn_hetero_mu += 1
                        tn_homo_mu += 1
                elif _true_ternary_value == -2:
                    if _ternary_value in [0, -1, -3]:
                        fn += 1
                        tn_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 1:
                        tp += 1
                        fp_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 2:
                        tp += 1
                        tn_hetero_mu += 1
                        fp_homo_mu += 1
                    elif _ternary_value == -2:
                        tp += 1
                        tn_hetero_mu += 1
                        tn_homo_mu += 1
        else:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred]
                )
                if _ternary_value in [0, -1, -3]:
                    tn += 1
                    tn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _ternary_value == 1:
                    fp += 1
                    fp_hetero_mu += 1
                    tn_homo_mu += 1
                    t_homo_ref_as_hetero_mu += 1
                elif _ternary_value == 2:
                    fp += 1
                    tn_hetero_mu += 1
                    fp_homo_mu += 1

    for _true_loci_info_index in range(0, len(true_loci_info)):
        if true_loci_info[_true_loci_info_index] not in detected_true_loci_info:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[_true_loci_info_index, row.index_in_true]
                )
                if _true_ternary_value in [0, -1, -3]:
                    tn += 1
                    tn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _true_ternary_value == 1:
                    fn += 1
                    fn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _true_ternary_value == 2:
                    fn += 1
                    tn_hetero_mu += 1
                    fn_homo_mu += 1
                    t_homo_mu_as_homo_ref += 1

    return (
        tp,
        fp,
        tn,
        fn,
        tp_hetero_mu,
        fp_hetero_mu,
        tn_hetero_mu,
        fn_hetero_mu,
        tp_homo_mu,
        fp_homo_mu,
        tn_homo_mu,
        fn_homo_mu,
        t_homo_ref_as_hetero_mu,
        t_homo_mu_as_homo_ref,
        t_homo_mu_as_hetero_mu,
    )


def _get_properties_del(
    loci_info: list[tuple[str, int]],
    input_data: None | list[tuple[str, int]],
    true_loci_info: list[tuple[str, int]],
    cell_names: pd.DataFrame,
    ternary_value: np.ndarray,
    true_ternary_value: np.ndarray,
):
    detected_true_loci_info = []
    tp = fp = tn = fn = tp_all = fp_all = tn_all = fn_all = tp_alt_left = (
        fp_alt_left
    ) = tn_alt_left = fn_alt_left = tp_ref_left = fp_ref_left = tn_ref_left = (
        fn_ref_left
    ) = t_ref_left_as_homo_ref = t_ref_left_as_hetero_mu = t_alt_left_as_hetero_mu = (
        t_alt_left_as_homo_mu
    ) = 0

    for _loci_info_index in range(0, len(loci_info)):
        if loci_info[_loci_info_index] in true_loci_info:
            detected_true_loci_info.append(loci_info[_loci_info_index])
            _true_loci_info_index = true_loci_info.index(loci_info[_loci_info_index])

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = true_ternary_value[
                    _true_loci_info_index, row.index_in_true
                ]
                _ternary_value = ternary_value[_loci_info_index, row.index_in_pred]

                if _true_ternary_value == -3:
                    if _ternary_value == -3:
                        tp += 1
                        tp_all += 1
                        tn_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -2:
                        tp += 1
                        fn_all += 1
                        fp_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -1:
                        tp += 1
                        fn_all += 1
                        tn_alt_left += 1
                        fp_ref_left += 1
                    else:
                        fn += 1
                        fn_all += 1
                        tn_alt_left += 1
                        tn_ref_left += 1
                elif _true_ternary_value == -2:
                    if _ternary_value == -3:
                        tp += 1
                        fp_all += 1
                        fn_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -2:
                        tp += 1
                        tn_all += 1
                        tp_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -1:
                        tp += 1
                        tn_all += 1
                        fn_alt_left += 1
                        fp_ref_left += 1
                    else:
                        fn += 1
                        tn_all += 1
                        fn_alt_left += 1
                        tn_ref_left += 1
                        if _ternary_value in [2, 3, 8, 9, 10]:
                            t_alt_left_as_homo_mu += 1
                        elif _ternary_value in [1, 5, 6, 7]:
                            t_alt_left_as_hetero_mu += 1
                elif _true_ternary_value == -1:
                    if _ternary_value == -3:
                        tp += 1
                        fp_all += 1
                        tn_alt_left += 1
                        fn_ref_left += 1
                    elif _ternary_value == -2:
                        tp += 1
                        tn_all += 1
                        fp_alt_left += 1
                        fn_ref_left += 1
                    elif _ternary_value == -1:
                        tp += 1
                        tn_all += 1
                        tn_alt_left += 1
                        tp_ref_left += 1
                    else:
                        fn += 1
                        tn_all += 1
                        tn_alt_left += 1
                        fn_ref_left += 1
                        if _ternary_value in [0, 4]:
                            t_ref_left_as_homo_ref += 1
                        elif _ternary_value in [1, 5, 6, 7]:
                            t_ref_left_as_hetero_mu += 1
                else:
                    if _ternary_value == -3:
                        fp += 1
                        fp_all += 1
                        tn_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -2:
                        fp += 1
                        tn_all += 1
                        fp_alt_left += 1
                        tn_ref_left += 1
                    elif _ternary_value == -1:
                        fp += 1
                        tn_all += 1
                        tn_alt_left += 1
                        fp_ref_left += 1
                    else:
                        tn += 1
                        tn_all += 1
                        tn_alt_left += 1
                        tn_ref_left += 1
        else:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _ternary_value = ternary_value[_loci_info_index, row.index_in_pred]
                if _ternary_value == -3:
                    fp += 1
                    fp_all += 1
                    tn_alt_left += 1
                    tn_ref_left += 1
                elif _ternary_value == -2:
                    fp += 1
                    tn_all += 1
                    fp_alt_left += 1
                    tn_ref_left += 1
                elif _ternary_value == -1:
                    fp += 1
                    tn_all += 1
                    tn_alt_left += 1
                    fp_ref_left += 1
                else:
                    tn += 1
                    tn_all += 1
                    tn_alt_left += 1
                    tn_ref_left += 1

    for _true_loci_info_index in range(
        0, len(true_loci_info if input_data is None else input_data)
    ):
        if input_data is None:
            __true_loci_info = true_loci_info[_true_loci_info_index]
        else:
            __true_loci_info = input_data[_true_loci_info_index]

        if __true_loci_info not in detected_true_loci_info:

            try:
                __true_loci_info_index = true_loci_info.index(__true_loci_info)
            except ValueError:
                continue

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = true_ternary_value[
                    __true_loci_info_index, row.index_in_true
                ]
                if _true_ternary_value == -3:
                    fn += 1
                    fn_all += 1
                    tn_alt_left += 1
                    tn_ref_left += 1
                elif _true_ternary_value == -2:
                    fn += 1
                    tn_all += 1
                    fn_alt_left += 1
                    tn_ref_left += 1
                elif _true_ternary_value == -1:
                    fn += 1
                    tn_all += 1
                    tn_alt_left += 1
                    fn_ref_left += 1
                    t_ref_left_as_homo_ref += 1
                else:
                    tn += 1
                    tn_all += 1
                    tn_alt_left += 1
                    tn_ref_left += 1

    return (
        tp,
        fp,
        tn,
        fn,
        tp_all,
        fp_all,
        tn_all,
        fn_all,
        tp_alt_left,
        fp_alt_left,
        tn_alt_left,
        fn_alt_left,
        tp_ref_left,
        fp_ref_left,
        tn_ref_left,
        fn_ref_left,
        t_ref_left_as_homo_ref,
        t_ref_left_as_hetero_mu,
        t_alt_left_as_hetero_mu,
        t_alt_left_as_homo_mu,
    )


def _get_properties_amp(
    loci_info: list[tuple[str, int]],
    input_data: None | list[tuple[str, int]],
    true_loci_info: list[tuple[str, int]],
    cell_names: pd.DataFrame,
    ternary_value: np.ndarray,
    true_ternary_value: np.ndarray,
):
    detected_true_loci_info = []
    tp = fp = tn = fn = 0

    for _loci_info_index in range(0, len(loci_info)):
        if loci_info[_loci_info_index] in true_loci_info:
            detected_true_loci_info.append(loci_info[_loci_info_index])
            _true_loci_info_index = true_loci_info.index(loci_info[_loci_info_index])

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[_true_loci_info_index, row.index_in_true],
                    2,
                )
                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred], 2
                )

                if _true_ternary_value == 0:
                    if _ternary_value == 0:
                        tn += 1
                    elif _ternary_value == 1:
                        fp += 1
                elif _true_ternary_value == 1:
                    if _ternary_value == 0:
                        fn += 1
                    elif _ternary_value == 1:
                        tp += 1
        else:
            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _ternary_value = _get_strict_ternary(
                    ternary_value[_loci_info_index, row.index_in_pred], 2
                )
                if _ternary_value == 0:
                    tn += 1
                elif _ternary_value == 1:
                    fp += 1

    for _true_loci_info_index in range(
        0, len(true_loci_info if input_data is None else input_data)
    ):
        if input_data is None:
            __true_loci_info = true_loci_info[_true_loci_info_index]
        else:
            __true_loci_info = input_data[_true_loci_info_index]

        if __true_loci_info not in detected_true_loci_info:

            try:
                __true_loci_info_index = true_loci_info.index(__true_loci_info)
            except ValueError:
                continue

            for _, row in cell_names.iterrows():
                if row.is_doublet:
                    continue

                _true_ternary_value = _get_strict_ternary(
                    true_ternary_value[__true_loci_info_index, row.index_in_true], 2
                )
                if _true_ternary_value == 0:
                    tn += 1
                elif _true_ternary_value == 1:
                    fn += 1

    return tp, fp, tn, fn


def _get_properties(
    true_cell_names: str,
    true_loci_info: str,
    true_ternary: str,
    cell_names: str,
    loci_info: str,
    input_data: None | str,
    ternary: str,
):
    _true_cell_names = _get_cell_names(true_cell_names, True)
    _cell_names = _get_cell_names(cell_names, False)
    assert _true_cell_names is not None and _cell_names is not None

    if _true_cell_names.shape[0] != _cell_names.shape[0]:
        raise ValueError(
            "Error! Unmatched number of cells in "
            + true_cell_names
            + " and "
            + cell_names
        )

    joined_cell_names = pd.merge(
        _cell_names,
        _true_cell_names,
        how="inner",
        on="cell",
    )
    assert joined_cell_names.shape[0] == _cell_names.shape[0]

    _true_loci_info, _true_evo_events = _get_loci_info(true_loci_info, True)
    _loci_info, _ = _get_loci_info(loci_info)
    assert (
        _true_loci_info is not None
        and _true_evo_events is not None
        and _loci_info is not None
    )

    _input_data = None
    if input_data is not None:
        _input_data = _get_input_data(input_data)

    _true_ternary = _get_ternary(
        len(_true_loci_info), joined_cell_names.shape[0], true_ternary
    )
    _ternary = _get_ternary(len(_loci_info), joined_cell_names.shape[0], ternary)
    assert _true_ternary is not None and _ternary is not None

    return (
        *_get_properties_all(
            _loci_info,
            [
                _true_loci_info[i]
                for i in range(len(_true_loci_info))
                if _true_evo_events["mu"][i]
            ],
            joined_cell_names,
            _ternary,
            _true_ternary,
        ),
        *_get_properties_mu(
            _loci_info,
            [
                _true_loci_info[i]
                for i in range(len(_true_loci_info))
                if _true_evo_events["mu"][i]
            ],
            joined_cell_names,
            _ternary,
            _true_ternary,
        ),
        len(_true_loci_info) * joined_cell_names.shape[0],
        *_get_properties_del(
            _loci_info,
            _input_data,
            [
                _true_loci_info[i]
                for i in range(len(_true_loci_info))
                if _true_evo_events["del"][i]
            ],
            joined_cell_names,
            _ternary,
            _true_ternary,
        ),
        *_get_properties_amp(
            _loci_info,
            _input_data,
            [
                _true_loci_info[i]
                for i in range(len(_true_loci_info))
                if _true_evo_events["amp"][i]
            ],
            joined_cell_names,
            _ternary,
            _true_ternary,
        ),
    )


def _get_metrics(
    true_positive: float,
    false_positive: float,
    true_negative: float,
    false_negative: float,
) -> tuple[float, float, float, float]:
    if true_positive + false_positive > 0.0:
        _precision = true_positive / (true_positive + false_positive)
    else:
        _precision = float("nan")

    if true_positive + false_negative > 0.0:
        _recall = true_positive / (true_positive + false_negative)
    else:
        _recall = float("nan")

    if _precision + _recall > 0.0:
        _f1_score = 2 * _precision * _recall / (_precision + _recall)
    else:
        _f1_score = float("nan")

    if true_negative + false_positive > 0.0:
        _fall_out = false_positive / (true_negative + false_positive)
    else:
        _fall_out = float("nan")

    return _precision, _recall, _f1_score, _fall_out


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
    fine_tune_types = np.repeat(
        np.nan,
        (
            len(snakemake.input["allFiles"])
            if snakemake.params["forSieve"]
            else len(snakemake.input["lociInfo"])
        ),
    )
elif len(snakemake.params["fineTuneType"]) == 1:
    fine_tune_types = np.repeat(
        snakemake.params["fineTuneType"],
        (
            len(snakemake.input["allFiles"])
            if snakemake.params["forSieve"]
            else len(snakemake.input["lociInfo"])
        ),
    )
else:
    fine_tune_types = []
    fine_tune_type_flag = True


tool = []
tool_setup = []
data_type = []
dataset = []
cell_names = []
loci_info = []
input_data = []
genotypes = []
ternary = []
genotype_probs = []
true_cell_names = []
true_loci_info = []
true_genotypes = []
true_ternary = []
TP = []
FP = []
TN = []
FN = []
recall = []
precision = []
fall_out = []
f1_score = []
TP_mu = []
FP_mu = []
TN_mu = []
FN_mu = []
recall_mu = []
precision_mu = []
fall_out_mu = []
f1_score_mu = []
TP_hetero_mu = []
FP_hetero_mu = []
TN_hetero_mu = []
FN_hetero_mu = []
recall_hetero_mu = []
precision_hetero_mu = []
fall_out_hetero_mu = []
f1_score_hetero_mu = []
TP_homo_mu = []
FP_homo_mu = []
TN_homo_mu = []
FN_homo_mu = []
recall_homo_mu = []
precision_homo_mu = []
fall_out_homo_mu = []
f1_score_homo_mu = []
true_homo_mu = []
true_homo_ref_as_hetero_mu = []
true_homo_mu_as_homo_ref = []
true_homo_mu_as_hetero_mu = []
TP_del = []
FP_del = []
TN_del = []
FN_del = []
recall_del = []
precision_del = []
fall_out_del = []
f1_score_del = []
TP_del_all = []
FP_del_all = []
TN_del_all = []
FN_del_all = []
recall_del_all = []
precision_del_all = []
fall_out_del_all = []
f1_score_del_all = []
TP_del_alt_left = []
FP_del_alt_left = []
TN_del_alt_left = []
FN_del_alt_left = []
recall_del_alt_left = []
precision_del_alt_left = []
fall_out_del_alt_left = []
f1_score_del_alt_left = []
TP_del_ref_left = []
FP_del_ref_left = []
TN_del_ref_left = []
FN_del_ref_left = []
recall_del_ref_left = []
precision_del_ref_left = []
fall_out_del_ref_left = []
f1_score_del_ref_left = []
true_del_ref_left_as_homo_ref = []
true_del_ref_left_as_hetero_mu = []
true_del_alt_left_as_hetero_mu = []
true_del_alt_left_as_homo_mu = []
TP_amp = []
FP_amp = []
TN_amp = []
FN_amp = []
recall_amp = []
precision_amp = []
fall_out_amp = []
f1_score_amp = []

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

            if snakemake.params["forSieve"]:
                __files = sorted(
                    [
                        j
                        for i in vec_get_content_files(
                            get_matched_files(snakemake.input["allFiles"], __pat)[1]
                        )
                        for j in i
                    ]
                )

                inferred_idx, __loci_info = get_matched_files(__files, r"\.loci_info$")
                cnt = np.sum(inferred_idx)

                if cnt == 0:
                    print(
                        f"Warning: no matches found for {'.loci_info$'} in {__files}."
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

                loci_info.extend(__loci_info)

                _, __cell_names = get_matched_files(__files, r"\.cell_names$")
                cell_names.extend(__cell_names)

                genotypes.extend(get_matched_files(__files, r"\.genotypes$")[1])

                _, __ternary = get_matched_files(__files, r"\.ternary$")
                ternary.extend(__ternary)

                genotype_probs.extend(get_matched_files(__files, r"\.probs$")[1])
            else:
                inferred_idx, __loci_info = get_matched_files(
                    snakemake.input["lociInfo"], __pat
                )
                cnt = np.sum(inferred_idx)

                if cnt == 0:
                    print(
                        f"Warning: no matches found for {__pat} in {snakemake.input['lociInfo']}."
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

                loci_info.extend(__loci_info)

                _, __cell_names = get_matched_files(snakemake.input["cellNames"], __pat)
                cell_names.extend(__cell_names)

                genotypes.extend(
                    get_matched_files(snakemake.input["genotypes"], __pat)[1]
                )

                _, __ternary = get_matched_files(snakemake.input["ternary"], __pat)
                ternary.extend(__ternary)

                genotype_probs.extend(
                    get_matched_files(snakemake.input["genotypeProbs"], __pat)[1]
                )

            if "inputData" in snakemake.input.keys():
                if search_for_tool_setup_flag:
                    __pat = "".join(
                        [
                            f"(?=.*{i})"
                            for i in vec_pattern(
                                [__tool_setup.split("-")[0], __data_type, dataset_name],
                                NON_WORDS_PATTERN,
                            )
                        ]
                    )

                __idx, __input_data = get_matched_files(
                    snakemake.input["inputData"], __pat
                )

                if np.sum(inferred_idx) != np.sum(__idx):
                    raise ValueError(
                        f'Files do not match:\n{",".join(__loci_info)}\nand\n{",".join(__input_data)}'
                    )

                input_data.extend(__input_data)

            if len(__loci_info) != len(__cell_names) != len(__ternary):
                raise ValueError(
                    f"lociInfo, cellNames, and ternary must have the same length ({len(__loci_info)} != {len(__cell_names)} != {len(__ternary)})."
                )

            __true_name_pat = "".join(
                [f"(?=.*{i})" for i in vec_pattern_2([dataset_name], NON_WORDS_PATTERN)]
            )
            true_cell_names_idx, __true_cell_names = get_matched_files(
                snakemake.input["trueCellNames"],
                __true_name_pat,
            )
            if np.sum(true_cell_names_idx) == 0 or np.sum(true_cell_names_idx) > 1:
                raise ValueError(
                    f"Pattern {__true_name_pat} should match exactly one object, but {__true_cell_names} found."
                )
            true_cell_names.extend(np.repeat(__true_cell_names, cnt))

            __true_pat = "".join(
                [
                    f"(?=.*{i})"
                    for i in vec_pattern([__data_type, dataset_name], NON_WORDS_PATTERN)
                ]
            )

            true_loci_info_idx, __true_loci_info = get_matched_files(
                snakemake.input["trueSNVSites"], __true_pat
            )
            if np.sum(true_loci_info_idx) == 0 or np.sum(true_loci_info_idx) > 1:
                raise ValueError(
                    f"Pattern {__true_pat} should match exactly one object, but {__true_loci_info} found."
                )
            true_loci_info.extend(np.repeat(__true_loci_info, cnt))

            true_genotypes_idx, __true_genotypes = get_matched_files(
                snakemake.input["trueSNVGenotypesNU"], __true_pat
            )
            if np.sum(true_genotypes_idx) == 0 or np.sum(true_genotypes_idx) > 1:
                raise ValueError(
                    f"Pattern {__true_pat} should match exactly one object, but {__true_genotypes} found."
                )
            true_genotypes.extend(np.repeat(__true_genotypes, cnt))

            true_ternary_idx, __true_ternary = get_matched_files(
                snakemake.input["trueSNVGenotypesTer"], __true_pat
            )
            if np.sum(true_ternary_idx) == 0 or np.sum(true_ternary_idx) > 1:
                raise ValueError(
                    f"Pattern {__true_pat} should match exactly one object, but {__true_ternary} found."
                )
            true_ternary.extend(np.repeat(__true_ternary, cnt))

            for idx in range(len(__loci_info)):
                (
                    tp,
                    fp,
                    tn,
                    fn,
                    tp_mu,
                    fp_mu,
                    tn_mu,
                    fn_mu,
                    tp_hetero_mu,
                    fp_hetero_mu,
                    tn_hetero_mu,
                    fn_hetero_mu,
                    tp_homo_mu,
                    fp_homo_mu,
                    tn_homo_mu,
                    fn_homo_mu,
                    t_homo_ref_as_hetero_mu,
                    t_homo_mu_as_homo_ref,
                    t_homo_mu_as_hetero_mu,
                    total_entries,
                    tp_del,
                    fp_del,
                    tn_del,
                    fn_del,
                    tp_del_all,
                    fp_del_all,
                    tn_del_all,
                    fn_del_all,
                    tp_del_alt_left,
                    fp_del_alt_left,
                    tn_del_alt_left,
                    fn_del_alt_left,
                    tp_del_ref_left,
                    fp_del_ref_left,
                    tn_del_ref_left,
                    fn_del_ref_left,
                    t_ref_left_as_homo_ref,
                    t_ref_left_as_hetero_mu,
                    t_alt_left_as_hetero_mu,
                    t_alt_left_as_homo_mu,
                    tp_amp,
                    fp_amp,
                    tn_amp,
                    fn_amp,
                ) = _get_properties(
                    true_cell_names[-1],
                    true_loci_info[-1],
                    true_ternary[-1],
                    __cell_names[idx],
                    __loci_info[idx],
                    (
                        __input_data[idx]
                        if "inputData" in snakemake.input.keys()
                        else None
                    ),
                    __ternary[idx],
                )

                TP.append(tp)
                FP.append(fp)
                TN.append(tn)
                FN.append(fn)
                TP_mu.append(tp_mu)
                FP_mu.append(fp_mu)
                TN_mu.append(tn_mu)
                FN_mu.append(fn_mu)
                TP_hetero_mu.append(tp_hetero_mu)
                FP_hetero_mu.append(fp_hetero_mu)
                TN_hetero_mu.append(tn_hetero_mu)
                FN_hetero_mu.append(fn_hetero_mu)
                TP_homo_mu.append(tp_homo_mu)
                FP_homo_mu.append(fp_homo_mu)
                TN_homo_mu.append(tn_homo_mu)
                FN_homo_mu.append(fn_homo_mu)
                true_homo_mu.append(tp_homo_mu + fn_homo_mu)
                true_homo_ref_as_hetero_mu.append(t_homo_ref_as_hetero_mu)
                true_homo_mu_as_homo_ref.append(t_homo_mu_as_homo_ref)
                true_homo_mu_as_hetero_mu.append(t_homo_mu_as_hetero_mu)
                TP_del.append(tp_del)
                FP_del.append(fp_del)
                TN_del.append(tn_del)
                FN_del.append(fn_del)
                TP_del_all.append(tp_del_all)
                FP_del_all.append(fp_del_all)
                TN_del_all.append(tn_del_all)
                FN_del_all.append(fn_del_all)
                TP_del_alt_left.append(tp_del_alt_left)
                FP_del_alt_left.append(fp_del_alt_left)
                TN_del_alt_left.append(tn_del_alt_left)
                FN_del_alt_left.append(fn_del_alt_left)
                TP_del_ref_left.append(tp_del_ref_left)
                FP_del_ref_left.append(fp_del_ref_left)
                TN_del_ref_left.append(tn_del_ref_left)
                FN_del_ref_left.append(fn_del_ref_left)
                true_del_ref_left_as_homo_ref.append(t_ref_left_as_homo_ref)
                true_del_ref_left_as_hetero_mu.append(t_ref_left_as_hetero_mu)
                true_del_alt_left_as_hetero_mu.append(t_alt_left_as_hetero_mu)
                true_del_alt_left_as_homo_mu.append(t_alt_left_as_homo_mu)
                TP_amp.append(tp_amp)
                FP_amp.append(fp_amp)
                TN_amp.append(tn_amp)
                FN_amp.append(fn_amp)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp), float(fp), float(tn), float(fn)
                )
                recall.append(_recall)
                precision.append(_precision)
                f1_score.append(_f1_score)
                fall_out.append(_fall_out)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_mu), float(fp_mu), float(tn_mu), float(fn_mu)
                )
                recall_mu.append(_recall)
                precision_mu.append(_precision)
                f1_score_mu.append(_f1_score)
                fall_out_mu.append(_fall_out)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_hetero_mu),
                    float(fp_hetero_mu),
                    float(tn_hetero_mu),
                    float(fn_hetero_mu),
                )
                recall_hetero_mu.append(_recall)
                precision_hetero_mu.append(_precision)
                f1_score_hetero_mu.append(_f1_score)
                fall_out_hetero_mu.append(_fall_out)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_homo_mu),
                    float(fp_homo_mu),
                    float(tn_homo_mu),
                    float(fn_homo_mu),
                )
                recall_homo_mu.append(_recall)
                precision_homo_mu.append(_precision)
                f1_score_homo_mu.append(_f1_score)
                fall_out_homo_mu.append(_fall_out)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_del), float(fp_del), float(tn_del), float(fn_del)
                )
                recall_del.append(_recall)
                precision_del.append(_precision)
                fall_out_del.append(_fall_out)
                f1_score_del.append(_f1_score)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_del_all),
                    float(fp_del_all),
                    float(tn_del_all),
                    float(fn_del_all),
                )
                recall_del_all.append(_recall)
                precision_del_all.append(_precision)
                fall_out_del_all.append(_fall_out)
                f1_score_del_all.append(_f1_score)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_del_alt_left),
                    float(fp_del_alt_left),
                    float(tn_del_alt_left),
                    float(fn_del_alt_left),
                )
                recall_del_alt_left.append(_recall)
                precision_del_alt_left.append(_precision)
                fall_out_del_alt_left.append(_fall_out)
                f1_score_del_alt_left.append(_f1_score)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_del_ref_left),
                    float(fp_del_ref_left),
                    float(tn_del_ref_left),
                    float(fn_del_ref_left),
                )
                recall_del_ref_left.append(_recall)
                precision_del_ref_left.append(_precision)
                fall_out_del_ref_left.append(_fall_out)
                f1_score_del_ref_left.append(_f1_score)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_amp), float(fp_amp), float(tn_amp), float(fn_amp)
                )
                recall_amp.append(_recall)
                precision_amp.append(_precision)
                f1_score_amp.append(_f1_score)
                fall_out_amp.append(_fall_out)

            if fine_tune_type_flag:
                fine_tune_type_idx = which_patterns(
                    __loci_info,
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


variant_collection = pd.DataFrame(
    {
        "cell_num": snakemake.params["cellNum"],
        "coverage_mean": snakemake.params["covMean"],
        "coverage_variance": snakemake.params["covVariance"],
        "eff_seq_err_rate": snakemake.params["effSeqErrRate"],
        "ado_rate": snakemake.params["adoRate"],
        "mutation_rate": snakemake.params["mutationRate"],
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
        "true_cell_names": pd.Series(true_cell_names),
        "cell_names": pd.Series(cell_names),
        "true_loci_info": pd.Series(true_loci_info),
        "loci_info": pd.Series(loci_info),
        "true_genotypes": pd.Series(true_genotypes),
        "genotypes": pd.Series(genotypes),
        "true_ternary": pd.Series(true_ternary),
        "ternary": pd.Series(ternary),
        "genotype_probs": pd.Series(genotype_probs),
        "true_positive": pd.Series(TP),
        "false_positive": pd.Series(FP),
        "true_negative": pd.Series(TN),
        "false_negative": pd.Series(FN),
        "recall": pd.Series(recall),
        "precision": pd.Series(precision),
        "f1_score": pd.Series(f1_score),
        "fall_out": pd.Series(fall_out),
        "true_positive_mu": pd.Series(TP_mu),
        "false_positive_mu": pd.Series(FP_mu),
        "true_negative_mu": pd.Series(TN_mu),
        "false_negative_mu": pd.Series(FN_mu),
        "recall_mu": pd.Series(recall_mu),
        "precision_mu": pd.Series(precision_mu),
        "f1_score_mu": pd.Series(f1_score_mu),
        "fall_out_mu": pd.Series(fall_out_mu),
        "true_positive_hetero_mu": pd.Series(TP_hetero_mu),
        "false_positive_hetero_mu": pd.Series(FP_hetero_mu),
        "true_negative_hetero_mu": pd.Series(TN_hetero_mu),
        "false_negative_hetero_mu": pd.Series(FN_hetero_mu),
        "recall_hetero_mu": pd.Series(recall_hetero_mu),
        "precision_hetero_mu": pd.Series(precision_hetero_mu),
        "f1_score_hetero_mu": pd.Series(f1_score_hetero_mu),
        "fall_out_hetero_mu": pd.Series(fall_out_hetero_mu),
        "true_positive_homo_mu": pd.Series(TP_homo_mu),
        "false_positive_homo_mu": pd.Series(FP_homo_mu),
        "true_negative_homo_mu": pd.Series(TN_homo_mu),
        "false_negative_homo_mu": pd.Series(FN_homo_mu),
        "recall_homo_mu": pd.Series(recall_homo_mu),
        "precision_homo_mu": pd.Series(precision_homo_mu),
        "f1_score_homo_mu": pd.Series(f1_score_homo_mu),
        "fall_out_homo_mu": pd.Series(fall_out_homo_mu),
        "true_homo_mu": pd.Series(true_homo_mu),
        "true_homo_ref_as_hetero_mu": pd.Series(true_homo_ref_as_hetero_mu),
        "true_homo_mu_as_homo_ref": pd.Series(true_homo_mu_as_homo_ref),
        "true_homo_mu_as_hetero_mu": pd.Series(true_homo_mu_as_hetero_mu),
        "true_positive_del": pd.Series(TP_del),
        "false_positive_del": pd.Series(FP_del),
        "true_negative_del": pd.Series(TN_del),
        "false_negative_del": pd.Series(FN_del),
        "recall_del": pd.Series(recall_del),
        "precision_del": pd.Series(precision_del),
        "fall_out_del": pd.Series(fall_out_del),
        "f1_score_del": pd.Series(f1_score_del),
        "true_positive_all_del": pd.Series(TP_del_all),
        "false_positive_all_del": pd.Series(FP_del_all),
        "true_negative_all_del": pd.Series(TN_del_all),
        "false_negative_all_del": pd.Series(FN_del_all),
        "recall_all_del": pd.Series(recall_del_all),
        "precision_all_del": pd.Series(precision_del_all),
        "fall_out_all_del": pd.Series(fall_out_del_all),
        "f1_score_all_del": pd.Series(f1_score_del_all),
        "true_positive_del_alt_left": pd.Series(TP_del_alt_left),
        "false_positive_del_alt_left": pd.Series(FP_del_alt_left),
        "true_negative_del_alt_left": pd.Series(TN_del_alt_left),
        "false_negative_del_alt_left": pd.Series(FN_del_alt_left),
        "recall_del_alt_left": pd.Series(recall_del_alt_left),
        "precision_del_alt_left": pd.Series(precision_del_alt_left),
        "fall_out_del_alt_left": pd.Series(fall_out_del_alt_left),
        "f1_score_del_alt_left": pd.Series(f1_score_del_alt_left),
        "true_positive_del_ref_left": pd.Series(TP_del_ref_left),
        "false_positive_del_ref_left": pd.Series(FP_del_ref_left),
        "true_negative_del_ref_left": pd.Series(TN_del_ref_left),
        "false_negative_del_ref_left": pd.Series(FN_del_ref_left),
        "recall_del_ref_left": pd.Series(recall_del_ref_left),
        "precision_del_ref_left": pd.Series(precision_del_ref_left),
        "fall_out_del_ref_left": pd.Series(fall_out_del_ref_left),
        "f1_score_del_ref_left": pd.Series(f1_score_del_ref_left),
        "true_del_ref_left_as_homo_ref": pd.Series(true_del_ref_left_as_homo_ref),
        "true_del_ref_left_as_hetero_mu": pd.Series(true_del_ref_left_as_hetero_mu),
        "true_del_alt_left_as_hetero_mu": pd.Series(true_del_alt_left_as_hetero_mu),
        "true_del_alt_left_as_homo_mu": pd.Series(true_del_alt_left_as_homo_mu),
        "true_positive_amp": pd.Series(TP_amp),
        "false_positive_amp": pd.Series(FP_amp),
        "true_negative_amp": pd.Series(TN_amp),
        "false_negative_amp": pd.Series(FN_amp),
        "recall_amp": pd.Series(recall_amp),
        "precision_amp": pd.Series(precision_amp),
        "fall_out_amp": pd.Series(fall_out_amp),
        "f1_score_amp": pd.Series(f1_score_amp),
    }
)

variant_collection.to_csv(
    snakemake.output[0],
    sep="\t",
    na_rep="NA",
    index=False,
)
