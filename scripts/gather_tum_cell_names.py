import os
import copy

import pandas as pd


def get_tum_cell_names(f: str) -> list[str]:
    assert os.path.isfile(filename)

    ret: pd.DataFrame = pd.read_csv(f)

    return ret["cell"].values


cell_names: list[str] = []
for filename in snakemake.input:
    _cell_names: list[str] = get_tum_cell_names(filename)

    if len(cell_names) > 0:
        assert len(cell_names) == len(_cell_names)
        for idx, val in enumerate(cell_names):
            assert val == _cell_names[idx]
    else:
        cell_names = copy.deepcopy(_cell_names)

with open(
    snakemake.output[0],
    "w",
    encoding="utf-8",
) as fh:
    for i in cell_names:
        fh.write(f"{i}\n")
