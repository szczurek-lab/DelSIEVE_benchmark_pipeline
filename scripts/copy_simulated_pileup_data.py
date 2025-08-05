import os
import re
import subprocess

for k, v in snakemake.input.items():
    assert k in snakemake.output.keys()

    if k == "mpileupDir":
        subprocess.call(
            [
                "mkdir",
                "-p",
                snakemake.output[k],
            ],
        )

        for i in os.listdir(v):
            if (
                re.fullmatch(
                    r".*\.pileup$",
                    i,
                    flags=re.I,
                )
                is not None
            ):
                subprocess.call(
                    [
                        "cp",
                        "-f",
                        os.path.join(
                            v,
                            i,
                        ),
                        snakemake.output[k],
                    ],
                )
    else:
        subprocess.call(
            [
                "cp",
                "-f",
                v,
                snakemake.output[k],
            ],
        )
