import os
import sys
import numpy as np
import pandas as pd
from warnings import warn

configfile: "config.yaml"

if not workflow.use_conda:
    sys.stderr.write("\nYou are not using conda. Pass '--use-conda' flag to snakemake.\n")
    sys.exit(1)

class Study:
    def __init__(self, file, outdir):
        self._file = file
        self._outdir = outdir
        out = self.study_file_io(file)
        self._filenames = out[0]
        self._replicates = out[1]
        self._groups = out[2]
        self.make_dirs()

    @property
    def file(self) -> str:
        return self._file

    @property
    def outdir(self) -> str:
        return self._outdir

    @property
    def filenames(self) -> list[str]:
        return self._filenames

    @property
    def replicates(self) -> list[str]:
        return self._replicates

    @property
    def groups(self) -> list[str]:
        return self._groups

    def study_file_io(self, file: str) -> (list[str], list[str], list[str]):
        if not os.path.exists(file):
            raise FileNotFoundError(
                f"Study File '{file}' not found. Be sure it is correctly marked in 'config.yaml' "
                f"or PROVIDE_STUDY_FILE is set to True."
            )
        sdf = pd.read_csv(file)
        if sdf.shape[1] < 3:
            raise ValueError(
                "Study File contains less than 3 columns. It must have 3 column in the following order: "
                "\n<filename.cram>, <replicate ID>, <Group ID>"
            )
        elif sdf.shape[1] > 3:
            warn("Study File contains more than 3 columns, additional columns will be ignored.")
                
        sdf.columns  = ["filename", "rep_id", "grp_id"]
        print("are u my type? ", type(sdf.iloc[:, 1]))
            
        return sdf.iloc[:, 0].tolist(), sdf.iloc[:, 1].tolist(), sdf.iloc[:, 2].tolist()

    def make_dirs(self) -> None:
        os.makedirs(self.outdir, exist_ok=True)
        for group, rep in zip(self.groups, self.replicates):
            os.makedirs(os.path.join(self.outdir, group), exist_ok=True)
            os.makedirs(os.path.join(self.outdir,group, rep),exist_ok=True)

        return

    def assign_checkfiles(self) -> list[str]:
        return [os.path.join(
                self.outdir, rep, f"{group}_{rep}_check.txt"
            ) for group, rep in zip(self.groups, self.replicates)
        ]

    def assign_realign_bams(self) -> list[str]:
        return [
            os.path.join(
                self.outdir, group, rep, f"{group}_{rep}_realigned.bam"
            ) for group, rep in zip(self.groups, self.replicates)
        ]

    def assign_vcf(self) -> list[str]:
        return [
            os.path.join(
                self.outdir, group, rep, f"{group}_{rep}_repeats.vcf"
            ) for group, rep in zip(self.groups, self.replicates)
        ]

    def assign_json(self) -> list[str]:
        return [
            os.path.join(
                self.outdir, group, rep, f"{group}_{rep}_repeats.json"
            ) for group, rep in zip(self.groups, self.replicates)
        ]


if config["STUDY_CONFIG"]["PROVIDE_STUDY_FILE"]:
    STUDY = Study(config["STUDY_CONFIG"]["STUDY_FILE_PATH"], config["OUT_DIR"])
    rule_all = []
else:  # TODO: prepose w/ rule to generate study file if set gets too big
    rule_all = [os.path.join(os.path.join(config["OUT_DIR"], config["STUDY_CONFIG"]["STUDY_FILE_PATH"]))]

rule_all.extend([
    expand(
        os.path.join(config["OUT_DIR"], "{groups}", "{reps}", "{groups}_{reps}_realigned.bam"),
        zip,
        groups=STUDY.groups,
        reps=STUDY.replicates
    ),
    expand(
        os.path.join(config["OUT_DIR"], "{groups}", "{reps}", "{groups}_{reps}_repeats.vcf"),
        zip,
        groups=STUDY.groups,
        reps=STUDY.replicates
    ),
    expand(
        os.path.join(config["OUT_DIR"],"{groups}","{reps}","{groups}_{reps}_repeats.json"),
        zip,
        groups=STUDY.groups,
        reps=STUDY.replicates
    )
])

rule all:
    input: rule_all


rule run_expansion_hunter:
    input:
        os.path.join(
            config["STUDY_DIR"],
            "{groups}",
            "{reps}",
            "high_cov_alignment",
            "{reps}.alt_bwamem_GRCh38DH.20150917.{groups}.high_coverage.cram"
        ),
    output: 
        realign=os.path.join(config["OUT_DIR"], "{groups}", "{reps}", "{groups}_{reps}_realigned.bam"),
        vcf=os.path.join(config["OUT_DIR"], "{groups}", "{reps}", "{groups}_{reps}_repeats.vcf"),
        json=os.path.join(config["OUT_DIR"], "{groups}", "{reps}", "{groups}_{reps}_repeats.json")
    params:
        out_dir=config["OUT_DIR"],
        var_cat=config["VAR_CAT"],
        groups="{groups}",
        reps="{reps}",
        exp_hunt_exec=config["EXP_HUNT_LOC"],
        samp_dir=config["STUDY_DIR"],
        ref_genome=config["REFERENCE_GENOME"],
        home_dir=config["HOME_DIR"]
    threads: 4
    resources:
        mem_mb = 500

    shell:
        """
        start=$SECONDS
        {params.exp_hunt_exec} \
            --reads  "{input}" \
            --reference "{params.ref_genome}" \
            --variant "{params.var_cat}" \
            --output-prefix "{params.reps}" \
            --analysis-mode streaming
            
        mv "{params.home_dir}/{params.reps}_realigned.bam" \
            "{params.out_dir}/{params.groups}/{params.reps}/{params.groups}_{params.reps}_realigned.bam"
        mv "{params.home_dir}/{params.reps}.vcf" \
            "{params.out_dir}/{params.groups}/{params.reps}/{params.groups}_{params.reps}_repeats.vcf"
        mv "{params.home_dir}/{params.reps}.json" \
            "{params.out_dir}/{params.groups}/{params.reps}/{params.groups}_{params.reps}_repeats.json"
        duration=$(( SECONDS - start ))
        echo "TIME"
        echo $duration
        """
