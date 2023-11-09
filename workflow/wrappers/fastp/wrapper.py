__author__ = "John Sundh"
__email__ = "john.sundh@nbis.se"
__license__ = "MIT"


import os
from snakemake.shell import shell
params = snakemake.params
log = snakemake.log
threads = snakemake.threads
output = snakemake.output
input = snakemake.input
dedup = snakemake.params.get("dedup", False)

if dedup:
    dedup_str = "--dedup"
else:
    dedup_str = ""

shell(
    """
    exec &>{log.shell}
    cat {input.R1} > {params.tmpR1}
    cat {input.R2} > {params.tmpR2}
    fastp --thread {threads} {dedup_str} -y -Y {params.complexity_threshold} -l {params.min_length} \
        -i {params.tmpR1} -I {params.tmpR2} -o {params.outR1} \
        -O {params.outR2} -j {log.json} {params.settings} > {log.log} 2>&1
    mv {params.outR1} {output.R1}
    mv {params.outR2} {output.R2}
    rm {params.tmpR1} {params.tmpR2}
    """
)