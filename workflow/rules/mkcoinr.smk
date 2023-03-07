rule coinr_db:
    input:
        expand("resources/coinr/{name}.sintax.fasta",
            name=config["coinr"]["dbname"],
        ),

rule download_coinr_src:
    output:
        format_db="resources/coinr/format_db.pl",
        mkdb="resources/coinr/mkdb.pm",
    log:
        "logs/coinr/download_coinr_src.log"
    retries: 2
    params:
        format_db_url="https://raw.githubusercontent.com/meglecz/mkCOInr/main/scripts/format_db.pl",
        mkdb_url="https://raw.githubusercontent.com/meglecz/mkCOInr/main/scripts/mkdb.pm",
    shell:
        """
        exec &>{log}
        curl -L -o {output.format_db} {params.format_db_url}
        curl -L -o {output.mkdb} {params.mkdb_url}
        """

rule download_coinr:
    output:
        tar=temp("resources/coinr/coinr.tar.gz")
    log:
        "logs/coinr/download.log"
    retries: 2
    params:
        url="https://zenodo.org/record/6555985/files/COInr_2022_05_06.tar.gz?download=1",
        tmpdir="$TMPDIR/coinr",
    shell:
        """
        mkdir -p {params.tmpdir}
        curl -L -o {params.tmpdir}/coinr.tar.gz {params.url} 2>{log}
        mv {params.tmpdir}/coinr.tar.gz {output.tar}
        rm -rf {params.tmpdir}
        """

rule extract_coinr:
    input:
        rules.download_coinr.output.tar,
    output:
        tsv="resources/coinr/COInr.tsv.gz",
        tax="resources/coinr/taxonomy.tsv.gz"
    log:
        "logs/coinr/extract.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.tsv),
        tmpdir="$TMPDIR/coinr"
    shell:
        """
        mkdir -p {params.tmpdir}
        tar -C {params.tmpdir} -xvf {input[0]} 2>{log}
        gzip {params.tmpdir}/*/*.tsv
        mv {params.tmpdir}/*/*.tsv.gz  {params.outdir}
        rm -rf {params.tmpdir} 
        """

rule format_coinr:
    output:
        fas="resources/coinr/{name}.qiime.fasta.gz",
        tax="resources/coinr/{name}.taxonomy.tsv.gz"
    input:
        tsv=rules.extract_coinr.output.tsv,
        tax=rules.extract_coinr.output.tax,
        format_db=rules.download_coinr_src.output.format_db,
        mkdb=rules.download_coinr_src.output.mkdb
    log:
        "logs/coinr/format_coinr_{name}.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.fas),
        tmpdir="$TMPDIR/format_coinr",
        out=lambda wildcards, output: os.path.basename(output.fas),
        name="{name}"
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.tax} > {params.tmpdir}/tax.tsv
        gunzip -c {input.tsv} > {params.tmpdir}/seqs.tsv
        perl {input.format_db} -tsv {params.tmpdir}/seqs.tsv \
            -taxonomy {params.tmpdir}/tax.tsv -outdir {params.tmpdir} \
            -out {params.name} -outfmt qiime > {log} 2>&1
        gzip {params.tmpdir}/{params.name}_trainseq.fasta 
        mv {params.tmpdir}/{params.name}_trainseq.fasta.gz {output.fas}
        cat {params.tmpdir}/{params.name}_taxon.txt | tr ' ' '\t' | gzip -c > {params.tmpdir}/tax.tsv.gz 
        mv {params.tmpdir}/tax.tsv.gz {output.tax}
        rm -rf {params.tmpdir}
        """

rule coinr2sintax:
    output:
        fas="resources/coinr/{name}.sintax.fasta"
    input:
        fas=rules.format_coinr.output.fas,
        tax=rules.format_coinr.output.tax,
    log:
        "logs/coinr/coinr2sintax.{name}.log"
    params:
        tmpdir="$TMPDIR/coinr2sintax",
        src=srcdir("../scripts/coinr2sintax.py")
    shell:
        """
        mkdir -p {params.tmpdir}
        gunzip -c {input.fas} > {params.tmpdir}/seqs.fasta
        gunzip -c {input.tax} > {params.tmpdir}/taxa.tsv
        python {params.src} {params.tmpdir}/seqs.fasta \
            {params.tmpdir}/taxa.tsv > {params.tmpdir}/sintax.fasta 2>{log} 
        mv {params.tmpdir}/sintax.fasta {output.fas}
        rm -rf {params.tmpdir}
        """

