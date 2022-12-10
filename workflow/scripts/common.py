def all_output(wildcards):
    output = []
    output.extend(
        expand(
            "{out_folder}/results/mappings/{map_name}/{mapper}/{sample}_{seqtype}.fastq.gz",
            out_folder=config["out_folder"],
            map_name=mappings.keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
            seqtype=["fwd", "rev", "unp"],
        )
    )
    # extend with fish mapping
    output.extend(
        expand(
            "{out_folder}/results/genome_mappings/{ref}/counts/{sample}_{t}_counts.txt",
            sample=samples.keys(),
            t=["taxid","species"],
            ref=config["mappings"]["genomes"].keys(),
            out_folder=config["out_folder"],
        )
    )
    output.extend(
        expand(
            "{out_folder}/results/multiqc/multiqc.html",
            out_folder=config["out_folder"],
        )
    )
    output.extend(
        expand(
            "{out_folder}/results/mappings/{map_name}/map_qc.html",
            out_folder=config["out_folder"],
            map_name=mappings.keys(),
        )
    )
    #TODO: Fail explicitly here if sintax database does not exist
    if os.path.exists(config["sintax"]["db"]):
        output.extend(
            expand(
                "{out_folder}/results/mappings/{map_name}/{mapper}/{sample}.sintax.{suff}",
                out_folder=config["out_folder"],
                map_name=mappings.keys(),
                mapper=config["mappers"],
                sample=samples.keys(),
                suff=["parsed.krona.txt", "parsed.tsv"],
            )
        )
        output.extend(
            expand(
                "{out_folder}/results/mappings/{map_name}/{mapper}/krona/krona.html",
                out_folder=config["out_folder"],
                map_name=mappings.keys(),
                mapper=config["mappers"],
            )
        )
    return output


def krona_input_string(wildcards):
    input_string = []
    for sample in samples.keys():
        text = f"{wildcards.out_folder}/results/mappings/{wildcards.map_name}/{wildcards.mapper}/{sample}.sintax.parsed.krona.txt,{sample}"
        input_string.append(text)
    return " ".join(input_string)
