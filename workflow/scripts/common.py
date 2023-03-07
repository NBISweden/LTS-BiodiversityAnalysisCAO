def samples_to_genomemap(config, samples):
    lib_types = config["lib_type_to_map"]
    return [x for x in samples.keys() if samples[x]["lib_type"] in lib_types]

def all_output(wildcards):
    output = []
    output.extend(
        expand(
            "{results_dir}/mappings/{map_name}/{mapper}/{sample}_{seqtype}.fastq.gz",
            results_dir=config["results_dir"],
            map_name=mappings.keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
            seqtype=["fwd", "rev", "unp"],
        )
    )
    # extend with fish mapping
    if "mappings" in config.keys():
        output.extend(
            expand(
                "{results_dir}/genome_mappings/{ref}/{f}",
                results_dir=config["results_dir"],
                f=["summary_raw_counts.csv", "summary_size_adjusted.csv", "lib_sizes.csv"],
                ref=config["mappings"]["genomes"].keys(),
            )
        )
        output.extend(
            expand(
                "{results_dir}/mappings/{map_name}/map_qc.html",
                map_name=mappings.keys(),
                results_dir=config["results_dir"],
            )
        )
    output.extend(
        expand("{results_dir}/multiqc/multiqc.html", results_dir=config["results_dir"])
    )
    # TODO: Fail explicitly here if sintax database does not exist
    if os.path.exists(config["sintax"]["db"]):
        output.extend(
            expand(
                "{results_dir}/mappings/{map_name}/{mapper}/{sample}.sintax.parsed.tsv",
                results_dir=config["results_dir"],
                map_name=mappings.keys(),
                mapper=config["mappers"],
                sample=samples.keys(),
            )
        )
        output.extend(
            expand(
                "{results_dir}/mappings/{map_name}/{mapper}/krona/krona.html",
                results_dir=config["results_dir"],
                map_name=mappings.keys(),
                mapper=config["mappers"],
            )
        )
        output.extend(
            expand(
                "{results_dir}/mappings/counts/sintax.{map_name}.{mapper}.tsv",
                results_dir=config["results_dir"],
                map_name=mappings.keys(),
                mapper=config["mappers"],
            )
        )
    return output


def mem_allowed(wildcards, threads):
    return max(threads * 6400, 6400)


def krona_input_string(wildcards):
    input_string = []
    for sample in samples.keys():
        text = f"results/mappings/{wildcards.map_name}/{wildcards.mapper}/{sample}.sintax.parsed.krona.txt,{sample}"
        input_string.append(text)
    return " ".join(input_string)
