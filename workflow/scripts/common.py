def samples_to_genomemap(config, samples):
    """
    Return samples that should be mapped to genomes

    :param config: config dict
    :param samples: samples dict
    :return: list of samples

    lib_types is a list of library types that should be mapped to genomes, currently only "DNA" is supported
    a list of samples is returned that have a library type specified in the sample list matching lib_types
    """
    lib_types = config["lib_type_to_map"]
    return [x for x in samples.keys() if samples[x]["lib_type"] in lib_types]


def all_output(wildcards):
    """
    Return all output files for given wildcards
    """
    output = []
    if not "mappings" in config.keys():
        return output
    ### Marker gene mappings ###
    # fastq files for reads mapped to the marker gene databases
    output.extend(
        expand(
            "{results_dir}/mappings/{map_name}/{mapper}/{sample}_{seqtype}.fastq.gz",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
            seqtype=["fwd", "rev", "unp"],
        )
    )
    # mapping statistics report for marker gene databases
    output.extend(
        expand(
            "{results_dir}/mappings/{map_name}/map_qc.html",
            map_name=config["mappings"]["marker_genes"].keys(),
            results_dir=config["results_dir"],
        )
    )
    # sintax parsed output for each sample and marker gene database
    output.extend(
        expand(
            "{results_dir}/mappings/{map_name}/{mapper}/{sample}.sintax.parsed.tsv",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
        )
    )
    # krona report for each marker gene database
    output.extend(
        expand(
            "{results_dir}/mappings/{map_name}/{mapper}/krona/krona.html",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
        )
    )
    output.extend(
        expand(
            "{results_dir}/mappings/counts/sintax.{map_name}.{mapper}.tsv",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
        )
    )
    ### Genome mappings ###
    output.extend(
        expand(
            "{results_dir}/genome_mappings/{ref}/{f}",
            results_dir=config["results_dir"],
            f=["summary_raw_counts.csv", "summary_size_adjusted.csv", "lib_sizes.csv"],
            ref=config["mappings"]["genomes"].keys(),
        )
    )
    # QC report for preprocessed samples
    output.extend(
        expand("{results_dir}/multiqc/multiqc.html", results_dir=config["results_dir"])
    )
    return output


def sintax_krona_input(wildcards):
    if "mappings" in config.keys():
        files = expand(
            "{results_dir}/mappings/{map_name}/{mapper}/krona/krona.html",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
        )
    else:
        files = []
    return files


def sintax_map_input(wildcards):
    if "mappings" in config.keys():
        files = expand(
            "{results_dir}/mappings/counts/sintax.{map_name}.{mapper}.tsv",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
        )
    else:
        files = []
    return files


def coi_map_input(wildcards):
    if "mappings" in config.keys():
        files = expand(
            "{results_dir}/mappings/{map_name}/{mapper}/{sample}.krakenuniq.filtered.fastq.gz",
            results_dir=config["results_dir"],
            map_name=config["mappings"]["marker_genes"].keys(),
            mapper=config["mappers"],
            sample=samples.keys(),
        )
    else:
        files = []
    return files


def genome_count_input(wildcards):
    if "mappings" in config.keys():
        files = expand(
            "{results_dir}/genome_mappings/{ref}/summary_{f}.csv",
            results_dir=config["results_dir"],
            f=["raw_counts", "size_adjusted"],
            t=["taxid", "species"],
            ref=config["mappings"]["genomes"].keys(),
        )
    else:
        files = []
    return files


def genome_map_input(wildcards):
    if "mappings" in config.keys():
        files = expand(
            "{results_dir}/genome_mappings/{ref}/bowtie2/{sample}/{sample}.reads.fa",
            results_dir=config["results_dir"],
            sample=samples_to_genomemap(config, samples),
            ref=config["mappings"]["genomes"].keys(),
        )
    else:
        files = []
    return files


def mem_allowed(wildcards, threads):
    """
    Return memory allowed for given wildcards and threads

    :param wildcards: wildcards object
    :param threads: number of threads
    :return: memory allowed

    The memory allowed is calculated by assuming each thread has 6.4 GB of RAM
    and multiplying this by the number of threads. The function returns at least 6.4 GB.
    """
    return max(threads * 6400, 6400)


def krona_input_string(wildcards):
    """
    Return input string for krona report

    :param wildcards: wildcards object
    :return: input string for krona report

    This function puts together the input string for the krona report. It iterates samples in the samples dict
    and adds the sample name and the path to the krona report for the sample to the input string.
    """
    input_string = []
    for sample in samples.keys():
        text = f"{wildcards.results_dir}/mappings/{wildcards.map_name}/{wildcards.mapper}/{sample}.sintax.parsed.krona.txt,{sample}"
        input_string.append(text)
    return " ".join(input_string)
