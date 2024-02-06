
rule download_gtdb:
    input:
        config["project_directory"]
    output:
        "{proj_dir}/data/refs/gtdb-species.fasta"
        "{proj_dir}/data/refs/gtdb-taxonomy.fasta"
    conda:
        "../etc/environment.yml"
    params:
        "{proj_dir}/data"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        scripts/download_gtdb.py -o {params}
        """

rule download_genbank:
    output:
        outdir = "{proj_dir}/data/ref",
        file = "{proj_dir}/data/ref/genbank-16s-{taxa}.fasta.gz"
    params:
        tax = "{taxa}"
        batch = config["batch_size"]
    conda:
        "../etc/environment.yml"
    threads: 32
    resources:
        mem_mb = 32000
    shell:
        """
        scripts/download_genbank.sh {output.outdir} {params.tax} {params.batch} {output.file}
        """

rule download_genbank:
    input:
    output:
    params:
    threads:
    resources:
    conda:
    shell:
    """
    
    """

rule extract_16s:
    input:
    output:
    params:
    threads:
    resources:
    conda:
    shell:
    """
    
    """


rule format_gtdb:
    input:
        "{proj_dir}/data/refs/gtdb-species.fasta"
        "{proj_dir}/data/refs/gtdb-taxonomy.fasta"
    output:
        "{proj_dir}/data/refs/gtdb-species.fasta"
        "{proj_dir}/data/refs/gtdb-taxonomy.fasta"
    conda:
        "../etc/ncbi-download-env.yml"
    params:
        "{proj_dir}/data"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        scripts/download_gtdb.py -o {params}
        """

