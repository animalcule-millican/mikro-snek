rule asv_fasta:
    input:
        "{proj_dir}/data/output/seqtab_nochim.RData"
    output:
        fasta = "{proj_dir}/data/output/asv.fasta",
        r_file = "{proj_dir}/data/output/asv-fasta.RData"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/write_fasta_file.r {input} {output.fasta} {output.r_file}
        """

rule align:
    input:
        "{proj_dir}/data/output/seqtab_nochim.RData"
    output:
        fasta = "{proj_dir}/data/output/asv.afa",
        r_file = "{proj_dir}/data/output/asv-align.RData"
    threads: 32
    resources:
        mem_mb = 64000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/align.r {input} {output.fasta} {output.r_file}
        """

rule fast_tree:
    input:
        "{proj_dir}/data/output/asv.aln"
    output:
        "{proj_dir}/data/output/asv.tree"
    threads: 32
    resources:
        mem_mb = 64000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        fasttree {input} {output}
        """