rule assign_taxonomy:
    input:
        "{proj_dir}/data/output/seqtab_nochim.RData"
    output:
        taxonomy = "{proj_dir}/data/output/taxonomy.RData",
        both = "{proj_dir}/data/output/taxonomy_seqtab.RData"
    params:
        tax = config["taxonomy_file"]
    threads: 32
    resources:
        mem_mb = 64000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/assign_taxonomy.r {input} {output.taxonomy} {output.both} {params.tax}
        """

rule add_species:
    input:
        "{proj_dir}/data/output/taxonomy.RData"
    output:
        taxonomy = "{proj_dir}/data/output/species_taxonomy.RData",
        both = "{proj_dir}/data/output/sp_mult_taxonomy.RData"
    params:
        tax = config["species_file"]
    threads: 32
    resources:
        mem_mb = 64000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/assign_taxonomy.r {input} {output.taxonomy} {output.both} {params.tax}
        """