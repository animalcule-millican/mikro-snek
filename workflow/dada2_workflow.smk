from evalfasta.evalfasta import *

seq_count, seq_len = get_seq_metrics(config["sequence_files"])
num_chunks = get_chunky(int(config["split_chunks"]))

rule filter_trim:
    input:
        path = config["sequence_files"]
    output:
        F = "{proj_dir}/data/output/filt_trm_fwd.RData",
        R = "{proj_dir}/data/output/filt_trm_rev.RData"
    params:
        trunc_len = seq_len
    threads: 32
    resources:
        mem_mb = 32000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/filter_trim.r {input.path} {output.F} {output.R} {params.trunc_len}
        """

rule learn_errors:
    input:
        input_file = "{proj_dir}/data/output/filt_trm_{orientation}.RData"
    output:
        output_file = "{proj_dir}/data/output/err_{orientation}.RData"
    params:
        orientation = "{orientation}"
    threads: 32
    resources:
        mem_mb = 32000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/learn_errors.r {input.input_file} {output.output_file} {params.orientation}
        """

rule split_derep:
    input:
        input_file = "{proj_dir}/data/output/err_{orientation}.RData"
    output:
        "{proj_dir}/data/output/derep_{orientation}_{chunks}.RData"
    params:
        orientation = "{orientation}",
        seq_count = seq_count,
        save_dir = "{proj_dir}/data/output",
        chunks = config["split_chunks"]
    threads: 1
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    resources:
        mem_mb = 10000
    shell:
        """
        scripts/split-derep.r {input.input_file} {params.save_dir} {params.seq_count} {params.chunks} {params.orientation}
        """

rule dada:
    input:
        "{proj_dir}/data/output/derep_{orientation}_{chunks}.RData"
    output:
        "{proj_dir}/data/output/dada_{orientation}_{chunks}.RData"
    params:
        "{orientation}"
    threads: 32
    resources:
        mem_mb = 32000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/dada.r {input} {output} {params}
        """

rule merge:
    input:
        expand("{{proj_dir}}/data/output/dada_{orientation}_{{chunks}}.RData", orientation = config["orientation"])
    output:
        output_file = "{proj_dir}/data/output/merged_{chunks}.RData",
    threads: 1
    params:
        fwd = "{proj_dir}/data/output/dada_fwd_{chunks}.RData",
        rev = "{proj_dir}/data/output/dada_rev_{chunks}.RData"
    resources:
        mem_mb = 10000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/merge.r {params.fwd} {params.rev} {output.output_file} {input}
        """

rule join_seqtab:
    input:
        expand("{proj_dir}/data/output/merged_{chunks}.RData", orientation = '', chunks=num_chunks, proj_dir=config["project_directory"])
    output:
        "{proj_dir}/data/output/seqtab.RData"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/join_seqtabs.r {input} {output}
        """

rule remove_bimera:
    input:
        "{proj_dir}/data/output/seqtab.RData"
    output:
        "{proj_dir}/data/output/seqtab_nochim.RData"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "{proj_dir}/etc/mikro-snek-env.yml"
    shell:
        """
        scripts/remove_bimera.r {input} {output}
        """