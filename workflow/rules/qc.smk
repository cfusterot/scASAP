import glob
import os

def fastqc_input(wc):
    f=units.loc[(wc.sample,wc.unit)]['fq' + wc.read]
    return f

rule fastqc:
    input:
        fastqc_input
    output:
        html="{}/{{sample}}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.html".format(OUTDIR),
        zip="{}/{{sample}}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.zip".format(OUTDIR)
# the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        lambda wc: "-t {}".format(get_resource("fastqc","threads")) 
        #expand("--quiet --outdir '{{OUTDIR}}/qc/{sample}'", sample=SAMPLES) 
    resources:
        mem_mb=get_resource("fastqc", "mem_mb"),
        walltime=get_resource("fastqc", "walltime")
    log:
        "{}/{{sample}}/fastqc_{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/{{sample}}/fastqc_{{unit}}.r{{read}}.bmk".format(LOGDIR)
    threads: 
        threads=get_resource("fastqc", "threads")
    wrapper:
        "v1.23.1/bio/fastqc"

rule fastq_screen_indexes:
    output:
        "resources/FastQ_Screen_Genomes/fastq_screen.conf"
    log:
        "{}/fastq_screen_idx.log".format(LOGDIR)
    benchmark:
        "{}/fastq_screen_idx.bmk".format(LOGDIR)
    conda:
        "envs/fastq_screen.yaml"
    threads: get_resource("fastq_screen", "threads")
    params:
        outdir=config["fastq_screen"]["index_dir"]
    resources:
        mem_mb=get_resource("fastq_screen", "mem_mb"),
        walltime=get_resource("fastq_screen", "walltime")
    shell:
        """
        fastq_screen --threads {threads} --get_genomes --outdir{params.outdir} /
        &> {log}
        """

rule fastq_screen:
    input:
        lambda wc: units.loc[(wc.sample,wc.unit)]['fq' + wc.read],
        conf="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["fastq_screen"]["index_dir"])
    output:
        txt="{}/{{sample}}/qc/fastq_screen/{{sample}}.{{unit}}.r{{read}}.fastq_screen.txt".format(OUTDIR),
        png="{}/{{sample}}/qc/fastq_screen/{{sample}}.{{unit}}.r{{read}}.fastq_screen.png".format(OUTDIR)
    log:
        "{}/{{sample}}/fastq_screen.{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/{{sample}}/fastq_screen.{{unit}}.r{{read}}.bmk".format(LOGDIR)
    threads: get_resource("fastq_screen","threads")
    resources:
        mem_mb=get_resource("fastq_screen","mem_mb"),
        walltime=get_resource("fastq_screen","walltime")
    params:
        fastq_screen_config="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["fastq_screen"]["index_dir"]),
        subset=100000,
        aligner='bowtie2'
    wrapper:
        "v1.23.4/bio/fastq_screen"

def multiqc_input(wc):
    f=expand("{OUTDIR}/{{sample}}/qc/fastqc/{unit.sample}.{unit.unit}.r{read}_fastqc.zip", unit=units.itertuples(), read=('1','2','3'), OUTDIR=OUTDIR)
    try: 
        if config["parameters"]["fastq_screen"]["enabled"]:
            f +=expand("{OUTDIR}/{{sample}}/qc/fastq_screen/{unit.sample}.{unit.unit}.r{read}.fastq_screen.txt",
unit=unit.itertuples(), OUTDIR=OUTDIR, read=('1', '2', '3'))
    except KeyError:
        print("FASTQ_SCREEN disabled by config file. Skipping...")
    return f

rule multiqc:
    input:
        multiqc_input
    output:
        report(f"{OUTDIR}/{{sample}}/qc/multiqc_report.html", category="1_QC")
    params: 
        config["multiqc"]
    benchmark:
        "{}/{{sample}}/multiqc.bmk".format(LOGDIR)
    log:
        "{}/{{sample}}/multiqc.log".format(LOGDIR)
    threads: get_resource("multiqc", "threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        walltime=get_resource("multiqc","walltime")
    wrapper:
        "v1.23.3/bio/multiqc"
