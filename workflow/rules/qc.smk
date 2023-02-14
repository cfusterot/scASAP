import glob
import os

FASTQDIR=config["cellranger"]["path_to_outs"]
SAMPLES=samples["sample"]

def fastqc_input(wc):
    f=os.path.join(FASTQDIR, units.loc[(wc.sample,wc.unit)]['fq' + wc.read])
    return f

rule fastqc:
    input:
        fastqc_input
    output:
        html="{}/qc/fastqc/{{sample}}/{{sample}}.{{unit}}.r{{read}}_fastqc.html".format(OUTDIR),
        zip="{}/qc/fastqc/{{sample}}/{{sample}}.{{unit}}.r{{read}}_fastqc.zip".format(OUTDIR)
# the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        lambda wc: "-t {}".format(get_resource("fastqc","threads")) 
        #expand("--quiet --outdir '{{OUTDIR}}/qc/{sample}'", sample=SAMPLES) 
    resources:
        mem_mb=get_resource("fastqc", "mem_mb"),
        walltime=get_resource("fastqc", "walltime")
    log:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.bmk".format(LOGDIR)
    threads: 
        threads=get_resource("fastqc", "threads")
    wrapper:
        "v1.23.1/bio/fastqc"

def multiqc_input(wc):
    f=expand("{OUTDIR}/qc/fastqc/{{sample}}/{unit.sample}.{unit.unit}.r{read}_fastqc.zip", unit=units.itertuples(), read=('1','2','3'), OUTDIR=OUTDIR)
    return f

rule multiqc:
    input:
        multiqc_input
    output:
        report(f"{OUTDIR}/qc/{{sample}}/multiqc_report.html", category="1_QC")
    params: 
        config["multiqc"]
    benchmark:
        "{}/multiqc/{{sample}}/multiqc.bmk".format(LOGDIR)
    log:
        "{}/multiqc/{{sample}}/multiqc.log".format(LOGDIR)
    threads: get_resource("multiqc", "threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        walltime=get_resource("multiqc","walltime")
    wrapper:
        "v1.0.0/bio/multiqc"
