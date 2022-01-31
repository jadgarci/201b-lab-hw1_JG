rule all:
    input:
        "SRR2584857_1.ecoli-rel606.sens.vcf",
        "SRR2584857_1.ecoli-rel606.spec.vcf",
        "SRR2584403_1.ecoli-rel606.sens.vcf",
        "SRR2584403_1.ecoli-rel606.spec.vcf",
        "SRR2584404_1.ecoli-rel606.sens.vcf",
        "SRR2584404_1.ecoli-rel606.spec.vcf",
        "SRR2584405_1.ecoli-rel606.sens.vcf",
        "SRR2584405_1.ecoli-rel606.spec.vcf"


#this rule was changed by rule download reads
#rule download_data:
#    conda: "env-wget.yml"
#    output: "SRR2584857_1.fastq.gz"
#    shell: """
#        wget https://osf.io/4rdza/download -O {output}
#    """


#####------######
##blocks to download multiple samples in single rule
# list sample names & download URLs.
sample_links = { "SRR2584403_1": "https://osf.io/6jx7z/download",
                 "SRR2584404_1": "https://osf.io/s24ky/download",
                 "SRR2584405_1": "https://osf.io/7qek6/download",
                 "SRR2584857_1": "https://osf.io/aksmc/download"}

# the sample names are dictionary keys in sample_links. extract them to a list we can use below
SAMPLES=sample_links.keys()

rule download_all:
    input:
        expand("{sample}.fq.gz", sample=SAMPLES)

# rule to download each individual file specified in sample_links
rule download_reads:
    output: "{sample}.fastq.gz" 
    params:
        # dynamically generate the download link directly from the dictionary
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell: """
        curl -L {params.download_link} -o {output}
        """

######-------######

rule download_genome:
    conda: "env-wget.yml"
    output: "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O {output}"

rule map_reads:
    conda: "env-minimap.yml"
    input: ref="ecoli-rel606.fa.gz", reads="{sample}.fastq.gz"
    output: "{sample}.ecoli-rel606.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.sam"
    output: "{sample}.ecoli-rel606.bam"
    shell: """
        samtools view -b -F 4 {input} > {output}
     """

rule sort_bam:
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.bam"
    output: "{sample}.ecoli-rel606.bam.sorted"
    shell: """
        samtools sort {input} > {output}
    """

rule gunzip_fa:
    input:
        "ecoli-rel606.fa.gz"
    output:
        "ecoli-rel606.fa"
    shell: """  
        gunzip -c {input} > {output}
    """

rule call_variants:
    conda: "env-bcftools.yml"
    input:
        ref="ecoli-rel606.fa",
        bamsort="{sample}.ecoli-rel606.bam.sorted"
    output:
        pileup="{sample}.ecoli-rel606.pileup",
        bcf="{sample}.ecoli-rel606.bcf",
        vcf="{sample}.ecoli-rel606.sens.vcf"
    shell: """
        bcftools mpileup -Ou -f {input.ref} {input.bamsort} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """

rule call_variants_spec:
    conda: "env-bcftools.yml"
    input:
        "{sample}.ecoli-rel606.sens.vcf"
    output:
        "{sample}.ecoli-rel606.spec.vcf"
    shell: """
        bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' {input} > {output}
    """
