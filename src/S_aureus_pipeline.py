from ruffus import *
import matplotlib
matplotlib.use("Agg")
from pylab import *
import seaborn as sns
import pandas as pd
import click
from Bio import SeqIO
from subprocess import Popen
from glob import glob
from os import path, system, makedirs, rename, environ
from shutil import copy


@click.command()
# - - - - - - Standard Options - - - - - - - -
@click.option(
    "-fqgzd",
    help="Fastq.gz containing folder",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-outd",
    help="Output directory",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-fd",
    help="Forward mate represetation",
    type=str,
    default='_1',
    show_default=True)
@click.option(
    "-rv",
    help="Reverse mate represetation",
    type=str,
    default='_2',
    show_default=True)
@click.option(
    "-sfx",
    help="Fastq Prefix",
    type=str,
    default="fastq.gz",
    show_default=True)
@click.option(
    "-ncr",
    help="Number of processes to be used",
    type=int,
    default=1,
    show_default=True)
# - - - - - - SRST2 Options - - - - - - - -
# srst2 was installed with venv py27
# symbolic link was generate to bin folder to make it in the path
# ln -s envs/py27/bin/srst2 bin/srst2
# ln -s envs/py27/bin/getmlst.py bin/getmlst.py
# and exports
# export SRST2_SAMTOOLS="envs/py27/bin/samtools"
# export SRST2_BOWTIE2="envs/py27/bin/bowtie2"
# export SRST2_BOWTIE2_BUILD="envs/py27/bin/bowtie2-build"
def run(fqgzd, outd, fd, rv, sfx, ncr):
    """
    NOTE: This is still under development. Use on your own risk.\n
    This pipeline is for analysing short read sequences from S. aureus.
    This pipeline uses following tools to perform different task under differnt
    functions as following.\n
    Tools           Tasks                   Functions\n
    -------------------------------------------------------------------------\n
    FastQC          Reads quality check     fastqc\n
    TrimeGalore     Removing bad quality bases and reads  trimglore\n
    Velvet          Constructing Conntigs   velvet\n
    Prokka          Annotating contigs      prokka\n
    Kraken          Species check           kraken\n
    SRST2           mlst and antibiotic     srst2\n

    Notes for SRST2\n
    ---------------\n
    SRSR2 related tools should be added in the paths as follows\n\n
    export SRST2_SAMTOOLS="/.anmol/anaconda3/envs/py27/bin/samtools"\n
    export SRST2_BOWTIE2="/.anmol/anaconda3/envs/py27/bin/bowtie2"\n
    export SRST2_BOWTIE2_BUILD="/.anmol/anaconda3/envs/py27/bin/bowtie2-build"\n

    Addition notes\n
    --------------\n
    data bases should be added in system path as follows\n\n
    export KRAKENDB="/.anmol/databases/minikraken_20171019_8GB"\n
    export ARDB="/.anmol/databases/ARmeta-genes.fa"\n
    """

    # Option checks
    if not fqgzd:
        exit("Fastqgz files containing folder not given. Exiting . . . .")
    if not path.isdir(fqgzd):
        exit("Given Fastqgz files folder path is not a diectory path. "
             "Exiting . . . .")
    #          "Exiting . . . .")
    if not outd:
        exit("Output directory path not given. Exiting . . . .")
    else:
        makedirs(outd, mode=0o777, exist_ok=True)

    # From system path
    kdb = environ.get('KRAKENDB')
    ar_fasta = environ.get('ARDB')
    mlstdb = environ.get('MLSTDB')
    # Ariba mlst
    if path.exists("aribaSaureusmlst") and path.isdir("aribaSaureusmlst"):
        pass
    elif not path.exists("aribaSaureusmlst"):
        system("""ariba pubmlstget "Staphylococcus aureus" aribaSaureusmlst""")

    files = glob("%s/*%s.%s" % (fqgzd, fd, sfx))
    # print("%s/*_%s.%s" % (fqgzd, fd, sfx), files)
    filepaires = [[f, f.replace("%s.%s" % (fd, sfx), "%s.%s" % (rv, sfx))]
                  for f in files]
    # print(filepaires)
    # Create MLST stuffs here

    @follows(mkdir("%s/FastQC" % outd))
    @transform(filepaires, formatter(".+/(?P<filebase>\w+)%s.%s" % (fd, sfx)),
               [
                   "%s/FastQC/{filebase[0]}%s_fastqc.html" % (outd, fd),
                   "%s/FastQC/{filebase[0]}%s_fastqc.html" % (outd, rv)
               ])
    def fastqc(inputfiles, outputfiles):
        for fl in inputfiles:
            fastqc = ['fastqc', fl, "-o", "%s/FastQC" % outd]
            Popen(fastqc).communicate()

    @follows(fastqc, mkdir("%s/TrimGalore" % outd))
    @transform(filepaires, formatter(".+/(?P<filebase>\w+)%s.%s" % (fd, sfx)),
               [
                   "%s/TrimGalore/{filebase[0]}_1.fq.gz" % (outd),
                   "%s/TrimGalore/{filebase[0]}_2.fq.gz" % (outd)
               ])
    def trim_galore(inputfiles, outputfiles):
        "Trimming bad quality bases from the sequences."
        fb = path.split(inputfiles[0])[1].split(".", 1)[0].rsplit('_', 1)[0]
        trim_galore = [
            "trim_galore", "--paired", "-o",
            "%s/TrimGalore" % outd, inputfiles[0], inputfiles[1]
        ]
        Popen(trim_galore).communicate()
        rename("%s/TrimGalore/%s%s_val_1.fq.gz" % (outd, fb, fd),
               "%s/TrimGalore/%s_1.fq.gz" % (outd, fb))
        rename("%s/TrimGalore/%s%s_val_2.fq.gz" % (outd, fb, rv),
               "%s/TrimGalore/%s_2.fq.gz" % (outd, fb))
        # pass

    @follows(trim_galore, mkdir("%s/Velvet" % outd))
    @transform(
        trim_galore, formatter(".+/(?P<filebase>\w+)_1.fq.gz"),
        "%s/Velvet/{filebase[0]}/contigs.fa" % outd)  # Needs some fixing here
    def velvet(inputfiles, outputfile):
        print(inputfiles, outputfile, "anmolkiran")
        vh = "'-shortPaired -fastq.gz -separate %s %s'" % (inputfiles[0],
                                                           inputfiles[1])
        fl_bs = path.split(inputfiles[0])[1].split("_1.fq.gz")[0]
        vel_cmd = [
            "VelvetOptimiser.pl", "-s", "51", "-e", "253", "-p", fl_bs, "-d",
            "%s/Velvet/%s" % (outd, fl_bs), "-f", vh
        ]
        system(" ".join(vel_cmd))

    @follows(trim_galore, mkdir("%s/Kraken" % outd))
    @transform(
        trim_galore, formatter(".+/(?P<filebase>\w+)_1.fq.gz"),
        "%s/Kraken/{filebase[0]}.kraken" % outd)  # Needs some fixing here
    def kraken(inputfiles, outputfile):
        fl_bs = path.split(inputfiles[0])[1].split("_1.fq.gz")[0]
        outfile = "%s/Kraken/%s" % (outd, fl_bs)
        kraken_cmd1 = [
            "kraken", "--threads", "31", "--fastq-input", "--gzip-compressed",
            "-db", kdb, "--output", outfile, "--paired", inputfiles[0],
            inputfiles[1]
        ]
        # print(" ".join(kraken_cmd1))
        kraken_cmd2 = [
            "kraken-report", "--db",
            "/home/anmol/binbag/minikraken_20171019_8GB", outfile
        ]
        Popen(kraken_cmd1).communicate()
        with open(outputfile, "w") as of:
            Popen(kraken_cmd2, stdout=of).communicate()

    @follows(velvet, mkdir("%s/Prokka" % outd))
    @transform(velvet, formatter(".+/contigs.fa"), [
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.err" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.faa" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.ffn" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.fna" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.fsa" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.gbk" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.gff" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.log" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.sqn" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.tbl" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.tsv" % outd,
        "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.txt" % outd,
    ])  # Needs some fixing here
    def prokka(inputfile, outputfiles):
        pkk_cmd = [
            "prokka", inputfile, "--outdir",
            path.split(outputfiles[0])[0], "--force", "--prefix",
            path.split(outputfiles[0])[1].split('.err')[0]
        ]
        Popen(pkk_cmd).communicate()

    @follows(trim_galore, mkdir("%s/MLST" % outd))
    @merge(trim_galore,
           "%s/MLST/srst2__mlst__Staphylococcus_aureus__results.txt" % outd
           )  # Needs some
    def mlst(inputfiles, outputfile):
        """
        Conda trick\n
        -----------\n
        conda create -n py27 python=2.7 samtools=0.1.18 scipy=0.12.0
        numpy=1.7.1 bowtie2=2.2.8\n
        source activate py27\n
        pip install srst2\n
        source deactivate
        """
        srst2_cmd = [
            "srst2", "--output",
            "%s/MLST/srst2" % outd, "--input_pe",
            "%s/TrimGalore/*.fq.gz" % outd, "--mlst_db",
            "%s/S_aureus/Staphylococcus_aureus.fasta" % mlstdb,
            "--mlst_definitions",
            "%s/S_aureus/saureus.txt" % mlstdb, "--mlst_delimiter", "'_'",
            "--threads", "30", "--gene_db", ar_fasta
        ]
        system(" ".join(srst2_cmd))

    @merge(velvet, [
        "%s/Velvet/stats.tsv" % outd,
    ])  # Needs some # "%s/Velvet/stats.png" % outd
    def velvet_stats(inputfiles, outputfiles):
        # TODO : Plan what to add as figure and in the table
        samp_details = {"samp_id": [], "genome_size": [], "mean_coverage": []}
        for fl in inputfiles:
            samp_id = fl.split("/")[-2]
            genome_size = 0
            reads = 0
            for rec in SeqIO.parse(fl, "fasta"):
                frag_depth = float(rec.id.rsplit("_", 1)[1])
                frag_size = len(rec.seq)
                genome_size += frag_size
                reads += frag_depth * frag_size
            samp_details["samp_id"].append(samp_id)
            samp_details["genome_size"].append(genome_size)
            samp_details["mean_coverage"].append(reads / genome_size)
        samp_details = pd.DataFrame.from_dict(samp_details)
        samp_details[["samp_id", "genome_size", "mean_coverage"]].to_csv(
            outputfiles[0], index=False, sep="\t")

        # TODO: Some plotting

    @follows(trim_galore, mkdir("%s/ARIBA_MLST" % outd))
    @transform(trim_galore, formatter(".+/(?P<filebase>\w+)_1.fq.gz"), [
        "%s/ARIBA_MLST/{filebase[0]}/mlst_report.tsv" % outd,
        "%s/ARIBA_MLST/{filebase[0]}/mlst_report.details.tsv" % outd,
    ])  # Needs some
    def ariba_mlst(inputfiles, outputfiles):
        ariba_mlst_cmd = [
            "ariba", "run", "aribaSaureusmlst/ref_db", inputfiles[0],
            inputfiles[1],
            path.split(outputfiles[0])[0]
        ]
        system(" ".join(ariba_mlst_cmd))

    @follows(velvet, mkdir("%s/VirSorter" % outd))
    @transform(velvet, formatter(".+/contigs.fa"),
               "%s/VirSorter/{subdir[0][0]}/{subdir[0][0]}.fa" % outd)
    def copy_velvet_contigs(inputfile, outputfile):
        """FCopy Velvet Files"""
        makedirs(path.split(outputfile)[0], mode=0o777, exist_ok=True)
        copy(inputfile, outputfile, follow_symlinks=True)

    @transform(
        copy_velvet_contigs, formatter(".+/(?P<filebase>\w+).fa"),
        "%s/VirSorter/{subdir[0][0]}/VIRSorter_global-phage-signal.csv" % outd)
    def virSorter(inputfile, outputfile):
        fd, fl = path.split(path.abspath(inputfile))
        vir_cmd = [
            "docker", "run", "--user", "1000:1000", "--cpus", "1", "-v",
            "/.anmol/docker_images/virsorter-data:/data", "-v",
            "%s:/wdir" % fd, "-w", "/wdir", "--rm",
            "discoenv/virsorter:v1.0.3", "--db", "2", "--fna",
            "/wdir/%s" % fl
        ]
        Popen(vir_cmd).communicate()
        # pass

    #
    # def ariba_merge():
    #     pass

    #     # Write your tricks here
    #     pass
    #
    # def emm():
    #     pass

    # TODO: Group them based on the organims. Create Separate folder for the each organism
    # TODO: Sequence typing. It may be multiple. Fetch details from kraken
    # TODO: Emm typing in some cases. Fetch details from kraken
    # TODO: Clustering
    # TODO : Phylogeny. Split paralogs before you doPhylogeny ask for what kind of phylogeny
    # TODO: Add sepration based on Kraken. Add a file satating organims name
    # TODO: Might need to push craken up in the queue

    # pass
    pipeline_run(verbose=1)


if __name__ == '__main__':
    run()

# CMD : python S_aureus_pipeline.py -fqgzd PoolC -outd Anmol -fd _R1 -rv _R2
