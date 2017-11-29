from ruffus import *
import pandas as pd
import click
from subprocess import Popen
from glob import glob
from os import path, system, makedirs


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
@click.option(
    "-kdb",
    help="Kraken DB folder path",
    type=str,
    default=None,
    show_default=True)
# - - - - - - SRST2 Options - - - - - - - -
@click.option(
    "-py27",
    help="Python 2.7 exe path",
    type=str,
    default="python2.7",
    show_default=True)
@click.option(
    "-srst2", help="Srst2 Path", type=str, default="srst2", show_default=True)
@click.option(
    "-smtl",
    help="Samtools ... Path",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-bwt2",
    help="Bowtie2 ... Path",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-bwt2b",
    help="Bowtie2 build... Path",
    type=str,
    default=None,
    show_default=True)
def run(fqgzd, outd, fd, rv, sfx, ncr, py27, srst2, smtl, bwt2, bwt2b):
    """
    NOTE: This is still under development. Use on your own risk.\n
    This pipeline is for analysing short read sequences from S. aureus.
    This pipeline uses following tools to perform different task under differnt
    functions as following.\n
    Tools           Tasks                   Functions\n
    -------------------------------------------------------------------------\n
    FastQC\n
    TrimeGalore\n
    Velvet\n
    Prokka\n
    Kraken\n

    Notes for SRST2\n
    ---------------\n
    """

    # Option checks
    if not fqgzd:
        exit("Fastqgz files containing folder not given. Exiting . . . .")
    if not path.isdir(fqgzd):
        exit("Given Fastqgz files folder path is not a diectory path. "
             "Exiting . . . .")
    if not kdb:
        exit("Kraken db folder not given. Exiting . . . .")
    if not path.isdir(kdb):
        exit("Given Kraken db folder path is not a diectory path. "
             "Exiting . . . .")
    if not outd:
        exit("Output directory path not given. Exiting . . . .")
    else:
        makedirs(outd, mode=0o777, exist_ok=True)


    #
    files = glob("%s/*%s.%s" % (fqgzd, fd, sfx))
    # print("%s/*_%s.%s" % (fqgzd, fd, sfx), files)
    filepaires = [[f, f.replace("%s.%s" % (fd, sfx), "%s.%s" % (rv, sfx))]
                  for f in files]
    # print(filepaires)

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
                   "%s/TrimGalore/{filebase[0]}%s_val_1.fq.gz" % (outd, fd),
                   "%s/TrimGalore/{filebase[0]}%s_val_2.fq.gz" % (outd, rv)
               ])
    def trim_galore(inputfiles, outputfiles):
        "Trimming bad quality bases from the sequences."
        trim_galore = [
            "trim_galore", "--paired", "-o",
            "%s/TrimGalore" % outd, inputfiles[0], inputfiles[1]
        ]
        Popen(trim_galore).communicate()
        # pass

    #
    #
    @follows(trim_galore, mkdir("%s/Velvet" % outd))
    @transform(
        trim_galore, formatter(".+/(?P<filebase>\w+)%s_val_1.fq.gz" % fd),
        "%s/Velvet/{filebase[0]}/contigs.fa" % outd)  # Needs some fixing here
    def velvet(inputfiles, outputfile):
        vh = "'-shortPaired -fastq.gz -separate %s %s'" % (inputfiles[0],
                                                           inputfiles[1])
        fl_bs = path.split(inputfiles[0])[1].split("%s_val_1.fq.gz" % fd)[0]
        vel_cmd = [
            "VelvetOptimiser.pl", "-s", "51", "-e", "253", "-p", fl_bs, "-d",
            "%s/Velvet/%s" % (outd, fl_bs), "-f", vh
        ]
        system(" ".join(vel_cmd))

    @follows(trim_galore, mkdir("%s/Kraken" % outd))
    @transform(
        trim_galore, formatter(".+/(?P<filebase>\w+)%s_val_1.fq.gz" % fd),
        "%s/Kraken/{filebase[0]}.kraken" % outd)  # Needs some fixing here
    def kraken(inputfiles, outputfile):
        fl_bs = path.split(inputfiles[0])[1].split("%s_val_1.fq.gz" % fd)[0]
        outfile = "%s/Kraken/%s" % (outd, fl_bs)
        kraken_cmd1 = [
            "kraken", "--threads", "31", "--fastq-input", "--gzip-compressed",
            "-db", "/home/anmol/binbag/minikraken_20171019_8GB", "--output",
            outfile, "--paired", inputfiles[0], inputfiles[1]
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

    def mlst():
        """
        Conda trick\n
        -----------\n
        conda create -n py27 python=2.7 samtools=0.1.18 scipy=0.12.0
        numpy=1.7.1 bowtie2=2.2.8\n

        SRST2 Trick\n
        -----------\n


        """
        # Write your tricks here
        pass

    def emm():
        pass

    # TODO: Group them based on the organims. Create Separate folder for the each organism
    # TODO: Sequence typing. It may be multiple. Fetch details from kraken
    # TODO: Emm typing in some cases. Fetch details from kraken
    # TODO: Clustering
    # TODO : Phylogeny. Split paralogs before you doPhylogeny ask for what kind of phylogeny
    # TODO: Add sepration based on Kraken. Add a file satating organims name
    # TODO: Might need to push craken up in the queue

    # pass
    pipeline_run(verbose=9)


if __name__ == '__main__':
    run()
