from glob import glob
from os import environ, makedirs, path, rename, symlink, system, unlink
from shutil import copy, rmtree
from subprocess import Popen

import click
import matplotlib
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from pylab import *
from ruffus import *

matplotlib.use("Agg")


@click.command()
# - - - - - - Standard Options - - - - - - - -
@click.option(
    "-fqgzd",
    help="Fastq.gz containing folder",
    type=str,
    default="/data/Data/MRSA/AllSamples",
    show_default=True)
@click.option(
    "-outd",
    help="Output directory",
    type=str,
    default="/data/Charlotte/Results/MRSA",
    show_default=True)
@click.option(
    "-fd",
    help="Forward mate represetation",
    type=str,
    default='_R1',
    show_default=True)
@click.option(
    "-rv",
    help="Reverse mate represetation",
    type=str,
    default='_R2',
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
    #  # Ariba mlst
    #  if path.exists("aribaSaureusmlst") and path.isdir("aribaSaureusmlst"):
    #      pass
    #  elif not path.exists("aribaSaureusmlst"):
    #      system("""ariba pubmlstget "Staphylococcus aureus" aribaSaureusmlst""")
    #
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
    #      # pass
    #
    @follows(trim_galore, mkdir("%s/Kraken" % outd))
    @transform(
        trim_galore, formatter(".+/(?P<filebase>\w+)_1.fq.gz"),
        "%s/Kraken/{filebase[0]}.kraken" % outd)  # Needs some fixing here
    def kraken(inputfiles, outputfile):
        fl_bs = path.split(inputfiles[0])[1].split("_1.fq.gz")[0]
        outfile = "%s/Kraken/%s" % (outd, fl_bs)
        kraken_cmd1 = [
            "kraken", "--threads", "1", "--fastq-input", "--gzip-compressed",
            "--db", kdb, "--output", outfile, "--paired", inputfiles[0],
            inputfiles[1]
        ]
        print(" ".join(kraken_cmd1))
        kraken_cmd2 = [
            "kraken-report", "--db",
            kdb, outfile
        ]
        Popen(kraken_cmd1).communicate()
        with open(outputfile, "w") as of:
            Popen(kraken_cmd2, stdout=of).communicate()
    #
    #  @follows(velvet, mkdir("%s/Prokka" % outd))
    #  @transform(velvet, formatter(".+/contigs.fa"), [
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.err" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.faa" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.ffn" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.fna" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.fsa" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.gbk" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.gff" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.log" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.sqn" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.tbl" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.tsv" % outd,
    #      "%s/Prokka/{subdir[0][0]}/{subdir[0][0]}.txt" % outd,
    #  ])  # Needs some fixing here
    #  def prokka(inputfile, outputfiles):
    #      pkk_cmd = [
    #          "prokka", inputfile, "--outdir",
    #          path.split(outputfiles[0])[0], "--force", "--prefix",
    #          path.split(outputfiles[0])[1].split('.err')[0]
    #      ]
    #      Popen(pkk_cmd).communicate()
    #
    @merge(kraken, "%s/s_aureus.tsv"%outd)
    def s_aureus(inputfiles, outputfile):
        res = []
        for fl in inputfiles:
            fl_bs = path.split(fl)[1].split(".kraken")[0]
            data = pd.read_csv(fl, header=None, sep="\t")
            data = data.drop_duplicates(3)
            data = data.loc[data[3]=="S", [0,5]]
            data["samp"] = fl_bs
            res.append(data)
        res = pd.concat(res)
        res[5] = res[5].apply(str.strip)
        res = res[res[5]=="Staphylococcus aureus"]
        res = res.rename(columns={0:'read_pt', 5:'species'})
        res.to_csv(outputfile, index=False, header=True)


    @mkdir("%s/S_aureus_Samples" % outd)
    @subdivide(s_aureus, formatter(),"%s/S_aureus_Samples/*.fq.gz" % outd)
    def s_aureus_split(inputfile, outputfiles):
        data = pd.read_csv(inputfile, usecols=["samp"])
        for i_d in data["samp"]:
            try:
                symlink(f"{outd}/TrimGalore/{i_d}_1.fq.gz",
                       f"{outd}/S_aureus_Samples/{i_d}_1.fq.gz")
            except:
                pass
            try:
                symlink(f"{outd}/TrimGalore/{i_d}_2.fq.gz",
                        f"{outd}/S_aureus_Samples/{i_d}_2.fq.gz")
            except:
                pass

    @mkdir("%s/Velvet" % outd)
    @transform(
        s_aureus_split, formatter(".+/(?P<filebase>\w+)_1.fq.gz"),
        "%s/Velvet/{filebase[0]}.fa" % outd)  # Needs some fixing here
    def velvet(inputfile, outputfile):
        print(inputfile)
        vh = "'-shortPaired -fastq.gz -separate %s %s'" % (inputfile,
                                                           inputfile.replace("_1.fq","_2.fq"))
        fl_bs = path.split(inputfile)[1].split("_1.fq.gz")[0]
        vel_cmd = [
                "/usr/bin/perl", "/usr/bin/velvetoptimiser", "-t", "1",
            "-s", "51", "-e", "250", "-p", fl_bs, "-d",
            fl_bs, "-f", vh
        ]
        system(" ".join(vel_cmd))

        # To generate empty file in case velvet fails to generate contigs
        if not path.exists(f"{fl_bs}/contigs.fa"):
            try:
                makedirs(fl_bs)
            except:
                pass
            open(f"{fl_bs}/contigs.fa","w").close()

        # Moving generated contigs to to desired folder
        system(f"cp {fl_bs}/contigs.fa {outputfile}")
        # removing all intemediate files
        try:
            system(f"rm -rf {fl_bs}*")
        except:
            pass
        #  if not path.exists(f"{fl_bs}")
    #


    @mkdir("%s/MLST" % outd)
    @merge(s_aureus_split,
           "%s/MLST/srst2__mlst__Staphylococcus_aureus__results.txt" % outd
           )  # Needs some
    def mlst(inputfiles, outputfile):
        srst2_cmd = [
            "$SRST2", "--output",
            "%s/MLST/srst2" % outd, "--input_pe",
            "%s/S_aureus_Samples/*.fq.gz" % outd, "--mlst_db",
            "%s/S_aureus/Staphylococcus_aureus.fasta" % mlstdb,
            "--mlst_definitions",
            "%s/S_aureus/saureus.txt" % mlstdb, "--mlst_delimiter", "'_'",
            "--threads", "25", "--gene_db", ar_fasta
        ]
        #  print(" ".join(srst2_cmd))
        system(" ".join(srst2_cmd))

    @mkdir("%s/ariba_db" % outd)
    @originate("%s/ariba_db/out.card.fa" % outd)
    def ariba_db_download(outputfile):
        db_down_cmd = ["ariba", "getref", "card", outputfile.split(".fa")[0]]
        system(" ".join(db_down_cmd))

    @transform(ariba_db_download, suffix(".fa"),
            ".prepareref/01.filter.check_metadata.log")
    def ariba_prepref(inputfile, outputfile):
        print(inputfile, outputfile)
        prepref_cmd = ["ariba", "prepareref",
                        "-f", inputfile,
                        "-m", inputfile.replace(".fa",".tsv"),
                        path.split(outputfile)[0]]
        system(" ".join(prepref_cmd))

    @mkdir("%s/ariba_run" % outd)
    @transform(s_aureus_split, formatter(".+/(?P<filebase>\w+)_1.fq.gz"),
            add_inputs(ariba_prepref),"%s/ariba_run/{filebase[0]}/report.tsv" % outd)
    def ariba_run(inputfiles, outputfile):
        print(inputfiles, outputfile)
        arb_rn_cmd = ["ariba", "run", path.split(inputfiles[1])[0],
                        inputfiles[0], inputfiles[0].replace("_1.fq", "_2.fq"),
                        path.split(outputfile)[0]]
        system(" ".join(arb_rn_cmd))
        if not path.exists(outputfile):
            open(outputfile,"w").close()

    @mkdir(f"{outd}/ariba_summary")
    @merge(ariba_run, [f"{outd}/ariba_summary/s_aureus.summary.csv",
                        f"{outd}/ariba_summary/s_aureus.summary.phandango.csv",
                        f"{outd}/ariba_summary/s_aureus.summary.phandango.tre"])
    def ariba_summary(inputfiles, outputfiles):
        print(inputfiles,outputfiles)
        smmry_cmd = ["ariba", "summary", outputfiles[0].split(".csv")[0]]
        for inputfile in inputfiles:
            smmry_cmd.append(inputfile)

        system(" ".join(smmry_cmd))

    #
    #  @merge(velvet, [
    #      "%s/Velvet/stats.tsv" % outd,
    #  ])  # Needs some # "%s/Velvet/stats.png" % outd
    #  def velvet_stats(inputfiles, outputfiles):
    #      # TODO : Plan what to add as figure and in the table
    #      samp_details = {"samp_id": [], "genome_size": [], "mean_coverage": []}
    #      for fl in inputfiles:
    #          samp_id = fl.split("/")[-2]
    #          genome_size = 0
    #          reads = 0
    #          for rec in SeqIO.parse(fl, "fasta"):
    #              frag_depth = float(rec.id.rsplit("_", 1)[1])
    #              frag_size = len(rec.seq)
    #              genome_size += frag_size
    #              reads += frag_depth * frag_size
    #          samp_details["samp_id"].append(samp_id)
    #          samp_details["genome_size"].append(genome_size)
    #          samp_details["mean_coverage"].append(reads / genome_size)
    #      samp_details = pd.DataFrame.from_dict(samp_details)
    #      samp_details[["samp_id", "genome_size", "mean_coverage"]].to_csv(
    #          outputfiles[0], index=False, sep="\t")
    #
        # TODO: Some plotting

    #  @follows(trim_galore, mkdir("%s/ARIBA_MLST" % outd))
    #  @transform(trim_galore, formatter(".+/(?P<filebase>\w+)_1.fq.gz"), [
    #      "%s/ARIBA_MLST/{filebase[0]}/mlst_report.tsv" % outd,
    #      "%s/ARIBA_MLST/{filebase[0]}/mlst_report.details.tsv" % outd,
    #  ])  # Needs some
    #  def ariba_mlst(inputfiles, outputfiles):
    #      ariba_mlst_cmd = [
    #          "ariba", "run", "aribaSaureusmlst/ref_db", inputfiles[0],
    #          inputfiles[1],
    #          path.split(outputfiles[0])[0]
    #      ]
    #      system(" ".join(ariba_mlst_cmd))
    #
    #  @follows(velvet, mkdir("%s/VirSorter" % outd))
    #  @transform(velvet, formatter(".+/contigs.fa"),
    #             "%s/VirSorter/{subdir[0][0]}/{subdir[0][0]}.fa" % outd)
    #  def copy_velvet_contigs(inputfile, outputfile):
    #      """FCopy Velvet Files"""
    #      makedirs(path.split(outputfile)[0], mode=0o777, exist_ok=True)
    #      copy(inputfile, outputfile, follow_symlinks=True)
    #
    #  @transform(
    #      copy_velvet_contigs, formatter(".+/(?P<filebase>\w+).fa"),
    #      "%s/VirSorter/{subdir[0][0]}/VIRSorter_global-phage-signal.csv" % outd)
    #  def virSorter(inputfile, outputfile):
    #      fd, fl = path.split(path.abspath(inputfile))
    #      vir_cmd = [
    #          "docker", "run", "--user", "1000:1000", "--cpus", "1", "-v",
    #          "/.anmol/docker_images/virsorter-data:/data", "-v",
    #          "%s:/wdir" % fd, "-w", "/wdir", "--rm",
    #          "discoenv/virsorter:v1.0.3", "--db", "2", "--fna",
    #          "/wdir/%s" % fl
    #      ]
    #      Popen(vir_cmd).communicate()
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
    pipeline_run(verbose=1, multithread=25)


if __name__ == '__main__':
    run()

# CMD : python S_aureus_pipeline.py -fqgzd PoolC -outd Anmol -fd _R1 -rv _R2
