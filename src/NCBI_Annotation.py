#!/usr/bin/env python

import urllib.parse as up
from glob import glob
from os import makedirs, path
from time import sleep

import click
import numpy as np
import pandas as pd
import requests


""" For more details
> https://stackoverflow.com/questions/17509607/submitting-to-a-web-form-using-python
> https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/help.html
> https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/
"""


def result_parse():
    pass


@click.command()
@click.option(
    "-inf",
    help="Contigs fasta containing folder",
    type=str,
    default="/home/devil/Documents/Mo",
    show_default=True)
@click.option(
    "-otf", help="Output folder", type=str, default="./", show_default=True)
@click.option(
    "-org",
    help="Organism name",
    type=click.Choice(["Flu"]),
    default="Flu",
    show_default=True)
@click.option(
    "-fld",
    help="Failed id list file",
    type=str,
    default="Failed.txt",
    show_default=True)
def run(inf, otf, org, fld):
    """Script is to generate annotation based on available data on NCBI Server.
    Please do not use it for large number of dataset."""
    if org == "Flu":
        root = "https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/"

    job_ids = []
    job_ids_dict = {}
    if not path.exists(otf) or not path.isdir(otf):
        exit("Given output folder is not provided. Exiting")

    for fl in glob("%s/*" % inf):
        samp = path.split(fl)[1].rsplit('.', 1)[0]
        if path.isfile("%s/%s/__done__" % (otf, samp)):
            continue
        click.echo("Task submitted for %s" % samp)

        # TODO: add option that the task has beed performed
        sequences = open(fl).read()
        data = {"sequence": sequences, "cmd": 'submit'}
        data = up.urlencode(data)
        req_out = requests.post(root, data=data)
        job_ids_dict[req_out.json()["job_id"]] = samp
        job_ids.append(req_out.json()["job_id"])
        sleep(10)
    # TODO: Generate a log
    # TODO: 10 attempts for Pending files
    datax = pd.DataFrame({
        'task': list(job_ids_dict.keys()),
        'samp': list(job_ids_dict.values())
    })
    datax.to_csv("tasks.tab", index=False, sep="\t")
    del datax

    failed_list = open(fld, "w")

    while len(job_ids):
        for job_id in job_ids:
            if not path.exists("%s/%s" % (otf, job_ids_dict[job_id])):
                makedirs("%s/%s" % (otf, job_ids_dict[job_id]))
            data = {"job_id": job_id, "cmd": "status"}
            data = up.urlencode(data)
            req_out = requests.post(root, data=data)
            status = req_out.json()["status"]
            if status in [
                    "NotFound", "Canceled", "Failed", "Done", "Reading",
                    "ReadFailed", "Deleted"
            ]:
                job_ids.remove(job_id)
                if status == "Done":
                    for fmt in [
                            "tbl", "tbl_clean", "gbf", "sqn", "xml", "faa",
                            "ffn", "aln"
                    ]:
                        outputfile = open("%s/%s/%s.%s" %
                                          (otf, job_ids_dict[job_id],
                                           job_ids_dict[job_id], fmt), "w")
                        download = {
                            "job_id": job_id,
                            "cmd": "download",
                            "format": fmt
                        }
                        download = up.urlencode(download)
                        req_out = requests.post(root, data=download)
                        outputfile.write(req_out.text)  # json()
                        outputfile.close()
                        sleep(10)
                    open("%s/%s/__done__" % (otf, job_ids_dict[job_id]),
                         "w").close()
                    click.echo("Task finished for %s" % job_ids_dict[job_id])

                else:
                    failed_list.write("%s\n" % job_ids_dict[job_id])
            else:  # ["Pending", "Running", "Confirmed"]:
                click.echo(
                    "Task for %s added for next cycle" % job_ids_dict[job_id])
                sleep(10)
                continue
    failed_list.close()
    print("Success......")


if __name__ == '__main__':
    run()
