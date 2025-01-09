# cancer-kinetics-

Here we'll perform scRNAseq primary and secondary anaylsis on a lung tumor dataset from [https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6653?query=%20E-MTAB-6149%20](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6653?query=%20E-MTAB-6149%20). Here we'll use google compute and terra.bio or nextflow.io to host and perform primary analysis. Secondary analysis will peformed locally using the Seurat package.

## Set-Up

In a google VM, create a bucket and mount it to a directory.
```
export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
echo "deb https://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
sudo apt-get update
sudo apt-get install fuse gcsfuse
gcloud auth application-default login
```

Authenticate gfuse and mount bucket

```
gcsfuse bc-scseq-data "$HOME/mybucket"
```
Use wget to get files from E-MTAB-6653. Nohup and '&' is to run in the background after closing the shell. It helps to use -nd just to get files and not the massive chain of directories. -nc for 'no-clobber' to not re-download files. Make sure you have permission to write to buckets or else this will just silently fail. [https://cloud.google.com/filestore/docs/copying-data](https://cloud.google.com/filestore/docs/copying-data)

```
nohup bash -c 'wget -nd -nc "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/653/E-MTAB-6653/Files/*"' &
```

### STAR-SOLO in Terra.bio
In order to get the data into terra, we need to create a workspace, note the bucket name, change permissions on the source account to access the buckets and copy the data into the local bucket.Im using rsync here to conviently move data around and have a copy of everything.

```
gcloud storage rsync --no-clobber gs://bc-scseq-data gs://fc-970a7650-2cd9-4d45-bed4-0580965a6184/sc-data
```
Now we need to set up STAR-SOLO to run in terra. The workflow uses a sample file as such: 


|     Sample    |    Reference  |                        Location                         |  Assay  |
| ------------- |:-------------:| -------------------------------------------------------:| -------:|
| scrBT1429m    | GRCh38-2020-A | gs://fc-970a7650-2cd9-4d45-bed4-0580965a6184/sc-data/s0 | tenX_5p |


After we're done we can pull the h5 files from the bucket to local. Here I've used a bash script to pull the files i want. I've saved the filenames to file via gcloud.

```
gcloud storage ls -r  gs://fc-970a7650-2cd9-4d45-bed4-0580965a6184/star-out | grep 'h5' > paths.txt  
```

And then download from them:
```
#!/bin/bash 

while IFS= read -r line
do
    prefix=`echo $line | pcregrep -o1 '(BT\\d+|scrBT\\d+)' -`
    filter="$(echo $line | pcregrep -o1 '(raw|filtered)' -)"
    filename="${prefix}_${filter}"
    filename+=".h5"
    echo $filename
    gcloud storage cp "${line}" "${filename}"
done < paths.txt
```

Which is executed via:
```
bash pull_h5.sh
```


### STAR-SOLO using NF-CORE

Alternatively we can use the nfcore pipelines in conjuction with google gcp. This gives quite a bit more control and less of hassle with moving stuff to different buckets.

First we make a local directory on a machine with nextflow installed. Follow instructions here for setting up glcoud: [https://seqera.io/blog/nextflow-with-gbatch/](https://seqera.io/blog/nextflow-with-gbatch/).

We are going to use this NF pipeline: [https://nf-co.re/scrnaseq/2.7.1](https://nf-co.re/scrnaseq/2.7.1)

Next we make a sample sheet, similar to the one above:

|     Sample    |    fastq1                                                   |                        fastq2                               |  expected_cells  |
| ------------- |:-----------------------------------------------------------:| -----------------------------------------------------------:| ----------------:|
| scrBT1429m    |gs://bc-scseq-data/fastqs/scrBT1429m_S0_L001_R1_001.fastq.gz |gs://bc-scseq-data/fastqs/scrBT1429m_S0_L001_R2_001.fastq.gz |       4000       |

Notice we cant mix assay types and references in this one. Not really an issue for our use case.

The nextflow.config file will look like this:
```
workDir = 'gs://bc-scseq-data/scrnaseq-out'

process {
executor = 'google-batch'
container = 'nextflow/rnaseq-nf'
errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
maxRetries = 5
}

google {
project = 'operating-land-440016-q9'
location = 'us-central1'
batch.spot = true
}

File: nextflow.config
workDir = 'gs://bc-scseq-data/scrnaseq-out'

process {
executor = 'google-batch'
container = 'nfcore/scrnaseq'
errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
maxRetries = 5
}
```

Lastly, we should save the command in a small script for reproducibility:

```
File: pipe.sh
#!/bin/bash

nextflow run nf-core/scrnaseq  -profile docker -c nextflow.config  --outdir ./pipe-out/ \
--input samplesheet.csv \
--fasta gs://bc-scseq-data/genomes/hg38.fa \
--gtf gs://bc-scseq-data/genomes/hg38.knownGene.gtf \
--protocol 10XV2 \
--aligner star
```

## Post-Processing

QC and secondary analysis will be performed in a Rmd and saved here in [clustering.md](https://github.com/Benjamin-R-Clark/cancer-kinetics/blob/main/scrnaseq/clustering.md).
