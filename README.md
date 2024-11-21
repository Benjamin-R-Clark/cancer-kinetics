# cancer-kinetics-

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
In order to get the data into terra, we need to create a workspace, note the bucket name, change permissions on the source account to access the buckets and copy the data into the local bucket.Im using rsync here to conviently move data around and have a copy of everything.

```
gcloud storage rsync --no-clobber gs://bc-scseq-data gs://fc-970a7650-2cd9-4d45-bed4-0580965a6184/sc-data
```
