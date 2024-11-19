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
use wget to get files from E-MTAB-6653. Nohup and '&' is to run in the background after closing the chell. It helps to use -nd just to get files and not the massive chain of directories. -nc for 'no-clobber' to not re-download files.

```
nohup bash -c 'wget -nd -nc "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/653/E-MTAB-6653/Files/*"' &
```
