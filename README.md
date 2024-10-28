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
Install ftp

```
sudo apt-get install ftp
```

Use the handy ftp script from ArrayExpress, here we're using E-MTAB-6653. 

```
sh download.sh &
```
