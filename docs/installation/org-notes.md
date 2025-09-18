
The following is an incomplete list of various locations where `cmuts` has been run, and the specific requirements for each one.

## Stanford Sherlock Cluster

Both building and running will require the following, which covers all dependencies.

```
ml load hdf5/1.14.4
ml load biology samtools/1.16.1
ml load cmake/3.31.4
```

## HHMI Janelia Cluster

Running will require the following two commands.

```
ml load samtools
LD_LIBRARY_PATH=/misc/local/samtools-1.22.1/lib:$LD_LIBRARY_PATH
```

In addition, building `cmuts` will require

```
ml load cmake/4.0.2
export PKG_CONFIG_PATH=/misc/local/samtools-1.22.1/lib/pkgconfig:$PKG_CONFIG_PATH
```
