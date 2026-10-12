# Cluster notes

An incomplete list of the clusters `cmuts` has been built and run on. Since `cmuts plot` requires a browser to view the report, it should be run locally.

## Stanford Sherlock

### Pre-Built

Only `cmuts align` has dependencies. minimap2 and samtools are loaded as modules:

```sh
ml load biology minimap2/2.30 samtools/1.16.1
```

For paired-end data, `vsearch` must be installed separately; see the [download instructions](https://github.com/torognes/vsearch#download-and-install). Sherlock needs the Linux x86_64 build whose name ends in `-ubi7`.

### Source Builds

The `cmuts align` dependencies listed above are unchanged. Running also needs HDF5 and htslib, the latter provided by the samtools module:

```sh
ml load hdf5/1.12.0
ml load biology samtools/1.16.1
```

In addition, building needs a later GCC:

```sh
ml load gcc/12.4.0
```

## HHMI Janelia

Running needs htslib on the path and its libraries alongside:

```sh
ml load samtools
export LD_LIBRARY_PATH=/misc/local/samtools-1.22.1/lib:$LD_LIBRARY_PATH
```

Building needs htslib where `pkg-config` will find it:

```sh
export PKG_CONFIG_PATH=/misc/local/samtools-1.22.1/lib/pkgconfig:$PKG_CONFIG_PATH
```
