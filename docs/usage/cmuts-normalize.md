

## Outputs

The output of `cmuts normalize` will be an HDF5 file with the following structure:

```
/
├── ROI
├── SNR
├── error
├── heatmap
├── norm
├── reactivity
└── reads
```

If a value was passed to `--group`, then the above will be contained in a group with that name.

**ROI** (Region Of Interest): The region of the RNA excluding the designated primers.

**SNR** (Signal-To-Noise): An estimate of the base-averaged SNR based on the computed error.

**error**: The standard error of the mean.

**heatmap**: The prevalence of each mutation, insertion, and deletion type.

**norm**: The normalization value used.

**reactivity**: The reactivity profiles.

**reads**: The number of reads used to compute each reactivity profile.
