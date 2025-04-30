# midas2smol

A code for unpacking MIDAS files directly to SMOL tree files (a minimal-filesize data format used in my [sort codes](https://github.com/e-j-w/GreasySortCodes)). Various bugs are known, user beware.

Based almost entirely on [grif-replay](https://github.com/GRIFFINCollaboration/grif-replay) by Chris Pearson and Adam Garnsworthy.

## Usage

Compile the code using `make`.


```
midas2smol midas_file output_SMOL_tree
```
- `midas_file` should be the path to the first subrun (eg. `run29623_000.mid`)
- `output_SMOL_tree` is the tree file to save to disk, which should have the `.smol` file extension

The data format for output SMOL trees is specfied in [smol-format.h](smol-format.h).

## Known issues
- Number of sorted events is not fully consistent each time the code is run.
- Timing is seemingly not consistent with GRSISort, though seems ok for the most part.