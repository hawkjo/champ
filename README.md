# CHAMP: Chip-Hybridized Association Mapping Platform

### Installation

CHAMP has only been run on Ubuntu. You'll need a few dependencies first:
```
sudo apt install -y git sextractor samtools bowtie2 virtualenv python-dev zlib1g-dev
git clone https://gitlab.com/jimrybarski/ngs_project.git
```

Optionally, you can install into a virtual environment (recommended):

```
virtualenv env
. env/bin/activate
```

Now install Python packages:

```
pip install numpy && pip install scipy && pip install -r ngs_project/requirements.txt

```

### Typical Pipeline

#### Mapping Reads

When a new chip is received, it needs to be analyzed once to determine which reads are fiducial markers (that is,
clusters of phiX genomic DNA) and which are to be used in the experiment. Sometimes this binary classification is not
enough - for example, you may have several different targets of interest on the same chip. In this case, you just pass
in additional BAM files for each target, in addition to the mandatory one for phiX.

In the example below, our chip has two targets that we care about, plus phiX. It may have some other reads that belong
to none of the groups, and those will get thrown into a catch-all file called "unclassified". Reads in this file will
still be used for later calculations.

`champ map SA16032/all_fastqs SA16032/mapped_reads ~/my_bamfiles/phix ~/my_bamfiles/target_1 ~/my_bamfiles/target_2`

#### New Image Data

When a new experiment is run and the image files are uploaded to the server, you'll need to run "champ init" to
associate some metadata about the experiment with the image files. There are several mandatory pieces of information and
 some optional ones:

`IMAGE_DIRECTORY` the directory that contains all of the TIF files (and will contain the HDF5 files)

`CHIP_NAME` the unique ID of the chip. This will show up in some graphs, but is mostly for recordkeeping.

`MAPPED_READS` the path to the directory that contains the text files produced by the `champ map` command.

`ALIGNMENT_CHANNEL` the name of the color channel that phiX is visible in. We actually recommend (require?) that a
subset of phiX be labeled in channels that proteins are visible in, to help with alignment of low concentrations, but
this refers specifically to the channel where 100% of phiX clusters are visible.

`--microns-per-pixel` the size of a side of one pixel, in microns

`--chip` the type of chip, either `miseq` or `hiseq`

`--ports-on-right` the stage adapters we use orient the two input ports on the chip to the left or right, depending on
which microscope we use. This needs to be known for alignment. By default, we assume they are on the left unless this
flag is passed.

`--flipud` invert all images through the horizontal axis (i.e. run numpy.flipud() on all images). This was added to handle
a quirk with the way MicroManager saves images.

`--fliplr` invert all images through the vertical axis (i.e. run numpy.fliplr() on all images). We have not needed this,
but if your images don't align you may try passing it in.

#### Preprocessing




####
