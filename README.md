# CHAMP: Chip-Hybridized Affinity Mapping Platform

This software was used for analyses described in the manuscript:

### Massively parallel biophysical analysis of a CRISPR-Cas complex on repurposed next generation sequencing chips
##### Cheulhee Jung, John Hawkins, Stephen K. Jones Jr, Yibei Xiao, James Rybarski, Kaylee Dillard, Jeffrey Hussmann, Mashelfatema A Saifuddin, Cagri A. Savran,  Andrew Ellington, Ailong Ke, William H. Press, and Ilya J. Finkelstein
##### Submitted 2016


### Installation

CHAMP has only been run on Ubuntu. You'll need a few dependencies first:

```
sudo apt install -y build-essential git sextractor samtools bowtie2 virtualenv python-dev zlib1g-dev
git clone https://github.com/hawkjo/champ.git
```

Optionally, you can install into a virtual environment (recommended):

```
cd champ
virtualenv env
. env/bin/activate
```

Now install Python packages and CHAMP:

```
pip install numpy==1.11.1 && pip install scipy==0.18.0 && pip install -r requirements.txt && python setup.py install

```

### Typical Pipeline

#### Mapping Reads

When a new chip is received, it needs to be analyzed once to determine which reads are fiducial markers (that is,
clusters of phiX genomic DNA) and which are to be used in the experiment. Target sequences are kept in a file in 
YAML format. Short read alignment files produced by Bowtie2 are used to classify genomic DNA. In the example below,
the phiX files are in a directory called `phix_bowtie` and the prefix `phix` is provided since the files all begin with
that (this is what Bowtie requires). `min-len` and `max-len` refer to the minimum and maximum length of the sequences
of interest (note that for CRISPR systems, this length includes the PAM).

`champ map SA16032/all_fastqs SA16032/read_names --target-sequence-file targets.yml --phix-bowtie phix_bowtie/phix --min-len 24 --max-len 46`

#### Setting Up a New Analysis

When a new experiment is run and the image files are uploaded to the server, you'll need to run `champ init` to
associate some metadata about the experiment with the image files. There are several mandatory pieces of information and
 some optional ones. This creates a file `champ.yml` that holds this metadata, and which is used during the alignment process
 to checkpoint progress.

`IMAGE_DIRECTORY` the directory that contains all of the TIF files (and will contain the HDF5 files)

`READ_NAMES_DIRECTORY` the path to the directory that contains the text files produced by the `champ map` command.

`ALIGNMENT_CHANNEL` the name of the color channel that phiX is visible in. We actually recommend (require?) that a
subset of phiX be labeled in channels that proteins are visible in, to help with alignment of low concentrations, but
this refers specifically to the channel where 100% of phiX clusters are visible.

`--perfect-target-name` the key used in the dictionary in the target YAML file that identifies your target sequence

`--alternate-perfect-reads` the path to a text file of read names that should be treated as the perfect target reads. Usually this is for experimenting with reads in case you're not sure what the protein might bind to

`--alternate-good-reads` just like alternate perfect reads above, except it is assumed that it contains some reads that will not bind as well

`--alternate-fiducial-reads` use read names in a given text file instead of phiX for the rough alignment step

`--microns-per-pixel` the size of a side of one pixel, in microns

`--chip` the type of chip, either `miseq` or `hiseq`

`--ports-on-right` the stage adapters we use orient the two input ports on the chip to the left or right, depending on
which microscope we use. This needs to be known for alignment. By default, we assume they are on the left unless this
flag is passed.

`--flipud` invert all images through the horizontal axis (i.e. run numpy.flipud() on all images). This was added to handle
a quirk with the way MicroManager saves images.

`--fliplr` invert all images through the vertical axis (i.e. run numpy.fliplr() on all images). If your images don't 
align you may try passing it in.

`-v -vv -vvv` set the verbosity level (-vvv is debug mode).

#### Generating HDF5 files

CHAMP uses HDF5 files with a specific format. If you're using MicroManager and generating OME-TIFF files, CHAMP's
built-in conversion tool (`champ h5`) will work out of the box. If your raw image files aren't formatted and named 
exactly as necessary, you'll need to generate the HDF5s yourself.

#### Aligning Images

CHAMP will attempt to align as many images as possible. The output will be the coordinates of each FASTQ read within 
an image, saved in text files in the `results` directory, along with a file containing the alignment parameters.

`IMAGE_DIRECTORY` the directory that contains all of the HDF5 image files

`--rotation-adjustment` rotational adjustment to apply to read coordinates before attempting alignment. Can be negative!
Even misalignment by a degree can prevent the rough alignment from working. If your alignments don't work, try a range 
of values from -5 to 5 degrees in 0.5 degree increments.

`--min-hits` the minimum number of exclusive hits required for a precision alignment to be considered valid

`--snr` the minimum signal-to-noise ratio (relative to random alignments) to consider a rough alignment valid. We have
found that 1.4 to be ideal under most scenarios.

`--make-pdfs` produce some diagnostic PDFs to examine the quality of the alignment

`--fiducial-only` only align the channel with the fiducial markers. 

`-v -vv -vvv` set the verbosity level (-vvv is debug mode).

#### Analyzing Results

Analyses of sequence specificity are performed using the Jupyter notebooks provided in the `notebooks` directory. The 
intended workflow is to copy them from this repo into each new experiment directory, edit the few variables as needed
at the top of each notebook, and run them.
