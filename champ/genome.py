from champ.readmap import FastqFiles
from champ.error import fail
import itertools
import os
import subprocess


def build_genomic_bamfile(fastq_directory, bowtie_directory_and_prefix='.local/champ/human-genome'):
    """ Aligns all gzipped FASTQ files in a given directory to a reference genome,
    creating a Bamfile that can be read by pysam, among other tools. This must be performed before
    other CHAMP genomic analyses can be run. """
    fastq_filenames = [os.path.join(fastq_directory, directory) for directory in os.listdir(fastq_directory)]
    fastq_files = FastqFiles(fastq_filenames)
    forward_paired = []
    reverse_paired = []
    unpaired = []
    # Trim Illumina adapters from all the sequences and store them temporarily in fastq.gz files
    for n, (forward_fastq, reverse_fastq) in enumerate(fastq_files.paired):
        forward_paired_file = 'genomic_1_trimmed_paired%d.fastq.gz' % n
        forward_paired.append(forward_paired_file)

        reverse_paired_file = 'genomic_2_trimmed_paired%d.fastq.gz' % n
        reverse_paired.append(reverse_paired_file)

        forward_unpaired_file = 'genomic_1_trimmed_unpaired%d.fastq.gz' % n
        reverse_unpaired_file = 'genomic_2_trimmed_unpaired%d.fastq.gz' % n
        unpaired.append(forward_unpaired_file)
        unpaired.append(reverse_unpaired_file)

        # Currently we hardcode the Trimmomatic binary and the adapters file.
        # These will either be bundled with subsequent versions, found dynamically, or given by the user
        result = subprocess.check_call(['/usr/bin/TrimmomaticPE', '-threads', '20', '-phred64',
                                        forward_fastq, reverse_fastq,
                                        forward_paired_file,
                                        forward_unpaired_file,
                                        reverse_paired_file,
                                        reverse_unpaired_file,
                                        'ILLUMINACLIP:/shared/trim/TruSeq-All.fa:2:30:10',
                                        'LEADING:3',
                                        'TRAILING:3',
                                        'MINLEN:25'])
        if result != 0:
            fail("Could not trim adapters from FASTQ files")

    forward_paired_filenames = ','.join(forward_paired)
    reverse_paired_filenames = ','.join(reverse_paired)
    unpaired_filenames = ','.join(unpaired)

    # Align the trimmed sequences with the reference genome. Again, we've hardcoded the path to the
    # Bowtie files. These are installed with the prepare-genomic-files.sh script, which needs to be bundled
    # with CHAMP or done in Python, or something.
    home = os.path.expanduser("~")
    result = subprocess.check_call(['/usr/bin/bowtie2', '-p', '16',
                                    '--no-unal', '-x', os.path.join(home, bowtie_directory_and_prefix),
                                    '-1', forward_paired_filenames,
                                    '-2', reverse_paired_filenames,
                                    '-U', unpaired_filenames,
                                    '-S', 'genomic.sam'])

    if result != 0:
        fail("Could not build samfile.")

    # Convert the Samfile to a Bamfile
    try:
        samtools_sort = subprocess.Popen(['/usr/bin/samtools', 'sort',
                                          '-', os.path.join(fastq_directory, 'final')],
                                         stdin=subprocess.PIPE)
        samtools_view = subprocess.Popen(['/usr/bin/samtools', 'view',
                                          '-bS', 'genomic.sam'],
                                         stdout=samtools_sort.stdin)
        samtools_sort.communicate()
        samtools_view.wait()
    except:
        fail("Problem with samtools.")

    # Delete the temporary files we created
    for filename in itertools.chain(forward_paired, reverse_paired, unpaired, ('genomic.sam',)):
        try:
            os.unlink(filename)
        except:
            continue

    print("Done aligning FASTQ reads to reference genome!")
