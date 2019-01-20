from champ.readmap import FastqFiles
from champ.error import fail
import distutils
import subprocess
import os
from Bio import SeqIO
from pysam import Samfile
import progressbar
import itertools
import numpy as np
import h5py

SUCCESS = 0
QUALITY_THRESHOLD = 20
MAXIMUM_REALISTIC_DNA_LENGTH = 1000


def build_genomic_bamfile(fastq_directory, output_filename, bowtie_directory_and_prefix='.local/champ/human-genome', reference_genome_fastq_filename=None):
    """ Aligns all gzipped FASTQ files in a given directory to a reference genome,
    creating a Bamfile that can be read by pysam, among other tools. This must be performed before
    other CHAMP genomic analyses can be run. """
    if reference_genome_fastq_filename is None:
        reference_genome_fastq_filename = os.path.join(os.path.expanduser("~"), '.local', 'champ', 'human-genome.fna')

    if not output_filename:
        output_filename = 'genomic'

    samtools = distutils.spawn.find_executable("samtools")
    if samtools is None:
        fail("samtools is not installed.")
    trimmomatic = distutils.spawn.find_executable("TrimmomaticPE")
    if trimmomatic is None:
        fail("trimmomatic is not installed.")
    bowtie2 = distutils.spawn.find_executable("bowtie2")
    if bowtie2 is None:
        fail("bowtie2 is not installed.")

    fastq_filenames = [os.path.join(fastq_directory, filename) for filename in os.listdir(fastq_directory)]
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
        result = subprocess.check_call([trimmomatic, '-threads', '20', '-phred64',
                                        forward_fastq, reverse_fastq,
                                        forward_paired_file,
                                        forward_unpaired_file,
                                        reverse_paired_file,
                                        reverse_unpaired_file,
                                        'ILLUMINACLIP:/shared/trim/TruSeq3-PE.fa:2:30:10',
                                        'LEADING:3',
                                        'TRAILING:3',
                                        'MINLEN:25'])
        if result != SUCCESS:
            fail("Could not trim adapters from FASTQ files")

    forward_paired_filenames = ','.join(forward_paired)
    reverse_paired_filenames = ','.join(reverse_paired)
    unpaired_filenames = ','.join(unpaired)

    # Align the trimmed sequences with the reference genome. Again, we've hardcoded the path to the
    # Bowtie files. These are installed with the prepare-genomic-files.sh script, which needs to be bundled
    # with CHAMP or done in Python, or something.
    home = os.path.expanduser("~")
    result = subprocess.check_call([bowtie2, '-p', '16',
                                    '--no-unal', '-x', os.path.join(home, bowtie_directory_and_prefix),
                                    '-1', forward_paired_filenames,
                                    '-2', reverse_paired_filenames,
                                    '-U', unpaired_filenames,
                                    '-S', os.path.join(fastq_directory, '{output_filename}.sam'.format(output_filename=output_filename))])

    if result != SUCCESS:
        fail("Could not build samfile.")

    # Convert the Samfile to a Bamfile
    try:
        samtools_sort = subprocess.Popen([samtools, 'sort',
                                          '-', os.path.join(fastq_directory, '{output_filename}'.format(output_filename=output_filename))],
                                         stdin=subprocess.PIPE)
        samtools_view = subprocess.Popen([samtools, 'view',
                                          '-bS', os.path.join(fastq_directory, '{output_filename}.sam'.format(output_filename=output_filename))],
                                         stdout=samtools_sort.stdin)
        samtools_sort.communicate()
        samtools_view.wait()
        bamfile_path = os.path.join(fastq_directory, '{output_filename}.bam'.format(output_filename=output_filename))
        result = subprocess.check_call([samtools, 'index', bamfile_path])
        if result != SUCCESS:
            fail("Could not index bamfile.")
        else:
            print("Building quality read name sequences.")
            read_name_sequences = get_quality_paired_end_read_sequences(bamfile_path, reference_genome_fastq_filename)
            save_quality_read_name_sequences(read_name_sequences, os.path.join(fastq_directory, 'quality-read-name-sequences.h5'))
            print("Created quality read name sequences.")
    except Exception as e:
        print(e)
        fail("Problem with samtools.")

    # Delete the temporary files we created
    for filename in itertools.chain(forward_paired,
                                    reverse_paired,
                                    unpaired,
                                    (os.path.join(fastq_directory,
                                                  '{output_filename}.sam'.format(output_filename=output_filename)),)):
        try:
            os.unlink(filename)
        except:
            continue

    print("Done aligning FASTQ reads to reference genome!")


def get_quality_alignment_start_and_end(alignment):
    if alignment.is_qcfail or alignment.mapq < 20:
        return None
    if alignment.is_proper_pair:
        # This read name appears twice in our data - once in the forward and once in the reverse direction
        if alignment.is_reverse:
            # For paired end reads, the forward and reverse reads are symmetric, so to avoid double counting the
            # bases between the two reads, we only count the bases once (with the forward read)
            return None
        if alignment.template_length <= alignment.reference_length or alignment.reference_length == 0:
            # the combined length of the first and second read is no greater than just one read - either the reads
            # perfectly overlapped or something weird is happening. We discard this read to be safe.
            return None
        start = alignment.reference_start
        end = start + alignment.template_length
    # elif alignment.reference_length == alignment.query_length and alignment.reference_length > 0:
    elif alignment.reference_length > 0:
        # This is an unpaired read, which is rare but not unheard of. We require that the read align perfectly to the
        # reference sequence, so we aren't dealing with indels. This might be a bit conservative - we'll come back to
        # this decision if the read counts are terrible
        start = alignment.reference_start
        end = start + alignment.reference_length
        assert start < end
    elif alignment.reference_length == 0:
        return None
    elif alignment.reference_length < 0:
        return None
    else:
        # The read is unpaired and the alignment was sketchy, so we can't trust this read
        return None
    if abs(end - start) > MAXIMUM_REALISTIC_DNA_LENGTH:
        return None
    return start, end


def get_quality_paired_end_read_sequences(bamfile, fasta_filename=None):
    """
    Associate genomic sequences with read names when those reads are paired and they pass all quality checks.
    We do this separately from the process where gene sequences are stored so that we can look at effects that might
    arise on physical clusters of DNA.

    """
    if fasta_filename is None:
        fasta_filename = os.path.join(os.path.expanduser("~"), '.local', 'champ', 'human-genome.fna')

    read_name_positions = {}
    try:
        with Samfile(bamfile) as samfile:
            contigs = list(reversed(sorted(samfile.references)))
            with progressbar.ProgressBar(max_value=len(contigs)) as pbar:
                for contig in pbar(contigs):
                    read_name_positions[contig] = {}
                    for alignment in samfile.fetch(contig):
                        start_end = get_quality_alignment_start_and_end(alignment)
                        if start_end is None:
                            continue
                        start, end = start_end
                        read_name_positions[contig][alignment.query_name] = (start, end)
        read_name_sequences = []
        with open(fasta_filename) as f:
            for record in SeqIO.parse(f, 'fasta'):
                contig = record.id
                for read_name, (start, end) in read_name_positions[contig].items():
                    read_name_sequences.append((read_name, record.seq[start:end]))
        return read_name_sequences
    except IOError:
        raise ValueError("Could not open either %s or %s. Does it exist and is it valid?" % (bamfile, fasta_filename))


def save_quality_read_name_sequences(read_name_sequences, hdf5_filename):
    string_dt = h5py.special_dtype(vlen=str)
    read_name_sequences_dt = np.dtype([('name', string_dt),
                                       ('sequence', string_dt)])
    with h5py.File(hdf5_filename, 'w') as h5:
        dataset = h5.create_dataset('/sequences', (len(read_name_sequences),), dtype=read_name_sequences_dt)
        dataset[...] = read_name_sequences
