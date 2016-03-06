import fastq
from model.tile import load_tile_manager


def main(arguments):

    read_data = fastq.load_classified_reads(arguments.sorted_reads_directory,
                                            arguments.alignment_reads,
                                            arguments.random_selection,
                                            arguments.ignore_side_1)
    tile_manager = load_tile_manager()