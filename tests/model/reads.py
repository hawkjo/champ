from model.reads import FastqFiles
import unittest


class FastqFilesTests(unittest.TestCase):
    def setUp(self):
        filenames = ['WOO-005_S1_L001_R1_001.fastq.gz',
                     'CC_S15_L001_I1_001.fastq.gz',
                     'WOO-005_S1_L001_R2_001.fastq.gz',
                     'temp_ideas.txt',
                     '.gitignore',
                     '.',
                     '..',
                     'Lulz_S16_L001_R1_001.fastq.gz',
                     'Jim11135-3_S22_L001_R2_001.fastq.gz',
                     'BOT-00_R1_0124124-WOO.fastq.gz',
                     'BOT-00_R2_0124124-WOO.fastq.gz']
        self.fq = FastqFiles(filenames)

    def test_single(self):
        expected = ['Lulz_S16_L001_R1_001.fastq.gz']
        actual = list(self.fq.single_files)
        self.assertListEqual(expected, actual)

    def test_paired(self):
        expected = [('BOT-00_R1_0124124-WOO.fastq.gz', 'BOT-00_R2_0124124-WOO.fastq.gz'),
                    ('WOO-005_S1_L001_R1_001.fastq.gz', 'WOO-005_S1_L001_R2_001.fastq.gz')]
        actual = list(sorted(self.fq.paired_files))
        self.assertEqual(expected, actual)
