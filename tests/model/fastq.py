from chimp.model.fastq import FastqFiles
import unittest


class FastqFilesTests(unittest.TestCase):
    def setUp(self):
        filenames = {'WOO-005_S1_L001_R1_001.fastq.gz': 123,
                     'CC_S15_L001_I1_001.fastq.gz': 123,
                     'WOO-005_S1_L001_R2_001.fastq.gz': 123,
                     'UGH_R1_I1_BADTIMES': 123,
                     'UGH_R2_I1_BADTIMES': 123,
                     'temp_ideas.txt': 123,
                     '.gitignore': 123,
                     '.': 0,
                     '..': 0,
                     'Lulz_S16_L001_R1_001.fastq.gz': 123,
                     'Jim11135-3_S22_L001_R2_001.fastq.gz': 123,
                     'BOT-00_R1_0124124-WOO.fastq.gz': 123,
                     'BOT-00_R2_0124124-WOO.fastq.gz': 123}
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
