# JIM NOTES: There are non-existent imports which prevent this script from being run.
# I'm guessing it's not used anymore and can be deleted.


# import sys
# import os
# import glob
# import local_config
# import pysam
# from subprocess import check_call
# from mypprint import pprint
#
#
# def get_relevant_dirs(projectname):
#     project_dir = os.path.join(local_config.data_dir,
#                                'from_fourierseq',
#                                projectname)
#     fq_dir = os.path.join(project_dir, 'all_fastqs')
#     phiX_sam_dir = os.path.join(project_dir, 'phiX_mappings')
#
#     for d in [project_dir, fq_dir]:
#         assert os.path.isdir(d), 'Directory not found: %s' % d
#
#     if os.path.isdir(phiX_sam_dir):
#         ans = raw_input('%s already exists. Continue? (y/[n])' % phiX_sam_dir)
#         if ans.lower() != 'y':
#             sys.exit('Aborted.')
#     else:
#         os.makedirs(phiX_sam_dir)
#
#     return project_dir, fq_dir, phiX_sam_dir
#
#
# def process_fnames(fq_dir):
#     fnames = os.listdir(fq_dir)
#     paired_fpaths, single_fpaths = set(), set()
#     for fname in fnames:
#         if 'fastq' not in fname:
#             continue
#         if '_I1_' in fname or '_I2_' in fname or '_I1.' in fname or '_I2.' in fname:
#             continue
#
#         fpath = os.path.join(fq_dir, fname)
#         if '_R1_' in fname or'_R1.' in fname:
#             fname2 = fname.replace('_R1_', '_R2_').replace('_R1.', '_R2.')
#             fpath2 = os.path.join(fq_dir, fname2)
#             if os.path.isfile(fpath2):
#                 paired_fpaths.add((fpath, fpath2))
#             else:
#                 single_fpaths.add(fpath)
#         elif '_R2_' in fname or '_R2.' in fname:
#             fname1 = fname.replace('_R2_', '_R1_').replace('_R2.', '_R1.')
#             fpath1 = os.path.join(fq_dir, fname1)
#             if os.path.isfile(fpath1):
#                 paired_fpaths.add((fpath1, fpath))
#             else:
#                 sys.exit('R2 file without R1: %s' % fname)
#         else:
#             single_fpaths.add(fpath)
#
#     paired_fpaths = list(sorted(paired_fpaths))
#     single_fpaths = list(sorted(single_fpaths))
#     return paired_fpaths, single_fpaths
#
#
# bowtie_cmd_starter = [
#                       'bowtie2',
#                       '--local',
#                       '-p 15',
#                       '--no-unal',
#                       '-x /home/hawkjo/genomes/phix/bowtie2/phiX_double',
#                      ]
#
#
# def bowtie_and_samtools(bowtie_cmd, sam_fpath):
#     assert sam_fpath.endswith('.sam')
#     base_fpath = sam_fpath[:-4]
#     bam_fpath = base_fpath + '.bam'
#     fpaths = {'base': base_fpath, 'sam': sam_fpath, 'bam': bam_fpath}
#
#     print 'Bowtie command:'
#     pprint(bowtie_cmd)
#     check_call(' '.join(bowtie_cmd), shell=True)
#
#     bam_cmd = 'samtools view -bS %(sam)s | samtools sort - %(base)s' % fpaths
#     print bam_cmd
#     check_call(bam_cmd, shell=True)
#
#     idx_cmd = 'samtools index %(bam)s' % fpaths
#     print idx_cmd
#     check_call(idx_cmd, shell=True)
#
#     os.unlink(sam_fpath)
#
#
# def get_output_fpaths(f, phiX_sam_dir):
#     sam_fname = os.path.split(f)[-1].split('.')[0].replace('_R1', '') + '.sam'
#     error_fname = sam_fname[:-4] + '_error.txt'
#     return [os.path.join(phiX_sam_dir, f) for f in [sam_fname, error_fname]]
#
#
# def paired_call(f1, f2, phiX_sam_dir):
#     assert '_R1' in os.path.split(f1)[-1], f1
#     sam_fpath, error_fpath = get_output_fpaths(f1, phiX_sam_dir)
#     bowtie_cmd = bowtie_cmd_starter + [
#                                        '-1 ' + f1,
#                                        '-2 ' + f2,
#                                        '-S ' + sam_fpath,
#                                        '2>&1 | tee ' + error_fpath,
#                                       ]
#     bowtie_and_samtools(bowtie_cmd, sam_fpath)
#
#
# def single_call(f, phiX_sam_dir):
#     sam_fpath, error_fpath = get_output_fpaths(f, phiX_sam_dir)
#     bowtie_cmd = bowtie_cmd_starter + [
#                                        '-U ' + f,
#                                        '-S ' + sam_fpath,
#                                        '2>&1 | tee ' + error_fpath,
#                                       ]
#     bowtie_and_samtools(bowtie_cmd, sam_fpath)
#
#
# def map_project(projectname):
#     project_dir, fq_dir, phiX_sam_dir = get_relevant_dirs(projectname)
#
#     paired_fpaths, single_fpaths = process_fnames(fq_dir)
#     paired = bool(paired_fpaths)
#     single = bool(single_fpaths)
#     assert paired != single, 'Not (single xor paired) reads in %s' % fq_dir
#
#     if paired:
#         for f1, f2 in paired_fpaths:
#             paired_call(f1, f2, phiX_sam_dir)
#     else:
#         for f in single_fpaths:
#             single_call(f, phiX_sam_dir)
#
#     bam_fpaths = glob.glob(os.path.join(phiX_sam_dir, '*.bam'))
#     read_names = set(read.qname for fpath in bam_fpaths for read in pysam.Samfile(fpath))
#     read_name_fpath = os.path.join(phiX_sam_dir, 'phiX_read_names.txt')
#     with open(read_name_fpath, 'w') as out:
#         out.write('\n'.join(read_names))
#
#
# if __name__ == '__main__':
#     if len(sys.argv) != 2:
#         sys.exit('Usage: %s <ProjectName>' % sys.argv[0])
#
#     map_project(sys.argv[1])
