from subprocess import check_output
import subprocess
import luigi
import os
import random
import sys

min_transcript_length = 200
fpkm_threshold = 0.5

EXTRACT_NOVEL_ISOFORM = '/home/lxgui/scripts/extract_novel_isoforms.py'

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    return output

class stringtie_merge(luigi.Task):
    GtfList = luigi.Parameter()
    RefGtf = luigi.Parameter()
    OutputDir = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('{self.OutputDir}/stringtie.merged.gtf'.format(**locals()))

    def run(self):
        tmp = run_cmd(["stringtie",
                        "--merge",
                        "-G",
                        self.RefGtf,
                        "-m",
                        str(min_transcript_length),
                        "-F",
                        str(fpkm_threshold),
                        self.GtfList])
        with self.output().open('w') as merged_gtf:
            merged_gtf.write(tmp)

class compare_with_ref(luigi.Task):

    GtfList = luigi.Parameter()
    RefGtf = luigi.Parameter()
    OutputDir = luigi.Parameter()

    def requires(self):
        return stringtie_merge(GtfList = self.GtfList, RefGtf = self.RefGtf, OutputDir = self.OutputDir)

    def output(self):
        return luigi.LocalTarget('{self.OutputDir}/cuffcompare.log'.format(**locals()))

    def run(self):
        tmp = run_cmd(['cuffcompare',
                        '-r',
                        self.RefGtf,
                        '-R',
                        '{self.OutputDir}/stringtie.merged.gtf'.format(**locals()),
                        '-o',
                        '{self.OutputDir}/ref'.format(**locals())])
        with self.output().open('w') as compare_log:
            compare_log.write('=====  cuffcompare with refference =======')

class extract_novel_isoforms(luigi.Task):

    GtfList = luigi.Parameter()
    RefGtf = luigi.Parameter()
    OutputDir = luigi.Parameter()

    def requires(self):
        return compare_with_ref(GtfList = self.GtfList, RefGtf = self.RefGtf, OutputDir = self.OutputDir)

    def run(self):
        tmp = run_cmd(['python',
                        EXTRACT_NOVEL_ISOFORM,
                        self.RefGtf,
                        '{self.OutputDir}/ref.combined.gtf'.format(**locals())])
        with self.output().open('w') as novel_inf:
            novel_inf.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutputDir}/novel.transcripts.gtf'.format(**locals()))

class sort_gtf(luigi.Task):

    GtfList = luigi.Parameter()
    RefGtf = luigi.Parameter()
    OutputDir = luigi.Parameter()

    def requires(self):
        return extract_novel_isoforms(GtfList = self.GtfList, RefGtf = self.RefGtf, OutputDir = self.OutputDir)

    def run(self):
        tmp = run_cmd(['sort',
                        '-k1,1',
                        '-k4,4g',
                        '{self.OutputDir}/novel.transcripts.gtf'.format(**locals())])
        with self.output().open('w') as sorted_gtf:
            sorted_gtf.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutputDir}/sorted.novel.transcripts.gtf'.format(**locals()))

if __name__=='__main__':
    luigi.run()
