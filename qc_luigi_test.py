import subprocess
import luigi
import os
import sys

FASTQC_SUMMERY = '/home/lxgui/scripts/fastqc_get_data_info_luigi.py'
sample_file = 'sample.list'

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output

class run_qc(luigi.Task):
    '''
    run fastqc
    '''
    sample = luigi.Parameter()
    CleanDir = luigi.Parameter()
    OutDir = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/{self.sample}.log'.format(**locals()))

    def run(self):
        tmp = run_cmd(['fastqc',
                        '{self.CleanDir}/{self.sample}_1_clean.fq.gz'.format(**locals()),
                        '{self.CleanDir}/{self.sample}_2_clean.fq.gz'.format(**locals()),
                        '-o',
                        '{self.OutDir}'.format(**locals())])
        with self.output().open('w') as qc_logs:
            qc_logs.write(tmp)     

class qc_report(luigi.Task) :
    '''
    generate qc report

    '''
    sample_list = [each.strip() for each in open(sample_file)]
    OutDir = luigi.Parameter()
    CleanDir = luigi.Parameter()

    def requires(self):
        return [run_qc(sample = sample, CleanDir = self.CleanDir, OutDir = self.OutDir) for sample in self.sample_list]

    def run(self):
        tmp = run_cmd(['python',
                        FASTQC_SUMMERY,
                        '{}'.format(sample_file),
                        '{self.OutDir}'.format(**locals())])
        with self.output().open('w') as qc_summary:
            qc_summary.write(tmp)

    def output(self):
        return luigi.LocalTarget('{self.OutDir}/qc.report.txt'.format(**locals()))

if __name__ == '__main__':
    luigi.run()