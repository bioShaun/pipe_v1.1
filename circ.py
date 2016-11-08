#!/usr/bin/env python
# -*- coding: utf-8 -*- #

# circ.py --- circRNA analysis library

# Author: kangmingming AT novogene DOT com
# Time-stamp: <2016-01-21 18:44:51 kangmingming>

# TODO:


'''circRNA analysis library.'''

import sys, os, re, argparse, glob
sys.dont_write_bytecode = True
from os.path import join as jp
from os.path import basename
from os.path import dirname
from time import localtime, strftime
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
#import misopy.cluster_utils as cluster_utils

LIB_CIRC = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, LIB_CIRC)
#import circ_var_defs as var


CIRC_PY_HEADER = '''\
#!/usr/bin/env python
# -*- coding: utf-8 -*- #

import os, sys, re
import misopy.cluster_utils as cluster_utils
sys.dont_write_bytecode = True

LIB_CIRC = '%s'
sys.path.insert(0, LIB_CIRC)
import circ
import circ_var_defs as var

''' % LIB_CIRC


def parse_fasta(f):
    retval = {}
    for record in SeqIO.parse(f, 'fasta'):
        retval[record.id] = str(record.seq)
    return retval


def write_str_to_file(string, fn, append=False):
    dirname = os.path.dirname(fn)
    if not os.path.isdir (dirname):
        circ_mkdir_unix (dirname)
    fh = open(fn, 'a' if append is True else 'w')
    fh.write(string)
    fh.close()


def write_obj_to_file(obj, fn, append=False):
    fh = open(fn, 'a' if append is True else 'w')
    if type(obj) is str:
        fh.write(obj)
    elif type(obj) is list:
        for item in obj:
            fh.write('%s\n' % item)
    elif type(obj) is dict:
        for key, val in obj.iteritems():
            fh.write('%s\t%s\n' % (key, val))
    else:
        raise TypeError('invalid type for %s' % obj)
    fh.close()


def circ_get_time():
    return strftime("[%Y-%m-%d %H:%M:%S]", localtime())


def rel2abs(path):
    return os.path.abspath(path)


def iswritable(path):
    return os.access(path, os.W_OK)


def log(msg, fn=None, add_blank_line_before=False, add_blank_line_after=False):
    msg = '%s [lncRNA] %s\n' % (circ_get_time(), msg)
    if add_blank_line_before is True:
        msg = '\n%s' % msg
    if add_blank_line_after is True:
        msg = '%s\n' % msg
    # print msg
    print >> sys.stderr, msg
    if fn:
        write_str_to_file(msg, fn, append=True)


def log_call_proc(msg, fn=None):
    msg = '%s [lncRNA] [calling process] %s\n' % (circ_get_time(), msg)
    # print msg
    print >> sys.stderr, msg
    if fn:
        write_str_to_file(msg, fn, append=True)


def error(msg, line=None):
    print >> sys.stderr, msg
    sys.exit(1)


def circ_call_process(cmd, allow_error=False, msg=None):
    log_call_proc(cmd)
    exit_status = os.system(cmd)
    if not allow_error and exit_status != 0:
        error_msg = cmd if msg is None else msg
        error('Error when running process: %s' % error_msg)


def circ_run_script(script, file_type, allow_error=False, msg=None):
    if file_type == 'sh':
        cmd = 'sh %s' % script
    elif file_type == 'R':
        cmd = '%s --vanilla %s' % (var.prog_Rscript, script)
    elif file_type == 'python':
        cmd = 'python %s' % script
    else:
        raise ValueError('invalid file_type: %s' % file_type)
    # run
    if os.path.isfile(script):
        circ_call_process(cmd, allow_error, msg)
    else:
        raise ValueError('no such file: %s' % script)


def circ_mkdir_unix(path):
    cmd = 'mkdir -p %s' % path
    if not os.path.isdir(path):
        circ_call_process(cmd)

def circ_mkdir(path):
    # cmd = 'mkdir -p %s' % path
    if not os.path.isdir(path):
        os.mkdir(path)
        # circ_call_process(cmd)


def circ_safe_ln(path1, path2):
    if os.path.exists(path2):
        time_stamp = strftime("%Y-%m-%d_%H:%M:%S", localtime())
        path2_bak = path2 + '.bak.%s' % time_stamp
        circ_call_process('mv -f %s %s' % (path2, path2_bak))
    cmd = 'ln -fs %s %s' % (path1, path2)
    circ_call_process(cmd)


def quote_str(obj):
    if type(obj) is str:
        obj = "'%s'" % obj
    else:
        obj
    return obj


def q(obj):
    return quote_str(obj)


def get_queue ():
    queue = os.popen("qselect -U $USER | grep -Ev 'mem|joyce' | cut -d@ -f1 | sort -u | tr '\\n' ' ' | sed -r 's/(\\S+)/-q \\1/g'").readlines()[0].strip()
    return queue


def qsub_cmd (mem, p, fn, error_dir=None, interpreter=None):
    """
    Keyword Arguments:
    mem --
    p   --
    fn  --
    """
    if error_dir is None:
        error_dir = os.path.dirname(fn)
    elif not os.path.isdir(error_dir):
        error('no such directory: %s' % error_dir)
    # determine which interpreter to use
    interpreter_switch = var.interpreter_switch
    # fn_suffix = re.sub ('^.+\.(.+)$',  '\g<1>',  basename (fn))
    fn_suffix = basename (fn).split('.')[-1]
    if interpreter_switch.has_key (fn_suffix):
        interpreter_auto = interpreter_switch[fn_suffix]
    else:
        interpreter_auto = 'bash'
    interpreter_auto = '$(which %s)' % interpreter_auto
    if interpreter is None:
        interpreter = interpreter_auto
    else:
        interpreter = '$(which %s)' % interpreter
    queue = get_queue ()
    return 'qsub -S {interp} {queue} -V -cwd -l vf={vf},p={p} -e {error_dir} -o {out_dir} {fn}'.format (
        interp=interpreter,
        queue=queue,
        vf=mem,
        p=p,
        error_dir=error_dir,
        out_dir=error_dir,
        fn=fn,
        )
    # else:
    #     return 'qsub -S /bin/bash -q all.q -q rna.q -V -cwd -l vf=%s,p=%s -e %s -o %s -S %s %s' % (mem, p, error_dir, error_dir, interpreter, fn)

# def launch_job(mem, p, fn, error_dir=None):
#     cmd = qsub_cmd(mem, p, fn, error_dir)
#     log_call_proc(cmd)
#     job_id = cluster_utils.launch_job(cmd, cmd_name='qsub')
#     return job_id

def launch_jobs(mem, p, fns, error_dir=None, interpreter=None, wait=False):
    if fns is None:
        log('(no jobs to run)')
        return
    job_ids = []
    if type(fns) is str:
        fns = [fns]
    for fn in fns:
        if not os.path.isfile(fn):
            raise ValueError('no such file: %s' % fn)
        cmd = qsub_cmd(mem, p, fn, error_dir, interpreter)
        log_call_proc(cmd)
        job_id = cluster_utils.launch_job(cmd, cmd_name='qsub')
        job_ids.append(job_id)
    if wait:
        wait_jobs(job_ids)
    else:
        return job_ids


def wait_jobs(job_ids):
    """
    Keyword Arguments:
    job_ids --
    """
    if job_ids is None:
        return None
    elif type(job_ids) is int:
        log('waiting for job: %s' % job_ids)
        cluster_utils.wait_on_job(job_ids, cluster_cmd='qsub', delay=60)
    elif type(job_ids) is list:
        for _id in job_ids:
            if type(_id) is not int:
                raise TypeError('invalid job ID: %s' % _id)
        log('waiting for jobs: %s' % job_ids)
        cluster_utils.wait_on_jobs(job_ids, cluster_cmd='qsub', delay=10)
    else:
        raise TypeError('invalid job IDs: %s' % job_ids)


def parse_compare_venn(Vsamples, Vgroupname, Vcompare, Vvenn):
    samples = re.split(':|,', Vsamples)
    for sample in samples:
        if re.match('[0-9]+', sample):
            raise ValueError('Invalid sample name: %s\nSample name should start with [A-Za-z]' % sample)
    sample_grp = Vsamples.split(',')
    replicate = True
    for samples_within_grp in sample_grp:
        if len(samples_within_grp.split(':')) == 1:
            replicate = False
            break
    groupname = Vgroupname.split(',')
    # if len(sample_grp) != len(groupname):
    #     raise ValueError('samples group and groupname are not equal in length')
    Vcompare = Vcompare.strip (' ')
    if Vvenn is not None:
        Vvenn = Vvenn.strip (' ')
    compare_list = Vcompare.split(',')
    compare = {}
    for i in range(0, len(compare_list)):
        compare_name = compare_list[i]
        compare_each_grp_idx = [int(j)-1 for j in compare_name.split(':')]
        compare_each_grpname = [groupname[j] for j in compare_each_grp_idx]
        samples_within_grp = [sample_grp[j].split(':') for j in compare_each_grp_idx]
        compare[compare_name] = [compare_each_grpname, samples_within_grp]
    if Vvenn:
        venn_list = Vvenn.split(',')
        venn = {}
        for venn_each in venn_list:
            diff_grp = venn_each.split(':')
            for compare_grp in diff_grp:
                grp = [groupname[int(i)-1] for i in compare_grp.split('_')]
                venn.setdefault(venn_each, []).append(grp)
    else:
        venn = None
    return samples, groupname, compare, venn, replicate


def gtf_to_genepred(novel_dir, gtf):
    """
    Keyword Arguments:
    novel_dir --
    gtf       --
    """
    circ_mkdir(novel_dir)
    genepred_new = var.circ_annotation_genepred(novel_dir, gtf)
    genepred_old = '%s.old' % genepred_new
    cmd_gtfToGenePred = '%s -genePredExt -ignoreGroupsWithoutExons %s %s' % (var.prog_gtfToGenePred, gtf, genepred_old)
    circ_call_process(cmd_gtfToGenePred)
    cmd_genepred_add_col1 = """cat %s | perl -fwane 'print join ("\\t", @F[11, 0 .. 9]), "\\n"' >%s""" % (genepred_old, genepred_new)
    # cmd_genepred_add_col1 = """awk '{print $1"\t"$0}' %s >%s""" % (genepred_old, genepred_new)
    circ_call_process(cmd_genepred_add_col1)
    return genepred_new


def run_by_run_type(run_type, cmd, fn=None):
    """Run external program, or write cmd to file and return the filename, depends on `run_type'.
    Keyword Arguments:
    cmd      -- external program and args to run.
    fn       -- filename that the `cmd' will be write to.
    run_type -- the external program will be run (run_type="call"), or the program (and args) will be write to file `fn' and `fn' is returned (run_type="file").
    """
    if run_type == 'call':
        circ_call_process(cmd)
    elif run_type == 'file':
        if type(fn) is not str:
            raise ValueError('fn should be a valid filename')
        write_str_to_file(cmd, fn)
        return fn
    else:
        raise ValueError('invalid run_type: %s' % run_type)


def as_list(obj):
    """ Return None or a list
    Keyword Arguments:
    obj -- None, str, or list
    """
    if type(obj) is str:
        return [obj]
    elif type(obj) is list:
        return obj
    else:
        raise ValueError('invalid object type: %s' % obj)


def cp (f1, f2, r=False, pdf=False):
    if os.path.exists(f1):
        if r:
            circ_call_process ('cp -r %s %s' % (f1, f2))
        else:
            circ_call_process ('cp %s %s' % (f1, f2))
            if pdf:
                f1pdf = re.sub('png$', 'pdf', f1)
                f2pdf = re.sub('png$', 'pdf', f2)
                circ_call_process ('cp %s %s' % (f1pdf, f2pdf))
    else:
        log('no such file or directory: %s' % f1)


def glob_cp (regex, dest, r=False, pdf=False):
    fns = glob.glob(regex)
    if len(fns) > 0:
        cp (fns[0], dest, r, pdf)
    else:
        log('no file or directory matching "%s"' % regex)


class circ_argv_parser:
    def __init__(self):
        self.argv = dict()
        self.argv_copy = dict()
        self.argv_parser()

    def argv_parser(self):
        parser = argparse.ArgumentParser(description="circRNA analysis pipline")
        parser.add_argument('--project', help='Project name, which will be displayed in the report title', required=True)
        parser.add_argument('--project-dir', help='Project dir, default is "os.getcwd()"', required=False)
        parser.add_argument('--sample', help="Sample names, separated by comma", default=None, required=True)
        parser.add_argument('--group', help="Group names, separated by comma", default=None, required=True)
        parser.add_argument('--groupname', help="Group names, separated by comma", default=None, required=True)
        parser.add_argument('--compare', help="Diff-expression list, separated by comma (e.g. 1:2,1:3,...)", default=None, required=True)
        parser.add_argument('--venn', help="Venn diagram list, separated by comma. defult is all groups, 1:2_1:3,1:3_2:3,... ", default=None)
        parser.add_argument('--genome', help="Genome sequence in fasta format", required=True)
        parser.add_argument('--gtf', help="GTF annotation corresponds to the given genome", required=True)
        parser.add_argument('--go', help="GO annotation file", required=True)
        parser.add_argument('--kegg', help="KEGG species abbreviation (e.g. 'hsa' for Human)", required=True)
        parser.add_argument('--report-language', help='Chinese (cn) or English (en) for the report ("cn")', choices=['cn', 'en'], default='cn', required=False)
        parser.add_argument('--abbr', help="Three-letter abbreviation of the species, used for naming novel circRNAs", required=True)
        parser.add_argument('--program', help="Program for identifying circRNAs", choices=['find_circ', 'ciri', 'both'], required=True)
        parser.add_argument('--circ-db', help="circRNA annotation file, currently only the circBase annotation is supported", required=False)

        # qc args
        arg_fq_source = parser.add_mutually_exclusive_group()
        arg_fq_source.add_argument('--fq', help="Fastq dir", required=False)
        arg_fq_source.add_argument('--raw-dir', help="Raw data dir", required=False)
        arg_fq_source.add_argument('--lnc-qc-dir', help="lncRNA QC dir, used only if the lncRNA QC dir is not 'byebyed'", default=None, required=False)

        parser.add_argument('--library-type', help="Library type", choices=['yes', 'reverse', 'no'], required=False)
        parser.add_argument('--mapfile', help="Mapfile", required=True)

        parser.add_argument('--linear-free', help="Is linear RNAs removed from the library?", choices=['yes', 'no'], required=True)
        parser.add_argument('--mirna', help="Three-letter miRNA abbreviation, or miRNA sequence file in fasta format", default=False, required=True)
        parser.add_argument('--kingdom', help="Either 'animal' or 'plant'", choices=['animal', 'plant'], required=True)

        # optional args that skip QC
        arg_group2 = parser.add_mutually_exclusive_group()
        arg_group2.add_argument('--tophat-dir', help="Tophat output dirs (in the same order of argv[sample]), separated by comma (',')", required=False)
        arg_group2.add_argument('--fq-clean', help="The prefixes of clean fq files (in the same order of argv[sample], please omit '_[12].clean.fq.gz'), separated by comma (',')", required=False)
        # parser.add_argument('--skip-qc', help="if we should skip the QC step.", required=False)

        # # optional args that skip the given steps
        # parser.add_argument('--skip-step', help="steps that are skipped to run, separated by comma (',')", default=None)

        # # for debugging
        # parser.add_argument('--debug', help="print more debug info", action='store_true', default=False)
        parser.add_argument('--test', help="if present, all other args will not be used", action='store_true', default=False)

        # argv = parser.parse_args()
        argv = vars(parser.parse_args())
        if not argv['test']:
            # self.project = argv['project']
            if argv['project_dir']:
                if os.path.isdir(argv['project_dir']):
                    # self.project_dir = argv['project_dir']
                    pass
                else:
                    raise ValueError('invalid project dir: %s' % argv['project_dir'])
            else:
                # self.project_dir = os.getcwd()
                argv['project_dir'] = os.getcwd()
            if not iswritable(argv['project_dir']):
                error('project directory not writable: %s' % argv['project_dir'])
            if argv['library_type']:
                # self.library_type = argv['library_type']
                pass
            else:
                # self.library_type = 'ss'
                argv['library_type'] = 'reverse'

            if argv['fq'] and argv['lnc_qc_dir']:
                raise ValueError ('Please specify either --fq or --lnc-qc-dir')
            # determine the start point of analysis
            if argv['fq_clean']:
                fq_cleans = re.split(':|,', argv['fq_clean'])
                for fq_prefix in fq_cleans:
                    fq1 = fq_prefix + '_1.clean.fq.gz'
                    fq2 = fq_prefix + '_2.clean.fq.gz'
                    if not os.path.isfile(fq1):
                        error('no such file: %s' % fq1)
                    elif not os.path.isfile(fq2):
                        error('no such file: %s' % fq2)
            elif argv['tophat_dir']:
                tophat_dirs = re.split(':|,', argv['tophat_dir'])
                for tophat_dir_old in tophat_dirs:
                    unmapped_bam = os.path.join(tophat_dir_old, 'unmapped.bam')
                    if not os.path.isfile(unmapped_bam):
                        error('no such file: %s' % unmapped_bam)
            elif argv['lnc_qc_dir']:      # fq and other qc args are not necessary
                if not os.path.isdir (argv['lnc_qc_dir']):
                    raise IOError ('no such directory: %s' % argv['lnc_qc_dir'])
            else:
                argv['fq'], argv['fq_raw_list'] = self.check_argv_fq(argv)
            # self.genome = argv['genome']
            # self.mirna = argv['mirna']
            # self.kingdom = argv['kingdom']
            # self.sample = argv['sample']
            # if argv.has_key ('skip_step') and argv['skip_step']:
            #     for skipped in argv['skip_step'].split(','):
            #         if skipped not in var.CIRC_SKIP_STEP_LIST:
            #             raise ValueError('no such step to skip: %s\navailable steps are: %s' % (skipped, var.CIRC_SKIP_STEP_LIST))
            # if not argv['fq_clean'] and not argv['tophat_dir']:
            #     if argv['skip_step']:
            #         if 'tophat' not in argv['skip_step'].split(','):
            #             argv['skip_step'] += ',tophat'
            #     else:
            #         argv['skip_step'] = 'tophat'

            if not argv['mapfile']:
                raise ValueError ('mapfile is required for QC')
            else:
                fh = open (argv['mapfile'], 'r')
                for line in fh:
                    line = line.strip ()
                    if re.match ('\s*(#.+)*$', line):
                        continue
                    elif re.match ('\S+', line):
                        if len (re.split ('\s+', line)) != 2:
                            raise ValueError ('Invalid mapfile in line:\n%s' % line)
                fh.close ()

            # FIXME: why the expression `argv['program'] is 'find_circ'` will fail ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if argv['program'] == 'find_circ':
                idx = '%s.1.bt2' % argv['genome']
                if not os.path.isfile (idx):
                    raise IOError ('Bowtie2 index not found for genome: %s' % argv['genome'])
            elif argv['program'] == 'ciri':
                idx = '%s.pac' % argv['genome']
                if not os.path.isfile (idx):
                    raise IOError ('BWA index not found for genome: %s\nUse "bwa index %s" to build' % (argv['genome'], argv['genome']))
            elif argv['program'] == 'both':
                idx = '%s.1.bt2' % argv['genome']
                if not os.path.isfile (idx):
                    raise IOError ('Bowtie2 index not found for genome: %s' % argv['genome'])
                idx = '%s.pac' % argv['genome']
                if not os.path.isfile (idx):
                    raise IOError ('BWA index not found for genome: %s\nUse "bwa index %s" to build' % (argv['genome'], argv['genome']))
            else:
                raise ValueError ('Invalid program: %s' % argv['program'])

            tmp = {}
            tmp['_sample'], tmp['_group'], tmp['_groupname'], tmp['_compare'], tmp['_venn'] \
              = argv['sample'], argv['group'], argv['groupname'], argv['compare'], argv['venn']
            self.argv_copy = tmp
            argv['samples'], argv['groupname'], argv['compare'], argv['venn'], argv['replicate'] = parse_compare_venn(argv['sample'], argv['groupname'], argv['compare'], argv['venn'])
        if argv['circ_db']:
            if len(argv['circ_db']) == 3:
                if var.circ_db_switch.has_key (argv['circ_db']):
                    argv['circ_db'] = var.circ_db_switch[argv['circ_db']]
                else:
                    raise ValueError ('Invalid circ-db: %s' % argv['circ_db'])
            else:
                if not os.path.isfile(argv['circ_db']):
                    raise IOError ('no such file: %s' % argv['circ_db'])
        self.argv = argv

    def check_argv_fq(self, argv):
        "stolen subr."
        # argv = vars(argv)
        samples = re.split(':|,|;', argv['sample'])
        PROJ_DIR = argv['project_dir']
        ln_raw_data = var.prog_ln_raw_data
        if argv['fq']:
            fq = os.path.abspath(argv['fq'].strip())
        else:
            # os.system('mkdir -p raw_data')
            circ_mkdir(jp(PROJ_DIR, 'raw_data'))
            circ_call_process ('perl %s %s %s pe raw_data' % (ln_raw_data,
                                                              argv['mapfile'],
                                                              argv['raw_dir']))
            fq = PROJ_DIR + '/raw_data'
        # self.fq = fq
        fq_raw_list = []
        sys.stderr.write('checking fq...')
        for each in samples:
            fq_1 = fq + '/' + each + '_1.fq.gz'
            fq_2 = fq + '/' + each + '_2.fq.gz'
            if not os.path.isfile(fq_1):
                raise IOError ('no such file: %s' % fq_1)
            if not os.path.isfile(fq_2):
                raise IOError ('no such file: %s' % fq_2)
            fq_raw_list.append(fq_1)
        sys.stderr.write(' OK\n')
        return fq, fq_raw_list


class circ_qc:
    '''QC.'''
    def __init__(self, argv):
        "circRNA QC."
        self.argv = argv

    def write_qc (self):
        if self.argv['lnc_qc_dir'] is not None:
            self.write_qc_report_only ()
        else:
            self.write_qc_full ()

    def write_qc_report_only (self):
        v = self.argv
        lnc_qc_dir = re.sub ('/*$', '/', v['lnc_qc_dir'])
        circ_qc_dir = var.circ_qc_dir (v['project_dir'])
        if os.path.exists (circ_qc_dir):
            circ_call_process ('mv -f %s %s' % (circ_qc_dir, circ_qc_dir + '.bak'))
        circ_call_process ('ln -fs %s %s' % (lnc_qc_dir, circ_qc_dir))
        sample = ','.join (re.split(':|,', v['sample']))
        project = v['project']
        qc_report_dir = jp (v['project_dir'], 'QC_report')
        pycode = '''
sh {nest_circ_qc_report} \\
   -dir {circ_qc_dir} \\
   -sample {sample} \\
   -title {project} \\
   -results {qc_report_dir}
'''.format (nest_circ_qc_report=var.nest_circ_qc_report,
            circ_qc_dir=circ_qc_dir,
            sample=sample,
            project=project,
            qc_report_dir=qc_report_dir, )
        qc_fn = var.circ_qc_fn_sh (v['project_dir'])
        write_str_to_file (pycode, qc_fn)

    def write_qc_full(self):
        v = self.argv
        sample_list = re.split(':|,', v['sample'])
        qc_dir = var.circ_qc_dir(v['project_dir'])
        circ_mkdir(qc_dir)
        qc_fn = var.circ_qc_fn(v['project_dir'])
        pycode = str()
        pycode +='''
import os, sys
import os.path
import misopy.cluster_utils as cluster_utils

LIB_CIRC = '{8}'
sys.path.insert(0, LIB_CIRC)
import circ
import circ_var_defs as var

project = '{0}'
proj_dir = '{1}'
fq = '{2}'
ss = {12}
fa = '{3}'
gtf = '{4}'
sample = '{5}'
sample_list = {6}
mapfile = '{7}'
gtf2bed = {9}
circ_run_qc = {10}
circ_qc_report = {11}
'''.format(v['project'], v['project_dir'], v['fq'], v['genome'], v['gtf'], re.sub(':', ",", v['sample']), sample_list, v['mapfile'], LIB_CIRC, q(var.nest_gtf2bed), q(var.nest_circ_run_qc), q(var.nest_circ_qc_report), q(v['library_type']))

        pycode += '''
qc_dir = var.circ_qc_dir(proj_dir)
circ.circ_mkdir(qc_dir)
os.chdir(qc_dir)
exon_gtf = os.path.join(qc_dir, 'exon.gtf')
sorted_gtf = os.path.join(qc_dir, 'sotred.gtf')
sorted_bed = os.path.join(qc_dir, 'sotred.bed')
generate_qc_fn = os.path.join(qc_dir, 'generate_QC.sh')

f = open(generate_qc_fn, 'w')
f.write('awk \\'{if($3==\\"exon\\"){print $0}}\\' %s > %s\\n' % (gtf, exon_gtf))
f.write('msort -k mf1 -k nf4 %s > %s\\n' % (exon_gtf, sorted_gtf))
f.write('perl %s %s %s\\n' % (gtf2bed, sorted_gtf, sorted_bed))
f.write(\'\'\'
perl %s \\\\
    -fq %s \\\\
    -se-pe pe \\\\
    -n %s \\\\
    -o %s \\\\
    -spe %s \\\\
    -R %s \\\\
    -G %s \\\\
    -bed %s \\\\
    -mapfile %s
\'\'\' %(circ_run_qc, fq, sample, qc_dir, ss, fa, gtf, sorted_bed, mapfile))
f.close()
circ.circ_call_process('sh %s' % generate_qc_fn)

job_ids_qc=[]
for s in sample_list:
    sample_dir = os.path.join(qc_dir, s)
    os.chdir(sample_dir)
    qc_qsub = circ.qsub_cmd('8G', 4, s + '_QC.sh', sample_dir)
    job_ids_qc.append(cluster_utils.launch_job(qc_qsub, cmd_name = 'qsub'))
os.chdir(proj_dir)
cluster_utils.wait_on_jobs(job_ids_qc, cluster_cmd='qsub')

# generate QCreport
qc_report_dir = os.path.join(qc_dir, 'QC_report')
circ.circ_mkdir(qc_report_dir)
qc_report_fn = os.path.join(qc_dir, 'qc_report.sh')
qc_report_cmd = 'sh %s -dir %s -sample %s -title %s -results %s' % (circ_qc_report, qc_dir, sample, project, qc_report_dir)
circ.write_str_to_file(qc_report_cmd, qc_report_fn)
qc_report_qsub = circ.qsub_cmd('1G', 1, qc_report_fn, qc_dir)
cluster_utils.launch_job(qc_report_qsub, cmd_name='qsub')

'''
        write_str_to_file(pycode, qc_fn)


class novel_circ_by_circexplorer:
    '''wrapper class for all samples.'''

    def __init__(self, PROJ_DIR, samples, genome_index, gtf, circ_finder='CIRCexplorer', tophat_dir=None, fq_clean=None, skip_step=None):
        "Find circRNA by CIRCexplorer."
        self.PROJ_DIR = PROJ_DIR
        self.samples = samples
        self.genome_index = genome_index
        self.gtf = gtf
        self.tophat_dir = tophat_dir
        self.fq_clean = fq_clean
        self.skip_step = skip_step
        self.novel_dir = var.circ_novel_dir(PROJ_DIR)

    def generate_fns(self):
        samples = self.samples
        qsub_fns = []
        # if not (argv['skip_step'] and ('find-circrna' in argv['skip_step'])):
        #     flag_run_circexplorer = True
        #     genepred = gtf_to_genepred(self.novel_dir, gtf)
        for i in range(0, len(samples)):
            sample = samples[i]
            # log('[novel circ] [START] analysis for sample "%s"' % sample)
            py_cmd = '''\
%s
PROJ_DIR = %s
sample = %s
genome = %s
gtf = %s
my_tophat_dir = %s
my_fq_clean = %s
skip_step = %s


my_circ = circ.novel_circ(PROJ_DIR, sample, genome, gtf, tophat_dir=my_tophat_dir, fq_clean=my_fq_clean, skip_step=skip_step)
my_circ.run_tophat()
my_circ.run_bam2fq()
my_circ.run_tophat_fusion()
my_circ.find_circ()
my_circ.parse_circ_table()
my_circ.get_circ_seq()
my_circ.write_read_count()

''' % (CIRC_PY_HEADER,
       q(self.PROJ_DIR),
       q(sample),
       q(self.genome_index),
       q(self.gtf),
       q(None if self.tophat_dir is None else as_list(self.tophat_dir)[i]),
       q(None if self.fq_clean is None else as_list(self.fq_clean)[i]),
       q(self.skip_step))
            py_fn = var.qsub_tophat_circexplorer_py_fn(self.novel_dir, sample)
            write_str_to_file(py_cmd, py_fn)
            sh_cmd = 'python %s\n' % py_fn
            sh_fn = var.qsub_tophat_circexplorer_sh_fn(self.novel_dir, sample)
            write_str_to_file(sh_cmd, sh_fn)
            qsub_fns.append(sh_fn)
        return qsub_fns


class novel_circ:
    '''circRNA'''
    def __init__(self, PROJ_DIR, sample, genome_index, gtf, circ_finder='CIRCexplorer', tophat_dir=None, fq_clean=None, skip_step=None, run_type='call'):
        self.PROJ_DIR = PROJ_DIR
        self.qc_dir = var.circ_qc_dir(PROJ_DIR)
        circ_mkdir(self.qc_dir)
        circ_mkdir(os.path.join(self.qc_dir, sample)) # NOTE: tophat needs this dir
        self.novel_dir = var.circ_novel_dir(PROJ_DIR)
        circ_mkdir(self.novel_dir)
        self.circ_novel_fn_prefix = var.circ_novel_fn_prefix(self.novel_dir, sample) # NOTE: CIRCexplorer will add sffix '_circ_txt' in the output filename
        self.circ_novel_fn = var.circ_novel_fn(self.novel_dir, sample)
        self.sample = sample
        self.circ_seq_fn = var.circ_seq_fn(self.novel_dir, sample)
        self.circ_seq_len_fn = var.circ_seq_len_fn(self.novel_dir, sample)
        self.tophat_dir = var.circ_sample_tophat_dir(self.qc_dir, self.sample)
        if tophat_dir:                    # run 'mkdir' and 'ln -s'
            tophat_dir_old = tophat_dir
            unmapped_bam_old = os.path.join(tophat_dir_old, 'unmapped.bam')
            unmapped_bam_new = os.path.join(self.tophat_dir, 'unmapped.bam')
            if not os.path.isfile(unmapped_bam_old):
                raise ValueError('no such file: %s' % unmapped_bam_old)
            elif unmapped_bam_old != unmapped_bam_new:
                circ_mkdir(self.tophat_dir)
                circ_safe_ln(unmapped_bam_old, unmapped_bam_new)
        self.unmapped_bam = os.path.join(self.tophat_dir, 'unmapped.bam')
        self.bam2fq_sh_fn = var.circ_bam2fq_sh_fn(self.tophat_dir)
        self.genome_index = genome_index
        self.gtf = gtf
        self.circ_finder = circ_finder
        self.fq_clean = fq_clean
        if skip_step:
            self.skip_step = skip_step.split(',')
        else:
            self.skip_step = None
        self.run_type = run_type
        self.seq = None
        self.seq_len = None
        self.circ_table = None
        self.read_count = None
        self.read_count_fn = var.circ_read_count_fn(self.novel_dir, sample)

        # for qsub and wait jobs
        self.tophat_sh_fn = var.circ_tophat_sh_fn(self.qc_dir, sample)
        self.tophat_fusion_sh_fn = var.circ_tophat_fusion_sh_fn(self.qc_dir, sample)
        self.job_id_tophat = None
        self.job_id_bam2fq = None
        self.job_id_tophat_fusion = None

    def run_tophat(self):
        if self.skip_step and 'tophat' in self.skip_step:
            log('skipping step: tophat')
        else:
            return self.__run_tophat_really_run()

    def __run_tophat_really_run(self):
        prog_tophat = var.prog_tophat2
        args = var.prog_tophat2_circ_specific_args
        qc_dir = self.qc_dir
        tophat_dir = self.tophat_dir
        cmd_tophat = '{0} {1} -G {2} -o {3} {4} {5}'.format(prog_tophat, args, self.gtf, tophat_dir, self.genome_index, '{0}_1.clean.fq.gz  {0}_2.clean.fq.gz'.format(self.fq_clean))
        tophat_sh_fn = self.tophat_sh_fn
        run_by_run_type(self.run_type, cmd_tophat, tophat_sh_fn)
        # write_str_to_file(cmd_tophat, tophat_sh_fn)
        # self.job_id_tophat = launch_jobs('20G', 10, tophat_sh_fn)
        # wait_job(self.job_id_tophat)
        # return tophat_sh_fn

    def run_bam2fq(self):
        # bamToFastq
        qc_dir = self.qc_dir
        tophat_dir = self.tophat_dir
        unmapped_bam = self.unmapped_bam
        unmapped_fq = re.sub('\.bam$', '.fq', unmapped_bam)
        cmd_bam2fq = '%s -i %s -fq %s' % (var.prog_bam2fq, unmapped_bam, unmapped_fq)
        if self.skip_step and 'bam2fq' in self.skip_step:
            log('skipping step: bam2fq')
        else:
            run_by_run_type(self.run_type, cmd_bam2fq, self.bam2fq_sh_fn)
        # write_str_to_file(cmd_bam2fq, self.bam2fq_sh_fn)
        # self.job_id_bam2fq = launch_jobs('2G', 1, self.bam2fq_sh_fn)
        # return self.bam2fq_sh_fn

    def run_tophat_fusion(self):
        if self.skip_step and 'tophat-fusion' in self.skip_step:
            log('skipping step: tophat-fusion')
        else:
            return self.__run_tophat_fusion_really_run()

    def __run_tophat_fusion_really_run(self):
        # bamToFastq
        qc_dir = self.qc_dir
        tophat_dir = self.tophat_dir
        unmapped_bam = self.unmapped_bam
        unmapped_fq = unmapped_bam.replace('.bam', '.fq')

        # tophat-fusion
        tophat_fusion_dir = var.circ_tophat_fusion_dir(tophat_dir)
        cmd_tophat_fusion = '%s -o %s -p 10 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search %s %s' % (var.prog_tophat2, tophat_fusion_dir, self.genome_index, unmapped_fq)
        tophat_fusion_sh_fn = self.tophat_fusion_sh_fn
        run_by_run_type(self.run_type, cmd_tophat_fusion, tophat_fusion_sh_fn)
        # write_str_to_file(cmd_tophat_fusion, tophat_fusion_sh_fn)
        # self.job_id_tophat_fusion = launch_jobs('15G', 10, tophat_fusion_sh_fn)
        # wait_job(self.job_id_tophat_fusion)
        # return tophat_fusion_sh_fn

    def find_circ(self):
        '''This step starts from clean fq.'''
        circ_finder = self.circ_finder
        if self.skip_step and 'find-circrna' in self.skip_step:
            log('skipping step: find-circrna')
        else:
            if circ_finder == 'CIRCexplorer':
                return self.find_circ_by_circexplorer()
            else:
                error('other circRNA finders are not supported currently.')

    def find_circ_by_circexplorer(self):
        '''This step starts from clean fq.'''
        # CIRCexplorer
        tophat_dir = self.tophat_dir
        tophat_fusion_dir = var.circ_tophat_fusion_dir(tophat_dir)
        fusion_hit = os.path.join(tophat_fusion_dir, 'accepted_hits.bam')
        circ_novel_fn_prefix = self.circ_novel_fn_prefix
        genepred = var.circ_annotation_genepred(self.novel_dir, self.gtf)
        cmd_CIRCexplorer = '''
%s \\
    --fusion %s
    --genome %s
    --ref %s
    --output %s

perl -i'.bak' -wane '($f[12] >= 2) and print' %s_circ.txt
''' % (var.prog_CIRCexplorer, fusion_hit, self.genome_index, genepred, circ_novel_fn_prefix)
        # TODO: filter by readcount >= 2 and backup CIRCexplorer output

        circexplorer_sh_fn = var.circexplorer_sh_fn(self.novel_dir, self.sample)
        run_by_run_type(self.run_type, cmd_CIRCexplorer, circexplorer_sh_fn)
        # write_obj_to_file(cmd_CIRCexplorer, circexplorer_sh_fn)
        # circ_call_process(cmd_CIRCexplorer)
        # return circexplorer_sh_fn

    def parse_circ_table(self):
        f = open(self.circ_novel_fn, 'r')
        circ_table = {}
        read_count = {}
        for line in f:
            array = re.split("\s+", line)
            circ_id = '_'.join([array[i] for i in (0, 1, 2, 5)])
            circ_table[circ_id] = array
            read_count[circ_id] = array[12]
        f.close()
        self.read_count = read_count
        self.circ_table = circ_table

    # def get_circ_seq1(self):
    #     # FIXME: parse_fasta() slurps lots of RAM for all samples
    #     genome_dict = parse_fasta(self.genome_index)
    #     circ_table = self.circ_table
    #     seq = {}
    #     seq_len = {}
    #     for _id in circ_table.keys():
    #         _chr, start, end, strand = [circ_table[_id][i] for i in (0, 1, 2, 5)]
    #         circ_seq = genome_dict[_chr][int(start)-1:int(end)-1]
    #         if strand is '-':
    #             circ_seq = str(Seq(circ_seq).reverse_complement())
    #         seq[_id] = circ_seq
    #         seq_len[_id] = len(circ_seq)
    #     self.seq = seq
    #     self.seq_len = seq_len
    #     self.write_circ_seq()
    #     # write seq length
    #     f = open(self.circ_seq_len_fn, 'w')
    #     f.write('ID\tlength\n')
    #     for _id, seq_len in self.seq_len.items():
    #         f.write('%s\t%s\n' % (_id, seq_len))
    #     f.close()

    def write_circ_seq(self):
        f = open(self.circ_seq_fn, 'w')
        for _id, seq in self.seq.items():
            f.write('>%s\n%s\n' % (_id, seq))
        f.close()

    def plot_circ_seq_len(self):
        # f = open(self.circ_seq_len_fn, 'w')
        # f.write('ID\tlength\n')
        # for _id, seq_len in self.seq_len.items():
        #     f.write('%s\t%s\n' % (_id, seq_len))
        # f.close()

        # plot seq length distribution
        rscript = self.circ_seq_len_fn + '.R'
        out = self.circ_seq_len_fn + '.pdf'
        rcode = '''
dat <- read.table("{0}", header=TRUE)
pdf("{1}")
hist(dat[, 2], main="Length distribution", xlab="Length (nt)")
dev.off()
'''.format(self.circ_seq_len_fn, out)
        # write_str_to_file(rcode, rscript)

    def write_read_count(self):
        read_count_all_fn = var.circ_read_count_all_fn(self.novel_dir, self.sample)
        fh_all = open(read_count_all_fn, 'w')
        fh = open(self.read_count_fn, 'w')
        for circ_id, count in self.read_count.items():
            fh_all.write('%s\t%s\n' % (circ_id, count))
            if int(count) >= var.CIRC_NOVEL_COUNT_MIN:
                fh.write('%s\t%s\n' % (circ_id, count))
        fh_all.close()
        fh.close()


def circ_novel_plot_length(samples, novel_dir, sample_fns=None, prog='find_circ'):
    """ density plot of length of all novel circRNAs.
    Keyword Arguments:
    samples -- samples (type: list)
    novel_dir -- novel  dir
    sample_fns -- input filenames (data source), which does not need novel_dir
    """
    out_dir = var.circ_length_plot_dir(novel_dir)
    circ_mkdir(out_dir)
    if sample_fns == None:
        if novel_dir == None:
            raise ValueError ('missing novel_dir, and sample_fns not given')
        else:
            sample_fns = []
            for s in samples:
                fn_switch = {
                    'find_circ' : var.annotated_circ (novel_dir, s),
                    'ciri' : var.ciri_as_out_fn (novel_dir, s),
                    }
                sample_fns.append (fn_switch[prog])
    elif type(sample_fns) is not list:
        raise TypeError('sample_fns is not of type list: %s' % sample_fns)
    rscript = var.circ_length_plot_rscript(out_dir)
    length_plot_q75_fn = var.circ_length_plot_q75_fn(out_dir)
    length_plot_q100_fn = var.circ_length_plot_q100_fn(out_dir)
    length_plot_data_fn = var.circ_length_plot_data_fn(out_dir)
    rcode = '''
prog <- "%s"
samples <- "%s"
samples <- unlist(strsplit(samples, ","))
sample_fns <- "%s"
sample_fns <- unlist(strsplit(sample_fns, ","))

get_length_data <- function(samples, fns) {
    len <- vector()
    grp <- vector()
    for (i in 1:length(fns)) {
        fn <- fns[i]
        sample <- samples[i]
        if (prog == "find_circ") {
            dat <- read.table(fn, quote=NULL, header=TRUE)[, 6]
            len <- c(len, dat)
            grp <- c(grp, rep(sample, length(dat)))
        } else if (prog == "ciri") {
            dat <- read.csv(fn, sep="\\t")
            dat <- dat[, c(3,4)]
            len <- c(len, dat[, 2] - dat[, 1])
            grp <- c(grp, rep(sample, nrow(dat)))
        } else stop(paste0("invalid program name: ", prog))
    }
    return(data.frame(len=len, sample=grp))
}

dat <- get_length_data(samples, sample_fns)
out_data_fn <- "%s"
write.table(dat, out_data_fn, sep="\\t", quote=FALSE, row.name=FALSE)

suppressMessages(library("ggplot2"))
p <- ggplot(dat, aes(len, fill=sample))
p <- p + geom_bar(position="dodge")
p <- p + xlab("Length (nt)")
p <- p + ylab("Count")

qall <- quantile(dat$len)
p.q75 <- p + xlim(qall[1], qall[4])
p.q100 <- p + xlim(qall[1], qall[5])

suppressMessages(ggsave("%s", plot=p.q75, width=7, height=7))
suppressMessages(ggsave("%s", plot=p.q100, width=7, height=7))
''' % (prog,
            ','.join(samples),
            ','.join(sample_fns),
            length_plot_data_fn,
            length_plot_q75_fn,
            length_plot_q100_fn,
            )
    write_str_to_file(rcode, rscript)
    circ_run_script(rscript, 'R')


def circ_novel_plot_genomic_feature(samples, novel_dir, prog='find_circ'):
    """ exon/intron ratio plot.
    Keyword Arguments:
    samples   --
    novel_dir --
    """
    out_dir = var.circ_genomic_feature_plot_dir(novel_dir)
    circ_mkdir(out_dir)
    plot_fn = var.circ_genomic_feature_plot_fn(out_dir)
    data_fn = var.circ_genomic_feature_data_fn(out_dir)
    novel_fns = []
    header_switch = {'circexplorer' : 'FALSE', 'ciri' : 'TRUE', 'find_circ' : 'TRUE', }
    comment_switch = {'circexplorer' : '', 'ciri' : ', skip=1, sep="\\t"', 'find_circ' : '', }
    # index_switch: one column for filtering junction reads >= 2
    index_switch = {'circexplorer' : 'c(13, 14)', 'ciri' : 'c(5, 9)', 'find_circ' : 'c(10,9)', }
    header = header_switch[prog]
    comment = comment_switch[prog]
    index = index_switch[prog]
    for sample in samples:
        fn_switch = {
            'circexplorer' : var.circ_novel_fn (novel_dir, sample),
            'ciri' : var.ciri_out_fn (novel_dir, sample),
            'find_circ' : var.annotated_circ (novel_dir, sample),
            }
        novel_fn = fn_switch[prog]
        novel_fns.append(novel_fn)
    rscript = var.circ_genomic_feature_plot_rscript(out_dir)
    rcode = '''
suppressMessages(library(stringr))
prog <- "%s"
samples <- "%s"
samples <- unlist(strsplit(samples, ","))
novel_fns <- "%s"
novel_fns <- unlist(strsplit(novel_fns, ","))

get_genomic_feature <- function(samples, fns) {
    feature <- vector()
    grp <- vector()
    num <- vector()
    label <- vector()
    for (i in 1:length(fns)) {
        fn <- fns[i]
        sample <- samples[i]
        dat <- read.table(fn, quote=NULL, header=%s%s)[, %s]
        colnames(dat) <- c("read_count", "type")
        dat <- dat[dat$"read_count" >= 2, 2]
        if (prog == "circexplorer") {
            num_exon <- length(which(dat == "circRNA"))
            num_intron <- length(which(dat == "ciRNA"))

            # return data
            feature <- c(feature, "exon", "intron")
            grp <- c(grp, rep(sample, 2))
            num <- c(num, num_exon, num_intron)
            label <- c(label, num_intron, num_exon)
        } else if (prog == "find_circ") {
            num_exon <- sum(sapply(dat, function(x) str_count(x, "exon")))
            num_intron <- sum(sapply(dat, function(x) str_count(x, "intron")))
            num_intergenic <- length(which(dat == "--"))

            # return data
            feature <- c(feature, "exon", "intron", "intergenic")
            grp <- c(grp, rep(sample, 3))
            num <- c(num, num_exon, num_intron, num_intergenic)
            label <- c(label, num_intron, num_exon, num_intergenic)
        } else if (prog == "ciri") {
            num_exon <- length(which(dat == "exon"))
            num_intron <- length(which(dat == "intron"))
            num_intergenic <- length(which(dat == "intergenic region"))

            # return data
            feature <- c(feature, "exon", "intron", "intergenic")
            grp <- c(grp, rep(sample, 3))
            num <- c(num, num_exon, num_intron, num_intergenic)
            label <- c(label, num_intron, num_exon, num_intergenic)
        }
    }
    return(data.frame(feature=feature, sample=grp, num=num, label=label))
}

dat <- get_genomic_feature(samples, novel_fns)
out_data_fn <- "%s"
write.table(dat[, -4], out_data_fn, sep="\\t", quote=FALSE, row.name=FALSE)

suppressMessages(library("ggplot2"))
p <- ggplot(dat, aes(x=sample, y=num, fill=feature))
p <- p + geom_bar(stat="identity")
# p <- p + geom_text(aes(label=dat$label))
p <- p + xlab("Sample")
p <- p + ylab("Number")

SAMPLES_CHAR_LEN_MAX <- 40
SAMPLE_NUMBER_MAX <- 6
sample_number <- length(unique(dat$sample))
samples_char_len <- length(paste0(unique(dat$sample),  collapse = ''))
if (samples_char_len > SAMPLES_CHAR_LEN_MAX || sample_number > SAMPLE_NUMBER_MAX) {
    p <- p + coord_flip()
}

suppressMessages(ggsave("%s", plot=p, width=7, height=7))
''' % (prog, ','.join(samples), ','.join(novel_fns), header, comment, index, data_fn, plot_fn)
    write_str_to_file(rcode, rscript)
    circ_run_script(rscript, 'R')


def wrapper_circ_exp(PROJ_DIR, samples, groupname, compare, venn, replicate, argv_copy, *args):
    pycode = '''\
%s

PROJ_DIR = %s
samples = %s
groupname = %s
compare = %s
venn = %s
replicate = %s
argv_copy = %s

my_exp = circ.circ_exp(PROJ_DIR, samples, groupname, compare, venn, replicate, argv_copy)
my_exp.merge_exp()
my_exp.merge_source_gene()
my_exp.diff_exp()

''' % (CIRC_PY_HEADER, q(PROJ_DIR), q(samples), q(groupname), q(compare), q(venn), q(replicate), q(argv_copy))
    exp_dir = var.circ_exp_dir(PROJ_DIR)
    py_fn = var.qsub_exp_py_fn(exp_dir)
    write_str_to_file(pycode, py_fn)
    return py_fn


class circ_exp:
    '''expression'''
    def __init__(self, PROJ_DIR, samples, groupname, compare, venn, replicate, argv_copy, prog='circexplorer', linear_free=True):
        self.samples, self.groupname, self.compare, self.venn, self.replicate \
          = samples, groupname, compare, venn, replicate
        self.argv_copy = argv_copy
        self.prog = prog
        self.linear_free = linear_free
        self.PROJ_DIR = PROJ_DIR
        self.novel_dir = var.circ_novel_dir(PROJ_DIR)
        self.exp_dir = var.circ_exp_dir(PROJ_DIR)
        circ_mkdir (self.exp_dir)
        self.source_gene_fn = var.circ_source_gene_fn (self.exp_dir)
        self.circ_merged_exp_fn = var.circ_merged_exp_fn (self.exp_dir)
        self.diff_exp_fn = {}

    def merge_exp(self):
        # NOTE: ID is named based on the same circRNA sequence (chr, start, end, strand)
        my_id = {}
        merged = {}
        index_id = 0
        index_sample = 9
        index_njexp = 11
        if self.linear_free:
            index_jexp = 10                            # junction reads
        else:
            index_jexp = 10
            # index_exp = (10, 11) # junction reads + non-junction reads
        fh = open (var.io_novel_circ (self.novel_dir))
        fh.readline ()
        for line in fh:
            # TODO: when calculating expression level, there is no difference
            # for different `linear_free', but this may change later.
            line = line.strip ()
            # print line
            if self.linear_free:
                _id, s, jexp, njexp = [line.split ('\t')[i] for i in (index_id, index_sample, index_jexp, index_njexp)]
                my_id[_id] = 1
                len1 = len (s.split (','))
                for i in range (0, len1):
                    merged.setdefault(s.split (',')[i], {})[_id] = jexp.split (',')[i] # + njexp.split (',')[i]
            else:
                _id, s, jexp, njexp = [line.split ('\t')[i] for i in (index_id, index_sample, index_jexp, index_njexp)]
                my_id[_id] = 1
                len1 = len (s.split (','))
                for i in range (0, len1):
                    merged.setdefault(s.split (',')[i], {})[_id] = jexp.split (',')[i] # + njexp.split (',')[i]
        fh.close ()
        # for s in self.samples ():
        #     fh = open (fn, 'r')
        #     for line in fh:
        #         _id, _exp = [re.split('\t', line)[i] for i in (index_id, index_exp)]
        #         if _exp >= 2:
        #             my_id[_id] = 1
        #             merged.setdefault(sample, {})[_id] = _exp
        #     fh.close()
        out = open(self.circ_merged_exp_fn, 'w')
        header_line = 'circRNA_id\t%s\n' % '\t'.join(self.samples)
        out.write(header_line)
        for _id in sorted(my_id.keys()):
            out.write (_id)
            # for sample in sorted(merged.keys()):
            for s in self.samples:   # fix bug (2015-08-21)
                if merged.has_key (s):
                    # print s
                    if merged[s].has_key (_id):
                        out.write ('\t%s' % merged[s][_id])
                    else:
                        out.write ('\t0')
                else:
                    out.write ('\t0')
            out.write ('\n')
        out.close ()

    def merge_source_gene(self):
        """
        merge circRNAs and source genes from all samples.
        """
        h = {}
        source_gene_fn = self.source_gene_fn
        # index_switch_id = {'circexplorer' : (0, 1, 2, 5), 'ciri' : (1, 2, 3)}
        # index_switch_source_gene = {'circexplorer' : 14, 'ciri' : 9}
        index_id = 0
        index_source_gene = 7
        fh = open (var.io_novel_circ (self.novel_dir))
        fh.readline ()
        for line in fh:
            circ_id, source_gene = [line.split ('\t')[i] for i in (index_id, index_source_gene)]
            # FIXME: the first source gene is reserved, since the pair should be '1-to-1' mapping
            source_gene = re.sub ('"?([^;]+?);.*', '\g<1>', source_gene)
            if (source_gene[0] is not '-') and (len (source_gene) > 1):
                h[circ_id] = source_gene
        fh.close ()
        write_obj_to_file (h, source_gene_fn)
        # for s in self.samples:
        #     # fn_switch = {'circexplorer' : var.circ_novel_fn(self.novel_dir, s), 'ciri' : var.ciri_out_fn(self.novel_dir, s)}
        #     circ_novel_fn = var.annotated_circ (self.novel_dir, s)
        #     fh = open(circ_novel_fn, 'r')
        #     for line in fh:
        #         array = re.split("\s+", line)
        #         # circ_id = '_'.join([array[i] for i in index_id])
        #         circ_id = array[index_id]
        #         source_gene = array[index_source_gene]
        #         merged_source_gene[circ_id] = source_gene
        #     fh.close()
        # write_obj_to_file(h, source_gene_fn)

    def diff_exp(self):
        if self.replicate:                # TODO: merge nested code
            self.__diff_exp_rep()
        for compare_list in self.compare.keys():
            grp, samples_within_grp = self.compare[compare_list]
            y2, y1 = [len(aref) for aref in samples_within_grp]
            arg_y = 'c(rep(2, %d), rep(1, %d))' % (y2, y1)
            # print arg_y
            resp_type = 'Two class unpaired' if len(grp) == 2 else 'Multiclass'
            str_grp = '-'.join(grp)
            str_samples = ','.join(','.join(s) for s in samples_within_grp)
            diff_exp_fn = var.circ_diff_exp_fn(self.exp_dir, str_grp)
            self.diff_exp_fn[str_grp] = diff_exp_fn
        fn_merged = self.circ_merged_exp_fn
        func_switch = { True : self.__diff_exp_rep, False : self.__diff_exp_norep }
        # func_switch[self.replicate]()
        # main entry: run_DE_sRNA_v2.pl
        _argv = self.argv_copy
        _sample = _argv['_sample'].replace(':', ",")
        _group = _argv['_group']
        _groupname = _argv['_groupname']
        _compare = _argv['_compare']
        _venn = _argv['_venn']
        # if venn is sth. like '2:1', then the venn arg is not passed to
        # run_DE_sRNA_v2.pl
        if _venn:
            if re.search('_', _venn):
                arg_venn = '--venn $venn'
            else:
                arg_venn = ''
        else:
            arg_venn = ''
        shfile = var.circ_generate_diff_fn(self.exp_dir)

        shcode = '''
readcount={0}
sample={1}
group={2}
groupname={3}
compare={4}
venn={5}
outdir={6}
DE=$outdir/DE.sh

perl /TJPROJ1/RNA/kangmingming/dev/circrna/script/nest/run_DE_sRNA_v2.pl \\
    -r $readcount \\
    -s $sample \\
    --group $group \\
    --groupname $groupname \\
    -g $compare \\
    {7} \\
    -o $outdir \\
    > $DE
'''.format(fn_merged, _sample, _group, _groupname, _compare, _venn, self.exp_dir, arg_venn)
        write_str_to_file(shcode, shfile)
        circ_run_script(shfile, 'sh')
        circ_run_script(os.path.join(os.path.dirname(shfile), 'DE.sh'), 'sh')

    def __diff_exp_rep(self):
        """
        differential expression for samples without replicates.
        """
        fn_merged = self.circ_merged_exp_fn
        for compare_list in self.compare.keys():
            grp, samples_within_grp = self.compare[compare_list]

            y2, y1 = [len(aref) for aref in samples_within_grp]
            arg_y = 'c(rep(2, %d), rep(1, %d))' % (y2, y1)

            name2, name1 = [','.join(aref) for aref in samples_within_grp]

            resp_type = 'Two class unpaired' if len(grp) == 2 else 'Multiclass'
            str_grp = '-'.join(grp)
            str_samples = ','.join(','.join(s) for s in samples_within_grp)
            diff_exp_fn = var.circ_diff_exp_fn(self.exp_dir, str_grp)
            self.diff_exp_fn[str_grp] = diff_exp_fn
            diff_exp_up_fn = var.circ_diff_exp_up_fn(self.exp_dir, str_grp)
            diff_exp_down_fn = var.circ_diff_exp_down_fn(self.exp_dir, str_grp)

            # output file name
            vs = var.circ_str_vs
            grp_vs = vs.join(grp)
            circ_mkdir_unix(os.path.join(self.exp_dir, grp_vs))
            diff_each_table_all_fn = var.circ_diff_each_table_all_fn(self.exp_dir, grp_vs)
            diff_each_table_up_fn = var.circ_diff_each_table_up_fn(self.exp_dir, grp_vs)
            diff_each_table_down_fn = var.circ_diff_each_table_down_fn(self.exp_dir, grp_vs)
            diff_each_list_all_fn = var.circ_diff_each_list_all_fn(self.exp_dir, grp_vs)
            diff_each_list_up_fn = var.circ_diff_each_list_up_fn(self.exp_dir, grp_vs)
            diff_each_list_down_fn = var.circ_diff_each_list_down_fn(self.exp_dir, grp_vs)
            diff_each_volcano_pdf = var.circ_diff_each_volcano_pdf(self.exp_dir, grp_vs)
            diff_each_volcano_png = var.circ_diff_each_volcano_png(self.exp_dir, grp_vs)

            rscript = var.circ_diff_each_rscript(self.exp_dir, grp_vs)
# library("PoissonSeq")

# my_data_orig <- read.table("%s", header=TRUE, stringsAsFactors=FALSE)
# grp <- unlist(strsplit("%s", "-"))
# sample_names2 <- unlist(strsplit("%s", ","))
# sample_names1 <- unlist(strsplit("%s", ","))
# my_data <- subset(my_data_orig, select = c(%s))
# gname <- my_data_orig[, 1]
# n <- as.matrix(my_data)
# y <- %s
# # type <- "%%"
# # pair <- FALSE

# # dat <- list(n=n, y=y, type=type, pair=pair, gname=gname)
# # PS.Main(dat, para = para)

# # seq.depth <- PS.Est.Depth(n)

# library("samr")
# samfit <- SAMseq(x = n, y = y, resp.type = "%s", geneid=gname, genenames=gname)
# fc_table_up <- samfit$siggenes.table$genes.up
# fc_table_lo <- samfit$siggenes.table$genes.lo
# fc_table_all <- rbind(fc_table_up, fc_table_lo)
# write.table(fc_table_up, file = "%s", quote = FALSE, sep = "\\t", row.names = FALSE)
# write.table(fc_table_lo, file = "%s", quote = FALSE, sep = "\\t", row.names = FALSE)
# write.table(fc_table_all, file = "%s", quote = FALSE, sep = "\\t", row.names = FALSE)
            rcode = '''
## samr
my_data_orig <- read.table("%s", header=TRUE, stringsAsFactors=FALSE)
grp <- unlist(strsplit("%s", "-"))
sample_name2 <- unlist(strsplit("%s", ","))
sample_name1 <- unlist(strsplit("%s", ","))
my_data <- subset(my_data_orig, select = c(%s))
gname <- my_data_orig[, 1]
n <- as.matrix(my_data)
y <- %s

suppressMessages(library("samr"))
samfit <- SAMseq(x = n, y = y, resp.type = "%s", genenames=gname)
## pval <- samr.pvalues.from.perms(samfit$samr.obj$tt,  samfit$samr.obj$ttstar)
fc_table_up <- as.data.frame(samfit$siggenes.table$genes.up, stringsAsFactors=FALSE)
diff_up <- data.frame(
    gene_id = as.character(fc_table_up$`Gene ID`),
    log2_fold_change = log2(as.double(fc_table_up$`Fold Change`)),
    qvalue = as.double(fc_table_up$`q-value(%%)`) / 100,
    significant = (as.double(fc_table_up$`q-value(%%)`) / 100) < 0.20
    )

fc_table_lo <- as.data.frame(samfit$siggenes.table$genes.lo, stringsAsFactors=FALSE)
diff_down <- data.frame(
    gene_id = as.character(fc_table_lo$`Gene ID`),
    log2_fold_change = log2(as.double(fc_table_lo$`Fold Change`)),
    qvalue = as.double(fc_table_lo$`q-value(%%)`) / 100,
    significant = (as.double(fc_table_lo$`q-value(%%)`) / 100) < 0.20
    )

samnorm <- as.data.frame(samr.norm.data(n))
samnorm$circ_id <- gname

## calculate mean for each group
grp2 <- grp[1]
grp1 <- grp[2]
eval(parse(text=sprintf("samnorm$mean_%%s <- apply(samnorm[, unlist(strsplit(\\"%%s\\", \\",\\"))], 1, mean)", grp2, paste(sample_name2, collapse=","))))
eval(parse(text=sprintf("samnorm$mean_%%s <- apply(samnorm[, unlist(strsplit(\\"%%s\\", \\",\\"))], 1, mean)", grp1, paste(sample_name1, collapse=","))))

if (nrow(diff_up) > 0) {
    for (i in 1:nrow(diff_up)) {
        id <- diff_up[i, 1]        # matrix [i, 1] = Gene ID
        mean2 <- eval(parse(text=sprintf("subset(samnorm, samnorm$circ_id == \\"%%s\\", select=mean_%%s)", id, grp2)))
        mean1 <- eval(parse(text=sprintf("subset(samnorm, samnorm$circ_id == \\"%%s\\", select=mean_%%s)", id, grp1)))
        eval(parse(text=sprintf("diff_up[%%s, \\"mean_%%s\\"] <- %%s", i, grp2, mean2)))
        eval(parse(text=sprintf("diff_up[%%s, \\"mean_%%s\\"] <- %%s", i, grp1, mean1)))
    }
}

if (nrow(diff_down) > 0) {
    for (i in 1:nrow(diff_down)) {
        id <- diff_down[i, 1]        # matrix [i, 1] = Gene ID
        mean2 <- eval(parse(text=sprintf("subset(samnorm, samnorm$circ_id == \\"%%s\\", select=mean_%%s)", id, grp2)))
        mean1 <- eval(parse(text=sprintf("subset(samnorm, samnorm$circ_id == \\"%%s\\", select=mean_%%s)", id, grp1)))
        eval(parse(text=sprintf("diff_down[%%s, \\"mean_%%s\\"] <- %%s", i, grp2, mean2)))
        eval(parse(text=sprintf("diff_down[%%s, \\"mean_%%s\\"] <- %%s", i, grp1, mean1)))
    }
}

diff_all <- rbind(diff_up, diff_down)

columns <- c("gene_id", paste0("mean_", grp2), paste0("mean_", grp1), "log2_fold_change", "qvalue")

if (nrow(diff_up) > 0) {
    diff_up <- diff_up[, columns]
    write.table(diff_up, file="%s", sep="\\t", quote=F, row.name=F)
    write.table(diff_up[, 1], file="%s", sep="\\t", quote=F, row.name=F, col.name=F)
}

if (nrow(diff_down) > 0) {
    diff_down <- diff_down[, columns]
    write.table(diff_down, file="%s", sep="\\t", quote=F, row.name=F)
    write.table(diff_down[, 1], file="%s", sep="\\t", quote=F, row.name=F, col.name=F)
}

diff_all <- diff_all[, columns]
write.table(diff_all, file="%s", sep="\\t", quote=F, row.name=F)
write.table(diff_all[, 1], file="%s", sep="\\t", quote=F, row.name=F, col.name=F)

k <- diff_all
for(i in 1:nrow(k)) {
    if(k[i, "qvalue"] <= 0.2 & k[i, "log2_fold_change"] > 0)
        k[i, "type"] <- "up"
    else if (k[i, "qvalue"] <= 0.2 & k[i, "log2_fold_change"] < 0)
        k[i, "type"] <- "down"
    else
        k[i, "type"] <- "false"
}
diff_all <- k

suppressMessages(library("ggplot2"))

volcano_plot <- function(dat, titles) {
    num_up <- length(which(dat$`log2_fold_change` > 0))
    num_down <- length(which(dat$`log2_fold_change` < 0))
    p <- ggplot(dat, aes(x=log2_fold_change, y=-log10(qvalue), colour=type))
    p <- p + geom_point()
    p <- p + coord_cartesian(xlim=c(-5, 5))
    p <- p + geom_hline(yintercept=-log10(0.05), linetype=2)
    p <- p + labs(title=titles,
                  x=bquote(paste(log[2], "(fold change)", sep="")),
                  y=bquote(paste(-log[10],"(padj)",sep="")))

    ## key, custom colour!
    p <- p + scale_colour_manual(
        limits=c("up","down","false"),
        values=c("red", "green", "blue"),
        breaks=c("up", "down"),
        labels=c(paste("up: ", num_up), paste("down: ", num_down)))

    ggsave(filename="%s", plot=p, width=7, height=7)
    ggsave(filename="%s", type="cairo-png", plot=p, width=7, height=7)
}

volcano_plot(diff_all, paste0(grp, collapse = " vs. "))

''' % (fn_merged,
       str_grp,
       name2,
       name1,
       str_samples,
       arg_y,
       resp_type,
       diff_each_table_all_fn,
       diff_each_table_up_fn,
       diff_each_table_down_fn,
       diff_each_list_all_fn,
       diff_each_list_up_fn,
       diff_each_list_down_fn,
       diff_each_volcano_pdf,
       diff_each_volcano_png)
            write_str_to_file(rcode, rscript)
            # circ_run_script(rscript, 'R', allow_error=True)

    def __diff_exp_norep(self):
        """
        differential expression for samples without replicates.
        """
        for compare_list in self.compare.keys():
            grp, samples_within_grp = self.compare[compare_list]
            y2, y1 = [len(aref) for aref in samples_within_grp]
            arg_y = 'c(rep(2, %d), rep(1, %d))' % (y2, y1)
            # print arg_y
            resp_type = 'Two class unpaired' if len(grp) == 2 else 'Multiclass'
            str_grp = '-'.join(grp)
            str_samples = ','.join(','.join(s) for s in samples_within_grp)
            diff_exp_fn = var.circ_diff_exp_fn(self.exp_dir, str_grp)
            self.diff_exp_fn[str_grp] = diff_exp_fn
            diff_exp_up_fn = var.circ_diff_exp_up_fn(self.exp_dir, str_grp)
            diff_exp_down_fn = var.circ_diff_exp_down_fn(self.exp_dir, str_grp)
            rscript = var.circ_diff_rscript(self.exp_dir, str_grp)

            # TODO:

    def diff_exp_venn(self):
        venn = self.venn
        for venn_each in venn.keys():
            grp_list = venn[venn_each]
            venn_infix = '__'.join('_'.join(s) for s in grp_list)
            venn_fn = var.circ_diff_exp_venn_fn(self.exp_dir, venn_infix)
            venn_out_fn = var.circ_diff_exp_venn_out_fn(self.exp_dir, venn_infix)
            grp_list_all = []
            diff_exp_fn_all = []
            for grp in grp_list:
                str_grp = '-'.join(grp)
                # FIXME: 'MT-WT' != 'WT-MT'
                if self.diff_exp_fn.has_key(str_grp):
                    diff_exp_fn = self.diff_exp_fn[str_grp]
                else:
                    # NOTE: str_grp not updated
                    diff_exp_fn = self.diff_exp_fn['-'.join(grp[::-1])]
                if not os.path.isfile(diff_exp_fn):
                    raise TypeError('no such file: %s' % diff_exp_fn)
                grp_list_all.append(str_grp)
                diff_exp_fn_all.append(diff_exp_fn)
            if len(grp_list_all) != len(diff_exp_fn_all):
                raise ValueError('group and filename are not equal in length')
            self.__plot_diff_exp_venn(venn_fn, venn_out_fn, ','.join(grp_list_all), ','.join(diff_exp_fn_all))

    def __plot_diff_exp_venn(self, venn_fn, venn_out_fn, grp_list, diff_exp_fn):
        rcode = '''
suppressMessages(library(VennDiagram))
grp_list <- "%s"
grp_list <- unlist(strsplit(grp_list, ","))
colors <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
color <- colors[1:length(grp_list)]
diff_exp_fn <- "%s"
diff_exp_fn <- unlist(strsplit(diff_exp_fn, ","))

get_venn_data <- function(grp_list, diff_exp_fn) {
    venn_data <- list()
    for (i in 1:length(grp_list)) {
        fn <- diff_exp_fn[i]
        dat <- read.table(fn, header = TRUE, sep="\\t")[, 1]
        text <- sprintf("venn_data$`%%s` <- dat", grp_list[i])
        eval(parse(text = text))
    }
    return(venn_data)
}

venn_data <- get_venn_data(grp_list, diff_exp_fn)
venn.diagram(venn_data, filename = "%s", alpha = 0.5, fill = color, imagetype = "png")
''' % (grp_list, diff_exp_fn, venn_out_fn)
        write_str_to_file(rcode, venn_fn)
        # circ_run_script(venn_fn, 'R')

    def tissue_specific_exp(self):
        fn_merged = self.circ_merged_exp_fn
        fn_heatmap_rscript = var.tse_heatmap_fn(self.exp_dir)
        # TODO: TSE
        rcode = '''

'''

class mirna_target ():
    '''miRNA binding site.'''
    def __init__ (self, PROJ_DIR, samples, mirna, kingdom, genome):
        self.PROJ_DIR = PROJ_DIR
        self.samples = samples
        self.novel_dir = var.circ_novel_dir(PROJ_DIR)
        self.mirna_dir = var.circ_mirna_dir(PROJ_DIR)
        circ_mkdir(self.mirna_dir)
        self.novel_seq_merged_fn = var.novel_seq_merged_fn (self.mirna_dir)
        self.mirna = mirna
        self.mirna_mature_fn = var.mirna_mature_fn (self.mirna_dir)
        if len (mirna) == 3:
            self.mirna = self.generate_mirna_seq_fasta ()
        elif not os.path.isfile (self.mirna):
            raise ValueError ('no such file: %s' % self.mirna)
        elif self.mirna != self.mirna_mature_fn:
            circ_call_process ('ln -fs %s %s' % (self.mirna, self.mirna_mature_fn))
        self.mirna_len_max = get_seq_len_max (self.mirna)
        self.kingdom = kingdom
        self.genome = genome
        self.io_novel_circ = var.io_novel_circ (self.novel_dir)
        self.novel_seq_merged_fn = var.novel_seq_merged_fn (self.mirna_dir)
        # self.sample = instance_novel_circ.sample
        # self.seq = instance_novel_circ.seq
        # self.circ_seq_fn = instance_novel_circ.circ_seq_fn
        self.target_software = None
        # super(self.__class__, self).__init__(in_circ, )

    def get_circ_seq_spliced (self):
        print 'Max length of miRNA:', self.mirna_len_max
        get_circ_seq (self.io_novel_circ, self.genome, self.novel_seq_merged_fn, self.mirna_len_max)

    def generate_mirna_seq_fasta (self):
        # mature miRNA: /TJPROJ1/RNA/kangmingming/dev/circrna/data/mature-v21.fa
        cmd = "grep -A1 '^>%s-' %s | grep -v '^--' >%s" % (self.mirna, var.mirbase_mature_fn, self.mirna_mature_fn)
        circ_call_process (cmd)
        return self.mirna_mature_fn

    def find_paired_mirna(self):
        kingdom = self.kingdom
        # http://cbio.mskcc.org/microrna_data/manual.html
        if kingdom is 'animal':
            self.target_software = var.mirna_target_software_animal
            self.__run_target_software_animal()
        elif kingdom is 'plant':
            self.target_software = var.mirna_target_software_plant
            self.__run_target_software_plant()
        else:
            raise ValueError('kingdom "%s" is not supported' % kingdom)

    def __run_target_software_animal(self):
        prog = self.target_software
        args = var.mirna_target_software_animal_args
        # fn_out = var.mirna_target_list_fn(self.sample)

        # cmd = '%s %s %s %s -out %s' % (prog, self.mirna, self.circ_seq_fn, args, fn_out)
        # print cmd
        # circ_call_process(cmd)

        # TODO: merge circrna seq

        shfile = var.mirna_target_sh_fn(self.mirna_dir)
        shcode = '''
mature=%s
circrna=%s
out_dir=%s
perl /TJPROJ1/RNA/kangmingming/dev/circrna/script/nest/run_3UTR_annimal_target_v2.pl -m $mature -r $circrna -o $out_dir
''' % (self.mirna, self.novel_seq_merged_fn, self.mirna_dir)
        write_str_to_file(shcode, shfile)
        circ_run_script(shfile, 'sh')
        circ_run_script(os.path.join(os.path.dirname(shfile), 'runtarget.sh'), 'sh')


    def __run_target_software_plant(self):
        prog = self.target_software
        args = var.mirna_target_software_plant_args
        # fn_out = var.mirna_target_list_fn(self.mirna_dir, self.sample)
        # cmd = '%s -s %s -t -o %s %s' % (prog, self.mirna, self.circ_seq_fn, fn_out, args)
        # circ_call_process(cmd)

        shfile = var.mirna_target_sh_fn(self.mirna_dir)
        shcode = '''
mature=%s
circrna=%s
out_dir=%s
perl /TJPROJ1/RNA/kangmingming/dev/circrna/script/nest/run_gene_plant_target.pl -m $mature -r $circrna -o $out_dir

''' % (self.mirna, self.novel_seq_merged_fn, self.mirna_dir)
        write_str_to_file(shcode, shfile)
        circ_run_script(shfile, 'sh')
        circ_run_script(os.path.join(os.path.dirname(shfile), 'runtarget.sh'), 'sh')

    def plot_mirna_circrna_network (self):
        shcode = '''
{Rscript_bin} {script} {data} {out_prefix}
'''.format(Rscript_bin=var.prog_Rscript,
           script=var.mirna_circrna_network_plot,
           # TODO: 'binding_site.pair.txt' is not defined in my circ lib, so
           # there is no choice but to use the filename itself.
           data=jp (self.mirna_dir, 'binding_site.pair.txt'),
           out_prefix=var.mirna_circrna_network_out_prefix (self.mirna_dir),
           )
        circ_call_process (shcode)


def wrapper_mirna_target(PROJ_DIR, samples, mirna, kingdom, genome):
    """
    Keyword Arguments:
    PROJ_DIR --
    samples  --
    mirna    --
    kingdom  --
    """
    mirna_dir = var.circ_mirna_dir(PROJ_DIR)
    pycode = '''\
{header}
PROJ_DIR = {proj_dir}
sys.path.insert(0, PROJ_DIR)
exec('from circ_vars import *')

my_target = circ.mirna_target (PROJ_DIR, samples, mirna, kingdom, genome)
my_target.get_circ_seq_spliced ()
my_target.find_paired_mirna ()
my_target.plot_mirna_circrna_network ()
'''.format (header=CIRC_PY_HEADER,
            proj_dir=q(PROJ_DIR),
            )
    py_fn = var.qsub_mirna_py_fn(mirna_dir)
    write_str_to_file(pycode, py_fn)
    # sh_code = 'python %s' % py_fn
    # sh_fn = var.qsub_mirna_sh_fn(mirna_dir)
    # write_str_to_file(sh_code, sh_fn)
    return py_fn


class enrichment:
    '''GO/KEGG enrichment.'''
    # $indir in enrich_ref_new.pl: out dir of diff
    # change to "_vs_": /TJPROJ1/RNA/kangmingming/dev/circrna/script/run_DE_sRNA_v2.pl
    def __init__(self, proj_dir, genome, gtf, go, kegg, compare, run_type='call'):
        "docstring"
        self.enrich_dir = var.enrich_dir(proj_dir)
        circ_mkdir(self.enrich_dir)
        self.generate_enrich_fn = var.circ_generate_enrich_fn(self.enrich_dir)
        self.exp_dir = var.circ_exp_dir(proj_dir)
        self.source_gene_fn = var.circ_source_gene_fn(self.exp_dir)
        self.genome = genome
        self.gtf = gtf
        self.go = go
        self.kegg = kegg
        self.compare = compare

    def run_enrich(self):
        shfile = self.generate_enrich_fn
        compare_name_list = []
        job_fns = []
        for key, grp_info in self.compare.iteritems():
            compare_name_each = var.circ_str_vs.join(grp_info[0])
            compare_name_list.append(compare_name_each)
        compare_name = ','.join(compare_name_list)
        shcode = '''
exp_dir='{0}'
target_pair='{1}'
genome='{2}'
gtf='{3}'
go='{4}'
kegg='{5}'
compare='{6}'
out_dir='{7}'

perl /TJPROJ1/RNA/kangmingming/dev/circrna/script/nest/enrich_ref_new.pl \\
     -i $exp_dir \\
     -c $compare \\
     -t $target_pair \\
     --species $kegg \\
     --genome $genome \\
     --gtf $gtf \\
     --goann $go \\
     -o $out_dir
'''.format(self.exp_dir, self.source_gene_fn, self.genome, self.gtf, self.go, self.kegg, compare_name, self.enrich_dir)
        write_str_to_file(shcode, shfile)
        # NOTE: 'the generated run_enrich.sh' will qsub *its* 2nd-genereted script
        circ_run_script (shfile, 'sh')
        circ_run_script (jp (dirname (shfile), 'run_enrich.sh'), 'sh')
        fns = glob.glob (jp (dirname (shfile), '*%s*_runenrich.sh' % var.circ_str_vs))
        if len (fns) == 0:
            raise ValueError ('enrich script(s) not found')
        return fns


class novel_circ_by_ciri:
    '''CIRI.'''
    def __init__(self, PROJ_DIR, samples, genome, gtf, fq_clean=None, mem='8G', cpu=4):
        self.PROJ_DIR = PROJ_DIR
        self.samples = samples
        self.genome = genome
        self.gtf = gtf
        # self.fq = fq                      # fq_dir actually
        self.fq_clean = fq_clean
        self.mem = mem
        self.cpu = cpu
        self.qc_dir = var.circ_qc_dir(PROJ_DIR)
        circ_mkdir(self.qc_dir)
        # circ_mkdir(os.path.join(self.qc_dir, sample))
        self.novel_dir = var.circ_novel_dir(PROJ_DIR)
        circ_mkdir(self.novel_dir)
        self.db_stat_dir = var.circ_db_stat_dir(self.novel_dir)
        circ_mkdir(self.db_stat_dir)

    def run_bwa(self):
        fq1 = []
        fq2 = []
        if self.fq_clean:
            for fq in self.fq_clean.split(','):
                _fq1 = fq + '_1.clean.fq.gz'
                _fq2 = fq + '_2.clean.fq.gz'
                fq1.append(_fq1)
                fq2.append(_fq2)
        else:
            # assuming that the QC step is done and has dir qc_dir/'cleandata'
            fq_clean_dir = var.fq_clean_dir(self.qc_dir)
            if os.path.isdir(fq_clean_dir):
                for s in self.samples:
                    _fq1 = glob.glob (jp (fq_clean_dir, s + '_1.clean.fq*'))[0]
                    _fq2 = glob.glob (jp (fq_clean_dir, s + '_2.clean.fq*'))[0]
                    fq1.append(_fq1)
                    fq2.append(_fq2)
            else:
                raise ValueError('no such dir (fq-clean-dir): %s' % fq_clean_dir)

        fns = []
        qsub_all = []
        for i in range(0, len(self.samples)):
            s = self.samples[i]
            out_sam = var.bwa_out_sam(self.novel_dir, s)
            out_log = var.bwa_out_log(self.novel_dir, s)
            cmd_bwa = '''
%s \\
    %s \\
    %s \\
    %s \\
    %s \\
    1>%s \\
    2>%s''' % (var.prog_bwa,
               var.prog_bwa_args,
               self.genome,
               fq1[i],
               fq2[i],
               out_sam,
               out_log
               )
            fn = var.bwa_run_fn(self.novel_dir, s)
            fns.append(fn)
            write_str_to_file(cmd_bwa, fn)
            qsub_all.append(qsub_cmd(self.mem, self.cpu, fn))
        qsub_fn= var.qsub_bwa_sh_fn(self.novel_dir)
        write_obj_to_file(qsub_all, qsub_fn)
        return qsub_fn, fns

    def run_ciri(self):
        fns = []
        qsub_all = []
        for s in self.samples:
            cmd_ciri = '''
set -eu
echo "[$(date)] running CIRI"
perl {prog_ciri} \\
     -P \\
     -I {sam} \\
     -O {ciri_out} \\
     -F {genome} \\
     -A {gtf}
echo "[$(date)]    done CIRI"

echo "[$(date)] running CIRI-AS"
perl {prog_ciri_as} \\
     -S {sam} \\
     -C {ciri_out} \\
     -O {ciri_as_out} \\
     -F {genome} \\
     -A {gtf}
echo "[$(date)]    done CIRI-AS"
'''.format (prog_ciri=var.prog_ciri,
            sam=var.bwa_out_sam (self.novel_dir, s),
            ciri_out=var.ciri_out_fn (self.novel_dir, s),
            genome=self.genome,        # NOTE: needs bwa index
            gtf=self.gtf,
            prog_ciri_as=var.prog_ciri_as,
            ciri_as_out=var.ciri_out_fn (self.novel_dir, s),
            )
            fn = var.ciri_run_fn(self.novel_dir, s)
            fns.append(fn)
            write_str_to_file(cmd_ciri, fn)
            qsub_all.append(qsub_cmd('1G', 1, fn))
        qsub_fn = var.qsub_ciri_sh_fn(self.novel_dir)
        write_obj_to_file(qsub_all, qsub_fn)
        return qsub_fn, fns

# def filter_ciri_result(self, db):

def circ_db_stat(PROJ_DIR, samples, db, OVERLAP=2, program='find_circ', db_type='circbase'):
    '''filter'''
    # samples = self.samples
    novel_dir = var.circ_novel_dir(PROJ_DIR)
    db_stat_dir = var.circ_db_stat_dir(novel_dir)
    db_stat_bar = var.circ_db_stat_bar_fn_data(db_stat_dir)
    db_stat_bar_pdf = var.circ_db_stat_bar_fn_pdf(db_stat_dir)
    db_stat_pie = var.circ_db_stat_pie_fn_data(db_stat_dir)
    db_stat_pie_pdf = var.circ_db_stat_pie_fn_pdf(db_stat_dir)
    circ_mkdir(db_stat_dir)
    if db_type == 'circbase':
        indexes_db = (0, 1, 2, 4)
    else:
        raise ValueError('invalid db_type: %s' % db_type)
    if program == 'find_circ':
        indexes_sample = (1, 2, 3, 0, 9, 10)   # _chr, start, end, _id, jread, njread
    elif program == 'ciri':
        raise ValueError('ciri is not supported')
        # indexes_sample = (1, 2, 3, 4)     # _chr, start, end, _id
    elif program == 'circexplorer':
        raise ValueError('circexplorer is not supported')
        indexes_sample = (0, 1, 2, 12)
    else:
        raise ValueError('invalid program: %s' % db_type)

    data_db = {}
    db_circ_id = {}
    fh = open(db, 'r')
    out_bar = open(db_stat_bar, 'w')
    fh.readline()
    for line in fh:
        _chr, start, end, _id = [line.split()[i] for i in indexes_db]
        data_db.setdefault(_chr, []).append([int(start), int(end), _id])
        db_circ_id[_id] = 1
    fh.close()

    # intersection_total = 0
    common_circ_id_all = {}
    for s in samples:
        if program == 'ciri':
            circ_out = var.ciri_out_fn (novel_dir, s)
        elif program == 'circexplorer':
            circ_out = var.circ_novel_fn (novel_dir, s)
        elif program == 'find_circ':
            circ_out = var.find_circ_out_fn (novel_dir, s)

        sample_db = {}
        fh = open(circ_out, 'r')
        fh.readline()
        for line in fh:
            _chr, start, end, _id, jread, njread = [line.split()[i] for i in indexes_sample]
            if jread >= 2:
                sample_db.setdefault(_chr, []).append([s, int(start), int(end), _id, jread, njread])
        fh.close()
        out_sample = var.circ_db_stat_fn(db_stat_dir, s)
        intersection, common_circ_id = get_circ_intersection(sample_db, data_db, OVERLAP, out_sample)
        # intersection_total += intersection
        for item in common_circ_id:
            common_circ_id_all[item] = 1
        out_bar.write('%s\t%s\n' % (s, intersection))
    # print 'total: %s' % intersection_total
    out_bar.close()
    id_count_sample = len(common_circ_id_all.keys())
    id_count_db  = len(db_circ_id.keys())
    # circ_ratio_in_db = float(id_count_sample) / float(id_count_db)
    # print 'ratio in db: %s / %s = %s' % (id_count_sample, id_count_db, circ_ratio_in_db)
    write_str_to_file('%s\t%s\n' % (id_count_sample, id_count_db), db_stat_pie)
    db_stat_plot_bar(db_stat_bar, db_stat_bar_pdf)
    db_stat_plot_pie(db_stat_pie, db_stat_pie_pdf)


def db_stat_plot_pie(data, out):
    shfile = '%s.sh' % data
    shcode = '''
Rscript /TJPROJ1/RNA/kangmingming/dev/circrna/script/lib/circ_db_stat_pie.R \\
    --fn %s \\
    --out %s
''' % (data, out)
    write_str_to_file(shcode, shfile)
    circ_run_script(shfile, 'sh')


def db_stat_plot_bar(data, out):
    shfile = '%s.sh' % data
    shcode = '''
Rscript /TJPROJ1/RNA/kangmingming/dev/circrna/script/lib/circ_db_stat_bar.R \\
    --fn %s \\
    --out %s
''' % (data, out)
    write_str_to_file(shcode, shfile)
    circ_run_script(shfile, 'sh')


def get_circ_intersection(sample, db, OVERLAP, out):
    # count_sample = 0
    common_circ_id = []
    count = 0
    fh = open(out, 'w')
    fh.write('chr\tdb_circRNA_ID\tsample\tsample_circRNA_ID\tjunction_read_count\tnon_junction_read_count\tstart_db\tend_db\tstart_sample\tend_sample\toverlap_start(nt)\toverlap_end(nt)\n')
    for _chr in sample.keys():
        if db.has_key(_chr):
            ref = sample[_chr]
            for ref1 in ref:
                s, start, end, _id, jread, njread = ref1
                for ref2 in db[_chr]:
                    start_db, end_db, id_db = ref2
                    distance_start = start - start_db
                    distance_end = end - end_db
                    if (abs(distance_start) <= OVERLAP) and (abs(distance_end) <= OVERLAP):
                        count += 1
                        common_circ_id.append(_id)
                        fh.write('\t'.join (map (lambda x: str(x), [_chr, id_db, s, _id, jread, njread, start_db, end_db, start, end, distance_start, distance_end])) + '\n')

        else:
            continue
    fh.close()
    return count, common_circ_id


class novel_circ_by_find_circ:
    ''''''
    def __init__(self, PROJ_DIR, samples, genome, gtf, abbr, fq_clean=None, cpu=8):
        "identify circRNAs by find_circ.py "
        self.PROJ_DIR = PROJ_DIR
        self.samples= samples
        self.genome = genome
        self.gtf = gtf
        self.abbr = abbr
        self.fq_clean = fq_clean
        self.cpu = cpu
        self.qc_dir = var.circ_qc_dir (PROJ_DIR)
        self.novel_dir = var.circ_novel_dir (PROJ_DIR)
        self.sorted_gtf_bed = var.sorted_gtf_bed (self.novel_dir)

    def generate_fns(self):
        fq1 = []
        fq2 = []
        fns = []
        if self.fq_clean:
            for fq in self.fq_clean.split(','):
                _fq1 = glob.glob (fq + '_1.*fq.gz')[0]
                _fq2 = glob.glob (fq + '_2.*fq.gz')[0]
                fq1.append(_fq1)
                fq2.append(_fq2)
        else:
            # assuming that the QC step is done and has dir qc_dir/'cleandata'
            fq_clean_dir = var.fq_clean_dir(self.qc_dir)
            if os.path.isdir(fq_clean_dir):
                for s in self.samples:
                    _fq1 = glob.glob (jp (fq_clean_dir, s + '_1.clean.fq*'))[0]
                    if (not os.path.isfile (_fq1)):
                        raise IOError ('no such file: %s' % _fq1)
                    _fq2 = glob.glob (jp (fq_clean_dir, s + '_2.clean.fq*'))[0]
                    if (not os.path.isfile (_fq2)):
                        raise IOError ('no such file: %s' % _fq2)
                    fq1.append(_fq1)
                    fq2.append(_fq2)
            else:
                raise IOError ('no such dir (fq-clean-dir): %s' % fq_clean_dir)
        shcode = '''
set -eu

# split genome to single chr
genome={genome}
split_genome_dir={split_genome_dir}
mkdir -p $split_genome_dir
python /TJPROJ1/RNA/kangmingming/dev/circrna/script/lib/split-genome.py \\
    --genome $genome \\
    --outdir $split_genome_dir

wrapper_find_circ={wrapper_find_circ}
local_copy={local_copy}
[ -e $local_copy ] && mv -f $local_copy $local_copy.bak.$(date +"%Y-%m-%d_%T")
cp $wrapper_find_circ $local_copy

# merge exon and intron gtf and sort
wd={novel_dir}
gtf={gtf}
gtf_exon=$wd/exon.gtf
gtf_intron=$wd/intron.gtf
gtf_merged=$wd/merged.gtf
sorted_gtf_bed={sorted_gtf_bed}

echo "[$(date)] extracting exon and intron"
# cat $gtf | perl -wane 'if (/gene_type "(.+?)";/) {{ next if $1 ne "mRNA"; }} ($F[2] eq "exon") and print' >$gtf_exon
cat $gtf | awk '$3 == "exon"' >$gtf_exon
perl /TJPROJ1/RNA/kangmingming/dev/circrna/script/lib/get_intron_gtf.pl $gtf >$gtf_intron
cat $gtf_exon $gtf_intron >$gtf_merged

echo "[$(date)] running gtf2bed"
__gtf2bed=/TJPROJ1/RNA/kangmingming/dev/circrna/bin/bedops_linux_x86_64-v2.4.14/convert2bed
mypath=/TJPROJ1/RNA/kangmingming/dev/circrna/bin/bedops_linux_x86_64-v2.4.14
export PATH=$mypath:$PATH       # convert2bed will use the sort-bed binary
$__gtf2bed --input=gtf <$gtf_merged \\
    | perl -F'\\t' -wlape '$F[3]=~s/^"//; $_ = join "\\t", @F;' \\
           > $sorted_gtf_bed
echo "[$(date)] Done"
'''.format(split_genome_dir=var.split_genome_dir (self.novel_dir),
           genome=self.genome,
           wrapper_find_circ=var.wrapper_find_circ,
           local_copy=var.wrapper_find_circ_local_copy (self.novel_dir),
           novel_dir=self.novel_dir,
           gtf=self.gtf,
           sorted_gtf_bed = self.sorted_gtf_bed,
           )
        shfile = var.prepare_find_circ_sh (self.novel_dir)
        write_str_to_file (shcode, shfile)
        split_genome_fn = shfile
        for i in range(0, len(self.samples)):
            sample = self.samples[i]
            out = var.find_circ_out_dir (self.novel_dir, sample)
            circ_mkdir_unix(out)
            shfile = var.run_find_circ_pipeline_sh (self.novel_dir, sample)
            shcode = '''
wrapper_find_circ={wrapper_find_circ}
local_copy={local_copy}
[ -e $local_copy ] || local_copy=$wrapper_find_circ

sh $local_copy \\
     --fq1 {fq1} \\
     --fq2 {fq2} \\
     --sample {sample} \\
     --out {out} \\
     --genome {genome} \\
     --gtf {gtf} \\
     --abbr {abbr} \\
     --sorted-gtf-bed {sorted_gtf_bed} \\
     --split-genome-dir {split_genome_dir} \\
     --cpu {cpu}
'''.format(wrapper_find_circ=var.wrapper_find_circ,
           local_copy=var.wrapper_find_circ_local_copy (self.novel_dir),
           split_genome_dir=var.split_genome_dir (self.novel_dir),
           fq1=fq1[i],
           fq2=fq2[i],
           sample=sample,
           out=out,
           genome=self.genome,
           gtf=self.gtf,
           abbr = self.abbr,
           novel_dir=self.novel_dir,
           sorted_gtf_bed = self.sorted_gtf_bed,
           cpu=self.cpu,
           )
            write_str_to_file(shcode, shfile)
            fns.append(shfile)
        return split_genome_fn, fns

    def get_intron_inside_circ(self):
        """
        get intron(s) within circRNAs.
        """
        pass


def get_seq_len_max (fasta):
    maxlen = 0
    for record in SeqIO.parse (fasta, 'fasta'):
        seqlen = len (str (record.seq))
        if seqlen > maxlen:
            maxlen = seqlen
    return maxlen


def get_circ_seq (io_novel_circ, genome, out, mirna_len_max):
    sys.stderr.write ('generating spliced circRNA sequence... ')
    # FIXME: parse_fasta() slurps lots of RAM for all samples
    genome_dict = parse_fasta(genome)
    # print genome_dict.keys()[0:2]
    fh = open (io_novel_circ)
    fh.readline ()
    oh = open (out, 'w')
    CIRC_SEQ_LEN_FILTER_MAX = 100000
    # CIRC_JUNCTION_READ_MIN = 10
    for line in fh:
        _id, _chr, start, end, strand, feature = [line.split ('\t')[i] for i in (0, 1, 2, 3, 4, 8)]
        if abs (int (start) - int (end)) > CIRC_SEQ_LEN_FILTER_MAX:
            continue
        seq = ''
        if feature == '--':
            seq = genome_dict[_chr][int(start)-1:int(end)]
        else:
            features = filter (lambda x: x != '', feature.split (';'))
            sites = []
            for feat in features:
                each = map (lambda x: re.sub ('.+?:', '', x).split ('-'), filter (lambda x: x != '', feat.split (',')))
                for start_end in each:
                    sites.append (start_end)
            sites = sorted (sites, key=lambda x: x[0])
            # print (sites)
            # print ('*** [start, end] = ', [start, end])
            # add additional sequence where circRNA starts before feature start
            sites[0][0] = start
            # if int (sites[0][0]) >= int (start):
            #     sites[0][0] = start
            # else:
            #     raise ValueError ('feature start (%s) is less that the full length start site (%s)' % (sites[0][0], start))
            # add additional sequence where circRNA ends after feature end
            sites[len (sites) - 1][1] = end
            # if int (sites[len (sites) - 1][1]) <= int (end):
            #     sites[len (sites) - 1][1] = end
            # else:
            #     raise ValueError ('feature end (%s) is greater that the full length end site (%s)' % (sites[len (sites) - 1][1], end))
            for start, end in sites:
                seq_fragment = genome_dict[_chr][int(start)-1:int(end)]
                seq += seq_fragment
        if strand is '-':
            seq = str (Seq (seq).reverse_complement ())

        # get junction-flanking sequence
        seqlen = len (seq)
        OVERLAP_LEN_MIN = 1
        jseq_3p_start = seqlen - mirna_len_max + OVERLAP_LEN_MIN # inclusive
        jseq_5p_end = mirna_len_max - OVERLAP_LEN_MIN      # exclusive
        junction_seq = seq[jseq_3p_start:seqlen] + seq[0:jseq_5p_end]
        # NOTE: To view the length distribution of known miRNAs, use: grep -v
        # '^>' /TJPROJ1/RNA/kangmingming/dev/circrna/data/mature-v21.fa | perl
        # -wne 'BEGIN {$h={};} chomp; $len = length; $h->{$len} ++; END {print
        # "$_\t$h->{$_}\n" for sort {$a cmp $b} keys %$h}'
        oh.write ('>%s\n%s\n' % (_id, seq))
        oh.write ('>%s_junction_seq\n%s\n' % (_id, junction_seq))
    fh.close ()
    oh.close ()
    sys.stderr.write ('done\n')


def annotate_circrna (PROJ_DIR, samples):
    novel_dir = var.circ_novel_dir (PROJ_DIR)
    for s in samples:
        out = var.find_circ_out_dir (novel_dir, sample)
        novel = jp (out, 'novel_circ_%s.txt' % s)
        spliced_circ_info = jp (out, '%s_spliced_circRNA_info.txt' % s)


def io_novel_circ_table (PROJ_DIR, samples, abbr, prog='find_circ'):
    h = {}
    novel_dir = var.circ_novel_dir (PROJ_DIR)
    for s in samples:
        switch = {
            'circexplorer' : [var.circ_novel_fn (novel_dir, s), ()],
            # NOTE: CIRI_AS output contains only partial circRNAs that are alternative spliced
            'ciri' : [var.ciri_out_fn (novel_dir, s), (1, 2, 3, 4, 10, 9, 5, 7, 12, 13)],
            'find_circ' : [var.annotated_circ (novel_dir, s), range (1, 11)],
            }
        fn, index = switch[prog][0:2]
        # cp (fn, jp (des['3.1'], re.sub('.txt$', '.xls', os.path.basename(fn_switch[prog]))))
        fh = open (fn, 'r')
        fh.readline ()
        if prog is 'find_circ':
            for line in fh:
            # _chr, start, end, strand, count = [line.split('\t')[i] for i in (0, 1, 2, 5, 12)]
            # array = re.split ('\t', line)[1:]
                array = \
                  _chr, start, end, strand, flen, slen, gene, feat, jread, njread = \
                  [re.split ('\t', line.strip ())[i] for i in index]
                if jread < 2:
                    continue
                id1 = '_'.join ([_chr, start, end, strand])
                length = abs (int (start) - int (end))
                h.setdefault (id1, {})['info'] = array[0:8]
                h.setdefault (id1, {}).setdefault ('sample', []).append (s)
                h.setdefault (id1, {}).setdefault ('jread', []).append (jread)
                h.setdefault (id1, {}).setdefault ('njread', []).append (njread)
        elif prog is 'ciri':
            g = {}
            for line in fh:
                # sys.stderr.write (line)
                if len (line.split ('\t')) != 18:
                    continue
                array = \
                  _chr, start, end, strand, gene, feat, jread, njread, ce_start, ce_end = \
                  [re.split ('\t', line.strip ())[i] for i in index]
                id1 = '_'.join ([_chr, start, end, strand])
                g.setdefault (id1, {})['part1'] = array[0:8]
                g.setdefault (id1, {}).setdefault ('part2', []).append ([ce_start, ce_end])
            for id1 in g.keys ():
                part1 = g[id1]['part1']
                part2 = get_non_overlap_region (g[id1]['part2'])
                # test:
                # atest1 = sorted (g[id1]['part2'], key=lambda x: x[0])
                # if len (atest1) != len (part2):
                #     print atest1, '----\n', part2, '\n'
                flen = abs (int (part1[2]) - int (part1[1]))
                slen = 0
                features = []
                for ce_start, ce_end in part2:
                    slen += abs (int (ce_start) - int (ce_end))
                    features.append ('%s:%s-%s' % (feat, ce_start, ce_end))
                if slen > flen:
                    slen = flen
                features = ','.join (features)
                info = part1[0:4] + [str (flen), str (slen)] + [part1[4], features]
                jread = part1[6]
                njread = part1[7]
                if njread is "":
                    njread = "0"
                h.setdefault (id1, {})['info'] = info
                h.setdefault (id1, {}).setdefault ('sample', []).append (s)
                h.setdefault (id1, {}).setdefault ('jread', []).append (jread)
                h.setdefault (id1, {}).setdefault ('njread', []).append (njread)
        fh.close ()
    i = 1
    header = 'circRNA_id\tchr\tstart\tend\tstrand\tfull_length\tspliced_length\tgene_id\tfeature\tsamples\tjunction_read\tnon_junction_read\n'
    fh = open (var.io_novel_circ (novel_dir), 'w')
    fh.write (header)
    for id1 in sorted (h.keys ()):
        uid = '%s_circ_%07d' % (abbr, i)
        info = '\t'.join (h[id1]['info'])
        sample = ','.join (h[id1]['sample'])
        jread = ','.join (h[id1]['jread'])
        njread = ','.join (h[id1]['njread'])
        fh.write('%s\t%s\t%s\t%s\t%s\n' % (uid, info, sample, jread, njread))
        i += 1
    fh.close ()


def get_non_overlap_region (aref):
    '''aref: 2D-array of [[start1, end1], [start2, end2], ..., [startN, endN]]'''
    # href = { i : aref[i] for i in range (0, len (aref) - 1) }
    # assuming that all items in aref are not overlapped by each other
    aref = sorted (aref, key=lambda x: x[0])
    new = [aref.pop (0)]
    if len (aref) == 0:
        return new
    else:
        for i in range (0, len (aref)):
            start, end = aref.pop (0)
            atmp = []
            for l in new:
                atmp.extend (l)
            atmp.append (start)
            idx = sorted (atmp).index (start)
            idx_elt = idx % 2
            if idx_elt == 0:                # idx_start_next outside a [start,end]
                if idx == len (atmp) - 1:
                    new.insert (idx / 2, [start, end])
                    # itype = 'type1'
                    continue
                else:
                    # itype = 'type2'
                    new_start_next = atmp[idx + 1]
                    new_end_next = atmp[idx + 2]
                    # new_start_next2 = new_start_next + 2
                    if end < new_start_next:           # this is a new feature to be inserted
                        new.insert (idx / 2, [start, end])
                    elif end <= new_end_next:        # this feature overlaps with [start,end]
                        new[idx / 2] = [start, new_end_next]
                    # elif (end > new_end_next) and (end < new_start_next2):
                    elif (end > new_end_next):
                        new[idx / 2] = [start, end]
            elif idx_elt == 1:              # inside a [start,end]
                new_start_this = atmp[idx - 1]
                new_end_this = atmp[idx + 1]
                # new_start_next = new[idx + 1]
                if end <= new_end_this:
                    pass
                # elif (end > new_end_this) and (end < new_start_next):
                elif (end > new_end_this):
                    # print 'idx =', idx
                    # print '[start,end]=', [start, end]
                    # print '[new_start,new_end]=', [new_start_this, new_end_this]
                    new[idx / 2] = [new_start_this, end]
                # elif end >= new_start_next:

            # print 'new[]:', new
    return new

    # |-----|-----|-----|
    #    A     B     C
    #
    # There are 3 + 2 + 1 = 6 combinations when the 'second' region compares
    # with the first pre-extracted region, represented as [start,end] below:
    # [A,A], [A,B], [A,C]
    # [B,B], [B,C]
    # [C,C]
    #
    # for j in new.keys ():
    #     n1, n2 = new[j]
    #     for i in href.keys ():
    #         start, end = href[i]
    #         if start >= n1 and end <= n2:   # href[i] within new[j]
    #             href.pop (i)
    #             continue
    #         elif start < n1 and end > n2:   # new[j] within href[i]
    #             new[j] = [start, end]
    #             href.pop (i)
    #             continue
    #         elif (start < n1 and end < n1) or (start > n2 and end > n2): # before or after
    #             new[j] =



# eval: (highlight-regexp "{[A-Za-z0-9_]+}" 'hi-green-b)

# Local Variables:
# eval: (highlight-regexp "var\\." 'hi-red-b)
# End:
