import sys
import os
from os.path import join as jp
from os.path import basename
from os.path import dirname
from time import localtime, strftime
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import json
import shutil

'''
Some functions for common needs

'''
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
        fh.write('%s\n' % obj)
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
    msg = '%s [onmath] %s\n' % (circ_get_time(), msg)
    if add_blank_line_before is True:
        msg = '\n%s' % msg
    if add_blank_line_after is True:
        msg = '%s\n' % msg
    # print msg
    print >> sys.stderr, msg
    if fn:
        write_str_to_file(msg, fn, append=True)


def log_call_proc(msg, fn=None):
    msg = '%s [onmath] [calling process] %s\n' % (circ_get_time(), msg)
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



def table_to_dict(fn,key_row = 1,val_row = 2,header = True,sep='\t'):
    fn_dict = {}
    with open(fn,'r') as fn_info :
        for n,eachline in enumerate(fn_info) :
            each_info = eachline.strip().split(sep)
            if n == 0 and header:
                continue
            else :
                fn_dict[each_info[key_row-1]] = each_info[val_row-1]
    return fn_dict

def Median(a):
    a = sorted(a)
    l = len(a)
    l2 = l / 2
    return a[l2] if l % 2 else (a[l2] + a[l2 - 1]) / 2.0

def add_dict_value(mydict,mykey,myvalue):
    if mykey not in mydict :
        mydict[mykey] = [myvalue]
    elif myvalue not in mydict[mykey] :
        mydict[mykey].append(myvalue)
    else :
        pass

def invert_dict(d):  
    return dict((v,k) for k,v in d.iteritems()) 

def multi_process_shell_script(sh_list,fn,multi_num = 12) :
    with open(fn,'w') as fn_info :
        fn_info.write('#! /bin/bash\n')
        fn_info.write('echo "####Process begins####"\n')
        fn_info.write('date\n')
        for n,each in enumerate(sh_list) :
            if (n+1) % multi_num == 0 :
                fn_info.write('wait\n')
            else :
                fn_info.write('sh %s &\n' % each)
        fn_info.write('date\n')
        fn_info.write('echo "####Process ends####"')


def file_to_list(fn):
    outlist = []
    with open(fn,'r') as fn_info :
        for eachline in fn_info :
            eachline = eachline.strip()
            outlist.append(eachline)
    return outlist

def th_launch_job_cmd(script) :
    cmd = 'yhrun -n 1 -N 1 -c 24 -x cn[9531] -p work %s &' % script
    return cmd

def write_obj_to_json(obj,fn):
    with open(fn,'w') as fn_info :
        json.dump(obj,fn_info)

def load_fn_to_obj(fn) :    
    with open(fn) as fn_info :
        obj = json.load(fn_info)
    return obj

def merge_files(filenames, output_file):
    outfh = open(output_file, "w")
    for filename in filenames:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()

def main():
    pass

if __name__ == '__main__' :
    main()
