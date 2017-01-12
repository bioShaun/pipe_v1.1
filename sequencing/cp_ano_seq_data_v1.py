"""Usage:
  cp_seq_data_v1.py [--cp_info_file=<file_path>] [--log_dir=<logs_directory>]
  cp_seq_data_v1.py -h | --help | --version
"""


from docopt import docopt
import sys
from os import path
from os import system
from os import listdir
from python_tools import write_obj_to_file

def shell_cmd_for_cp(from_dir, to_dir, log_dir):
    from_dir_name = path.basename(from_dir)
    cmd_list = ['echo start cp {0}\ndate'.format(from_dir_name)]
    to_dir = path.join(to_dir, from_dir_name)
    if not path.exists(to_dir):
        system('mkdir -p "{0}"'.format(to_dir))
    for each_file in listdir(from_dir):
        each_file_path = path.join(from_dir, each_file)
        if path.isfile(each_file_path):
            if each_file not in listdir(to_dir):
                ## cp md5, fq stat files
                cmd_list.append('echo "start cp {0}"'.format(each_file))
                cmd_list.append('cp "{0}" "{1}"'.format(each_file_path, to_dir))
                cmd_list.append('echo "finished cp {0}"'.format(each_file))
        else:
            ## cp clean and raw data
            fq_dir_list = listdir(each_file_path)
            to_fq_dir = path.join(to_dir, each_file)
            if each_file not in listdir(to_dir):
                system('mkdir -p "{0}"'.format(to_fq_dir))
            for each_fq_dir in fq_dir_list:
                if each_fq_dir not in listdir(to_fq_dir):
                    each_fq_dir_path = path.join(each_file_path, each_fq_dir)
                    cmd_list.append('echo "start cp {0}/{1}"'.format(each_file, each_fq_dir))
                    cmd_list.append('cp -r "{0}" "{1}"'.format(each_fq_dir_path, to_fq_dir))
                    cmd_list.append('echo "finished cp {0}/{1}"'.format(each_file, each_fq_dir))
    cmd_list.append('tree -s "{0}" > "{1}/{2}.ori"'.format(from_dir, log_dir, from_dir_name))
    cmd_list.append('tree -s "{0}" > "{1}/{2}.bak"'.format(to_dir, log_dir, from_dir_name))
    cmd_list.append('diff "{0}/{1}.ori" "{0}/{1}.bak"'.format(log_dir, from_dir_name))
    cmd_list.append('echo "finished cp {0}"\ndate'.format(from_dir_name))
    return cmd_list

def nohup_sh_job(sh_script, log_dir):
    shell_name = path.basename(sh_script)
    log_file = path.join(log_dir, '{0}.log'.format(shell_name))
    system('nohup sh {0} > {1} 2>&1 &'.format(sh_script, log_file))

if __name__ == '__main__':
    arguments = docopt(__doc__, sys.argv[1:], version='v1')
    cp_info_file = arguments['--cp_info_file']
    log_dir = arguments['--log_dir']
    if not cp_info_file or not log_dir:
        print(__doc__)
        sys.exit(1)

    ## get cp info dict
    cp_info_dict = {}
    with open(cp_info_file) as cp_info_file_inf:
        for eachline in cp_info_file_inf:
            eachline_inf = eachline.strip().split('\t')
            cp_info_dict[eachline_inf[0]] = eachline_inf[1]

    ## creat cp shell script for each directory and lauch it
    cp_info_file_name = path.basename(cp_info_file)
    shell_script = '{0}.sh'.format(cp_info_file_name)
    cmd_list = []
    for each in cp_info_dict:
        from_dir = each
        to_dir = cp_info_dict[each]
        cmd_list.extend(shell_cmd_for_cp(from_dir, to_dir, log_dir))
    write_obj_to_file(cmd_list, shell_script)
    nohup_sh_job(shell_script, log_dir)

