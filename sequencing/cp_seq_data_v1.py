"""Usage:
  cp_seq_data_v1.py [--cp_info_file=<file_path>] [--log_dir=<logs_directory>]
  cp_seq_data_v1.py -h | --help | --version
"""


from docopt import docopt
import sys
from os import path
from os import system
from python_tools import write_obj_to_file

def shell_cmd_for_cp(from_dir, to_dir, log_dir):
    cmd_list = ['echo start\ndate']
    from_dir_name = path.basename(from_dir)
    to_dir = path.join(to_dir, from_dir_name)
    cmd_list.append('cp -r "{0}" "{1}"'.format(from_dir, to_dir))
    cmd_list.append('tree -s "{0}" > "{1}/{2}.ori"'.format(from_dir, log_dir, from_dir_name))
    cmd_list.append('tree -s "{0}" > "{1}/{2}.bak"'.format(to_dir, log_dir, from_dir_name))
    cmd_list.append('diff "{0}/{1}.ori" "{0}/{1}.bak"'.format(log_dir, from_dir_name))
    cmd_list.append('echo finished\ndate')
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
    for each in cp_info_dict:
        from_dir = each
        from_dir_name = path.basename(each)
        shell_script = '{0}.sh'.format(from_dir_name)
        to_dir = cp_info_dict[each]
        cmd_list = shell_cmd_for_cp(from_dir, to_dir, log_dir)
        write_obj_to_file(cmd_list, shell_script)
        nohup_sh_job(shell_script, log_dir)

