#! /usr/bin/python

import sys
import os
import python_tools

if len(sys.argv) != 3:
    print 'python '+sys.argv[0]+' input_file output_dir'
    sys.exit(0)

input_file = sys.argv[1]
name = "lxgui"
port = "158"
output_dir = sys.argv[2]

port_name = '%s@10.128.76.%s' % (name,port)

cmd = 'scp -r "%s" %s:%s' % (input_file,port_name,output_dir)
python_tools.circ_call_process(cmd)
