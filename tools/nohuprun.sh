#!/bin/bash
for args in $@
do
        cmd_file=$args
        log_file=${cmd_file/.sh/.log}".$(date +"%Y%m%d%H%M%S")"
        CurrentPath=$(pwd)
        chmod +x "$CurrentPath/$cmd_file"
        nohup sh "$CurrentPath/$cmd_file" > "$CurrentPath/$log_file" 2>&1&
done
