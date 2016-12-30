#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Date:   December 7, 2016 9:33:35 AM
# @Filename: param_parser.py
# @Last modified by:   dileep
# @Last modified time: December 7, 2016 9:56:27 AM

import yaml

def yaml_parser(yaml_file):
    with open(yaml_file, 'r') as fid:
        record = yaml.load(fid)
    return record
