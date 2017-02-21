#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Date:   December 7, 2016 9:05:15 AM
# @Filename: renderer.py
# @Last modified by:   dileep
# @Last modified time: February 2, 2017 9:39:08 AM

import jinja2
import os
import sys
from param_parser import yaml_parser

def render_module(module_name, module_params):
    os.chdir('modules/')
    module_loader = jinja2.FileSystemLoader(os.getcwd())
    jenv = jinja2.Environment(loader=module_loader)
    module_template = jenv.get_template(module_name+'.j2')
    # params = dict()
    # params[module_name] = module_params
    rendered_module = module_template.render(module_params)
    os.chdir('..')
    return rendered_module

def get_header(record_dict):
    rendered_file = ''
    rendered_file += render_module('header', record_dict['header'])
    rendered_file += '\n\n'
    return rendered_file

def get_modules_needed(rendered_file, record_dict):
    module_list = [dirname for dirname in os.listdir('modules') if '.j2' in dirname]
    module_list = [mod.split('.')[0] for mod in module_list]
    modules_needed = record_dict['Modules']
    for module in modules_needed:
        if module in module_list:
            print(modules_needed[0])
            rendered_file += render_module(module, record_dict[module])
            rendered_file += '\n\n'
        else:
            print("Requested module {0} not available".format(module))
            sys.exit("Error")
    return rendered_file

def update_channels(rendered_file, record_dict):
    modules = record_dict["Modules"]
    channel_set = {"reads", "genome", "annotation"}
    channel_dict = {key:[] for key in channel_set}
    for module in modules:
        meta_file = module + '.yml'
        meta_data = yaml_parser(meta_file)
        if set(meta_data["input"]["default"]) < channel_set:
            for in_channel in meta_data["input"]["default"]:
                channel_dict[in_channel].append(module)
        else:
            for missing_channel in set(meta_data["input"]["default"]) - channel_set:
                channel_dict[meta_data["input"]["alternate"][missing_channel]].append(module)
    return channel_dict

def main(yaml_file):
    records = yaml_parser(yaml_file)
    rendered_file = get_header(records)
    rendered_file = get_modules_needed(rendered_file, records)
    output_file = 'pipeline.nf'
    with open(output_file, 'w') as fid:
        fid.write(rendered_file)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: renderer.py <param.yml>")
    # tools = ['quality', 'trimming', 'indexing', 'aligning', 'quantification', 'report']
    main(sys.argv[1])
