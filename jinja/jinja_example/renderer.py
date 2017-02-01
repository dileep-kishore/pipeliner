#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Date:   November 30, 2016 2:02:24 PM
# @Filename: renderer.py
# @Last modified by:   dileep
# @Last modified time: December 7, 2016 9:36:32 AM

import jinja2
import os

if __name__ == '__main__':
    loader = jinja2.FileSystemLoader(os.getcwd())
    jenv = jinja2.Environment(loader=loader)
    template = jenv.get_template('jinja.txt')
    param = dict()
    param["name"] = 'Jinja Test'
    param["subtitle"] = 20161130
    param["my_list"] = ['first', 'second', 'third']
    rendered_filetxt = template.render(param)
    print(rendered_filetxt)
    with open('rendered_file.txt', 'w') as fid:
        fid.write(rendered_filetxt)
