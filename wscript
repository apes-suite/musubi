#!/usr/bin/env python3
# encoding: utf-8
# 2023 Harald Klimach <harald.klimach@dlr.de>

APPNAME = 'musubi'

top = '.'
out = 'build'

def options(opt):
    '''Building options provided by Musubi.
       Remember, all options can be displayed with waf --help.'''
    opt.recurse('bin')
    opt.recurse('aotus')
    opt.recurse('tem')
    opt.recurse('mus')


def configure(conf):
    '''Project configuration'''
    import os
    conf.recurse('aotus', 'subconf')
    conf.recurse('bin', 'preconfigure')
    # Use default.coco as coco settings file by default
    if not conf.options.coco_set:
        conf.options.coco_set = 'default.coco'
    conf.recurse('tem')
    conf.recurse('mus')
    conf.recurse('bin', 'postconfigure')


def build(bld):
    '''Build the Musubi project'''
    from revision_module import fill_revision_string
    bld.recurse('bin')
    fill_revision_string(bld, subdir='mus')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('aotus')
        bld.recurse('tem')
    else:
        bld.load('coco')
        bld(rule='cp ${SRC} ${TGT}',
            source = bld.path.find_node(['tem', 'source', 'arrayMacros.inc']),
            target = bld.path.find_or_declare('arrayMacros.inc'))
    bld.recurse('mus')

#clean build directory and coco completely to create the build from scratch
def cleanall(ctx):
    from waflib import Options
    Options.commands = ['distclean'] + Options.commands
    ctx.exec_command('rm coco')
