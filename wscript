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
    # Initialize the coco preprocessing tool
    conf.load('coco')
    conf.env['COCOSET'] = 'default.coco'
    if not conf.options.coco_reports:
      # Make coco silent, if not explicitly asked for reports:
      if conf.env.COCOFLAGS:
        conf.env.COCOFLAGS.insert(0, '-s')
        conf.env.COCOFLAGS.append('-ad')
      else:
        conf.env.COCOFLAGS = ['-s', '-ad']
    conf.recurse('tem')
    conf.recurse('mus')
    conf.recurse('bin', 'postconfigure')


def build(bld):
    '''Build the Musubi project'''
    from revision_module import fill_revision_string
    bld.recurse('bin')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('aotus')
    fill_revision_string(bld, subdir='mus')
    bld(rule='cp ${SRC} ${TGT}', source=bld.env.COCOSET, target='coco.set')
    bld.recurse('tem')
    bld.recurse('mus')
