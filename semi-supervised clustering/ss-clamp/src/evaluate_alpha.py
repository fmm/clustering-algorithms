#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import subprocess

prog = "./a.out"

flag = { 
  '[INFO]' : [],
  '[DATABASE]' : [],
  '[CLASS_VARIABLE]' : [],
  '[INITIALIZATION]' : [],
  '[MAXIMUM_ITERATION]' : [],
  '[CLUSTERS]' : [],
  '[PROTOTYPES]' : [],
  '[RELEVANCE_V]' : [],
  '[ALPHA]' : [],
  '[TIME_LIMIT]' : [],
  '[INPUT]' : []
}

def load_config(path):
  f = open(path,"r")
  key = "[]"
  for line in f:
    token = line.rstrip('\n')
    if len(token) > 0:
      if token[0] == '[':
        key = token
        flag[key] = []
      else:
        flag[key].append(token)
  f.close()

def generate_config(path):
  f = open(path,'w')
  for key in flag:
    f.write("{0}\n".format(key))
    for value in flag[key]:
      f.write("{0}\n".format(value))
  f.close()

def get_pwc_name(path,label,repeat):
  return path + "-" + "l" + str(label) + "r" + str(repeat) + ".pwc"

# run once before the main execution
def prepare_restrictions(path):
  for label in xrange(10,100+1,10):
    for repeat in xrange(1,30+1):
      npath = get_pwc_name(path,label,repeat)
      cmd = []
      cmd.append(prog)
      cmd.append("config.txt")
      cmd.append("-g {0}".format(label))
      cmd.append("-o {0}".format(npath))
      print cmd
      subprocess.call(cmd)

# TODO
path = "input/iris/iris"

if False:
  load_config(path + ".txt")
  generate_config("config.txt")
  prepare_restrictions(path)
  exit(0)

for label in xrange(10,100+1,10):
  for repeat in xrange(1,30+1):
    load_config(path + ".txt")
    npath = get_pwc_name(path,label,repeat)
    flag["[INPUT]"].append(npath)
    generate_config("config.txt")
    subprocess.call(["cat","config.txt"])
    cmd = []
    cmd.append(prog)
    cmd.append("config.txt")
    cmd.append("-m prl")
    subprocess.call(cmd)
