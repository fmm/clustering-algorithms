#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import subprocess, shlex
import time

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

def get_pwc_token(label,repeat):
  return "l" + str(label) + "r" + str(repeat)

def get_pwc_name(path,label,repeat):
  return path + "-" + get_pwc_token(label,repeat) + ".pwc"

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
      time.sleep(2)
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
    #extract seed from previous execution
    if False:
      cmd = "./get_seed.sh result_iris " + get_pwc_token(label,repeat)
      output = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE).communicate()
      seed = int(output[0])
      flag["[SEED]"] = [seed]
    generate_config("config.txt")
    subprocess.call(["cat","config.txt"])
    cmd = []
    cmd.append(prog)
    cmd.append("config.txt")
    cmd.append("-m prl")
    subprocess.call(cmd)
