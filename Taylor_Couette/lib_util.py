#!/usr/bin/env python3
#
# lib_util.py
#  Yannick Ponty
#
#
########################################################   
###########   Library needed  ##########################

import glob
import os
from collections import OrderedDict
import multiprocessing

#########################################################################
def all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)
####################################################################################################
class MyException(Exception):
    pass
####################################################################################################
def OptionValue(pattern,**kwargs):
    try: 
        value=kwargs[pattern]
        return value
    except KeyError:
        raise MyException(pattern+' Option is not defined !')
###################################################################
def OptionValueD(pattern, default,warn=True, **kwargs):
    try: 
        value=kwargs[pattern]
        return value
    except:
        if warn :
            print('Warning '+pattern+' option is not defined ! default='+str(default))
        return default
#########################################################################
###############################################################################
def map_char(dir, pattern_filter, bool_filter, ext='.cub'):  
    olddir= os.getcwd()
    if dir : os.chdir(dir)
    basedir0= '*_*'+ext
    listfiles = glob.glob(basedir0)
    Char_output= [s.split('_')[0] for s in listfiles]                              
    map_char_output=list(OrderedDict.fromkeys(Char_output))
    if '' in map_char_output :
        map_char_output.remove('')
    map_char_output.sort()
    if (pattern_filter) :
        if (bool_filter): iterator=filter(lambda c : c.find(pattern_filter)> -1, map_char_output)
        if (not bool_filter): iterator=filter(lambda c : c.find(pattern_filter)==-1, map_char_output)
        map_char_output_filtered=list(iterator)
    else :
         map_char_output_filtered=map_char_output
    #print('Field Char Outputs=',map_char_output_filtered)
    os.chdir(olddir)
    return map_char_output_filtered
#################################################################################
def map_step(dir, stepinf, stepsup, ext='.cub'):
    olddir= os.getcwd()
    if dir : os.chdir(dir)
    basedir= '*_*'+ext
    start = stepinf   #start step
    overstep= stepsup #start step
    temp = glob.glob(basedir)
    temp_r = [int((f.split('_')[1]).split(ext)[0]) for f in temp] 
    temp_r2=list(OrderedDict.fromkeys(temp_r))
    temp_r2.sort()                                  
    map_step_output = []
    for t in temp_r2:
        if t>start  and t<overstep:
            map_step_output.append(t)
    #print('Step Outputs =',map_step_output)
    os.chdir(olddir)
    return  map_step_output
#################################################################################
def map_files(dir, nostep, ext='.cub'):
    olddir= os.getcwd()
    if dir : os.chdir(dir)
    stepstring='*_'
    if nostep : stepstring=''
    basedir0= stepstring+'*'+ext
    temp_all = glob.glob(basedir0)
    Char_output= [s.split(ext)[0] for s in temp_all]
    Char_output.sort()
    os.chdir(olddir)
    return  Char_output
#################################################################################
def cubby_varname_step(input_name):
    name_list=input_name.split('_')
    step=int(name_list[1])
    varname=name_list[0]
    return [varname,step]
########################################################################################
#   Static counter 
######################################################################################
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate
################################################################################################    
@static_vars(counter=-1)
def count():
    count.counter += 1
    return count.counter
###############################################################################################
#######   cpu number on the node. 
#################################################################################################
def cpu_number(options_ref):
    nprocs_max = multiprocessing.cpu_count()
    print('nprocs max=',nprocs_max)
    if options_ref.cpu:
        nprocs=min(nprocs_max,options.cpu)
    else :
        nprocs=nprocs_max
    if options_ref.verbose : print('nprocs=',nprocs)
    return nprocs
###################################################################################################