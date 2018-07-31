#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 15:35:10 2018

@author: mathiasberger
"""

from json import dump, load
from pyomo.environ import Var, Param
from os import mkdir, listdir
from os.path import isdir

class PostProcess(object):
    
    def __init__(self, results_path, data, model,):
        
        self.results_path = results_path
        self.data = data
        self.model = model
        
    def save_results(self):
        if not isdir("./results"):
            mkdir("./results")

        var = [v for v in self.model.component_objects(ctype=Var, active=True, descend_into=True)]
        #param = [p for p in model.component_objects(ctype=Param, active=True, descend_into=True)]
        indexed_var = []
        for v in var:
            if not isinstance(v._index, set):
                indexed_var.append(v)
    
        for v in indexed_var:
            var.remove(v)
            dict_tmp = dict(zip(self.data.time, [v[t].value for t in self.data.time]))
            with open("./results/"+str(v)+".json", "w") as f:
                dump(dict_tmp, f)
            f.close()
            
    def load_results(self):
            files = [f for f in listdir(self.results_path) if f.endswith('.json')]
            res = []
            for f in files:
                res.append(read_json(self.results_path+"str"))
            return res
                
                
def read_json(json_object_path):
        if isinstance(json_object_path, str):
            with open(json_object_path, "r") as in_file:
                dict_object = load(in_file)
                keys = list(dict_object.keys())
                keys_mapping, dict_out = dict(zip([keys[i] for i in range(len(keys))], keys)), {}
                sorted_keys = sorted(list(keys_mapping.keys()))
                for i in range(len(keys_mapping.keys())):
                    dict_out.update({sorted_keys[i]: dict_object[keys_mapping[sorted_keys[i]]]})
            in_file.close()
            return dict_out
        else:
            raise TypeError("Input parameter must be of type str, type %s passed." % type(json_object_path))

#def plot_json(json_object_path, name_str="", show=True, out=False):
#    read_json = Plotman.read_json
#    dict_to_plot = read_json(json_object_path)
#    time_index, nodes = sorted(dict_to_plot.keys()), len(dict_to_plot[dict_to_plot.keys()[0]])
#    fig = plt.figure(figsize=(16, 9))
#    for n in range(1, nodes + 1):
#        if name_str == "":
#            plt.plot(time_index, [dict_to_plot[i][n - 1] for i in time_index])
#        else:
#            plt.plot(time_index, [dict_to_plot[i][n - 1] for i in time_index],
#                     label="%s at node %d" % (name_str, n))
#    if name_str != "":
#        plt.xlabel("Simulation Time [s]")
#        plt.ylabel("%s" % str)
#        plt.title("%s evolution over time" % name_str)
#    if show:
#        fig.show()
#    if out:
#        return fig
