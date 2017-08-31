#!/usr/bin/env python
# -*- coding: utf-8 -*-



import sys, os, re, string, time
from parsePPlacerOut import parsePPlacerOut
from optparse import OptionParser

class formatPPlacer (object):
    def __init__(self, ctrf="", refTree="",  inCSV = ""):
        self.ctrf = ctrf
        self.refTree = refTree
        self.inCSV = inCSV
        self.cf_dico = {}
        self.outpars = {}


    def define_argfromcmdline(self):
        x = 'Parse pplacer output files'
        parser = OptionParser(description = x)
        # parser.add_option ('-i', '--inparser', dest= 'inparser', action = "store",   type = "string", help = "outfileof the first parser")
        parser.add_option ('-f', '--ctrf', dest= 'ctrf', action = "store", default='control file', type = "string", help = "control file to defile by user")
        parser.add_option("-r", "--refTree", dest = "refTree", action = "store",   type = "string", help = ".jplace format file")
        parser.add_option("-c", "--csv", dest = "csvFile", action = "store",   type = "string", help = ".csv format file")

        options, args = parser.parse_args()
        self.ctrf = options.ctrf
        # self.inparser = options.inparser
    	self.refTree = options.refTree
        self.inCSV = options.csvFile
        self.i = parsePPlacerOut(self.refTree, self.inCSV)
        self.i.run()

    def read_controlF (self):

        try:
            with open (self.ctrf, 'r') as cf:
                for line in cf:
                    tmp = line.strip().split(':')
                    grp_name = str(tmp[0]).strip()
                    node_numbers = tmp[-1].strip().split(';')
                    # self.cf_dico[grp_name]=(node_number)
                    for j in node_numbers:
                        j=j.strip()
                        if not grp_name in self.cf_dico : self.cf_dico[grp_name] = []
                        self.cf_dico[grp_name]= self.cf_dico[grp_name]+ [j] + self.i.nodes[j].childs

        except IOError:
            pass



    def write_controlF(self):

        base= os.path.splitext(os.path.basename(self.inCSV))[0]
        code_sp = base.split('_')[0]
        with open ("%s.CTF.csv" %(code_sp),'w') as f:
            f.write ("%s\t%s\t%s\n" %("Cluster", "Specie", "number"))
            for grp_name in self.cf_dico :
                l_nodes = set(self.cf_dico[grp_name])
                total=0
                for node in l_nodes :
                    number = len(self.i.nodes[node].species)
                    total += number
                f.write("%s\t%s\t%s\n" %(grp_name, code_sp, total))


    # def read_input (self):
    #     for j in self.i.node


    def run (self):
        self.read_controlF
        self.i

if __name__ == "__main__":
    p = formatPPlacer()
    p.define_argfromcmdline()
    t = time.time()
    p.read_controlF()
    p.write_controlF()
    # p.run()
    print("the execution time", time.time() - t)
