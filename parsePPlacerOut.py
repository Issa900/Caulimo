#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Parsing of pplacer output files
    input : .jplace | for the correspondance between taxa and nodes number
            .csv | to count the number of sequence placement for each alignment as query

'''
import sys, os, re, string, time
from optparse import OptionParser

class node(object):
	def __init__(self):

		self.name =""
		self.path =()
		self.feuille = False
		self.parent =""
		self.species=[]
		self.de = float(0)
		self.deFromRoot = float(0)
		self.childs=[]

	def __repr__(self):
		# return ("node({}), type of node({}), path to root({}), number of sequence({}), evolutionary distance({})".format(str(self.node), self.feuille, ";".join(self.path), str(self.lght), str(self.distEvol))
		# return "toto"
		lght= str(len(self.species))
		distEvol = str(self.deFromRoot)
		if self.feuille :
			feuille = "Leaf***"
		else:
			feuille = "internal node"
		path = ";".join(self.path)
		return ("%s\n" %(",".join([self.name,feuille,path,lght,distEvol])))

class parsePPlacerOut(object):
	def __init__(self, refTree ="", inCSV = ""):
		self.refTree = refTree
		self.inCSV = inCSV
		self.nodes={}
		self.count = 0

	def checkOptions(self):
		if not os.path.exists(self.refTree) :
			raise Exception ("ERROR: file not found %s" %(self.refTree))
		if not os.path.exists(self.inCSV) :
			raise Exception ("ERROR: file not found %s" %(self.inCSV))

	def setAttributesFromCmdLine(self):
		description = "Read jplace file and return for each node the corresponding number"
		parser = OptionParser(description = description)
		parser.add_option("-r", "--refTree", dest = "refTree", action = "store",   type = "string", help = ".jplace format file")
		parser.add_option("-c", "--csv", dest = "csvFile", action = "store",   type = "string", help = ".csv format file")
		options, args = parser.parse_args()
		self.refTree = options.refTree
		self.inCSV = options.csvFile

	def sT(self, sub):
		stringToParse="".join(sub)
		position = stringToParse.find('(')
		subT =""
		subTL=""
		if (position==-1):
			subT=stringToParse
		else :
			subTL=stringToParse[:position]
			subT=stringToParse[position+1:]
		return subTL, subT

	def write_gvz (self):

		with open ("%s.gvz" %(self.refTree),'w') as graphfile :
			graphfile.write('digraph G {\n\trankdir=LR;splines=polyline ;\n')
			for name in self.nodes :
				if self.nodes[name].feuille :
					graphfile.write("\""+ name + "\" [shape = polygon, sides=9, label=\"" + name +": " + self.nodes[name].name +" "+ str(self.nodes[name].de) + "\", color=\"darkseagreen1\",  style=filled ]; \n")
				else:
					graphfile.write("\""+ name + "\" [shape = \"none\", label=\"" + name +": " + self.nodes[name].name +" "+ str(self.nodes[name].de) + "\", color=\"darkslategrey\"]; \n")
			for name in self.nodes :
				# print ("node", name, ":" , self.nodes[name].species, " path is ", self.nodes[name].path)
				if self.nodes[name].parent > -1 :
					graphfile.write("%s -> %s ;\n"%(self.nodes[name].parent,name))
			graphfile.write('}')

	def path_tab(self):
		for name in self.nodes :
			self.nodes[name].de

	def saveNode(self,nb,info,parent,token=""):
			i=node()
			i.name = info.rsplit(':',1)[0]
			if parent >-1 and info :
				i.de = float(info.rsplit(':',1)[-1])
				i.deFromRoot = self.nodes[parent].deFromRoot + i.de
			i.parent = parent
			i.feuille = not(token == ')')
			if parent >-1 :
				i.path = self.nodes[parent].path + [parent]
				self.nodes[parent].childs=self.nodes[parent].childs + [nb]
			else :
				i.path =[]

			self.nodes[nb]=i



	def read_tree (self, na, sub):

		self.count +=1
		#print ("TOTO Iteration", self.count, 'noeud', na, 'suite', sub)
		node =""
		info=""
		st=""
		sB=""
		nblevel = 0
		accolade_r= sub.rfind('}')
		accolade_l= sub.rfind('{')
		info_l= max(sub[:accolade_l].rfind(')'),sub[:accolade_l].rfind(','),sub[:accolade_l].rfind('('),sub[:accolade_l].rfind('}'))

		if accolade_r>-1:
			#un noeud existe extraire noeud et info
			node= sub[accolade_l+1:accolade_r]
			#info
			info=sub[info_l+1:accolade_l]
			token = sub[info_l]
			self.saveNode(node,info,na,token)

		if info_l>-1:
			# une suite existe connaitre quel est le premier token trouve
			#s'il existe des freres
			if token ==")" : #trouver la fermante correspondante pour faire le sous-arbre
				nblevel=1
				newToken = token
				newST=sub[:info_l]
				sT_l = sub.find("(")
				while nblevel> 0 :
					#print ("TOTO ",nblevel," - ", sT_l,"newST", newST)
					sT_l= max(newST.rfind(')'),newST.rfind('('))
					if sub[sT_l] == ')': # on augmente le niveau de profondeur du sous arbre
						nblevel+=1
					else:
						nblevel-=1
					newtoken = sub[sT_l]
					newST=sub[:sT_l]
				sT=sub[sT_l+1:info_l]
				afterST=sub[:sT_l] # ce qui suit n'appartient pas au sous arbre donc a un meme ancetre
				#print ("TOTO sT",sT, "afterST", afterST)
				self.read_tree (node, sT)
				self.read_tree (na, afterST)
			elif token =="," : #ce qui suit est un frere donc a le meme ancestre
				sB_l= sub.rfind(",")
				sB=sub[sB_l+1:info_l]
				afterSB = sub[:sB_l] # ce qui suit est un sous arbre du grand-parent
				#print ("TOTO sB",sB, "afterSB", afterSB)
				self.read_tree (na, sB)
				self.read_tree (na, afterSB)

	def read_csv (self):
		with open (self.inCSV, 'r') as csv :
			next(csv)
			for line in csv:
				line = line.strip().split(",")
				node_position = str(line[3])
				self.nodes[node_position].species.append(line[1])

	def write_outf (self):

		with open ("%s.outf.csv" %(self.inCSV),'w') as f:
			f.write ("%s\n" %(",".join(["node","type of node","path to root","number of sequence","evolutionary distance"])))
			for node in self.nodes :

				lght= str(len(self.nodes[node].species))
				distEvol = str(self.nodes[node].deFromRoot)
				if self.nodes[node].feuille :
				    feuille = "Leaf***"
				else:
				    feuille = "internal node"
				path = ";".join(self.nodes[node].path)
				f.write ("%s\n" %(",".join([node,feuille,path,lght,distEvol])))

	def run (self):
		self.checkOptions()
		line = open (self.refTree, 'r').readlines()[1]
		line = line.rstrip().strip(';",')
		self.read_tree ( -1, line)
		self.read_csv()
		self.write_outf()
		self.write_gvz()

if __name__ == "__main__":
	i=parsePPlacerOut()
	i.setAttributesFromCmdLine()
	deb= time.time()
	i.run()
	print ("execution", time.time()-deb)
