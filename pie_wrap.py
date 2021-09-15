#!/usr/bin/env python3

import argparse
from ete3 import Tree
from glob import glob
import os

def getRootName():
  species_tree = Tree(args.species_tree)
  names = [a.name for a in species_tree.children[0].get_leaves()]
  print(names)
  return names

def handleTrees(rootTaxa, rootSp):
  new_path = os.path.join(args.tree_directory, "..", "PhyParts")
  new_tree_directory = os.path.join(new_path, "testd")

  if not os.path.isdir(new_path):
    os.mkdir(new_path)
  if not os.path.isdir(new_tree_directory):
    os.mkdir(new_tree_directory)

  raw_trees = []

  if "dna.tree" in args.species_tree:
    raw_trees = glob(os.path.join(args.tree_directory, "*NT*.treefile"))
  else:
    raw_trees = glob(os.path.join(args.tree_directory, "*AA*.treefile"))

  num_trees = len(raw_trees)
  for treename in raw_trees:
    basename = treename.split("/")[-1]
    tree = Tree(treename)
    try:
      if len(rootTaxa) == 1:
        tree.set_outgroup(rootTaxa[0])
        tree.write(format=9, outfile=os.path.join(new_tree_directory, basename + "_REROOTED"))
      else:
        #re-root to taxon not in rootTaxa, then root to rootTaxa
        #this keeps things from breaking if the common ancestor of rootTaxa is the root node of the tree
        for taxon in tree.get_leaves():
          taxon = taxon.name
          if taxon not in rootTaxa:
            tree.set_outgroup(taxon)
            break

        common_taxa = [leaf.name for leaf in tree.get_leaves() if leaf.name in rootTaxa]
        common_node = tree.get_common_ancestor(common_taxa)
        tree.set_outgroup(common_node)
        tree.write(format=9, outfile=os.path.join(new_tree_directory, basename + "_REROOTED"))
    except:
      print("{} does not have the outgroup".format(treename))
      num_trees -= 1

  #new_spTree = ""
  #if do_astral:
  spBasename = args.species_tree.split("/")[-1]
  spTree = Tree(args.species_tree)

  if rootSp:
    if len(rootTaxa) == 1:
      spTree.set_outgroup(rootTaxa[0])
    else:
      common_node = spTree.get_common_ancestor(rootTaxa)
      spTree.set_outgroup(common_node)

  new_spTree = os.path.join(new_path, spBasename + "_REROOTED")
  spTree.write(format=9, outfile=new_spTree)

  return new_spTree, num_trees

def main():
  rootTaxa = []
  if not args.outgroup:
    rootTaxa = getRootName()
  else:
    rootTaxa = args.outgroup

  #print(root)
  new_spTree = ""
  num_trees = 0

  if not args.outgroup:
    new_spTree, num_trees = handleTrees(rootTaxa, False)
  else:
    new_spTree, num_trees = handleTrees(rootTaxa, True)

  #print(new_spTree)

  new_path = os.path.join(args.tree_directory, "..", "PhyParts")
  out_path = os.path.join(new_path, "out")
  script_directory = "/".join(os.path.realpath(__file__).split("/")[:-1])
#  print(script_directory)

  phyparts_location = os.path.join(script_directory, "dependencies", "phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar")
  phyparts_command = "java -jar {} -a 1 -v -d {} -m {} -o {}".format(phyparts_location, os.path.join(new_path, "testd", ""), new_spTree, out_path)
  os.system(phyparts_command)
  #print(phyparts_command)
  phypartspiecharts_location = os.path.join(script_directory, "dependencies", "phypartspiecharts.py")
  phypartspiecharts_command = "python {} {} {} {} --svg_name {}".format(phypartspiecharts_location, new_spTree, out_path, num_trees, os.path.join(new_path, "pies.svg"))
  os.system(phypartspiecharts_command)
  #print(phypartspiecharts_command)

  print("Trees used in PhyParts analysis:", num_trees)

parser = argparse.ArgumentParser(description="Helper script for PhyParts and PhyPartsPieCharts")
parser.add_argument('-t', '--tree_directory', help = 'The path to the directory with the single locus gene trees')
parser.add_argument('-p', '--species_tree', help = 'The tree you consider to be the species tree')
parser.add_argument('-g', '--outgroup', nargs="+", help = 'The sample you consider the outgroup for rooting purposes. If left out, trees will automatically be rooted to have the same root as your species tree.')
#parser.add_argument('-p', '--phyparts_path', help = 'The path to the phyparts .jar')
#parser.add_argument('-c', '--phypartspiecharts_path', help = 'The path to phypartspiecharts.py')
args = parser.parse_args()

if __name__=='__main__':
  main()
