import argparse
from ete3 import Tree
from glob import glob
import os

def handleTrees():
  new_path = os.path.join(args.tree_directory, "..", "PhyParts")
  new_tree_directory = os.path.join(new_path, "testd")

  if not os.path.isdir(new_path):
    os.mkdir(new_path)
  if not os.path.isdir(new_tree_directory):
    os.mkdir(new_tree_directory)

  raw_trees = glob(os.path.join(args.tree_directory, "*AA*.treefile"))
  num_trees = len(raw_trees)
  for treename in raw_trees:
    basename = treename.split("/")[-1]
    tree = Tree(treename)
    try:
      tree.set_outgroup(args.outgroup)
      tree.write(format=9, outfile=os.path.join(new_tree_directory, basename + "_REROOTED"))
    except:
      print("{} does not have the outgroup".format(treename))
      num_trees -= 1

  spBasename = args.species_tree.split("/")[-1]
  spTree = Tree(args.species_tree)
  spTree.set_outgroup(args.outgroup)
  new_spTree = os.path.join(new_path, spBasename + "_REROOTED")
  spTree.write(format=9, outfile=new_spTree)
  return new_spTree, num_trees

def main():


  new_spTree, num_trees = handleTrees()
  new_path = os.path.join(args.tree_directory, "..", "PhyParts")
  out_path = os.path.join(new_path, "out")
  phyparts_command = "java -jar {} -a 1 -v -d {} -m {} -o {}".format(args.phyparts_path, os.path.join(new_path, "testd", ""), new_spTree, out_path)
  os.system(phyparts_command)
  #print(phyparts_command)
  phypartspiecharts_command = "python {} {} {} {} --svg_name {}".format(args.phypartspiecharts_path, new_spTree, out_path, num_trees, os.path.join(new_path, "pies.svg"))
  os.system(phypartspiecharts_command)
  #print(phypartspiecharts_command)

  print(num_trees)

parser = argparse.ArgumentParser(description="Helper script for PhyParts and PhyPartsPieCharts")
parser.add_argument('-t', '--tree_directory', help = 'The path to the directory with the single locus gene trees')
parser.add_argument('-s', '--species_tree', help = 'The tree you consider to be the species tree')
parser.add_argument('-g', '--outgroup', help = 'The sample you consider the outgroup for rooting purposes')
parser.add_argument('-p', '--phyparts_path', help = 'The path to the phyparts .jar')
parser.add_argument('-c', '--phypartspiecharts_path', help = 'The path to phypartspiecharts.py')
args = parser.parse_args()

if __name__=='__main__':
  main()
