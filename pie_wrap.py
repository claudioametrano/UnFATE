import argparse
from ete3 import Tree
from glob import glob
import os

def getRootName():
  species_tree = Tree(args.species_tree)
  root_node = species_tree.get_tree_root()
  root_sample = ""
  if root_node.children[0].name != "":
    root_sample = root_node.children[0].name
  elif root_node.children[1].name != "":
    root_sample = root_node.children[1].name
  else:
    print("Couldn't find a root, using args.outgroup.")
  return root_sample

def handleTrees(root_sample, do_astral):
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
      tree.set_outgroup(root_sample)
      tree.write(format=9, outfile=os.path.join(new_tree_directory, basename + "_REROOTED"))
    except:
      print("{} does not have the outgroup".format(treename))
      num_trees -= 1

  #new_spTree = ""
  #if do_astral:
  spBasename = args.species_tree.split("/")[-1]
  spTree = Tree(args.species_tree)

  if do_astral:
    spTree.set_outgroup(root_sample)

  new_spTree = os.path.join(new_path, spBasename + "_REROOTED")
  spTree.write(format=9, outfile=new_spTree)

  return new_spTree, num_trees

def main():

  root = getRootName()
  #print(root)
  new_spTree = ""
  num_trees = 0
  if not args.outgroup:
    new_spTree, num_trees = handleTrees(root, False)
  elif args.outgroup == root:
    new_spTree, num_trees = handleTrees(args.outgroup, False)
  elif args.outgroup != root:
    new_spTree, num_trees = handleTrees(args.outgroup, True)

  #print(new_spTree)

  new_path = os.path.join(args.tree_directory, "..", "PhyParts")
  out_path = os.path.join(new_path, "out")
  script_directory = "/".join(os.path.realpath(__file__).split("/")[:-1])
#  print(script_directory)

  phyparts_location = os.path.join(script_directory, "phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar")
  phyparts_command = "java -jar {} -a 1 -v -d {} -m {} -o {}".format(phyparts_location, os.path.join(new_path, "testd", ""), new_spTree, out_path)
  os.system(phyparts_command)
  #print(phyparts_command)
  phypartspiecharts_location = os.path.join(script_directory, "phypartspiecharts.py")
  phypartspiecharts_command = "python {} {} {} {} --svg_name {}".format(phypartspiecharts_location, new_spTree, out_path, num_trees, os.path.join(new_path, "pies.svg"))
  os.system(phypartspiecharts_command)
  #print(phypartspiecharts_command)

  print("Trees used in PhyParts analysis:", num_trees)

parser = argparse.ArgumentParser(description="Helper script for PhyParts and PhyPartsPieCharts")
parser.add_argument('-t', '--tree_directory', help = 'The path to the directory with the single locus gene trees')
parser.add_argument('-s', '--species_tree', help = 'The tree you consider to be the species tree')
parser.add_argument('-g', '--outgroup', help = 'If left out, trees will automatically be rooted to have the same root as your species tree (if it is rooted). The sample you consider the outgroup for rooting purposes')
#parser.add_argument('-p', '--phyparts_path', help = 'The path to the phyparts .jar')
#parser.add_argument('-c', '--phypartspiecharts_path', help = 'The path to phypartspiecharts.py')
args = parser.parse_args()

if __name__=='__main__':
  main()
