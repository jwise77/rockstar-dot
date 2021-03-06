import os, sys
import pydot
import numpy as np

if len(sys.argv) < 3:
    print ("usage: %s tree_file halo_id" % sys.argv[0])
    sys.exit()

min_mass = 1e4

halo_id = int(sys.argv[-1])
fname = sys.argv[-2]
if not os.path.exists(fname):
    print ("Cannot find %s" % fname)
    sys.exit()
lines = open(fname, "r").readlines()

# Search for tree
start_line = None
rank = -halo_id
num = 0
if halo_id <= 0:
    for l in lines:
        if l.startswith("#tree"):
            if num == -halo_id:
                halo_id = int(l.split()[1])
                break
            num += 1
    print ("Halo %d (by mass) has ID %d" % (rank, halo_id))

header = "#tree %d\n" % halo_id
for i,l in enumerate(lines):
    if l.startswith("#Omega_M"):
        h = float(l.split("=")[3])
    if l == header:
        start_line = i+1
    elif l.startswith("#tree") and start_line != None:
        end_line = i
        break

# Load data from strings
nleaves = end_line - start_line
nfields = len(lines[start_line].split())
data = np.empty((nleaves, nfields))
for i in range(start_line, end_line):
    data[i-start_line,:] = np.array(lines[i].split(), dtype='float')
del lines

# Create the pydot graph object
graph = pydot.Dot("galaxy", graph_type="digraph")
halo_shape = "rect"
z_shape = "plaintext"

# Get scale factors and make subgrids to align the halos in a given output
auniq = np.unique(data[:,0])
auniq.sort()
subgs = []
for a in auniq:
    subgs.append(pydot.Subgraph('', rank='same'))
    graph.add_subgraph(subgs[-1])

aeps = 1e-6
global_mmp = halo_id
nump = 0
for i in range(nleaves):
    mmp = int(data[i,14]) == 1
    hid = int(data[i,1])
    desc = int(data[i,3])
    # The actual mass is (mass [Msun/h]) = mass/h [Msun]
    mass = data[i,10]/h
    if mass < min_mass: continue
    color = "red" if mmp else "black"
    lvl = np.where(np.abs((data[i,0] - auniq)/data[i,0]) < aeps)[0][0]
    next_lvl = np.where(np.abs((data[i,2] - auniq)/data[i,0]) < aeps)[0]

    halo_str = "C%d_H%d" % (lvl, hid)
    # Don't do this for halos without a descendant (usually the last redshift)
    if next_lvl.size != 0: 
        prog_str = "C%d_H%d" % (next_lvl, desc)
        edge = pydot.Edge(halo_str, prog_str, color=color)
        graph.add_edge(edge)

    if desc == global_mmp and mmp:
        global_mmp = hid
    all_desc = np.where(data[:,3].astype('int64') == hid)[0]
    color = "grey"
    if all_desc.size > 0:
        major_merger = np.any(np.logical_and(data[all_desc,10] > 0.3*data[i,10],  # Mass ratio > 0.3
                                             data[all_desc,14] < 1))
        if major_merger:
            color = "lightblue"
    style = "filled" if hid == global_mmp else "solid"
    label = "Halo %d\\n%.3g" % (int(data[i,1]), mass)
    node = pydot.Node(halo_str, label=label, style=style, shape=halo_shape, color=color)
    graph.add_node(node)
    # Add this node to the correct level subgraph
    subgs[lvl].add_node(pydot.Node(halo_str))
    nump += 1

# Add subgraphs
print ("Plotting %d of %d halos" % (nump, nleaves))
for i in range(len(auniq)):
    # If there are no halos at this redshift, don't add the subgraph
    if len(subgs[i].get_node_list()) == 0: continue
    redshift = 1.0/auniq[i] - 1.0
    node = pydot.Node("%1.5e" % (redshift), shape=z_shape,
                      label = "z=%0.3f" % (redshift))
    subgs[i].add_node(node)

graph.write("mt-halo%d.pdf" % halo_id, format="pdf")
