# Powered by Python 2.7

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import tlp

# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined 
# to run the script on the current graph

viewLabel = graph.getStringProperty("viewLabel")
viewColor = graph.getColorProperty("viewColor")
viewSize = graph.getSizeProperty("viewSize")
Locus = graph.getStringProperty("Locus")
viewNegative = graph.getBooleanProperty("Negative")
viewPositive = graph.getBooleanProperty("Positive") 
viewLayout = graph.getLayoutProperty("viewLayout")
viewMetric = graph.getDoubleProperty("viewMetric")
viewSize = graph.getSizeProperty("viewSize")
blue = tlp.Color(0,0,255)
red = tlp.Color(255,0,0)
green = tlp.Color(0,255,0)

 
def labels(graph):
  for n in graph.getNodes():
    viewLabel[n] = Locus[n]
    viewColor[n] = blue
    viewSize[n]=(2.0,2.0)

def dessinerModeleForce(graph, layout):
	parametre = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
#	parametre["Unit edge length"] = 10
#	parametre["Edge Length Property"] = viewMetric
#	parametre["Node Size"] = viewSize 
	graph.applyLayoutAlgorithm("FM^3 (OGDF)", layout, parametre)


def edges(graph):
	for e in graph.getEdges():
		#print viewNegative[e], viewPositive[e]
		if viewNegative[e] == True and viewPositive[e] == False:
			viewColor[e] = red
		if viewNegative[e] == False and viewPositive[e] == True:
			viewColor[e] = green

	
def constructTree(tree, current_cluster, root):
	# root = racine de mon arbre
	# cluster_courant = Gene Interaction ? 
	# node = noeuds a ajouter dans l'arbre 
	# ajouter un sommet n dans l'arbre
	# relier le root au sommet n 
	# appeler la fonction sur chaque sous cluster
	
		if current_cluster.numberOfSubGraphs() == 0: #feuilles
			for n in current_cluster.getNodes(): #retourne une liste de noeuds, il faut donc la parcourir
				tree.addNode(n)
				tree.addEdge(n, root)
			return 
			
		for subgraph in current_cluster.getSubGraphs(): #branches
			node = tree.addNode()
			tree.addEdge(node, root) 
			constructTree(tree, subgraph, node)

def radial(graph): 
	params = tlp.getDefaultPluginParameters('Tree Radial', graph)
	resultLayout = graph.getLayoutProperty('resultLayout')
	success = graph.applyLayoutAlgorithm('Tree Radial', resultLayout, params)

def main(graph): 
   labels(graph)
   dessinerModeleForce(graph, viewLayout)
   edges(graph)
   tree = graph.addSubGraph("Hierarchic Tree")
   root = tree.addNode()
   current_cluster = graph.getSubGraph("Genes interactions")
   constructTree(tree, current_cluster, root)
   #radial(tree)
   

for n in graph.getNodes():
  print n
