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




###############################################################
###############################################################

##### Script Python by COTTAIS Deborah and SCHRIEKE Hans ######
################ Master in Bioinformatics #####################
###################### 2018-2019 ##############################
### Tulip project : Visual analysis of gene expression data ###

###############################################################
###############################################################



from tulip import *

blue = tlp.Color(0,0,255)
red = tlp.Color(255,0,0)
green = tlp.Color(0,255,0)
orchid=tlp.Color(218,112,214)
fuchsia=tlp.Color(255,0,255)



###############################################################
###################### Part I #################################
########## Pre-treatment & first visualization ################
###############################################################

def colorLabels(graph, color, size, label, locus):

	## This function allows to put colors, size and the name of genes (locus) for all the nodes in the graph "Gene interactions".
	## It takes as arguments only the graph and the properties : viewColor, viewSize, viewLabel and Locus.

	for n in graph.getNodes():
		color[n] = blue
		size[n]=(2.0,2.0)
		label[n]=locus[n]


def drawForceModel(graph, layout):

	##  This function allows to visualize the graph "Gene interactions" as a graph (whith nodes and edges) and not as a circle.
	##  It takes as arguments the graph of interest and the property viewLayout of this graph.

	parameter = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
	graph.applyLayoutAlgorithm("FM^3 (OGDF)", layout, parameter)


def colorEdges(graph, negative, positive, color):

	## This function allows to create a scale of colors about the values of their negative or positive statutes for
	## each edge of the "Gene interactions" graph.
	## It takes as arguments the graph and the properties viewNegative, viewPositive and viewColor.

	for e in graph.getEdges():
		#print viewNegative[e], viewPositive[e]
		if negative[e] == True and positive[e] == False:
			color[e] = red
		if negative[e] == False and positive[e] == True:
			color[e] = green
		if negative[e] == True and positive[e] == True:
			color[e] = orchid
		if negative[e] == False and positive[e] == False:
			color[e] = fuchsia



###############################################################
###################### Part II ################################
######### Networked drawing of interactions network ###########
###############################################################

def constructTree(tree, cluster, root):

	## This function allows to build a subgraph (hierarchic tree) from the first graph (gene interactions) which represents
	## the interaction network partitioning. This new tree had as top the differents clusters and its leaves represented
	## all the genes of cluster.
	## This function takes as arguments the new graph ("hierarchic tree"), the first graph ("Gene interactions")
	## and the root of hierarchic tree (which allows to add a new node).

	if cluster.numberOfSubGraphs() != 0:
		for subgraph in cluster.getSubGraphs(): #branches
			node = tree.addNode()
			tree.addEdge(root, node)
			constructTree(tree, subgraph, node)
	if cluster.numberOfSubGraphs() == 0: #leaves
		for n in cluster.getNodes():
			tree.addNode(n)
			tree.addEdge(root, n)


def drawRadial(graph, layout):

	## This function allows to draw radially the graph hierarchic tree thanks to the nodes coordinates.
	## It takes as arguments the graph of interest ("hierarchic tree") and the layout.

	params = tlp.getDefaultPluginParameters('Tree Radial', graph)
	success = graph.applyLayoutAlgorithm('Tree Radial', layout, params)

def graphProperty(graph, prop, color):

	## This function allows to color the tops of hierachic tree according to their double property (the differents tp_*).
	## We decided to take 5 as threshold because the values of tp_* are between 0 and 15.
	## It takes as arguments the hierarchic tree, the proprerties tp_* and viewColor. We can chose in the main which property that we want.

	for n in graph.getNodes():
		if prop[n] <= 5:
			color[n]=red
		if prop[n] >6:
			color[n]=green


def rootPathNodes(graph, nodeList, node):

	## This function allows to calculate a path between a node and the root of the graph
	## So it takes as arguments the hierarchic tree, a list which will contain the path nodes, and the node of interest.

	if graph.getInNodes(node)!=0:
		nodeList.append(node)
		for n in graph.getInNodes(node):
			rootPathNodes(graph, nodeList, n)


def shortestPath(tree, nodeSource, nodeTarget):

	## This function allows to calculate the shortest path between two nodes
	## This path is a list of the nodes between the first and the second selected node
	## This function uses the rootPathNodes function to create the list of source nodes which store all the nodes between the source node and the root
	## and  a list of target nodes which store all nodes between the target node and the root
	## From these lists, we can determine the common nodes shared by source and target paths ; we take only the first, it's the common ancestor
	## We then compare the source list and the target list with the common ancestor in order to store the source nodes and the target nodes before the commun ancestor in different lists
	## Finally, we store the list of source nodes, the commun ancestor and the list of target nodes (reversed) in a new list which is the shortest path
	## The function returns the shortest path list

	sourceList = []
	targetList = []
	rootPathNodes(tree, sourceList, nodeSource)
	rootPathNodes(tree, targetList, nodeTarget)

	commonAncestor = []
	process=True
	for i in range(len(sourceList)):
		for j in range(len(targetList)):
			if process == True :
				if sourceList[i] == targetList[j]:
					commonAncestor.append(sourceList[i])
					process = False


	newNodeList = []
	newTargetList = []
	pathNodesList = []
	for s in range(len(sourceList)):
		if sourceList[s] == commonAncestor[0]:
			break
		else:
			newNodeList.append(sourceList[s])
	for t in range(len(targetList)):
		if targetList[t] == commonAncestor[0]:
			break
		else:
			newTargetList.append(targetList[t])

	newTargetList.reverse()
	pathNodesList = pathNodesList+newNodeList
	pathNodesList.append(commonAncestor[0])
	pathNodesList = pathNodesList+newTargetList
	del pathNodesList[0],pathNodesList[-1]

	return pathNodesList


def constructBundles(graph, tree, layout, treeLayout, shape):

	## This function assings control points to the edges of graph for draw bundles on the tree
	## First, it uses the shortestPath function to calculate the shortest path between each source nodes and target nodes of each edge of the graph
	## Then, each node coordinates of this shortestPath list is store in a list of coordinates
	## This coordinates are assings to the edge value and applied on the graph layout property
	## Finally we applies the cubic b-spline filter to set the shape of edges, it allows to fit the bundles

	for edge in graph.getEdges():
		coordControlPoints = []
		eachShortestPath = shortestPath(tree, graph.source(edge), graph.target(edge))

		for node in eachShortestPath:
			coordControlPoints.append(treeLayout[node])
		layout.setEdgeValue(edge, coordControlPoints)

	for edge in graph.getEdges():
		shape[edge] = tlp.EdgeShape.CubicBSplineCurve




###############################################################
###################### Part III ###############################
############### Small multiples construction  #################
###############################################################

def constructSmallMultiples(tree, cluster, tps, color, layout):

	## This function allows to create the 17 small sub multiples which represents
	## the different levels of gene expression (the 17 clusters). It store them
	## in the subgraph small multiples. We recover the tp_* property and store it
	## in view metric property in order to build the different 17 clusters.
	## A color mapping is used on the purpose of differentiating and visualizating
	## gene expression levels.
	## This function takes as arguments "hierarchic tree", "Gene interactions", the values
	## contained in tp_* property, viewColor and  viewLayout.

	smallMultiples=graph.addSubGraph("Small Multiples")
	for tpNumber in range(1,len(tps)+1):
		smallSubMultiples=smallMultiples.addSubGraph("tp"+str(tpNumber))
		tlp.copyToGraph(smallSubMultiples,cluster)
		tpValue=smallSubMultiples.getDoubleProperty("tp"+str(tpNumber)+" s")
		for n in smallSubMultiples.getNodes():
			tpNodeValue=tpValue.getNodeDoubleValue(n)
			valueMetric={"viewMetric":tpNodeValue}
			smallSubMultiples.setNodePropertiesValues(n,valueMetric)
			colorSmallSubMultiples(smallSubMultiples, color)
		makeGrid(smallSubMultiples, smallMultiples, tpNumber, 5, layout)


def colorSmallSubMultiples(graph, color):

	## This function allows to create a color mapping. In fact, each sub small multiples
	## created (thanks to the function from above) corresponds to a different level
	## of gene expression. In order to visualize this expression, we need a color mapping.
	## This function takes so as arguments the sub small multiples and viewColor.

	colorScale = tlp.ColorScale([tlp.Color.Red,tlp.Color.Red,tlp.Color(192,0,0),tlp.Color(192,0,0),tlp.Color(128,0,0),tlp.Color(128,0,0),tlp.Color(64,0,0),tlp.Color(0,0,0),
															tlp.Color(0,0,0),tlp.Color(0,51,0),tlp.Color(0,51,0),tlp.Color(0,102,0),tlp.Color(0,153,0),tlp.Color(0,153,0),tlp.Color(0,204,0),tlp.Color(0,204,0),
															tlp.Color.Green,tlp.Color.Green,tlp.Color.Green])
	parameter = tlp.getDefaultPluginParameters("Color Mapping", graph)
	parameter["color scale"] = colorScale
	graph.applyColorAlgorithm("Color Mapping", color, parameter)


def translate(x, y, graph, layout):

	## This function allows the translation of a graph coordinates ; from the layout of graph
	## It takes a "x" value for horizontal move and "a" y value for vertical move

	layout.translate(tlp.Vec3f(x,y,0), graph)


def makeGrid(subGraph, graph, tpNumber, numberOfColumn, layout):

	##
	##

	x = -15000
	y = 0
	if tpNumber == numberOfColumn or tpNumber==numberOfColumn*2 or tpNumber==numberOfColumn*3:
		x = x + 15000*numberOfColumn
		y = y + 15000
	translate(x, y, graph, layout)



###############################################################
###################### Main ###################################
###############################################################

def main(graph):

	##### Variables from Tulip API #####

	label = graph.getStringProperty("viewLabel")
	color = graph.getColorProperty("viewColor")
	locus = graph.getStringProperty("Locus")
	size = graph.getSizeProperty("viewSize")
	layout = graph.getLayoutProperty("viewLayout")
	negative = graph.getBooleanProperty("Negative")
	positive = graph.getBooleanProperty("Positive")
	treeLayout = graph.getLayoutProperty("viewLayout")
	shape = graph.getIntegerProperty("viewShape")

	tp1_s = graph.getDoubleProperty("tp1 s")
	tp10_s = graph.getDoubleProperty("tp10 s")
	tp11_s = graph.getDoubleProperty("tp11 s")
	tp12_s = graph.getDoubleProperty("tp12 s")
	tp13_s = graph.getDoubleProperty("tp13 s")
	tp14_s = graph.getDoubleProperty("tp14 s")
	tp15_s = graph.getDoubleProperty("tp15 s")
	tp16_s = graph.getDoubleProperty("tp16 s")
	tp17_s = graph.getDoubleProperty("tp17 s")
	tp2_s = graph.getDoubleProperty("tp2 s")
	tp3_s = graph.getDoubleProperty("tp3 s")
	tp4_s = graph.getDoubleProperty("tp4 s")
	tp5_s = graph.getDoubleProperty("tp5 s")
	tp6_s = graph.getDoubleProperty("tp6 s")
	tp7_s = graph.getDoubleProperty("tp7 s")
	tp8_s = graph.getDoubleProperty("tp8 s")
	tp9_s = graph.getDoubleProperty("tp9 s")


	#### Functions call ####

	### Part I ###
	colorLabels(graph, color, size, label, locus)
	drawForceModel(graph, layout)
	colorEdges(graph, negative, positive, color)

	### Part II ###
	tree = graph.addSubGraph("Hierarchic Tree")
	root = tree.addNode()
	currentCluster = graph.getSubGraph("Genes interactions")
	constructTree(tree, currentCluster, root)
	drawRadial(tree, treeLayout)
	graphProperty(tree, tp4_s, color)
	constructBundles(currentCluster, tree, layout, treeLayout, shape)

	### Part III ###
	tps = [tp1_s, tp2_s, tp3_s, tp4_s, tp5_s, tp6_s, tp7_s, tp8_s, tp9_s, tp10_s, tp11_s, tp12_s, tp13_s, tp14_s, tp15_s, tp16_s, tp17_s]
	constructSmallMultiples(tree, currentCluster, tps, color, treeLayout)
