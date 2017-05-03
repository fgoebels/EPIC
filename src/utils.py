from __future__ import division

import CalculateCoElutionScores as CS
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import GoldStandard as GS
import copy
import json

def bench_clf(scoreCalc, train, eval, clf, outDir, verbose=False, format = "pdf"):
	_, data_train, targets_train = scoreCalc.toSklearnData(train)
	_, data_eval, targets_eval = scoreCalc.toSklearnData(eval)
	clf.fit(data_train, targets_train)
	precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.eval(data_eval, targets_eval)
	plotCurves([("", curve_roc)], outDir + ".roc." + format, "False Positive rate", "True Positive Rate")
	recall_vals, precision_vals, threshold = curve_pr
	plotCurves([("", (precision_vals, recall_vals))], outDir + ".pr." + format, "Recall", "Precision")

	threshold = np.append(threshold, 1)
	plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], outDir + ".cutoff." + format, "Cutoff", "Evaluation metric score")
	if verbose:
		rownames = ["Precision", "Recall", "F-Measure", "AUC PR", "AUC ROC"]
		val_scores = [precision, recall, fmeasure, auc_pr, auc_roc]
		for i in range(len(rownames)):
			print rownames[i]
			print val_scores[i]

# @author: Florian Goebels
# makes precision recall plot for mutliple rp ccurves
# @Param:
#	curves list of tuples with (name, precision, recall) which should be plotted
#	outF Pdf file location for the created plot
def plotCurves(curves, outF, xlab, ylab):
	plt.clf()
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	cols = ['b', 'r', 'c', 'm', 'y', 'k']
	for (name, curve) in curves:
		x, y = curve[0:2]
		if name != "":
			plt.plot(x, y, label=name, color = cols.pop())
		else:
			plt.plot(x, y, color=cols.pop())
	art = []
	if len(curves)>1:
		lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1),  ncol = 5, fontsize=8)
		art.append(lgd)
	plt.savefig(outF, additional_artists=art, bbox_inches="tight")
	plt.close()

# @author Florian Goebels
def predictInteractions(scoreCalc, clf, gs, verbose= False):
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(gs)
	clf.fit(data_train, targets_train)
	num_features = data_train.shape[1]
	def getPredictions(scores, edges, clf):
		out = []
		pred_prob = clf.predict_proba(scores)
		pred_class = clf.predict(scores)
		for i, prediction in enumerate(pred_class):
			if prediction == 1:
				out.append("%s\t%f" % (edges[i], pred_prob[i]))	#Alternative code that also print label:out.append("%s\t%f\t%i" % (edges[i], pred_prob[i], prediction))
		return out
	out = []
	tmpscores = np.zeros((100000, num_features))
	edges = [""]*100000
	k = 0
	chunk_num=1
	scoreCalc.open()
	print "to predict: %i" % scoreCalc.to_predict
	for line in range(scoreCalc.to_predict):
		if k % 100000==0 and k != 0:
			out.extend(getPredictions(tmpscores[0:k, :], edges[0:k], clf))
			tmpscores = np.zeros((100000, num_features))
			edges = [""] * 100000
			if verbose:
				print "Completed chunk %i" % chunk_num
				chunk_num += 1
			k = 0
		edge, edge_scores = scoreCalc.get_next()
		if edge == "" or edge_scores == []: continue
		edge_scores = edge_scores.reshape(1, -1)
		edges[k] = edge
		tmpscores[k,0:(edge_scores.shape)[1]] = edge_scores
		k += 1
	scoreCalc.close()
	out.extend(getPredictions(tmpscores[0:k,:], edges[0:k], clf))
	return out

def get_FA_data(anno_source, file=""):
	functionalData = ""
	if anno_source == "GM":

		genemania = CS.Genemania("6239")
		functionalData = genemania.getScoreCalc()

	elif anno_source == "STRING":

		string = CS.STRING("6239")
		functionalData = string.getScoreCalc()

	elif anno_source == "FILE":

		# the supplied functional evidence data needs to have the correct header row...
		externaldata = CS.ExternalEvidence(file)
		externaldata.readFile()
		functionalData = externaldata.getScoreCalc()

	else:
		print "EPIC only support GeneMane, STRING, and flat file input please use the followign tags for anno_source GM, STRING, FILE. Returning empty string object."
	return functionalData

def make_predictions(score_calc, mode, clf, gs, fun_anno="", verbose = False):
	def get_edges_from_network(network):
		out = {}
		for edge in network:
			edge, score = edge.rsplit("\t", 1)
			out[edge] = score
		return out

	networks = []
	# predicts using experiment only
	if mode == "exp" or mode == "BR": networks.append(predictInteractions(score_calc, clf, gs, verbose))

	#predicts using fun_anno only
	if mode == "fa"or mode == "BR":
		if fun_anno=="":
			# TODO make illigal argument error
			print "if using only functional annotation for prediction functional annotation (fun_anno param != "") must not be empty"
			sys.exit()
		networks.append(predictInteractions(fun_anno, clf, gs, verbose))

	#predict using both functional annotation and exp
	if mode == "comb"or mode == "BR":
		tmp_score_calc = copy.deepcopy(score_calc)
		tmp_score_calc.add_fun_anno(fun_anno)
		networks.append(predictInteractions(tmp_score_calc, clf, gs, verbose))

	# return error when no networks is predicted
	if len(networks) == 0:
		print "Error no networks predicted"
		sys.exit()
	# return finised network when only one network is predicted, which happens in any mode expect final
	elif len(networks) ==1:
		return networks[0]
	# use bias reduced method to merge experimental, functional annotation, and combined network
	else:
		exp = get_edges_from_network(networks[0])
		fa = get_edges_from_network(networks[1])
		merged = get_edges_from_network(networks[2])
		br_edges = set(exp.keys()) | (set(merged.keys()) - set(fa.keys()))
		br_network = []
		for edge in br_edges:
			if edge in exp: score = exp[edge]
			else: score = merged[edge]
			br_network.append("%s\t%s" % (edge, score))
		return br_network

def predict_clusters(predF, outF):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	clustering_CMD = "java -jar %s/cluster_one-1.0.jar %s > %s" % (dir_path, predF, outF)
	os.system(clustering_CMD)

def load_data(data, scores, orthmap=""):

	if type(data) is list:
		paths = data
	else:
		paths = [os.path.join(data,fn) for fn in next(os.walk(data))[2]]

	elutionDatas = []
	elutionProts = set([])
	for elutionFile in paths:
		if elutionFile.rsplit(os.sep, 1)[-1].startswith("."): continue
		elutionFile = elutionFile.rstrip()
		elutionData = CS.ElutionData(elutionFile)
		if orthmap !="":
			if orthmap != False:
				mapper = GS.Inparanoid("", inparanoid_cutoff=1)
				mapper.readTable(orthmap, direction=0)
				elutionData.orthmap(mapper)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())
		for score in scores:
			score.init(elutionData)
	return elutionProts, elutionDatas

def create_goldstandard(target_taxid, valprots):
	def create_gs_set(cluster_obj, orthmap, name, valprots):
		gs =  GS.Goldstandard_from_Complexes(name)
		gs.make_reference_data(cluster_obj, orthmap, found_prots=valprots)
		return gs

	if target_taxid !="9606":
		orthmap = GS.Inparanoid(taxid=target_taxid)
		reference_clusters = [GS.Intact_clusters(True), GS.CORUM(True), GS.QuickGO("9606", True), GS.QuickGO(target_taxid, False)]
		#reference_clusters = [GS.Intact_clusters(True)]
	else:
		reference_clusters = [GS.Intact_clusters(False), GS.CORUM(False), GS.QuickGO("9606", False)]
		orthmap = ""
	all_gs =  create_gs_set(reference_clusters, orthmap, "Training", valprots)
	return all_gs


def clustering_evaluation(eval_comp, pred_comp, prefix, verbose= True):
	head = "\t".join(["%s %s" % (prefix, h) for h in ["mmr", "overlapp", "simcoe", "mean_simcoe_overlap", "sensetivity", "ppv", "accuracy", "sep"]])
	cluster_scores = "\t".join(map(str, pred_comp.clus_eval(eval_comp)))
	if verbose:
		tmp_head = head.split("\t")
		tmp_scores = cluster_scores.split("\t")
		for i in range(len(tmp_head)):
			print "%s\t%s" % (tmp_head[i], tmp_scores[i])
	return cluster_scores, head

def clusters_to_json(clusters, network, frac_names, eData):
	graph = {}
	for line in network:
		edge, score = line.rsplit("\t", 1)
		graph[edge] = score

	cy_elements = []
	nodes = [] # used to send entwork to cytoscape
	edges = [] # used to send entwork to cytoscape

	net_nodes = set([])
	for complex in clusters.complexes:
		prots = list(clusters.complexes[complex])
		for i in range(len(prots)):
			protA = prots[i]
			for j in range(i+1, len(prots)):
				protB = prots[j]
				edge = "\t".join(sorted([protA, protB]))
				score = 0.5
				if edge in graph: score = graph[edge]
				nodeA = "%s_%s" % (protA, str(complex))
				nodeB = "%s_%s" % (protB, str(complex))
				net_nodes.add(nodeA)
				net_nodes.add(nodeB)
				edge = {
					'group': 'edges',
					'data': {
						'source': nodeA,
						'target': nodeB,
						'score': float(score),
					}
				}
				cy_elements.append(edge)
				edges.append(edge)

	for gene in net_nodes:
		name, cluster_id = gene.split("_")
		node = {
			'group': 'nodes',
			'data': {
				'id': str(gene),
				'name': name,
				'cluster_id': cluster_id
			}}
		for i in range(len(frac_names)):
			score = 0
			if name in eData: score = eData[name][i]
			node['data'][frac_names[i]] = float(score)
		cy_elements.append(node)
		nodes.append(node)

	return json.dumps(cy_elements, default=lambda cy_elements: cy_elements.__dict__), edges, nodes


def json_to_cy_js(div_id, json_str):
	return """

	                $('#cy').show();
	                var cy = window.cy = cytoscape({
	                    container: document.getElementById('%s'),
	                    layout: { },
	                    elements: %s,
	                    style: [
	                     {
	                        selector: 'node',
	                        style: {
	                          'content': 'data(name)',
	                          'font-size': 12,
	                          'text-valign': 'center',
	                          'text-halign': 'center',
	                          'background-color': '#555',
	                          'text-outline-color': '#555',
	                          'text-outline-width': 1.75,
	                          'color': '#fff',
	                          'overlay-padding': 6,
	                          'z-index': 10
	                        }
	                      },
	                      {
	                        selector: 'edge',
	                        style: {
                              'content': 'data(score)',
	                          'line-color': 'black',
	                          'width': 'mapData(score, 0.5, 1, 0, 20)',
	                        }
	                      },
	                    ]
	                });

                    cy.elements().components().forEach( (eles, i, components) => {
                    let n = Math.floor( Math.sqrt( components.length ) );
                    let w = 2000; // width of bb for 1 cmp
                    let h = 2000; // height "

                    eles.makeLayout({
                        name: 'circle',
                        boundingBox: {
                        x1: w * (i %% n),
      x2: w * (1 + (i %% n)),  // this line fixed
      y1: Math.floor(i / n) * h,
      y2: (1 + Math.floor(i / n)) * h
                    }
                    }).run();
                    });
	            """ % (div_id, json_str)

def elutionDatas_to_treeview(eDatas, foundprots):
	out = {}
	colnums = {}
	header = []
	for eData in eDatas:
		name = eData.name
		colnum = eData.elutionMat.shape[1]
		colnums[name] = colnum
		prefix = "%s.F" % name
		header.extend(map(lambda x: "%s%s" % (prefix, x), range(1,colnum+1)))

	for prot in foundprots:
		out[prot] = []
		for eData in eDatas:
			scores = [0]*colnums[eData.name]
			if eData.hasProt(prot):
				scores = eData.getElution(prot)
			out[prot].extend(scores)
	return header, out


def prep_network_for_cy(nodes, edges):
	# Basic Setup
	PORT_NUMBER = 1234
	# IP = '192.168.1.1'
	IP = 'localhost'
	BASE = 'http://' + IP + ':' + str(PORT_NUMBER) + '/v1/'

	# Header for posting data to the server as JSON
	HEADERS = {'Content-Type': 'application/json'}
	network_cy = {
		'data': {'name': "EPIC clusters"},
		"elements": {"nodes": nodes, 'edges': edges}
	}
	return BASE, json.dumps(network_cy), HEADERS