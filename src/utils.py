from __future__ import division

import CalculateCoElutionScores as CS
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import GoldStandard as GS
import copy

def bench_clf(scoreCalc, train, eval, clf, outDir, verbose=False):
	_, data_train, targets_train = scoreCalc.toSklearnData(train)
	_, data_eval, targets_eval = scoreCalc.toSklearnData(eval)
	clf.fit(data_train, targets_train)
	precision, recall, fmeasure, auc_pr, auc_roc, curve_pr, curve_roc = clf.eval(data_eval, targets_eval)
	plotCurves([("", curve_pr)], outDir + ".pr.pdf", "Recall", "Precision")
	plotCurves([("", curve_roc)], outDir + ".roc.pdf", "False Positive rate", "True Positive Rate")
	recall_vals, precision_vals, threshold = curve_pr
	threshold = np.append(threshold, 1)
	plotCurves([("Precision", (precision_vals, threshold)), ("Recall", (recall_vals, threshold))], outDir + ".cutoff.pdf", "Cutoff", "Evaluation metric score")
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

# @author Florian Goebels
def predictInteractions(scoreCalc, clf, gs, verbose= False):
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(gs)
	clf.fit(data_train, targets_train)
	print data_train.shape
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
	for line in range(scoreCalc.to_predict):
		if k % 100000==0 and k != 0:
			tmpscores = np.zeros((100000, num_features))
			if verbose:
				print "Completed chunk %i" % chunk_num
				chunk_num += 1
			k = 0
		edge, edge_scores = scoreCalc.get_next()
		if edge == None: continue
		print edge
		edge_scores = edge_scores.reshape(1, -1)
		edges[k] = edge
		tmpscores[k,0:(edge_scores.shape)[1]] = edge_scores
		k += 1
	scoreCalc.close()
	out.extend(getPredictions(tmpscores[0:k,:], edges[0:k], clf))
	return out


def make_predictions(score_calc, mode, clf, gs, fun_anno=""):
	def get_edges_from_network(network):
		out = {}
		for edge in network:
			edge, score = edge.rsplit("\t", 1)
			out[edge] = score
		return out

	networks = []
	# predicts using experiment only
	if mode == "exp" or mode == "BR": networks.append(predictInteractions(score_calc, clf, gs))

	#predicts using fun_anno only
	if mode == "fa"or mode == "BR":
		if fun_anno=="":
			# TODO make illigal argument error
			print "if using only functional annotation for prediction functional annotation (fun_anno param != "") must not be empty"
			sys.exit()
		networks.append(predictInteractions(fun_anno, clf, gs))

	#predict using both functional annotation and exp
	if mode == "comb"or mode == "BR":
		tmp_score_calc = copy.deepcopy(score_calc)
		tmp_score_calc.add_fun_anno(fun_anno)
		networks.append(predictInteractions(tmp_score_calc, clf, gs))

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
		br_edges = set(exp.keys()) + (set(merged.keys()) - set(fa.keys()))
		br_network = []
		for edge in br_edges:
			if edge in exp: score = exp[edge]
			else: score = fa[edge]
			br_network.append("%s\t%s" % (edge, score))
		return br_network

def predict_clusters(predF, outF):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	clustering_CMD = "java -jar %s/cluster_one-1.0.jar %s > %s" % (dir_path, predF, outF)
	os.system(clustering_CMD)

def load_data(data_dir, scores, orthmap=""):
	paths = [os.path.join(data_dir,fn) for fn in next(os.walk(data_dir))[2]]
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