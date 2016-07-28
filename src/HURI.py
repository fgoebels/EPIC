#!/usr/bin/env python
from __future__ import division
import sys
import os
import re
import subprocess
import gzip
import math
import numpy as np

def predictFromFile(scoreTableF, referenceF, outF):
	reference = GS.Goldstandard_from_File(referenceF)
	scoreCalc = CalculateCoElutionScores()
	scoreCalc.initFromFile(scoreTableF, reference.goldstandard)
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(labels=set(["positive", "negative"]))
#	print data_train.shape
#	ids_pred, data_pred, targets_pred = scoreCalc.toSklearnData(labels=set(["?"]))
#	clf = CLF_Wrapper(data_train, targets_train)
#	preds = clf.predict(data_train)
	outFH = open(outF, "w")
	for i in range(len(ids_train)):
		print >> outFH, "%s\t%f" % ("\t".join(ids_train[i]), targets_train[i])
	outFH.close()


def main():
	print "fubar"


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
