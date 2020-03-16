#!/usr/bin/env python

import sys, os
if len(sys.argv) < 4:
    print "%s MACRO INPUT.root OUTPUT.root" % sys.argv[0]
    sys.exit(1)
macro = sys.argv[1]
inFile = sys.argv[2]
outFile = sys.argv[3]

from ROOT import *
delphesPath = "Delphes"
gSystem.AddIncludePath('-I"%s"' % delphesPath)
gSystem.AddDynamicPath(delphesPath)
gSystem.AddLinkedLibs('-L"%s"' % delphesPath)
gSystem.Load("libDelphes")

gROOT.ProcessLine(".L %s+" % macro)
gROOT.ProcessLine('%s("%s", "%s");' % (macro.rsplit('.')[0], inFile, outFile))

