# nurion4hep/DelphesTools
Setup Delphes to produce/preprocess input data for DL

## Prerequisites
- A Linux workstation configured with CVMFS, CMSSW

## Step1: Install packages
```
git clone https://github.com/nurion4hep/DelphesTools
./install.sh
```

## Step2: Run the Delphes
Example: extract prunedGenParticles+packedGenParticles from CMS MiniAOD and run the Delphes simulator.
```
cd Delphes
./DelphesCMSFWLite cards/delphes_card_CMS.tcl ../DELPHES.root root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FED04F50-A3BE-E611-A893-0025900E3502.root
cd ..
```
You will have DELPHES.root 

## Step3: Reduce root file / project on a MxM "image"
```
./run.py delphes2Image.C DELPHES.root FLAT.root
```
Then you will have FLAT.root which contains a flat TTree "Events".
Branch names are similar to the CMS NanoAOD, hoping to be analyzed with almost same analysis macro.

## Step4: Convert flat trees to hdf5
TBA
