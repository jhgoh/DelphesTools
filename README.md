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

## Step3: project on a MxM "image"
We use NERSC's script to convert Delphes root files to hdf5 with image projections.
- File names should be in a form of SAMPLENAME-SUFFIX.root
- SAMPLENAME should be in the cross section table, config/DelphesXSec
  - RPV10\_1400\_850
  - QCDBkg\_JZ3\_160\_400
  - QCDBkg\_JZ4\_400\_800
  - QCDBkg\_JZ5\_800\_1300
  - QCDBkg\_JZ6\_1300\_1800
  - QCDBkg\_JZ7\_1800\_2500
  - QCDBkg\_JZ8\_2500\_3200
  - QCDBkg\_JZ9\_3200\_3900
  - QCDBkg\_JZ10\_3900\_4600
  - QCDBkg\_JZ11\_4600\_5300
  - QCDBkg\_JZ12\_5300\_7000
- List of files in a txt file

Example:
```
mv DELPHES.root RPV10_1400_850-xxxx.root
echo ../RPV10_1400_850-xxxx.root > fileList.txt ## NOTE the relative path
```

```
git clone https://github.com/eracah/atlas_dl
cd atlas_dl
./scripts/prepare_data.py --input-type delphes --output-h5 ../../data.h5 --bins 64 ../fileList.txt 
```
