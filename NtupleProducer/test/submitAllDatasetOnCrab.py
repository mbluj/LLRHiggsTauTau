import os
import sys
import re

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

###################################################################
#### Parameters to be changed for each production

#PROCESS = ["HHBACKGROUNDS"] # select blocks of datasets to be processed
#tag = "produzione_MC_3Ago2015" # a folder with cra3_[tag] is created, ad appended in many places
#PROCESS = ["2015RUNBDATA"]
#tag = "produzione_DATA_3Ago2015"

#datasetsFile = "datasets.txt" # name of file containing datasets

PROCESS = ["PROD_PARTIAL"]
tag = "llrNtuples_partial_5Ago_resub"
datasetsFile = "datasets_Enriched.txt"

#PROCESS = ["HHBACKGROUNDS_RES"]
#tag = "produzione_DATA_3Ago2015_resubTTJets_enrich"

#datasetsFile = "datasets.txt" # name of file containing datasets

EnrichedToNtuples = False # do not create ntuples on CRAB because it is very slow

###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
sectionBeginEnd = "==="

# check if file with dataset exist
if not os.path.isfile(datasetsFile):
    print "File %s not found!!!" % datasetsFile
    sys.exit()

#check if directory exists
crabJobsFolder = "crab3_" + tag
if os.path.isdir(crabJobsFolder):
    print "Folder %s already exists, please change tag name or delete it" % crabJobsFolder
    sys.exit()

# grep all datasets names, skip lines with # as a comment
# block between === * === are "sections" to be processed

currSection = ""
dtsetToLaunch = []

# READ INPUT FILE
with open(datasetsFile) as fIn:
    for line in fIn:
        line = line.strip() # remove newline at the end and leding/trailing whitespaces
        
        if not line: #skip empty lines
            continue

        if comment in line:
            continue
        
        #print line        
        words = line.split()
        if len(words) >= 3:
            if words[0] == sectionBeginEnd and words[2] == sectionBeginEnd: 
                currSection = words[1]
        else:
            if currSection in PROCESS:
                dtsetToLaunch.append(line)

# CREATE CRAB JOBS
os.system ("voms-proxy-init -voms cms")

for name in PROCESS: crabJobsFolder + "_" + name
print crabJobsFolder
os.system ("mkdir %s" % crabJobsFolder)

counter = 1 # appended to the request name to avoid overlaps between datasets with same name e.g. /DoubleEG/Run2015B-17Jul2015-v1/MINIAOD vs /DoubleEG/Run2015B-PromptReco-v1/MINIAOD
outlog = open ((crabJobsFolder + "/submissionLog.txt"), "w")
for dtset in dtsetToLaunch:
    dtsetNames = dtset
    if '/MINIAODSIM' in dtset:
        dtsetNames = dtset.replace('/MINIAODSIM', "")
    elif '/MINIAOD' in dtset:
        dtsetNames = dtset.replace('/MINIAOD', "")
    dtsetNames = dtsetNames.replace('/', "__")
    dtsetNames = dtsetNames.strip("__") # remove leading and trailing double __ 
    #dtSetName = dtsetNames[1]
    command = "crab submit -c crab3_template.py"
    command += " General.requestName=%s" % (dtsetNames + "_" + tag + "_" + str(counter))
    command += " General.workArea=%s" % crabJobsFolder
    command += " Data.inputDataset=%s" % dtset
    command += " Data.outLFNDirBase=/store/user/lcadamur/HHNtuples/%s" % tag # FIXME: stesso folder per datasets che hanno lo stesso nome?
    command += " Data.publishDataName=%s" % tag
    if (EnrichedToNtuples): command += " Data.inputDBS=phys03" # if I published the dataset need to switch from global (default)
    if (EnrichedToNtuples): command += " JobType.psetName=ntuplizer.py" # run a different python config for enriched
    if (EnrichedToNtuples): command += " Data.publication=False" # cannot publish flat root ntuples
    print command ,  "\n"
    os.system (command)
    outlog.write(command + "\n\n")
    counter = counter + 1