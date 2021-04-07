#voms-proxy-init --voms cms --valid 100:00

# do we submit or just generate submit scripts
dosubmit=false
mode='MINIAOD'
nversion='Phase2_Rechit'
# make the directory where we'll submit from
thesubdir="./gitignore/${nversion}"
mkdir -p ${thesubdir}

# copy necessary files into submit directory
if [ ${mode} = reco ]
then
 msubmitconfig="run_mc_reco113X.py"
elif [ $mode} = RAW ]
then 
 msubmitconfig="run_mc_RAW113X.py"
else
 msubmitconfig="run_mc_113X.py"
fi

# copy cmsRun configuration to submit directory
cp ${msubmitconfig}  ${thesubdir}


# sample names to run over
samples=( \
"ZMM"                   \
#"WJetsToLNu"                   \
)
# loop over mc samples
for samplename in ${samples[@]}
do
 datasetname="/RelValZMM_14/CMSSW_11_3_0_pre2-PU25ns_113X_mcRun4_realistic_v2_2026D49PU200-v2/GEN-SIM-RECO"
 #datasetname="/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"

 submitname="submit_${samplename}"
 submitfile="${thesubdir}/${submitname}.py"

 # set variables for submitting this specific sample
 WORKAREA="'crabsubmits_${nversion}'"
 CMSRUNCONFIG="'${msubmitconfig}'"
 UPERJOB="30"
 SPLITTING="'FileBased'"

 NUNITS="-1"
 REQUESTNAME="'${samplename}'"
 DATASET="'${datasetname}'"
 STORESITE="'T3_US_FNALLPC'"
 OUTLFNBASE="'/store/group/lpchbb/LLDJntuples/${nversion}'"
 MAXMEM="5000"

 # copy and then fill template for crab submits
 cp ./crab_template.py             "${submitfile}"
 sed -i "s@WORKAREA@${WORKAREA}@g"         "${submitfile}"
 sed -i "s@CMSRUNCONFIG@${CMSRUNCONFIG}@g" "${submitfile}"
 sed -i "s@NUNITS@${NUNITS}@g"             "${submitfile}"
 sed -i "s@UPERJOB@${UPERJOB}@g"           "${submitfile}"
 sed -i "s@SPLITTING@${SPLITTING}@g"       "${submitfile}"
 sed -i "s@REQUESTNAME@${REQUESTNAME}@g"   "${submitfile}"
 sed -i "s@DATASET@${DATASET}@g"           "${submitfile}"
 #sed -i "s@LUMIMASK@${LUMIMASK}@g"         "${submitfile}"
 sed -i "s@STORESITE@${STORESITE}@g"       "${submitfile}"
 sed -i "s@OUTLFNBASE@${OUTLFNBASE}@g"     "${submitfile}"
 sed -i "s@MAXMEM@${MAXMEM}@g"             "${submitfile}"

 # submit the jobs
 if [ ${dosubmit} = true ]
 then
  pushd ${thesubdir} > /dev/null
  python ${submitfile}
  popd > /dev/null
 fi

done
