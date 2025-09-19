export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt ROOT"
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt matplotlib";
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt lxml";
source setup.sh

