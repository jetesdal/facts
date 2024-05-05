#!/bin/bash

if [ -z "$WORKDIR" ]; then  
   WORKDIR=/scratch/`whoami`/test.`date +%s`
fi
mkdir -p $WORKDIR

if [ -z "$OUTPUTDIR" ]; then  
   OUTPUTDIR=/scratch/`whoami`/test.`date +%s`/output
fi
mkdir -p $OUTPUTDIR
BASEDIR=`pwd`

#EXPERIMENT STEP:  sealevel_step 


# - Pipeline onemodule.onemodule.tlm.sterodynamics:


PIPELINEDIR=$WORKDIR/onemodule.onemodule.tlm.sterodynamics
mkdir -p $PIPELINEDIR

cd $BASEDIR
cp test/onemodule/location.lst $PIPELINEDIR

# ---- Stage preprocess:

cd $BASEDIR
cp ./modules/tlm/sterodynamics/readMarzeion.py ./modules/tlm/sterodynamics/read_CSIRO.py ./modules/tlm/sterodynamics/IncludeCMIP6ZOSModels.py ./modules-data/tlm_sterodynamics_preprocess_data.tgz ./modules-data/tlm_sterodynamics_cmip6_data.tgz ./modules-data/ipccar6_climate_data.tgz ./modules/tlm/sterodynamics/IncludeCMIP6Models.py ./modules/tlm/sterodynamics/SmoothZOSTOGA.py ./modules/tlm/sterodynamics/read_locationfile.py ./modules/tlm/sterodynamics/Import2lmData.py ./modules/tlm/sterodynamics/tlm_sterodynamics_preprocess.py ./modules/tlm/sterodynamics/Smooth.py $PIPELINEDIR
cd $PIPELINEDIR
. $RP_PILOT_SANDBOX/env/rp_named_env.rp.sh || true
tar -xvf tlm_sterodynamics_preprocess_data.tgz 2> /dev/null; rm tlm_sterodynamics_preprocess_data.tgz
tar -xvf tlm_sterodynamics_cmip6_data.tgz 2> /dev/null; rm tlm_sterodynamics_cmip6_data.tgz
tar -xvf ipccar6_climate_data.tgz 2> /dev/null; rm ipccar6_climate_data.tgz
pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml h5py
python3 tlm_sterodynamics_preprocess.py --scenario ssp585 --baseyear 2005 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --pipeline_id onemodule.onemodule.tlm.sterodynamics


# ---- Stage fit:

cd $BASEDIR
cp ./modules/tlm/sterodynamics/tlm_sterodynamics_fit.py $PIPELINEDIR
cd $PIPELINEDIR
. $RP_PILOT_SANDBOX/env/rp_named_env.rp.sh || true
python3 tlm_sterodynamics_fit.py --pipeline_id onemodule.onemodule.tlm.sterodynamics


# ---- Stage project:

cd $BASEDIR
cp ./modules/tlm/sterodynamics/tlm_sterodynamics_project.py $PIPELINEDIR
cd $PIPELINEDIR
. $RP_PILOT_SANDBOX/env/rp_named_env.rp.sh || true
python3 tlm_sterodynamics_project.py --pipeline_id onemodule.onemodule.tlm.sterodynamics --nsamps 2000

cp onemodule.onemodule.tlm.sterodynamics_globalsl.nc $OUTPUTDIR

# ---- Stage postprocess:

cd $BASEDIR
cp ./modules/tlm/sterodynamics/tlm_sterodynamics_postprocess.py $PIPELINEDIR
cd $PIPELINEDIR
. $RP_PILOT_SANDBOX/env/rp_named_env.rp.sh || true
python3 tlm_sterodynamics_postprocess.py --nsamps 2000 --pipeline_id onemodule.onemodule.tlm.sterodynamics

cp onemodule.onemodule.tlm.sterodynamics_localsl.nc $OUTPUTDIR
