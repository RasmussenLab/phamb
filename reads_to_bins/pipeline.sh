#!/bin/bash
#
# 
#
# Usage:
# First copy "config.file" file to your working directory and edit it
# Prepare Fastq files manually by cleaning names, merging lanes and putting them in a 'fastq' directory. Or automatically specifying a Directory of fastq run-directories using prepare_fastq_files.sh.
# make sample_table.txt using create_sample_table.sh


readonly VERSION="2020"

### If sourced on computerome and an Error is encountered, the session will be terminated. Can be Outcommented for testing purposes.
#set -euo pipefail

if [ "$#" -gt 0 ]; then
    echo "ERROR: ... does not take arguments, but you provided: $@"
    exit 1
fi

### Read in the configuration file:
source ./pipeline.config.file

### Set up the log file (all output from the rest of this script goes to both screen and log file)
exec &> >(tee -a "${LOGFILE}")

### Setting source working directory
readonly SRC_DIR=reads_to_bins


echo "This is version ${VERSION}"
echo "Running script from: ${SRC_DIR}"
echo "Current time: `date`"
echo "Current working directory: `pwd`"
echo "SAMPLETABLE=${SAMPLETABLE}"
echo "TRIMMOMATIC_PARAMS=${QC_TRIMMOMATIC_PARAMS}"
echo "QC_CLEAN=${QC_CLEAN}"
echo "QC_HOST_DATABASE=${QC_HOST_DATABASE}"
echo "ASSEMBLY_METHOD=${ASSEMBLY_METHOD}"
echo "ASSEMBLY_GSIZE=${ASSEMBLY_GSIZE}"
echo "ONLY_ASSEMBLIES=${ONLY_ASSEMBLIES}"

### Create log dir for PBS logs
mkdir -p log/pbs

if [ "${ONLY_ASSEMBLIES}" == "true" ]; then
        echo "Skipping QC and Assembly"	
else
	## List all samples and (if doing assembly) check that FASTQ files exist
	NSAMPLES=0
	N_UNPAIRED_SAMPLES=0
	echo -n "Input files:"
	while read -r SAMPLEID FQ1 FQ2 PLACEHOLDER; do
	    echo -n "${SAMPLEID},"
	    [ "${DO_ASSEMBLY}" != "true" ] || [ -f "${FQ1}" ] || { echo -e "\nFwd read file ${FQ1} not found!" ; exit 1; }
	    if [ -z "${FQ2}" ]; then
		 N_UNPAIRED_SAMPLES=$((${N_UNPAIRED_SAMPLES} + 1))
	    else
		[ "${DO_ASSEMBLY}" != "true" ] || [ -f "${FQ2}" ] || { echo -e "\nRev read file ${FQ2} not found!" ; exit 1; }
	    fi
	    NSAMPLES=$((${NSAMPLES} + 1))
	done < "${SAMPLETABLE}"
	echo "total=${NSAMPLES}"
	echo "of which ${N_UNPAIRED_SAMPLES} are unpaired."


	### Run QC / Trimming of FastQ samples

	[[ $QC_PPN -gt 40 ]] || [[ $QC_PPN -lt 1 ]] && { echo "ERROR: 'QC_PPN' must be between 1 and 40." ; exit 1; }


	## Submit read mapping jobs
	if [ "${DO_QC}" = "true" ]; then
	    readonly MAXMEM=$((130*$QC_PPN/40))
	    readonly QC_RESOURCES="nodes=1:ppn=${QC_PPN}:thinnode,walltime=24:00:00,mem=${MAXMEM}gb"
	    readonly QC_LOGS="-o log/pbs/QC.pbs.o -e log/pbs/QC.pbs.e"
	    readonly QC_ARGS="${SAMPLETABLE} ${QC_HOST_DATABASE} ${QC_PPN} ${QC_CLEAN} ${QC_METHOD} ${QC_PREQC} ${QC_TRIMMOMATIC_PARAMS}"
	    QCJOBID=$(myqsub -d `pwd` -t ${RUN_FROM_NSAMPLES}-${RUN_NSAMPLES}%${NJOBS} ${QC_LOGS} -l ${QC_RESOURCES} -F "${QC_ARGS}" ${SRC_DIR}/QC.pbs  | cut -d . -f 1)
	    QCDEPENDENCY="-W depend=afterokarray:$QCJOBID"
	    echo "QC job submitted: ${QCJOBID}"
	else
	    echo "Skipping QC step"
	    QCDEPENDENCY=""
	fi


	### Running Assemblies 

	### Check NNodes input.
	[[ $ASSEMBLY_PPN -gt 40 ]] || [[ $ASSEMBLY_PPN -lt 1 ]] && { echo "ERROR: 'ASSEMBLY_PPN' must be between 1 and 40." ; exit 1; }

	if [ "${DO_ASSEMBLY}" = "true" ]; then

	    readonly ASSEMBLY_MAXMEM=$((130*$ASSEMBLY_PPN/40))
	    readonly ASSEMBLY_RESOURCES="nodes=1:ppn=${ASSEMBLY_PPN}:thinnode,walltime=24:00:00,mem=${ASSEMBLY_MAXMEM}gb"
	    readonly ASSEMBLY_LOGS="-o log/pbs/ASSEMBLY.pbs.o -e log/pbs/ASSEMBLY.pbs.e"
	    readonly ASSEMBLY_ARGS="${SAMPLETABLE} ${ASSEMBLY_GSIZE} ${ASSEMBLY_METHOD} ${ASSEMBLY_TMPDIR} ${ASSEMBLY_PPN} ${ASSEMBLY_MAXMEM} ${ASSEMBLY_ISOLATE} ${ASSEMBLY_MINIMUM_CONTIG_LENGTH}"
	    ASSEMBLYJOBID=$(myqsub -d `pwd` ${QCDEPENDENCY} -t ${RUN_FROM_NSAMPLES}-${RUN_NSAMPLES}%${NJOBS} ${ASSEMBLY_LOGS} -l ${ASSEMBLY_RESOURCES} -F "${ASSEMBLY_ARGS}" ${SRC_DIR}/ASSEMBLY.pbs  | cut -d . -f 1)
	    ASSEMBLYDEPENDENCY="-W depend=afterokarray:$ASSEMBLYJOBID"

	    echo "Assembly job submitted: ${ASSEMBLYJOBID}"
	else
	    echo "Skipping Assembly step"
	    ASSEMBLYDEPENDENCY=""
	fi

fi 

### Running Mapping 

[[ $MAP_PPN -gt 40 ]] || [[ $MAP_PPN -lt 1 ]] && { echo "ERROR: 'MAP_PPN' must be between 1 and 40." ; exit 1; }

if [ "${DO_MAP}" = "true" ]; then

    readonly MAP_MAXMEM=$((130*$MAP_PPN/40))
    readonly MAP_RESOURCES="nodes=1:ppn=${MAP_PPN}:thinnode,walltime=24:00:00,mem=${MAP_MAXMEM}gb"
    readonly MAP_LOGS="-o log/pbs/MAP.pbs.o -e log/pbs/MAP.pbs.e"
    readonly MAP_ARGS="${SAMPLETABLE} ${MAP_PPN} ${MAP_MAXMEM} ${MAP_BAI} ${MAP_METHOD} ${MAP_REFERENCE} ${MAP_TO_SAME_REF}"
    MAPJOBID=$(myqsub -d `pwd` ${ASSEMBLYDEPENDENCY} -t ${RUN_FROM_NSAMPLES}-${RUN_NSAMPLES}%${NJOBS} ${MAP_LOGS} -l ${MAP_RESOURCES} -F "${MAP_ARGS}" ${SRC_DIR}/MAP.pbs  | cut -d . -f 1)
    MAPDEPENDENCY="-W depend=afterokarray:$MAPJOBID"

    echo "Mapping job submitted: ${MAPJOBID}"
else
    echo "Skipping Annotation step"
    MAPDEPENDENCY=""
fi




### Running Annotation

### Check NNodes input.
[[ $ANNOTATION_PPN -gt 40 ]] || [[ $ANNOTATION_PPN -lt 1 ]] && { echo "ERROR: 'ANNOTATION_PPN' must be between 1 and 40." ; exit 1; }

if [ "${DO_ANNOTATION}" = "true" ]; then

    readonly ANNOTATION_MAXMEM=$((130*$ANNOTATION_PPN/40))
    readonly ANNOTATION_RESOURCES="nodes=1:ppn=${ANNOTATION_PPN}:thinnode,walltime=80:00:00,mem=${ANNOTATION_MAXMEM}gb"
    readonly ANNOTATION_LOGS="-o log/pbs/ANNOTATION.pbs.o -e log/pbs/ANNOTATION.pbs.e"
    readonly ANNOTATION_ARGS="${SAMPLETABLE} ${ANNOTATION_TMPDIR} ${ANNOTATION_PPN} ${ANNOTATION_CONTIG_SUFFIX}"
    ANNOTATIONJOBID=$(myqsub -d `pwd` ${ASSEMBLYDEPENDENCY} -t ${RUN_FROM_NSAMPLES}-${RUN_NSAMPLES}%${NJOBS} ${ANNOTATION_LOGS} -l ${ANNOTATION_RESOURCES} -F "${ANNOTATION_ARGS}" ${SRC_DIR}/ANNOTATION.2.pbs  | cut -d . -f 1)
    ANNOTATIONDEPENDENCY="-W depend=afterokarray:$ANNOTATIONJOBID"

    echo "Annotation job submitted: ${ANNOTATIONJOBID}"
else
    echo "Skipping Annotation step"
    ANNOTATIONDEPENDENCY=""
fi





echo "Done with Pipeline" `date`
