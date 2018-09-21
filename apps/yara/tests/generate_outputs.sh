#!/bin/bash
#
# (DREAM-)Yara output generation.

EXEC_PATH=~/workspace/development/builds/dream_yara/debug/bin
YARA_INDEXER=$EXEC_PATH/yara_indexer
YARA_MAPPER=$EXEC_PATH/yara_mapper

DREAM_YARA_INDEXER=$EXEC_PATH/dream_yara_indexer
DREAM_YARA_MAPPER=$EXEC_PATH/dream_yara_mapper
IBF_INDEXER=$EXEC_PATH/dream_yara_build_filter


# ============================================================
# Run Yara Indexer
# ============================================================

# Run with different organisms.
for organism in 64-viral; do
    ${YARA_INDEXER} input/$organism-genomes.fa -o gold/$organism-genomes
done

# ============================================================
# Run DREAM-Yara indexer and IBF
# ============================================================

# Run with different organisms.
for organism in 64-viral; do
    python manage_bins.py split $organism;
    mkdir gold/$organism-binned-indices/;
    ${DREAM_YARA_INDEXER} input/$organism-binned-genomes/* -o gold/$organism-binned-indices/;
    # do not creat IBF it is too big
    # ${IBF_INDEXER} -b 64 -t 4 -k 19 -nh 2 -bs 1 input/$organism-binned-genomes/ -o gold/$organism-binned-genomes.filter;
    python manage_bins.py clear $organism;
done

# ============================================================
# Run Single-End Yara Mapper
# ============================================================

# Run with different threads.
MAPPER_ARGS=( "-e 3 --threads 1" "-e 3 --threads 1 -sm record -s 10" "-e 3 --threads 1 -sm tag -s 10" "-e 3 --threads 1 -sm record -as -s 10" "-e 3 --threads 1 -sm tag -as -s 10" ) # "--threads 8")
MAPPER_SUFFIX=( "t1" "rec.t1" "tag.t1" "as.rec.t1" "as.tag.t1" ) # "t8")

# Run with different organismnisms.
for organism in 64-viral; do
    # Run with different arguments.
    for ((a=0; a<${#MAPPER_ARGS[@]}; a++)); do
        ${YARA_MAPPER} gold/$organism-genomes input/$organism-reads.fa -o gold/$organism-reads.${MAPPER_SUFFIX[$a]}.sam ${MAPPER_ARGS[$a]}
    done
done


# ============================================================
# Run Single-End DREAM-Yara Mapper
# ============================================================
# Run with different threads.

MAPPER_ARGS=("-e 3 --threads 1" "-e 3 --threads 1 -sm record -s 10" "-e 3 --threads 1 -sm tag -s 10") # "--threads 8")
MAPPER_SUFFIX=("t1-dis" "rec.t1-dis" "tag.t1-dis") # "t8")

# Run with different organismnisms.
for organism in 64-viral; do
    # Run with different arguments.
    for ((a=0; a<${#MAPPER_ARGS[@]}; a++)); do
        ${DREAM_YARA_MAPPER} -fi gold/$organism-binned-genomes.filter gold/$organism-binned-indices/ input/$organism-reads.fa -o gold/$organism-reads.${MAPPER_SUFFIX[$a]}.sam ${MAPPER_ARGS[$a]}
    done
done
