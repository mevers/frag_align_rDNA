#!/bin/bash

# Create DAGs
#snakemake --dag | dot -Tpng > workflow_dag.png
#snakemake --rulegraph | dot -Tpng > workflow_rulegraph.png

# Execute workflow
# To debug the DAG use `--debug-dag`
snakemake --use-conda --cores
