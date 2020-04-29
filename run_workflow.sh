#!/bin/bash

# Create DAGs
#snakemake --dag | dot -Tpng > workflow_dag.png
#snakemake --rulegraph | dot -Tpng > workflow_rulegraph.png

# Execute workflow
snakemake --use-conda --debug-dag
