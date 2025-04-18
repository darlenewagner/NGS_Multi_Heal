## Requires Biopython in a virtual environment
#### After git clone of NGS_Multi_Heal set up virtual environment in BiopythonBased folder:
##### HPC 'module load' shown

`module load python/3.7.4`

`cd NGS_Multi_Heal/`

`virtualenv -p /apps/x86_64/Python/3.7.4-GCCcore-8.3.0/bin/python3 BiopythonBased/`

##### Note: virtualenv command will vary according to python3 executable path

`cd BiopythonBased/`

`bin/pip install biopython`

`bin/python readShuffledFastq.py ../ExampleData/Viral/Pool-1_S1_adenovirus_B3_001.fastq`

