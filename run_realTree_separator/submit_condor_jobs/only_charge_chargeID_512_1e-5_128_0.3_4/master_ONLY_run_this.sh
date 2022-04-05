eval "$(/afs/desy.de/user/a/axelheim/miniconda3/bin/conda shell.bash hook)" 
conda activate baum
jupyter nbconvert --to python ~/private/baumbauen/notebooks/30.6_ah_simpleDecayTree_separator.ipynb 


condor_submit trainingGPU.submit


sleep 5


#rm ~/private/baumbauen/notebooks/30.6_ah_simpleDecayTree_separator.py
