#import logging
import sys
from pathlib import Path

#import click
from IPython.core import ultratb

import torch
import ignite as ig
#import baumbauen as bb
#import yaml
import optuna
import datetime


def objective(trial):   
    


    
    nhid = trial.suggest_categorical("n_hid", [128,256,512])
    n_blocks = trial.suggest_categorical("n_blocks", [2,4,6,8,10])
    drop  = 0.0 #trial.suggest_categorical("drop", [0.1,0.2,0.3])
    print("nhid:",nhid,"n_blocks:",n_blocks,"drop:",drop)
    


    return int(str(datetime.datetime.now())[-1])

study_name = "optunaParallelTrial"
#study_file = Path("/afs/desy.de/user/a/axelheim/private/HTCondor_interactive_optuna/" + study_name)

study = optuna.create_study(
    direction="maximize",
    #pruner=optuna.pruners.MedianPruner(n_warmup_steps=15),
    study_name=study_name#,
    #storage=study_file,
    #load_if_exists=True,
)

study.optimize(
    objective,
    n_trials=90,
    # timeout=configs['train']['optuna']['timeout'],
    # gc_after_trial=True,
)

print("Number of finished trials: ", len(study.trials), "\n")

trial = study.best_trial
print("Best trial:", trial.value, "\n")

print("  Params: ")
for key, value in trial.params.items():
    print("    {}: {} \n".format(key, value))


print(study.best_trial)
