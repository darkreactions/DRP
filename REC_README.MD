# DRP Recommender

## Relevant Folders

1. DRP/recommender - holds the main recommender files
2. DRP/ml_models - holds the files necessary to train new models, store models and interface with WEKA
3. DRP/plugins - holds python scripts that calculate the reaction descriptors through ChemAxon
4. DRP/models - database model files
5. DRP/tests - holds files for Django test suite
6. DRP/management - holds Django files for defining Django terminal commands

## Overview

# Note: Make sure you are working on the rec_engine branch in git and github.

# Note: Relevant DRP database and setup is working only on Musk, so make sure you ssh into Musk. TODO: get DRP recommender running on other machines

# Note: Musk is setup to automatically commit through Chris's github account. TODO: remove Chris's account and make it so everyone has to manually enter username + password to commit

To run the recommender and see output on Musk: type ```python3 manage.py recommend "Norquist Oxide Lab" --testing=1```
from the DRP/DRP directory. The recommender is currently only recommending 2 reactions since calculating the descriptors is costly.



