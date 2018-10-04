# DRP Recommender

## Relevant Folders

1. DRP/recommender - holds the main recommender files
2. DRP/ml_models - holds the files necessary to train new models, store models and interface with WEKA
3. DRP/plugins - holds python scripts that calculate the reaction descriptors through ChemAxon
4. DRP/models - database model files
5. DRP/tests - holds files for Django test suite
6. DRP/management - holds Django files for defining Django terminal commands

## Overview

### Note: Make sure you are working on the rec_engine branch in git and github.

### Note: Relevant DRP database and setup is working only on Musk, so make sure you ssh into Musk. TODO: get DRP recommender running on other machines

### Note: Musk is setup to automatically commit through Chris's github account. TODO: remove Chris's account and make it so everyone has to manually enter username + password to commit

### Note: Testing is turned off on Musk. TODO: change this

To run the recommender and see output on Musk: type ```python3 manage.py recommend "Norquist Oxide Lab" --testing=1```
from the DRP/DRP directory. The recommender is currently only recommending 2 reactions since calculating the descriptors is costly.

## Recommender Outline

When the recommender is run it goes through several stages.

1. The recommender gets the averages of some base descriptor values (need to determine which base descriptors are necessary) and passes those as a grid parameter dictionary into a ReactionRecommender object. (recommend.py)

2. ReactionRecommender instantiates a ReactionGenerator and generates reasonable compound amounts for each compound in the database. (ReactionRecommender.py and ReactionGenerator.py)

3. ReactionGenerator selects a subset of triples to add to the grid parameter search. (ReactionGenerator.py)

4. ReactionGenerator uses said triples to generate new reactions using the compound amounts and grid parameter values. (ReactionGenerator.py) This step first creates Reaction objects in the database then for each of the Reactions it runs the descriptor calculator to populate the database.

5. ReactionSieve takes the id's of these generated reactions and runs them through our model to create a prediction and then promptly removes any reactions that do not comform to the desired_desc_dict passed into ReactionRecommender.

6.* After obtaining a subset of plausible reactions form the sieve, the recommender scores these reactions based on mutual information gain (mutual_information.py) and scales these based on the average distance to k-nn.

7. The recommender takes the top k scoring reactions and saves them to the database as the recommended reactions and displays said reaction ids and scores.

8. One can access the details of recommendations from database.

## Database

The DRP uses an SQL database to store the reaction data, descriptor data, and statistical model data.

To access the database directly:

1. open a terminal and type ```mysql -u root -p``` and hit enter.
2. Next enter the mySQL password (you can find it within the DRP/DRP/settings.py file)
3. Once within mysql, you need to type ```use DRP_2``` where DRP_2 is the name of the sql database holding the relevant reacition data.

Now you can execute SQL commands to explore the database. A good starting place is to type ```show tables;``` to see all the tables in the database.

TODO: add an in-depth description of relevant tables




