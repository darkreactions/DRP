Dark Reaction Project README
===========================

######Last Updated by Casey Falk -- 1/13/15

1. **Setup and General Information**
  1. Setup of the Server/a Test-bed
  2. Important Directories
  3. Django Directory Architecture
  4. Creating a User
  5. Connecting a User to BitBucket with an SSH Key
  6. Using a Test Bed ON the DRP Server
  7. Using a Test Bed OFF of the DRP Server
  8. Accessing the GitHub Repo and notes on the Repo Structure
2. **Django Management Commands**
  1. generate_model
3. **Database**
  1. Using the ORM
  2. Getting a CSV File of a Database Table
  3. Independent Processes
  4. Editing the Database Schema
  5. Database Backups
4. **Accounts**
5. **Machine Learning Models**
  1. Creating a Model
  2. Testing a Model
  3. Getting/Setting the "Current" Model


Setup and General Information
=============================

**Accessing the Server**

The website can be accessed at [darkreactions.haverford.edu](http://darkreactions.haverford.edu) -- this domain is managed by Haverford. The server itself (named "drp") is in the KINSC Server Room and can be accessed by SSH while on campus.  Note that if you are off-campus and need to access drp, you will need to tunnel through another server on campus -- such as those hosted by FIG or a CS Lab Computer. That is, SSH there and THEN SSH into drp.

**Setup of the Server/a Test-bed**

For a complete, step-by-step guide on this, please check out "setup.py". I leave the installation process for that process and will focus on the actual architecture of the system here. Also note that you can use `git clone https://github.com/cfalk/DRP.git` to clone the repo to your machine if you have been granted read-access to the the private GitHub repo.


**Important Directories**

There are a few directories that are vital to hosting the site.  Notably, these are the directories of NGINX and UWSGI -- respectively "/etc/nginx/" and "/etc/uwsgi/". NGINX is responsible for serving static files and distributing HTTP requests; UWSGI acts as the gateway from NGINX to the Django files (located in the home directory of the "drp" user: "/home/drp/web/darkreactions.haverford.edu/app/").


**Django Directory Architecture**

The DRP Django Project is set up as follows inside the DRP directory:

1. **./** *(The DRP Project)*
  - README	*-- This file.*
  - setup.txt	*-- Instructions for setting up the server.*
  - HowToCommitToTheGitRepo.txt	*-- Instructions for pairing an SSH Key with BitBucket.*
  - DRP_nginx	*-- The NGINX configuration. Move to /etc/nginx/sites-enabled/*
  - DRP_uwsgi	*-- The UWSGI configuration. Move to /etc/uwsgi/apps-enabled/*
  - ChemAxonLicense.cxl	*-- A ChemAxon license file that should go in /<home>/.chemaxon/*
  - manage.py	*-- The Django management file -- leave it as it is.*

1. **DRP/** *(The DRP App)*
  - retrievalFunctions.py	*-- Contains functions for retrieving database entries.*
  - database_construction.py	*-- "                    "  adding database entries.*
  - models.py		*-- Contains all the Django models and some important accessors.*
  - settings.py	*-- Contains settings for the database, admin emails, and more.*
  - urls.py		*-- Responsible for mapping a URL call to a function call.*
  - data_config.py	*-- Sets the license file, some static directories, and more.*
  - [more files].py	*-- The logging, email, and view helper functions (etc.).*

  - research/		*-- Research files incorporated into the DRP scripts.*
  - management/	*-- Directory to add custom `python manage.py` commands.*
  - templatetags/	*-- Directory to add custom Django template tags.*
  - views/		*-- Directory of various views sorted into different files.*
  - compound_calculations/
  - migrations/		*-- South Migration Files -- don't modify manually.*
  - model_building/		*-- Scripts for building the model itself.*
  - recommendation/		*-- Scripts for calculating and storing recommendations.*
2. **logs** *-- (Directory for the log files of worker processes.)*
  - compound_calculations/ 	*-- The worker process logs for compound property calcs.*
  - seed_recommend/		* -- logs for seed recommendations.*
3. **templates** *(The HTML templates for Django Views.)*
4. **research** *(Any "research" scripts that are still being explored.)*
5. **static** *(Any static files for which Django can skip templating.)*
  - favicon.ico	*-- The "favicon" for the site (actually served by NGINX).*

  1. js/	*-- Any Javascript for any view belongs in this directory.*
       universal.js *-- Contains Javascript used on nearly every page.*
  2.  css/	*-- The CSS for any view should go here.*
       universal.css *-- Contains styling used on nearly every page.*
  3.  icons/	*-- Small "icon" images belong here -- such as the hover-button images.*
  4.  images/	*-- Any large images for the site should go here (eg: the logo).*
  5.  licenses/	*-- The current and past licenses.*
  6.  admin/	*-- CSS for the /admin/ page. Feel free to ignore.*
6. **tmp**		*-- Calculation files for compounds. Don't modify manually.*
7. **test_scripts** *-- Scripts that test site status (etc).*

There are many files that are not listed above in order to elucidate the "framework" of the DRP Django Project succinctly. Notably, there are many python files in the views and research directories that are not listed above -- but which should contain explanatory comments in the files themselves.


**Creating a User**

After logging in, you'll want to make another user. The command for
that is pretty simple: "sudo adduser <theNewUserName>". It should
prompt you then to make a password and for some other details
(which you can skip).

Next, give the account sudo access; to do so, all we need to do
is add the new user to the "sudo" group. Run the command:
"sudo adduser <theNewUserName> sudo". This change will take effect
the next time the new user logs in. Hurray!


**Connecting a User to BitBucket with an SSH Key**

To start, hop over to BitBucket and log in/make an account. Then, follow the instructions in ".../DRP/HowToCommitToTheGitRepo.txt" to create and pair your development user with the BitBucket Repo. Note that to push, you'll still need to be added to the BitBucket Organization.


**Using a Test Bed ON the DRP Server**

This should only be used as a last-resort when you cannot develop
locally -- as this can endanger the integrity and processes running
on the production server.

Another disclaimer: note that any test bed should be placed in the
home directory of your own user and should always use the test database.
Only connect to the actual database if you need the most up-to-date
database version -- and always make sure to back up the database before
you perform any script that may modify the database in any way.

With those warnings in place, setting up a test bed is simple. First,
"cd" to your home directory and check BitBucket for the repository URL
that you can use to copy the repository to your development directory with
the command: "git clone <theOverlyLongURLForTheGitRepo>".

There are a few settings you must change in ".../DRP/settings.py"
(notably the database password and the database name "DRP_db_test" is the
test database, whereas "DRP_db" is the actual database.). Remember, data
is gold.

Django provides a nifty test-server through the command:
`python manage.py runserver <ip>:<port>`. This "runserver" serves static
files and operates slowly -- but allows rapid prototyping and development.
Definitely use this in development rather than NGINX and UWSGI, as it will
remove many obstacles and allow you to work faster and more efficiently.
Check more out online: https://docs.djangoproject.com/en/dev/ref/django-admin/


**Using a Test Bed OFF of the DRP Server**

This is CERTAINLY the recommended development strategy. Note that the
same process as described above can be used to set up the Django Project.
However, you'll want to SCP a version of the database over to your
development box and set up MySQL appropriately (see the setup.txt file above).
Any of the backups in the DropBox backup folder should suffice.

If you do not set up ChemAxon or WEKA, the ConfigManager in `data_config.py` will
throw errors when you try to start the project. This error-checking can be disabled
by using `validate_config=False` instead of `True` in [data_config.py](https://github.com/cfalk/DRP/blob/master/DRP/data_config.py#L18).


**Accessing the GitHub Repo and Notes on the Repo Structure**

Firstly, you'll need a GitHub account and you'll need someone with access to
the repo to grant your account access (though if you can view this README without
access to the GitHub repo, you should tell someone). Then, you should be able to
use `git clone https://github.com/cfalk/DRP.git` to copy the repository to
your workstation. Note that if you have a branch set up using the old BitBucket,
you'll want to delete that section in your ".../.git/config" file and use
`git remote add <branch> https://github.com/cfalk/DRP.git` to transition smoothly
to the new repo.

Django has an important file, settings.py, which is not tracked by this git
repository for security reasons. If you add or remove items in this file
whilst developing with the code, please ensure that you update the
file settings_example.py in step.

Lastly, the repository utilises branches heavily. At present:
    
    - Release-1.0 is for bugfixing our next stable release version.
    - dev contains our development versions of code.
    - master has been frozen excepting documentation changes.

The branching scheme follows recommendations made [here](http://nvie.com/posts/a-successful-git-branching-model/).

Django Management Commands
=========================

Django management commands (that is, the commands that pass through
manage.py) can be called as `python manage.py <the command>` when
in the main DRP project directory. These commands are each stored in
a separate python file in the .../DRP/management/commands directory. To
add a new command, use the existing files as examples.

It also should be noted that using `python manage.py help` will detail
all of the available management commands (of which there are many).

**generate_model**

To generate a new model using the current pipeline, run:
`python manage.py generate_model '<Model_Title>' '<An explanatory description.>'`
The script encompasses the model construction, training, and testing
stages and ultimately creates a new ModelStats entry in the Django database.
This entry is complete with "stats" and a reference to the model name,
description, and unique ID. The estimated run time is ~20 minutes.


Database
=======

The DRP uses a MySQL database that should be accessed by running
MySQL queries using the "root" MySQL user on the "DRP_db" database.
However, that being said, these queries should be performed by
utilizing the [Django ORM](https://docs.djangoproject.com/en/dev/topics/db/)).


** Using the ORM **

The ORM can best be thought of as a generic wrapper that will
translate Model methods into database queries (for example, MySQL
wrappers). All database retrieval functions should be placed in
"retrievalFunctions.py" and any database construction functions
should be placed in "database_construction.py" so that they can
be easily imported into Python scripts. For example, to import
all of the retrieval and construction functions, just use:

    from DRP.retrievalFunctions import *
    from DRP.database_construction import *

And now, to query the database for all of the reactions (which
are stored as the Data model (aka Python "class" in models.py),
now just use:

    Data.objects.all()

For more complicated queries, check out the examples in
"retrievalFunctions.py" -- and likewise for information on how to
create new entries in the database, check out "database_construction.py".
The ORM is flexible and capable of any query you should be making -- and
allows us to abstract complicated MySQL queries out of our scripts.


** Getting a CSV File of a Database Table **

The best way to interact with the database is through the Python functions
designed to retrieve, "expand", and send the data. However, there may be
times when you wish to download a CSV; while the "simple" CSV (with the
original input features) can be downloaded through the site interface, the
"expanded" CSV is only accessible on the backend. For convenience, the
"fileFunctions.py" python file contains a `writeExpandedCSV(filename)`
function that will write a full, expanded CSV to the system. For even more
convenience, this function is wrapped by a management command:

> python manage.py writeCSV filename.csv

If you need to access the data in a CSV-like format through a library
such as d3.js, it is recommended that you use a Django view and
**stream** the data in CSV format directly where you need it.
For an example of using a view to download a CSV, check out
".../DRP/views/download.py".

Note that the same process can and should be used to get JSON-formatted
entries from the database. Find examples in ".../DRP/views/jsonViews.py".


**Independent Processes**

Often, you may want to run independent Python processes that access
the database. This requires the Python Path to be set to the DRP
project directory so that Django can include the appropriate files.
To do so, just include the following lines to in top of the script:

    import os, sys
    full_path = os.path.dirname(os.path.realpath(__file__))+"/"
    django_path = full_path[:full_path.rfind("/DRP/")]

    if django_path not in sys.path:
      sys.path = [django_path] + sys.path
      os.environ['DJANGO_SETTINGS_MODULE'] = 'DRP.settings'


**Editing the Database Schema **

Be wary about adding fields to the Django Models and adding new Models
in general; there are a few steps that must be completed. The DRP
uses a database migration tool called South which handles most
schema migrations. First, modify the Django Model (in models.py) and save
the file. Second, run the South command that detects model changes:
`python manage.py schemamigration DRP --auto`. Finally, perform
the migration: `python manage.py migrate DRP`.

Note that if this is your first migration, you may need to "fake" an initial
migration so South knows what fields to keep and what fields to change.

Functionally, South re-creates any database tables that are modified and
migrates the data in the old table over to the new table. There is more that
goes on in the back-end of South (such as versioning and reverse-update
control), but I leave that to the [South Documentation.](http://south.aeracode.org/).

***NOTE: South was EOL-ed since Django 1.7 will have [Migrations](https://docs.djangoproject.com/en/dev/topics/migrations/) built in. When we transition to 1.7, we'll most likely no longer need South.***


**Database Backups**

The database is dumped to a file and placed in a DropBox folder
every night. The DropBox folder can be found on the drp server at
"/home/drp/database_backups/Dropbox/". The DropBox account
uses the email "darkreactionsproject@gmail.com" and the standard
password. Any database dumps that are performed should be placed
in this folder to avoid clutter.


Accounts
========

There are a few accounts that are needed by the DRP to operate
successfully -- they range from database accessors to mail serves.

*Google:* darkreactionsproject@gmail.com -- Used to distribute automatic emails
from the DRP site, such as error messages, contact form completions,
and seed recommendation results.

*ChemSpider:* darkreactionsproject@gmail.com -- Used to allow developer access
to the ChemSpider database; we query ChemSpider to get reactant
images, molecular weights, and SMILES strings.

*ChemAxon:* darkreactionsproject@gmail.com -- Used to perform chemical calculations
on various reactants.

*BitBucket:* "Dark Reaction Project" Team -- The BitBucket "team" (aka "organization") that controls the DRP site repository. Each developer should have
a personal BitBucket account with read/write access to this team.

[Markdown Cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


Machine Learning Models
========

The machine learning models ("ML models") are the calculated models that predict a set of responses from a set of observations; for DRP, the ML models predict the "outcome" of a reaction based on the reaction's expanded feature-set. These models are *not* the same as the Django models that can be found and edited in the `.../DRP/models` directory; however, you *should* use the `ModelStats` Django model when interacting with the ML models.


**Creating a Model**

If you want to create a model, use the `generate_model` management command (refer to the section on the `generate_model` management command for details). Feature usage or modifications should be handled in the pre-/post-processing functions, **however**, do not edit such functions that already exist as this will cause confusion when looking at previous models. Instead, create a new pre-/post-processing function. Likewise, if you are making substantial changes to the splitting method, just create a new function. To "choose" which functions are used in the model-generation, change the respective import in [.../DRP/model_building/generate_models.py](https://github.com/cfalk/DRP/blob/master/DRP/model_building/generate_models.py#L152-L154) before running `generate_model`. Please document all of your model-generation functions, especially.


**Testing a Model**

It is recommended that you test the model using Django's Python Shell. This can be done by calling `python manage.py shell` in the main project directory. Then, import the `ModelStats` objects using `from DRP.models import ModelStats`; this allows us to access all of the ML model information stored in the database. You can then use Django's ORM to retrieve the ModelStats object for the ML model that you want to test. For example, to grab the latest ML model, you could just use `model = ModelStats.objects.last()`. Finally, pass the data you want to test against to the `_test_model` method of the `ModelStats` object: `model._test_model(test_data, all_data)`. Note that the `test_data` should have been pre-processed and post-processed into a matrix (a list of observations) with features in the same order as the model's stored headers; to check what headers were used in a model, print the result of `model.get_headers()`. To see the results of the test, call `model.summary()` and look for the "test" section.

Note that one current limitation of this implementation is that *all* of the data must be present (in `all_data`, above) for WEKA to know all possible values of a given feature; what we **should** do is store these feature-value sets in the `ModelStats` object itself -- though this has not been implemented yet. This makes retro-generation of a model slightly more involved. Another limitation is that only the latest test is stored in the database for any given model.

For a complete (though somewhat more complex than necessary) example, see [here](https://github.com/cfalk/DRP/blob/master/DRP/research/casey/temp_test_models.py).


**Getting/Setting the "Current" Model**

To get the current model, simply query for the latest "active" model: `ModelStats.objects.filter(active=True).order_by("end_time").last()`. In that query, the `order_by("end_time")` makes sure we are sorting the list in chronological order based on the time that model-generation was completed (the "end time"). Likewise, to "set" a model, just make it the latest active model.

Currently, if you want to set an old model to be the latest "active" model, you'd need to modify the "end_time" of that ModelStats object in the database. In the future, we should just add an additional feature to the `ModelStats` Django model so that we can just toggle any ML model "on" or "off" (with only one model being "on" at any time), but as of yet this has not been implemented.

To re-create the current model, you should only need to call `generate_model` using the default pre-/post-processors and the default split on the entire dataset. Make sure to double-check that no [research filters](https://github.com/cfalk/DRP/blob/master/DRP/model_building/generate_models.py#L11) are being applied to the dataset; specifically, make sure no data is being removed in the `research_data_filter` function in `.../DRP/model_building/generate_models.py` and that [all valid data](https://github.com/cfalk/DRP/blob/master/DRP/model_building/generate_models.py#L139-L142) is being used in `gen_model`.
