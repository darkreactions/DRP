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
2. **Django Management Commands**
  1. generate_model
3. **Database**
  1. Using the ORM
  2. Getting a CSV File of a Database Table
  3. Independent Processes
  4. Editing the Database Schema
  5. Database Backups
4. **Accounts**

Setup and General Information
=============================

**Accessing the Server**

The website can be accessed at [darkreactions.haverford.edu](http://darkreactions.haverford.edu) -- this domain is managed by Haverford. The server itself (named "drp") is in the KINSC Server Room and can be accessed by SSH while on campus.  Note that if you are off-campus and need to access drp, you will need to tunnel through another server on campus -- such as those hosted by FIG or a CS Lab Computer. That is, SSH there and THEN SSH into drp.

**Setup of the Server/a Test-bed**

For a complete, step-by-step guide on this, please check out "setup.py". I leave the installation process for that process and will focus on the actual architecture of the system here.


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

The best way to do this is to use the site interface, but you may need
to load the data as a CSV temporarily for use with a library such as D3.
To do so, use a Django view and stream the data directly in CSV format.
The benefits of this are reduced file-clutter, the most up-to-date version
of the database, and improved security (since we can enforce only users
that are authenticated can access the database). For an example of using
a view to download a CSV, check out ".../DRP/views/download.py".

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
