Dark Reaction Project README
===========================

######Last Updated by Casey Falk -- 1/13/15

1. **Setup and General Information**
  1. Accessing the Main Server
  2. Setup of the Server/a Test-bed
  3. Important Directories
  3. Django Directory Architecture
  4. Creating a User
  5. Connecting a User to BitBucket with an SSH Key
  6. Using a Test Bed ON the DRP Server
  7. Using a Test Bed OFF of the DRP Server
  8. Accessing the GitHub Repo and notes on the Repo Structure
2. **Django Management Commands**
  1. check\_hash\_collisions
  2. import\_data
  3. port\_database
  4. re\_save\_reactions
  5. run\_tests
3. **Accounts**
4. **Machine Learning Models**

Setup and General Information
=============================

**Accessing the Server**

The website can be accessed at [darkreactions.haverford.edu](http://darkreactions.haverford.edu) -- this domain is managed by Haverford. The server itself (named "drp") is in the KINSC Server Room and can be accessed by SSH while on campus.  Note that if you are off-campus and need to access drp, you will need to tunnel through another server on campus -- such as those hosted by FIG or a CS Lab Computer. That is, SSH there and THEN SSH into drp.

**Setup of the Server/a Test-bed**

For a complete, step-by-step guide on this, please check out "setup.md". I leave the installation process for that process and will focus on the actual architecture of the system here. Also note that you can use `git clone https://github.com/cfalk/DRP.git` to clone the repo to your machine if you have been granted read-access to the the private GitHub repo.

This README assumes more than a passing familiarity with the Django framework

**Django Directory Architecture**

The DRP Django Project is set up as follows inside the DRP directory:

1. **./** *(The DRP Project)*
  - README	*-- This file.*
  - setup.txt	*-- Instructions for setting up the server.*
  - HowToCommitToTheGitRepo.txt	*-- Instructions for pairing an SSH Key with BitBucket.*
  - DRP_nginx	*-- The NGINX configuration. Move to /etc/nginx/sites-enabled/*
  - DRP_uwsgi	*-- The UWSGI configuration. Move to /etc/uwsgi/apps-enabled/*
  - manage.py	*-- The Django management file -- leave it as it is.*

1. **DRP/** *(The DRP App)*
  - models/		*-- Contains all the Django models and some important accessors.*
  - settings_example.py	*-- Contains example settings for the database; true settings are placed in settings.py (see setup.md).*
  - urls/		*-- Responsible for mapping a URL call to a function call; see the django documentation for details.*
  - research/		*-- A folder for very experimental/highly unstable code.*
  - management/	*-- Directory to add custom `python manage.py` commands- see the django documentation.*
  - templatetags/	*-- Directory to add custom Django template tags- see the django documentation.*
  - views/		*-- Directory of various views sorted into different files.*
  - migrations/		*-- South Migration Files; see South documentation for details.*
  - model_building/		*-- Scripts for building the model itself.*
  - recommendation/		*-- Scripts for calculating and storing recommendations.*
2. **logs** *-- (Directory for the log files of worker processes.)*
3. **templates** *(The HTML templates for Django Views.)*
4. **research** *(Any "research" scripts that are still being explored.)*
5. **static** *(Any static files for which Django can skip templating.)*
  - favicon.ico	*-- The "favicon" for the site (actually served by NGINX).*

  1. js/	*-- Any Javascript for any view belongs in this directory.*
       universal.js *-- Contains Javascript used on nearly every page.*
  2.  css/	*-- The CSS for any view should go here.*
  3.  icons/	*-- Small "icon" images belong here -- such as the hover-button images.*
  4.  images/	*-- Any large images for the site should go here (eg: the logo).*
  6.  admin/	*-- CSS for the /admin/ page. Used in many django installations.*

There are many files that are not listed above in order to elucidate the "framework" of the DRP Django Project succinctly. Notably, there are many python files in the views and research directories that are not listed above -- but which should contain explanatory comments in the files themselves.

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

Lastly, the repository utilises branches heavily. The master branch
is reserved for releases. There are no working long-term support branches
at present. Persons editing this code should set up their own branches,
with one marked as stable, such that all tests pass in that branch.    

Porting to Django 1.8
=============================

Pull django1.6 tag:
git pull origin django1.6

Make sure your database (both main and testing) is up to date with migrations:
Set Testing = False in settings.py
./manage.py migrate
Set Testing = True in settings.py

Restart nginx and uwsgi
sudo service nginx restart && sudo service uwsgi restart

Run all tests and make sure you pass them all:
./manage.py run_tests

Pull django1.8 tag:
git pull origin django1.8

Update your settings:
Remove "south" from your installed apps

Install django 1.8:
pip install -U Django==1.8.9

Do initial migrations:
Delete all .pyc files in your migrations folder
rm DRP/migrations/*.pyc
./manage.py migrate --fake-initial


Restart nginx and uwsgi
sudo service nginx restart && sudo service uwsgi restart

Run all tests and make sure you pass them all:
./manage.py run_tests

Migrate your main database:
Set Testing to False in settings.py
./manage.py migrate --fake-initial

Uninstall south. You don't need it anymore:
pip uninstall south

Bake Geoffrey a cake because what just took you 20 minutes took him 12 hours to make sure everything worked properly


Django Management Commands
=========================

Django management commands (that is, the commands that pass through
manage.py) can be called as `python manage.py <the command>` when
in the main DRP project directory. These commands are each stored in
a separate python file in the .../DRP/management/commands directory. To
add a new command, use the existing files as examples.

It also should be noted that using `python manage.py help` will detail
all of the available management commands (of which there are many).

**check\_hash\_collisions**

This command checks the hashed values of the reactions calculated by
the DRP descriptor plugin for clashes. If there is a clash,
the lead developer should be notified ASAP.

**import_data**

Imports reaction data from the main Haverford Dark Reactions server.
Accepts one positional argument importing the corresponding number
of Performed Reactions and their associated data. Not available to
persons who are not administrators on the main server.

Note that at present this command does not import reaction descriptors.

**re_save_reactions**

Starts a batch parralell task to re-save each reaction, forcing
descriptor calculation. Useful in conjunction with the above
command. Note that descriptors will only be calculated
for plugins present in the code-base, which are correctly
set up.

**run_tests**

Runs all of the tests correctly imported in the test suite.

