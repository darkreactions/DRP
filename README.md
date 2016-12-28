Dark Reaction Project README
===========================

######Last Updated by Philip Adler 27 June 2016

General Information
=============================

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

Always make sure to check the latest version of the README on the master branch.

This repository contains the software for the [https://www.djangoproject.com/](Django)-based source code for the Dark Reactions Project Software. If you are looking to contribute to the chemistry aspects of the project, please visit the main project site at [http://darkreactions.haverford.edu](http://darkreactions.haverford.edu). If you are looking to contribute to the source code of the project __and are not a member of haverford college__, please fork this repository and issue a pull request with any changes or fixes you may have made. A list of known bugs can be found at [http://bugs.darkreactions.haverford.edu](our instance of Mantis Bug Tracker). Please note that you will to sign up for an account, and that the authentication credentials for the bug reporting page and the main project page are separate.

Setting up your own instance of the DRP
=============================

The following instructions are written to work with Ubuntu 14 and have (mostly) been tested. These instructions assume familiarity with Linux and a Command Line, and that you are using nginx as your webserver.

Install the necessary programs.

`sudo apt-get install python3-dev python3-pip mysql-server libmysqlclient-dev nginx uwsgi-plugin-python3 git weka graphviz memcached python-memcache mailutils python3-scipy python3-Pillow`

`sudo pip3 install numpy pygraphviz mysqlclient sqlparse`

Install Django. The current version of DRP is designed to work with Django 1.8

`sudo pip3 install django==1.8`

Install required pip python libraries
`sudo pip3 install chemspipy requests pep8 pep257 xxhash`

####For RDKit Descriptors

If you wish to use calculations from RDKit as a part of your install, the following is also necessary *these installation instructions are ubuntu 14.04 specific*.

`sudo apt-get install bison cmake flux build-essential sqlite3 libsqlite3-dev libboost-all-dev`

`sudo pip3 install cairocffi`

Then, in a directory *that is not* your main installation directory for git.

`git clone https://github.com/shadowadler/rdkit.git`

Change into the rdkit repository directory and then

`export RDBASE=$(pwd)`

`export LD_LIBRARY_PATH="$(pwd)/lib"`

`export PYTHONPATH="$(pwd)/lib"

`mkdir build`

`cd build`

`cmake -DRDK_BUILD_INCHI_SUPPORT=ON -D PYTHON_LIBRARY=/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4.so -D PYTHON_INCLUDE_DIR=/usr/include/python3.4/ -D PYTHON_EXECUTABLE=/usr/bin/python3.4 -DBOOST_ROOT=/usr/lib/x86_64-linux-gnu/ ..`

`make install`

`ctest`

If any of this generates an error, unless you are *very* familiar with compiling new code for linux operating systems, seek the assistance of Philip Adler via the rdkit repository provided in this document.

Otherwise:

`unset LD_LIBRARY_PATH`
`unset RDBASE`
`unset PYTHONPATH`

`cd ../rdkit`

`sudo ln -s "$(pwd)" /usr/lib/python3.4/rdkit

`cd ../lib`

`sudo cp -i *.so.2 /usr/lib`

`python3.4 -c "import rdkit.Chem"`

If that runs without error messages, congratulations, you have compiled and installed rdkit for use with DRP.

###Clone from the Git Repository into your directory of choice.

`git clone git@github.com:darkreactions/DRP`

###Release Versions

DRP is distributed in release versions. To use a specific version of the code, use the following template command:

`git checkout <version number>`

###Git Hooks

DRP comes distributed with a number of useful git hooks in the drp\_hooks directory in this repo. These warn you if the expected structure of the settings.py file changes, or if you need to run database migrations for the DRP Django application. They also remove orphaned .pyc files, which have been known to confuse the test suite historically. These should be added to your local git repository as per the git documentation.

####Server settings

In the DRP repository there is a file DRP\_uwsgi.ini and another DRP\_nginx. Both should be modified to suit your local server *after* having been placed in the relevant locations:

`/etc/uwsgi/apps-enabled/DRP_uwsgi.ini`
`/etc/nginx/sites-enabled/DRP_nginx`

It should be noted that the uwsgi.ini is backwards compatible with older version of this repo, but that an old DRP_uwsgi file will need replacing.

Both uwsgi and nginx must be restarted (in that order) for the server to work.

###Set up the settings.py file

In DRP/DRP, there is a file called 'settings\_example.py'. This must be copied to 'settings.py', and the settings therein set to the appropriate values for your server. At present, the available fields should be fairly self explanatory, though the following should be noted:

`ALLOWED_HOSTS` should be set to an iterable containing only the element '\*'.

To pass the unit tests, at least one `ADMIN_EMAILS` should be provided

To pass the unit tests, the EMAIL\_HOST\_USER and related settings should be set.  

###https

If you are setting up a publicly viewable instance of DRP, there are additional settings for these protocols present in the configuration file for nginx, and the Django settings_example.py file, which have been left commented out. Setting up https access varies greatly depending on your local server environment and organisation so there will be no further documentation here. Additional information for a simple method to set up https can be found at [https://letsencrypt.org](https://letsencrypt.org).

###Running tests

In order to run tests you must have the following environment variables set up in your shell session:

export PYTHONPATH=/path/to/DRP/

export DJANGO_SETTINGS_MODULE=DRP.settings

You must also have `TESTING` set to `True` in your `settings.py` file.

To run specific tests, provided the test is conformant to the template test (which they should be if you are writing new tests!), one can simply execute the test:

`path/to/DRP/test.py`

Else, one can run the entire test suite from the management script:

`./manage.py run_tests`

###On Development Servers

The `ALLOWED_HOSTS` option in 'settings.py' should be set to an iterable containing only the localhost ip address as a string.

In the '/etc/nginx/sites-enabled/DRP\_nginx' file, the host name that is being listened to should only be localhost.

##Servers with multiple web applications.

If you are only developing DRP on your server, the setup you have should be sufficient, however, people running other applications on their local development server should note the following.

If you are running the django testing server, this requires you to select a port which is unoccupied. By default, the nginx settings file listens for port 8000, which is also the default port of the django test server; you will need to configure one or the other so that this clash does not occur. The Django documentation addresses this for django, whilst in the DRP\_nginx file, the only change that needs to be made is to delete the line:

`listen		8000`

 #dnsmasq

For instances where you are hosting multiple development projects on your local server, it may be beneficial to install dnsmasq:

`sudo apt-get install dnsmasq`

dnsmasq is a powerful tool for rerouting and managing dns requests. This makes it extremely helpful in managing multiple local development projects.

Having installed dnsmasq, open the file `/etc/dnsmasq.conf` in your favourite text editor, and add the following line into the file:

`address=/loc/127.0.0.1`

Save the change, and then on the command line:

`sudo service dnsmasq restart`

In the DRP\_nginx file, change the `server_name` configuration to something like `darkreactions.loc`. It does not matter what this is set to, provided it:

a. is unique on your development server
b. ends in `.loc`

Don't forget to set the `SERVER_NAME` setting in your settings.py file to the same value!

Restart nginx:

`sudo service nginx restart`

When you open your browser and direct yourself to darkreactions.loc (or whatever you named the server), the dark reactions project should display.

###DRP Versions < 0.1

Versions of code prior to 0.1 are not compatible with versions above.

**Upgrading from versions prior to 0.7**

Version 0.6 of DRP was the last version to use Django 1.6, subsequent versions use Django 1.8. There are, therefore, necessary transition steps to be made. 

`git fetch --all`
`git checkout 0.7`

Make sure your database (both main and testing) is up to date with migrations:
`Set Testing = False in settings.py`
`./manage.py migrate`

Set Testing = True in settings.py

Restart nginx and uwsgi
`sudo service nginx restart && sudo service uwsgi restart`

Run all tests and make sure you pass them all:
`./manage.py run_tests`

git checkout 0.81

Remove "south" from your installed apps

Install django 1.8:
pip install -U Django==1.8.9

Delete all .pyc files in your migrations folder
`rm DRP/migrations/*.pyc`
`./manage.py migrate DRP --fake`
`./manage.py migrate --fake-initial`
`./manage.py migrate`

Repeat this process for your Main database.

Restart nginx and uwsgi
sudo service nginx restart && sudo service uwsgi restart

Run all tests and make sure you pass them all:
`./manage.py run_tests`

**Notes for Local (Haverford) Developers**

The website can be accessed at [darkreactions.haverford.edu](http://darkreactions.haverford.edu) -- this domain is managed by Haverford. The server itself (named "drp") is in the KINSC Server Room and can be accessed by SSH while on campus.  Note that if you are off-campus and need to access drp, you will need to tunnel through another server on campus -- such as those hosted by FIG or a CS Lab Computer. That is, SSH there and THEN SSH into drp.

**DRP Project Structure**

1. **./** *(The DRP Project)*
  - README	*-- This file.*
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
  - migrations/		*-- Migration Files*
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
use `git clone https://github.com/darkreactions/DRP.git` to copy the repository to
your workstation. 

Django has an important file, settings.py, which is not tracked by this git
repository for security reasons. If you add or remove items in this file
whilst developing with the code, please ensure that you update the
file settings_example.py in step.

Lastly, the repository utilises branches heavily. The master branch
is reserved for releases. There are no working long-term support branches
at present. Persons editing this code should set up their own branches,
with one marked as stable, such that all tests pass in that branch.    

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

**build\_model**

Builds a machine learning model. The most basic usage is
`python manage.py build_model -p 'reaction_temperature'`
this will build a model using only the reaction temperature
as a descriptor to predict the default outcome descriptor
(currently boolean crystallisation outcome) using the default
model (currently a Weka SVM using the PUK kernel and cross-validated using
a 4-fold split to analyze model performance).
More advanced usage should be well documented in the help text
`python manage.py build_model -h`

**import\_data**

Imports reaction data from the main Haverford Dark Reactions server.
Accepts one positional argument importing the corresponding number
of Performed Reactions and their associated data. Not available to
persons who are not administrators on the main server.

Note that at present this command does not import reaction descriptors.

**re\_save\_reactions**

Starts a batch parralell task to re-save each reaction, forcing
descriptor calculation. Useful in conjunction with the above
command. Note that descriptors will only be calculated
for plugins present in the code-base, which are correctly
set up.

**run_tests**

Runs all of the tests correctly imported in the test suite.
To only run some tests, one may enter a list of test modules
to run as positional arguments. The --failfast option causes
tests to halt on the first failure. Otherwise all tests will
be run and error details output at the end.
