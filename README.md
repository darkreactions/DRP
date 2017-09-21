Dark Reaction Project README
Last updated by Monique Byars 09/21/17
General Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
Always make sure to check the latest version of the README on the master branch.
This repository contains the software for the https://www.djangoproject.com/-based source code for the Dark Reactions Project Software. If you are looking to contribute to the chemistry aspects of the project, please visit the main project site at http://darkreactions.haverford.edu. If you are looking to contribute to the source code of the project and are not a member of haverford college, please fork this repository and issue a pull request with any changes or fixes you may have made. A list of known bugs can be found at [http://bugs.darkreactions.haverford.edu] (our instance of Mantis Bug Tracker). Please note that you will to sign up for an account, and that the authentication credentials for the bug reporting page and the main project page are separate.
Setting up your own instance of the DRP
The following setup methods are suitable for use with Ubuntu 16.04.
On a Virtual Machine
This software supports development using https://www.vagrantup.com, and has a Vagrantfile included in the root directory.
Vagrant can be installed using apt-get:
sudo apt-get install vagrant virtualbox
Then setting up a virtual machine should be as simple as issuing the command vagrant up from anywhere in the repository.
Manual Instructions - for a physical machine

The following instructions are written to work with Ubuntu 14 and 16 and have (mostly) been tested. These instructions assume familiarity with Linux and a Command Line, and that you are using nginx as your webserver. 
STEP 1: Installing the necessary programs
* (Required)
The DRP *
Mysql * 
Nginx *
Uwsgi  
Note: other wsgi solutions are available and may be compatible with DRP.
Chemspider*
Firstly, download or clone the DRP from our git repo: https://github.com/darkreactions/DRP
Then:
sudo apt-get install python3 python3-dev python3-pip mailutils mysql-server libmysqlclient-dev nginx uwsgi uwsgi-plugin-python3 python-rdkit git weka graphviz memcached python-memcache python3-scipy python3-pillow cmake libboost-all-dev python3-cffi graphviz-dev pkg-config pwgen dnsmasq
sudo pip3 install numpy pygraphviz mysqlclient
sudo -H pip3 install chemspipy requests pep8 pep257 xxhash sqlparse
Install Django*. The current version of DRP is designed to work with Django 1.8
sudo pip3 install django==1.8
Install required pip python libraries*.
sudo pip3 install chemspipy requests pep8 pep257 xxhash
STEP 2: Setting up MySQL and Chemspider
	Consult documentation for outstanding questions: https://dev.mysql.com/doc/
Login to your account: 
mysql -u root -p
The user is initialized as root and the password is blank in Ubuntu. This may vary depending on your setup. 
Once you login to mysql, and create two databases:
  CREATE DATABASE DRP CHARACTER SET utf8 COLLATE utf8_bin;
Expected output: Query OK, 1 row affected (0.00 sec)

   CREATE DATABASE DRP_test CHARACTER SET utf8 COLLATE utf8_bin;
Expected output: Query OK, 1 row affected (0.00 sec)
Exit out of MySQL my typing “exit.”
You will also need a chemspider token. Follow the directions on the website for webservices: 
“Some operations require a security token; to obtain a token please complete the registration process –when you are registered the Security Token is listed on the Profilepage. For web services which require a “Service Subscriber” role, then email us above to discuss upgrading your user account.”
http://www.chemspider.com/AboutServices.aspx
STEP 3: Setting up the settings.py file
In the DRP repository, copy the settings_example.py to settings.py.
cp settings_example.py settings.py
In the settings.py file, you will need to edit:
SERVER_NAME (full name of PC)
CHEMSPIDER\_TOKEN
MAIN\_SERVER
MAIN\_SERVER\_USER
MAIN\_SERVER\_PASS 
EMAIL\_HOST
EMAIL\_HOST\_USER
EMAIL\_HOST\_PASS
EMAIL\_IMAP\_HOST
ADMINS 


These have placeholders that indicate their use. 
You will also need to set the standard database settings as per the [https://docs.djangoproject.com/en/1.8/](django documentation).
ALLOWED_HOSTS should be empty.
To pass the unit tests, at least one ADMIN_EMAILS should be provided
STEP 4: Nginx and Uwsgi
Copy DRP_nginx from the DRP folder to the /etc/nginx/sites-available folder.
	cp DRP/DRP_nginx /etc/nginx/sites-available
Create a symlink between sites-available/DRP_nginx and sites-enabled/DRP_nginx
ln -s sites-available/DRP_nginx sites-enabled/DRP_nginx


Copy DRP_uwsgi.ini file to the /etc/uwsgi/apps-available folder. 
	cp DRP/DRP_uwsgi.ini /etc/uwsgi/apps-available
Create a symlink between apps-available/DRP_uwsgi.ini to apps-enabled/DRP_uwsgi.ini
ln -s /etc/uwsgi/apps-available/DRP_uwsgi.ini /etc/uwsgi/apps-enabled/DRP_uwsgi.ini

Note: For symlinks, you may need to create the directory before linking it.
Now, edit each fine. Replace the placeholders with the relevant values. Ex: log path. 
Note that the “location” paths in the nginx DRP_nginx file will need to be changed if you are not working on the server.  
Now restart both Uwsgi: 
sudo service uwsgi restart 
And nginx:
sudo service nginx restart
STEP 5: RDKIT Descriptors
	Option 1 (recommended):
We have a repo setup for a python3 specific build of rdkit which should not clash with other packages in the ubuntu repositories, however, we make no guarantees to that effect, and installation is at your own risk.
In the file /etc/apt/sources.list add the line: 
deb [trusted=yes] https://darkreactions.haverford.edu/software ./
If you need to change the permissions on the /etc/apt/sources.list using chmod to edit the file, be sure to change it back. 
Then execute:
sudo apt-get update
sudo apt-get install python3-rdkit
	Option 2 (manual):
sudo apt-get install bison cmake flux build-essential sqlite3 libsqlite3-dev libboost-all-dev
sudo pip3 install cairocffi
Then, in a directory that is not your main installation directory for DRP.
git clone https://github.com/darkreactions/rdkit.git
Change into the rdkit repository directory and then
export RDBASE=$(pwd)
export LD_LIBRARY_PATH="$(pwd)/lib"
export PYTHONPATH="$(pwd)/lib"
mkdir build
cd build
cmake -DRDK_BUILD_INCHI_SUPPORT=ON -D PYTHON_LIBRARY=/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4.so -D PYTHON_INCLUDE_DIR=/usr/include/python3.4/ -D PYTHON_EXECUTABLE=/usr/bin/python3.4 -DBOOST_ROOT=/usr/lib/x86_64-linux-gnu/ ..
make install
ctest
If any of this generates an error, unless you are very familiar with compiling new code for linux operating systems, seek the assistance of Philip Adler via the rdkit repository provided in this document.
Otherwise:
unset LD_LIBRARY_PATH unset RDBASE unset PYTHONPATH
cd ../rdkit
`sudo ln -s "$(pwd)" /usr/lib/python3.4/rdkit
cd ../lib
sudo cp -i *.so.2 /usr/lib
python3.4 -c "import rdkit.Chem"
If that runs without error messages, congratulations, you have compiled and installed rdkit for use with DRP.
STEP 6: Installing Chemaxon (optional)
Follow the documentation to install and validate your license for ChemAxon
https://www.chemaxon.com/download/marvin-suite/#marvin
STEP 7: Remaining loose ends
Make sure you run the command: 
	./manage.py migrate
(If you are experiencing errors with management commands, make sure you are in the right directory and try putting python3 or whichever python is relevant before the command.)

_________________________________________________________________
General Information
Logging in
When first registering through the website on a server, the user is default inactive. 
You must run:
	./manage.py createsuperuser
Enter your new username and password. 
Login to the website with /admin at the end of the URL
Go into the users and change the new user to active and whatever permissions you deem necessary. 
Release Versions
DRP is distributed in release versions. To use a specific version of the code, use the following template command:
git checkout <version number>
Git guidelines 
###Git Hooks
DRP comes distributed with a number of useful git hooks in the drp_hooks directory in this repo. These warn you if the expected structure of the settings.py file changes, or if you need to run database migrations for the DRP Django application. They also remove orphaned .pyc files, which have been known to confuse the test suite historically. These should be added to your local git repository as per the git documentation.
####Server settings
In the DRP repository there is a file DRP_uwsgi.ini and another DRP_nginx. Both should be modified to suit your local server after having been placed in the relevant locations:
/etc/uwsgi/apps-enabled/DRP_uwsgi.ini /etc/nginx/sites-enabled/DRP_nginx
It should be noted that the uwsgi.ini is backwards compatible with older version of this repo, but that an old DRP_uwsgi file will need replacing.
Both uwsgi and nginx must be restarted (in that order) for the server to work.
###https
If you are setting up a publicly viewable instance of DRP, there are additional settings for these protocols present in the configuration file for nginx, and the Django settings_example.py file, which have been left commented out. Setting up https access varies greatly depending on your local server environment and organisation so there will be no further documentation here. Additional information for a simple method to set up https can be found at https://letsencrypt.org.
###Running tests
In order to run tests you must have the following environment variables set up in your shell session:
export PYTHONPATH=/path/to/DRP/
export DJANGO_SETTINGS_MODULE=DRP.settings
You must also have TESTING set to True in your settings.py file.
To run specific tests, provided the test is conformant to the template test (which they should be if you are writing new tests!), one can simply execute the test:
path/to/DRP/test.py
Else, one can run the entire test suite from the management script:
./manage.py run_tests
###On Development Servers
The ALLOWED_HOSTS option in 'settings.py' should be set to an iterable containing only the localhost ip address as a string.
In the '/etc/nginx/sites-enabled/DRP_nginx' file, the host name that is being listened to should only be localhost.
##Servers with multiple web applications.
If you are only developing DRP on your server, the setup you have should be sufficient, however, people running other applications on their local development server should note the following.
If you are running the django testing server, this requires you to select a port which is unoccupied. By default, the nginx settings file listens for port 8000, which is also the default port of the django test server; you will need to configure one or the other so that this clash does not occur. The Django documentation addresses this for django, whilst in the DRP_nginx file, the only change that needs to be made is to delete the line:
listen 8000
#dnsmasq
For instances where you are hosting multiple development projects on your local server, it may be beneficial to install dnsmasq:
sudo apt-get install dnsmasq
dnsmasq is a powerful tool for rerouting and managing dns requests. This makes it extremely helpful in managing multiple local development projects.
Having installed dnsmasq, open the file /etc/dnsmasq.conf in your favourite text editor, and add the following line into the file:
address=/loc/127.0.0.1
Save the change, and then on the command line:
sudo service dnsmasq restart
In the DRP_nginx file, change the server_name configuration to something like darkreactions.loc. It does not matter what this is set to, provided it:
a. is unique on your development server b. ends in .loc
Don't forget to set the SERVER_NAME setting in your settings.py file to the same value!
Restart nginx:
sudo service nginx restart
When you open your browser and direct yourself to darkreactions.loc (or whatever you named the server), the dark reactions project should display.
###DRP Versions < 0.1
Versions of code prior to 0.1 are not compatible with versions above.
Upgrading from versions prior to 0.7
Version 0.6 of DRP was the last version to use Django 1.6, subsequent versions use Django 1.8. There are, therefore, necessary transition steps to be made.
git fetch --all git checkout 0.7
Make sure your database (both main and testing) is up to date with migrations: Set Testing = False in settings.py./manage.py migrate
Set Testing = True in settings.py
Restart nginx and uwsgi sudo service nginx restart && sudo service uwsgi restart
Run all tests and make sure you pass them all: ./manage.py run_tests
git checkout 0.81
Remove "south" from your installed apps
Install django 1.8: pip install -U Django==1.8.9
Delete all .pyc files in your migrations folder rm DRP/migrations/*.pyc ./manage.py migrate DRP --fake./manage.py migrate --fake-initial ./manage.py migrate
Repeat this process for your Main database.
Restart nginx and uwsgi sudo service nginx restart && sudo service uwsgi restart
Run all tests and make sure you pass them all: ./manage.py run_tests
Notes for Local (Haverford) Developers
The website can be accessed at darkreactions.haverford.edu -- this domain is managed by Haverford. The server itself (named "drp") is in the KINSC Server Room and can be accessed by SSH while on campus. Note that if you are off-campus and need to access drp, you will need to tunnel through another server on campus -- such as those hosted by FIG or a CS Lab Computer. That is, SSH there and THEN SSH into drp.
DRP Project Structure
./ (The DRP Project)
README	-- This file.
DRP_nginx	-- The NGINX configuration. Move to /etc/nginx/sites-enabled/
DRP_uwsgi	-- The UWSGI configuration. Move to /etc/uwsgi/apps-enabled/
manage.py	-- The Django management file -- leave it as it is.
DRP/ (The DRP App)
models/	-- Contains all the Django models and some important accessors.
settings_example.py	-- Contains example settings for the database; true settings are placed in settings.py (see setup.md).
urls/	-- Responsible for mapping a URL call to a function call; see the django documentation for details.
research/	-- A folder for very experimental/highly unstable code.
management/	-- Directory to add custom python manage.py commands- see the django documentation.
templatetags/	-- Directory to add custom Django template tags- see the django documentation.
views/	-- Directory of various views sorted into different files.
migrations/	-- Migration Files
templates (The HTML templates for Django Views.)
research (Any "research" scripts that are still being explored.)
static (Any static files for which Django can skip templating.)
favicon.ico	-- The "favicon" for the site (actually served by NGINX).
js/	-- Any Javascript for any view belongs in this directory. universal.js -- Contains Javascript used on nearly every page.
css/	-- The CSS for any view should go here.
icons/	-- Small "icon" images belong here -- such as the hover-button images.
images/	-- Any large images for the site should go here (eg: the logo).
admin/	-- CSS for the /admin/ page. Used in many django installations.
There are many files that are not listed above in order to elucidate the "framework" of the DRP Django Project succinctly. Notably, there are many python files in the views and research directories that are not listed above -- but which should contain explanatory comments in the files themselves.
Accessing the GitHub Repo and Notes on the Repo Structure
Firstly, you'll need a GitHub account and you'll need someone with access to the repo to grant your account access (though if you can view this README without access to the GitHub repo, you should tell someone). Then, you should be able to use git clone https://github.com/darkreactions/DRP.git to copy the repository to your workstation.
Django has an important file, settings.py, which is not tracked by this git repository for security reasons. If you add or remove items in this file whilst developing with the code, please ensure that you update the file settings_example.py in step.
Lastly, the repository utilises branches heavily. The master branch is reserved for releases. There are no working long-term support branches at present. Persons editing this code should set up their own branches, with one marked as stable, such that all tests pass in that branch.
Django Management Commands
Django management commands (that is, the commands that pass through manage.py) can be called as python manage.py <the command> when in the main DRP project directory. These commands are each stored in a separate python file in the .../DRP/management/commands directory. To add a new command, use the existing files as examples.
It also should be noted that using python manage.py help will detail all of the available management commands (of which there are many).
check_hash_collisions
This command checks the hashed values of the reactions calculated by the DRP descriptor plugin for clashes. If there is a clash, the lead developer should be notified ASAP.
build_model
Builds a machine learning model. The most basic usage is python manage.py build_model -p 'reaction_temperature'this will build a model using only the reaction temperature as a descriptor to predict the default outcome descriptor (currently boolean crystallisation outcome) using the default model (currently a Weka SVM using the PUK kernel and cross-validated using a 4-fold split to analyze model performance). More advanced usage should be well documented in the help textpython manage.py build_model -h
import_data
Imports reaction data from the main Haverford Dark Reactions server. Accepts one positional argument importing the corresponding number of Performed Reactions and their associated data. Not available to persons who are not administrators on the main server.
Note that at present this command does not import reaction descriptors.
re_save_reactions
Starts a batch parralell task to re-save each reaction, forcing descriptor calculation. Useful in conjunction with the above command. Note that descriptors will only be calculated for plugins present in the code-base, which are correctly set up.
run_tests
Runs all of the tests correctly imported in the test suite. To only run some tests, one may enter a list of test modules to run as positional arguments. The --failfast option causes tests to halt on the first failure. Otherwise all tests will be run and error details output at the end.



