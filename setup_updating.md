#Setup Page for the Dark Reaction Project
######Last updated by Philip Adler 2015-09-02

1. **Ubuntu**

The following instructions are written to work with Ubuntu 12+ and have (mostly) been tested on Ubuntu 13.

##Ubuntu development server

###Set up a virtual environment
Install a virtual environment.
`sudo pip install virtualenv`

Create a virtual environment.
`cd <project_folder>`
`virtualenv <venv_name>`

Activate the virtual environment.
`source <venv_name>/bin/activate`

Pip install packages as usual except DON'T USE SUDO.

To exit the virtual environment
`deactivate`

For more information on using virtual environments, check out the documentation (http://docs.python-guide.org/en/latest/dev/virtualenvs/)

###Install scipy
`sudo apt-get install python-dev python-scipy python-rdkit`

###Install other non-python things
`sudo apt-get install python-pip mysql-server libmysqlclient-dev nginx git weka graphviz memcached mailutils`


###Clone from the Git Repository into your directory of choice.

`git clone git@github.com:cfalk/DRP`

####For the time being
Switch your branch to the `phil_refactor` branch.


###Set up environment

Activate virtual environment. Go into the git repository and (NO SUDO)
`pip install -r requirements.txt`


### OLD REMOVE
                    Install the necessary programs.
                    `sudo apt-get install python-dev python-pip mysql-server python-mysqldb nginx uwsgi uwsgi-plugin-python python-rdkit git virtualenvwrapper weka graphviz memcached python-memcache mailutils`
                    
                    Note which version of Django gets installed.
                    `sudo pip install numpy scipy Django pygraphviz`
                    
                    Install required pip python libraries
                    `sudo pip install chemspipy requests pep8 pep257 xxhash South uwsgi`
                    
                    pip install uwsgi

### OLD REMOVE
                

####Server settings

In the DRP repository there is a file DRP_uwsgi.ini and another DRP_nginx. Both should be modified to suit your local server *after* having been placed in the relevant locations:

`/etc/uwsgi/apps-enabled/DRP_uwsgi.ini`
`/etc/nginx/sites-enabled/DRP_nginx`

It should be noted that the uwsgi.ini is backwards compatible with older version of this repo, but that an old DRP_uwsgi file will need replacing.

Both uwsgi and nginx must be restarted (in that order) for the server to work.

###Set up the settings.py file

In DRP/DRP, there is a file called 'settings_example.py'. This must be copied to 'settings.py', and the settings therein set to the appropriate values for your server. At present, the available fields should be fairly self explanatory, though the following should be noted:

`ALLOWED_HOSTS` should be set to an iterable containing only the element '*'.

To pass the unit tests, at least one `ADMIN_EMAILS` should be provided

To pass the unit tests, the EMAIL_HOST_USER and related settings should be set.  

### Working with south

South is our library of choice for dealing with database migrations. It is already "installed" by the default settings file.

To get DRP running you need to setup a database with the correct encoding, in mysql:

`CREATE DATABASE db_name DEFAULT CHARACTER SET utf8 DEFAULT COLLATE utf8_bin;`

Without this, all string comparisons done at a database level are case insensitive, which causes spurious query results.

Ensure that `db_name` matches with the name specified in your settings.py file.

Then, in the DRP directory, run `./manage.py syncdb`

This sets up the south tables.

You should create a "superuser" called "root" when prompted.

Then run `./manage.py migrate DRP`

This should set up the schema to the most recent version. Note that you may need to do this periodically as another person updates the model.

###Git hooks.

Inside DRP there is a folder called drp_hooks. Among these is a file called pre-push, and this should be moved in your .git/hooks directory for this project. This runs all of the tests
in the test suite before pushing your copy of the repository, and hence forces our code into compliance.

###Running tests

In order to run tests you must have the following environment variables set up in your shell session:

export PYTHONPATH=/path/to/DRP/

export DJANGO_SETTINGS_MODULE=DRP.settings

This also applies to the pre-push hook.

You must also have `TESTING` set to `True` in your `settings.py` file.

To run specific tests, provided the test is conformant to the template test (which they should be if you are writing new tests!), one can simply execute the test:

`path/to/DRP/test.py`

Else, one can run the entire test suite from the management script:

`./manage.py run_tests`

###For the time being...

If you are reading this setup instruction manual, you are working on the refactoring mission of the DRP. This comes with some rules:

1. Changes to the schema should be run past Phil, who will implement them; this makes dealing with the south migrations much easier.
2. You should make your own branch from the phil_refactor branch, and push that to the Git repo (`push -u origin <branchname>`)- Phil will merge periodically or on request.
3. Don't attempt to circumvent the pre-push tests.
4. Comment all the things.
5. If you change the pre-push hook, tell everyone.

There are additional developer notes in teh developer_notes.md file.

###On Development Servers

The `ALLOWED_HOSTS` option in 'settings.py' should be set to an iterable containing only the localhost ip address as a string.

In the '/etc/nginx/sites-enabled/DRP_nginx' file, the host name that is being listened to should only be localhost.

##Servers with multiple web applications.

If you are only developing DRP on your server, the setup you have should be sufficient, however, people running other applications on their local development server should note the following.

If you are running the django testing server, this requires you to select a port which is unoccupied. By default, the nginx settings file listens for port 8000, which is also the default port of the django test server; you will need to configure one or the other so that this clash does not occur. The Django documentation addresses this for django, whilst in the DRP_nginx file, the only change that needs to be made is to delete the line:

`listen		8000`

 #dnsmasq

For instances where you are hosting multiple development projects on your local server, it may be beneficial to install dnsmasq:

`sudo apt-get install dnsmasq`

dnsmasq is a powerful tool for rerouting and managing dns requests. This makes it extremely helpful in managing multiple local development projects.

Having installed dnsmasq, open the file `/etc/dnsmasq.conf` in your favourite text editor, and add the following line into the file:

`address=/loc/127.0.0.1`

Save the change, and then on the command line:

`sudo service dnsmasq restart`

In the DRP_nginx file, change the `server_name` configuration to something like `darkreactions.loc`. It does not matter what this is set to, provided it:

a. is unique on your development server
b. ends in `.loc`

Don't forget to set the `SERVER_NAME` setting in your settings.py file to the same value!

Restart nginx:

`sudo service nginx restart`

When you open your browser and direct yourself to darkreactions.loc (or whatever you named the server), the dark reactions project should display.

###Additional Notes for Production Servers

There are a number of cron-jobs that need to be set to ensure good functioning of a production grade server.

Firstly, there is the reaction hash checking cron-job. This is a database integrity check that cannot be done using django's inbuilt framework. On failure, it notifies your local administrators that there is a problem. It is advised that should this check fail, that your local administrator file a bug report with the development group.

The command for the cron-job is ./manage.py check_hash_collisions


