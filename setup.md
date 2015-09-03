#Setup Page for the Dark Reaction Project
######Last updated by Philip Adler 2015-09-02

1. **Ubuntu**
3. **Git Hooks**

##Ubuntu development server

###Set up environment

The following instructions are written to work with Ubuntu 12+ and have (mostly) been tested on Ubuntu 13.

Install the necessary programs.
`sudo apt-get install python-dev python-pip mysql-server python-mysqldb nginx uwsgi uwsgi-plugin-python python-rdkit git virtualenvwrapper weka graphviz memcached python-memcache mailutils`

Note which version of Django gets installed.
`sudo pip install numpy scipy Django pygraphviz`

Install required pip python libraries
`sudo pip install chemspipy requests`

`sudo easy_install South`

###Clone from the Git Repository into your directory of choice.

`git clone git@github.com:cfalk/DRP`

####For the time being
Switch your branch to the `phil_refactor` branch.

####Server settings

In the DRP repository there is a file DRP_uwsgi.ini and another DRP_nginx. Both should be modified to suit your local server *after* having been placed in the relevant locations:

`/etc/uwsgi/apps-enabled/DRP_uwsgi.ini`
`/etc/nginx/sites-enabled/DRP_nginx`

It should be noted that the uwsgi.ini is backwards compatible with older version of this repo, but that an old DRP_uwsgi file will need replacing.

Both uwsgi and nginx must be restarted (in that order) for the server to work.

###Set up the settings.py file

In DRP/DRP, there is a file called 'settings_example.py'. This must be copied to 'settings.py', and the settings therein set to the appropriate values for your server. At present, the available fields should be fairly self explanatory, though the following should be noted:

`ALLOWED_HOSTS` should be set to your local ip

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


###Note that everything from here on has been retained for future use, but does not currently apply in this development version

###Skip if you don't need model-building.
Install [ChemAxon](http://www.chemaxon.com/)'s JChem and Marvin. These do intense, high-level chem-calcs for us. Note that these DON'T need to be installed on a test-bed if you don't need to test the compound property calculators. Also note that ChemAxon requires an install key in order to run.

###Skip if you don't want backups.
On a production server, use a "cron" job to back up the database (the command "sudo scrontab -e" initiates crontab editor). Note that the password needs to change from what it is below. Ideally the password will be stored more discretely in one place on the server. Add the following line to dump the database to a unique file every day in a DropBox backup folder:
`30 3 * * * mysqldump -uroot -p SecurePassword DRP_db > /home/drp/database_backups/Dropbox/DRP_db__$(date +\%m_\%d_\%y).sql`

Likewise, set up nginx to start every hour in case it goes happens to go down:
`0 * * * * service nginx start`

And back up the models directory each day:
`0 3 * * * cp -rn /home/drp/web/darkreactions.haverford.edu/app/DRP/models/* /home/drp/model_backups/`


###Set up a virtual environment.
Install a virtualenv if you are on the production server. This helps prevent hackers from accessing file and system information outside of the directory of the projects. There are also other perks (which can be researched online).

The URL below might change depending on your username.
`git clone git@bitbucket.org:darkreactionproject/dark-reaction-site.git`

###Set up Nginx and uWSGI (Not necessary if you develop with Django's runserver).
Move the Nginx and uWSGI files to their appropriate directories. Note that in development, you can just the the Django runserver (command: "python manage.py runserver") and thus can skip anything related to nginx/uwsgi.

`sudo mv DRP_nginx aetc/nginx/sites-available/DRP_nginx`

`sudo mv DRP_uwsgi /etc/uwsgi/apps-available/DRP_uwsgi`

Set up "sym links" so that the files are "activated." This is good practice so that the actual configuration file itself is never discarded when a site is de-activated. Make sure to use absolute file locations.
`sudo ln -s /etc/nginx/sites-available/DRP_nginx /etc/nginx/sites-enabled/DRP_nginx`

`sudo ln -s /etc/uwsgi/apps-available/DRP_uwsgi /etc/uwsgi/apps-enabled/DRP_uwsgi`

Restart nginx and uwsgi so that they know that we changed files. If we don't, they will continue to distribute the obsolete "cached" versions. This should be done every time a change is made.
`sudo services nginx restart`

`sudo services uwsgi restart`

###Set up the MySQL database.
Copy over a version of the database from the production server or a backup. Choose the most recent mysqldump file to avoid data-loss.
` scp drp@darkreactions.haverford.edu:/home/drp/database_backups/Dropbox/DRP_db__03_30_14.sql ./database.sql`

Make sure the MySQL user accessing the database has a secure password.

Create an empty database for DRP using the MySQL client:
`mysql -uroot -p`

`mysql> CREATE DATABASE DRP_db;`

`exit`

Load the mysqldump into this empty DRP database. Note that "SecurePassword" should be more secure in production.
`mysql -u root --password=SecurePassword DRP_db < database.sql`

Remove any migration history that might exist. ALWAYS BE CAREFUL WHAT YOU DELETE.
`rm -r DRP/migrations/`

###Set up South for Database Migrations (see README.md for more).
Convert the database to be managed by South. The last "--fake" and "--delete-ghost-migrations" options specify that the database version that South has stored in its personal database tables should be ignored.
`python manage.py syncdb`

`python manage.py convert_to_south DRP`

`python manage.py schemamigration DRP --initial`

`python manage.py migrate DRP --fake --delete-ghost-migrations`

When you change a model changes or add a field, you must run the following two commands. They'll also need to be run on the production server. Be warned that if you change the Django models before South is loaded, it may blow up on you.
`python manage.py schemamigration --auto DRP`

`python manage.py migrate DRP`

###In Production
At this point, you should have a fully-functional production-bed or test-bed. Verify by going to wherever the webapp is hosted ([darkreactions.haverford.edu](http://darkreactions.haverford.edu)).

#In Development:
Verify by running `python manage.py runserver 0.0.0.0:8000` and going to port 8000 of the server hosting your workspace (ie: if you're on drp, go to [darkreactions.haverford.edu:8000](http://darkreactions.haverford.edu:8000)). Note that you may not be able to view ports other than port 80 while outside of campus.

