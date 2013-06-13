# Initialize App Engine and import the default settings (DB backend, etc.).
# If you want to use a different backend you have to remove all occurences
# of "djangoappengine" from this file.
from djangoappengine.settings_base import *

import os

# Activate django-dbindexer for the default database
DATABASES['native'] = DATABASES['default']
DATABASES['default'] = {'ENGINE': 'dbindexer', 'TARGET': 'native'}
AUTOLOAD_SITECONF = 'indexes'

STATIC_ROOT = (os.path.join(os.path.dirname(__file__)))
STATIC_URL = '/static/'
STATICFILES_DIRS = (os.path.join(os.path.dirname(__file__), 'static'),)

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)


SECRET_KEY = '=r-$b*8hglm+858&9t043hlm6-&6-3d3vfc4((7yd0dbrakhvi'

INSTALLED_APPS = (
	'django.contrib.admin',
	'django.contrib.contenttypes',
	'django.contrib.auth',
	'django.contrib.sessions',

	#Modified:
	'django.contrib.sites',
	'django.contrib.messages',
	'django.contrib.staticfiles',

	#Django-nonrel needs these: ###C
	'djangotoolbox',
	'autoload',
	'dbindexer',
	'formsite',

	# djangoappengine should come last, so it can override a few manage.py commands
	'djangoappengine',
)

MIDDLEWARE_CLASSES = (
	# This loads the index definitions, so it has to come first
	'autoload.middleware.AutoloadMiddleware',

	'django.middleware.common.CommonMiddleware',
	'django.contrib.sessions.middleware.SessionMiddleware',
	'django.contrib.auth.middleware.AuthenticationMiddleware',
)

TEMPLATE_CONTEXT_PROCESSORS = (
	'django.contrib.auth.context_processors.auth',
	'django.core.context_processors.request',
	'django.core.context_processors.media',
	'django.core.context_processors.static',
)

# This test runner captures stdout and associates tracebacks with their
# corresponding output. Helps a lot with print-debugging.
TEST_RUNNER = 'djangotoolbox.test.CapturingTestSuiteRunner'

ADMIN_MEDIA_PREFIX = '/media/admin/'
TEMPLATE_DIRS = (os.path.join(os.path.dirname(__file__), 'templates'),)

ROOT_URLCONF = 'urls'

#Apply a user-profile to users: ###C
AUTH_PROFILE_MODULE = 'Lab_Member'
