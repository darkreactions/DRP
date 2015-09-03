#Django settings for DRP project.
#TODO: Set / to being os.sep
import os
SITE_ID = 1

SERVER_NAME = ''
TESTING = False #tells the system that you are happy for tests, some of which run on the live database to be run. DO NOT SET TO TRUE IN PRODUCTION.
EXTERNAL_HTML_VALIDATOR = 'http://validator.w3.org/nu/' #The external html validator. You shouldn't need to change this.

CHEMSPIDER_TOKEN = ''

LOGIN_REDIRECT_URL = '/'

APP_DIR = (os.path.join(os.path.dirname(__file__)))
BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# URLs
STATIC_URL = '/static/'
STATIC_ROOT = APP_DIR + "/static_served/"
TEMPLATE_DIRS = (os.path.join(BASE_DIR, 'templates'),)

# Directories
STATIC_DIR = os.path.join(BASE_DIR, "static")
TMP_DIR = os.path.join(BASE_DIR, "tmp")
RESEARCH_DIR = os.path.join(BASE_DIR + "research")
LOG_DIR = os.path.join(BASE_DIR, "logs")
MODEL_DIR = os.path.join(BASE_DIR, "models")

CHEMAXON_DIR = {
  }
# {version:directory}
WEKA_PATH = {
  '3.6':'/usr/share/java/weka.jar' #default path on Ubuntu
}

if TESTING:
  MOL_DESCRIPTOR_PLUGINS=('DRP.plugins.moldescriptors.example',)
  RXN_DESCRIPTOR_PLUGINS=()
else:
  MOL_DESCRIPTOR_PLUGINS=('DRP.plugins.moldescriptors.example',)
  RXN_DESCRIPTOR_PLUGINS=()

STATICFILES_DIRS = (STATIC_DIR,)

#Changes to Default Django Behavior
LOGIN_URL = "/login.html"

#Email Settings
EMAIL_USE_TLS = True
EMAIL_HOST = ""
EMAIL_PORT = 587
EMAIL_HOST_USER = ""
EMAIL_HOST_PASSWORD = "" #TODO: Change me in production!
DEFAULT_FROM_EMAIL = EMAIL_HOST_USER
EMAIL_IMAP_HOST = '' #leave blank in production. Necessary for unit tests.
EMAIL_IMAP_INBOX= 'Inbox' #The inbox to use for email tests. Inbox for gmail
SKIP_EMAIL_TESTS=False #These tests are slow, so if you haven't tweaked this then skip the tests, but don't abuse this.


#Change to "False" to see standard errors:
DEBUG = True #TODO: Change me in production to "False"!
#DEBUG = False if TESTING else True #useful in dev environments
TEMPLATE_DEBUG = DEBUG

ALLOWED_HOSTS = []
ADMINS = (
     )

MANAGERS = ADMINS
#Emails of the Site Admins and Project Managers for the DRP.

ADMIN_EMAILS = tuple(admin[0] for admin in ADMINS)
MANAGER_EMAILS = tuple(manager[0] for manager in MANAGERS)

DATABASES = { #Production database.
    'default': {
        'ENGINE': '', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': '',                      # Or path to database file if using sqlite3.
        # Test_DRP_db should be used next time!
	'USER': '',
	#TODO: Change me in production!
        'PASSWORD': 'SecurePassword', ###Delete me in your commits!!!###################################
        'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
        'PORT': '3306',                      # Set to empty string for default.
    }
}

# ==== Another Useful Hack for Testing Environments - use this INSTEAD of the above block ===
#if TESTING:  
#  DATABASES = { #Production database.
#      'default': {
#          'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
#          'NAME': 'Test_DRP_db_2',                      # Or path to database file if using sqlite3.
#          # Test_DRP_db should be used next time!
#  	'USER': 'root',
#          'PASSWORD': '', ###Delete me in your commits!!!###################################
#          'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
#          'PORT': '3306',                      # Set to empty string for default.
#      }
#  }
#else:
#  DATABASES = { #Production database.
#      'default': {
#          'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
#          'NAME': 'Test_DRP_db_3',                      # Or path to database file if using sqlite3.
#          # Test_DRP_db should be used next time!
#  	'USER': 'root',
#  	#TODO: Change me in production!
#          'PASSWORD': '', ###Delete me in your commits!!!###################################
#          'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
#          'PORT': '3306',                      # Set to empty string for default.
#      }
#  }

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = "America/New_York"

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = False

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/var/www/example.com/media/"
MEDIA_ROOT = os.path.join(BASE_DIR, 'media_root')

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://example.com/media/", "http://media.example.com/"
MEDIA_URL = '/media/'


# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = ''

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'DRP.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'DRP.wsgi.application'

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    # Uncomment the next line to enable the admin:
     'django.contrib.admin',
    # Uncomment the next line to enable admin documentation:
     'django.contrib.admindocs',
     "DRP",
     "south"
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        }
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

TEMPLATE_CONTEXT_PROCESSORS = ("django.contrib.auth.context_processors.auth",
"django.core.context_processors.debug",
"django.core.context_processors.i18n",
"django.core.context_processors.media",
"django.core.context_processors.static",
"django.core.context_processors.request",
"django.core.context_processors.tz",
"django.contrib.messages.context_processors.messages",
'DRP.context_processors.testing')

#Set up Memcached caching:
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',
        'LOCATION': '0.0.0.0:11211',
    }
}

#Force users to log out when the browser is closed.
SESSION_EXPIRE_AT_BROWSER_CLOSE = True
LIBRARY_CHOICES=()
TOOL_CHOICES=()

EMPTY_LABEL='----'

LAB_GROUP_HASH_SALT = ''

#force temporary file creation for uploads (required for some views to work)
FILE_UPLOAD_HANDLERS=("django.core.files.uploadhandler.TemporaryFileUploadHandler",)
