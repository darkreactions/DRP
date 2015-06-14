#Django settings for DRP project.
import os
SITE_ID = 1

APP_DIR = (os.path.join(os.path.dirname(__file__)))
BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# URLs
STATIC_URL = '/static/'
LICENSE_URL = STATIC_URL + "/licenses/"
DYNAMIC_URL = '/dynamic/'
STATIC_ROOT = APP_DIR + "/static_served/"
TEMPLATE_DIRS = (os.path.join(BASE_DIR, 'templates'),)

# Directories
STATIC_DIR = BASE_DIR + "/static/"
TMP_DIR = BASE_DIR + "/tmp/"
RESEARCH_DIR = BASE_DIR + "/research/"
LOG_DIR = BASE_DIR + "/logs/"
MODEL_DIR = BASE_DIR + "/models/"
LICENSE_DIR = STATIC_DIR + "/licenses/"
LIBRARY_DIR = os.path.expanduser("~") #Defaults to "home"

STATICFILES_DIRS = (STATIC_DIR,)

#Changes to Default Django Behavior
LOGIN_URL = "/login/"

#Email Settings
EMAIL_USE_TLS = True
EMAIL_HOST = ""
EMAIL_PORT = 587
EMAIL_HOST_USER = ""
EMAIL_HOST_PASSWORD = "" #TODO: Change me in production!
DEFAULT_FROM_EMAIL = EMAIL_HOST_USER

#Emails of the Site Admins and Project Managers for the DRP.
ADMIN_EMAILS = [""]
MANAGER_EMAILS = [""]

#Change to "False" to see standard errors:
DEBUG = True #TODO: Change me in production to "False"!
TEMPLATE_DEBUG = DEBUG

ALLOWED_HOSTS = []
ADMINS = (
     )

MANAGERS = ADMINS

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
     "DRP",
     "south",
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

#Apply a user-profile to users:
AUTH_PROFILE_MODULE = 'DRP.lab_member'
ACCESS_CODE_LENGTH = 20 #Designates the max_length of user access_codes

#Set up Memcached caching:
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',
        'LOCATION': '0.0.0.0:11211',
    }
}

#Force users to log out when the browser is closed.
SESSION_EXPIRE_AT_BROWSER_CLOSE = True
