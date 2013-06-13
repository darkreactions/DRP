#!/bin/sh

mkdir build

pip install --download='./build' --no-install -r requirements.txt

unzip -q build/django-autoload-0.01.zip -d build
unzip -q build/django-dbindexer-0.3.zip -d build
unzip -q build/django-nonrel-1.4.5.zip -d build
unzip -q build/djangoappengine-1.0.zip -d build
unzip -q build/djangotoolbox-0.9.2.zip -d build

cp -r build/django-autoload/autoload ./autoload
cp -r build/django-dbindexer/dbindexer ./dbindexer
cp -r build/django-nonrel/django ./django
cp -r build/djangoappengine/djangoappengine ./djangoappengine
cp -r build/djangotoolbox/djangotoolbox ./djangotoolbox

rm -r ./djangoappengine/djangoappengine.egg-info
