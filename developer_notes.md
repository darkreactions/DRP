#Developer notes for the Dark Reaction Project

##Writing Tests

Tests for the DRP are stored in the `tests` folder (imaginitively).

These tests follow a layout as indicated by the `TemplateTests.py` file.

If you create a new tests file which matches that layout (and it should), then you must also update
the `__init__.py` file of that folder, again, following the format already therein.

This applies to any subfolders within the tests directory.

All test classes should inherit from the `DRPTestCase` class, although not necessarily directly;
there is a special folder of HttpTests, which test items which require a httprequest to be made.
There are additional decorators in the `HttpTests` folder, along with many helpful classes.

###Decorators

Inside the file `decorators.py` there are a variety of helpful class decorators which make it
easier to perform common operations like creating users as a part of your test.

###pep8 and pep257

pep8 and pep257 are the style documents for python. We need to work to make sure that our code
conforms to these standards as much as possible. Obviously, this is made harder by the
fact that we are dealing with a lot of legacy code.

As such, as you work on a file, whether you wrote it or not make sure to add it (in the appropriate format)
to the _pep8files variable in the /DRP/tests/fileTests.py file.

If a large number of files in the same directory are present in this variable, delete them and replace them
with a representation of the directory in question.

Then, fix the files as they inevitably fail the checking sequence.
