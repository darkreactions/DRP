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

