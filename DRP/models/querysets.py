'''A module for turning querysets into csv files or output'''
import csv
import numpy as np
from django.db import models
import abc
from collections import OrderedDict
from itertools import islice, chain


class QuerySetSet(object):
    """
    A set of querysets that quacks like a query set.
    For any method called on it, it calls the corresponding method on each of its querysets.
    Be careful that this makes sense for the given querysets.
    For some queryset methods, this doesn't make sense (like order_by) and more work is needed.
    Not all of these are implemented.
    Partially stolen from here: http://stackoverflow.com/questions/431628/how-to-combine-2-or-more-querysets-in-a-django-view
    and partially inspired from this: http://ramenlabs.com/2010/12/08/how-to-quack-like-a-queryset/
    """
    # queryset methods that return a queryset.
    # These are passed through to the underlying querysets unless defined otherwise
    # Does django indicate these in any way other than reading them from the docs?
    qs_methods = ['filter', 'exclude', 'annotate', 'order_by', 'reverse', 'distinct', 'values', 'values_list', 'dates', 'datetimes',
                  'none', 'all', 'select_related', 'prefetch_related', 'extra', 'defer', 'only', 'using', 'select_for_update', 'raw']

    def __init__(self, *args):
        self.querysets = args

    def __getattr__(self, name):
        """
        Deals with all names that are not defined explicitly
        """
        # Check to make sure this is actua
        if name not in models.query.QuerySet.__dict__:
            raise AttributeError("This is not a queryset method")

        # This passing through only works for methods that return a queryset
        if name not in QuerySetSet.qs_methods:
            raise AttributeError("This queryset method does not return a queryset and therefore cannot be passed through to underlying querysets.")

        def map_attr(*args, **kwargs):
            return QuerySetSet(*[getattr(qs, name)(*args, **kwargs) for qs in self.querysets])

        return map_attr

    def order_by(self, *args):
        """
        Not implemented because it's a pain and we don't (yet) need it.
        See here for an example of how to do so: http://ramenlabs.com/2010/12/08/how-to-quack-like-a-queryset/
        """
        raise NotImplementedError("No order_by implemented for querysetset. File a feature request if you need this feature.")

    def count(self):
        """
        Performs a .count() for all subquerysets and returns the number of
        records as an integer.
        """
        return sum(qs.count() for qs in self.querysets)

    def _clone(self):
        "Returns a clone of this queryset chain"
        return self.__class__(*self.querysets)

    def _all(self):
        "Iterates records in all subquerysets"
        return chain(*self.querysets)

    def __getitem__(self, ndx):
        """
        Retrieves an item or slice from the chained set of results from all
        subquerysets.
        """
        if type(ndx) is slice:
            return list(islice(self._all(), ndx.start, ndx.stop, ndx.step or 1))
        else:
            return islice(self._all(), ndx, ndx+1).next()

    def exists(self):
        return any(qs.exists() for qs in self.querysets)


class CsvQuerySet(models.query.QuerySet):
    '''This queryset permits the output of the data from a model as a csv'''

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def expandedCsvHeaders(self):
        '''For some classes, like Compounds and Reactions, there needs to be the option to send additional data (descriptors)
        as part of the csv. This method permits that expansion, and defaults to return the non-expanded headers.'''
        return self.csvHeaders

    @abc.abstractproperty
    def csvHeaders(self):
        '''The basic headers to be used for the model. Note that the implementation on the CsvQuerySet class is extremely basic,
        and will fail if any field holds a relationship, and will not include automagically generated fields.'''
        return [field.name for field in self.model._meta.fields]

    def toCsv(self, writeable, expanded=False, whitelistHeaders=None, missing="?"):#TODO:figure out most sensible default for missing values
        '''Writes the csv data to the writeable (file, or for Django a HttpResponse) object. Expanded outputs any expanded
            information that the corresponding methods provide- this requires the model being called to have a property
            'expandedValues', which should be a dictionary like object of values, using fieldNames as keys as output
            by fetchExpandedHeaders.
        '''

        if expanded:
            headers = self.expandedCsvHeaders
        else:
            headers = self.csvHeaders

        if whitelistHeaders is not None:
            headers = [h for h in headers if h in whitelistHeaders]

        writer = csv.DictWriter(writeable, fieldnames=headers, restval=missing)

        writer.writeheader()
        for row in self.rows(expanded):
            writer.writerow({k:row.get(k, missing) for k in row.keys() if k in headers})

    def rows(self, expanded):
        for item in self:
            yield {field.name:getattr(item, field.name) for field in self.model._meta.fields}

class ArffQuerySet(models.query.QuerySet):
    '''This queryset class permits data from a model to be output as a .arff file'''

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def expandedArffHeaders(self):
        '''returns expanded headers, designed to be overridden by classes that need it'''
        return self.headers

    @abc.abstractproperty
    def arffHeaders(self):
        '''the basic headers to be used for the mode. Note that this imlementation is extremely basic, though not as much so as
        the csv file query set. This will manage foreignkey relations (make sure to define __unicode__ on your models!),
        but won't handle automagic fields or fields to many objects. It will silently ignore
        fields that it does not know how to manage.'''
        headers=OrderedDict()
        for field in self.model._meta.fields:
            if isinstance(field, models.IntegerField) or isinstance(field, models.FloatField) or isinstance(field, models.DecimalField):
                headers[field.name] = '@attribute {} numeric'.format(field.name)
            elif isinstance(field, models.CharField) or isinstance(field, models.TextField):
                headers[field.name] = '@attribute {} string'.format(field.name)
            elif isinstance(field, models.DateTimeField):
                headers[field.name] = '@attribute {} date "yyyy-MM-dd HH:mm:ss"'.format(field.name)
            elif isinstance(field, models.DateField):
                headers[field.name] = '@attribute {} date "yyyy-MM-dd"'.format(field.name)
            elif isinstance(field, models.BooleanField):
                headers[field.name] = '@attribute {} {{True, False}}'.format(field.name)
            elif isinstance(field, models.ForeignKey):
                value_set = {"\"{}\"".format(choice[1]) for choice in field.get_choices()[1:]}
                headers[field.name] = '@attribute {} {{{}}}'.format(field.name, ','.join(value_set))
        return headers

    def toArff(self, writeable, expanded=False, relationName='relation', whitelistHeaders=None, missing="?"):
        '''outputs to an arff file-like object'''
        writeable.write('%arff file generated by the Dark Reactions Project provided by Haverford College\n')
        writeable.write('\n@relation {}\n'.format(relationName))
        if expanded:
            headers = self.expandedArffHeaders
        else:
            headers= self.arffHeaders

        if whitelistHeaders is not None:
            headers = OrderedDict(((k, v) for k, v in headers.items() if k in whitelistHeaders))

        writeable.write('\n'.join(headers.values()))

        writeable.write('\n\n@data\n')
        for row in self.rows(expanded):
            writeable.write(','.join(('"'+str(row.get(key))+'"' if (row.get(key) is not None) else missing) for key in headers.keys()))
            writeable.write('\n')

    def toNPArray(self, expanded=False, whitelistHeaders=None, missing=np.nan):
        '''returns a numpy array'''

        matrix = []
        if expanded:
            headers = self.expandedArffHeaders
        else:
            headers= self.arffHeaders

        if whitelistHeaders is not None:
            headers = OrderedDict(((k, v) for k, v in headers.items() if k in whitelistHeaders))

        for row in self.rows(expanded):
            matrix.append([(row.get(key) if (row.get(key) is not None) else missing) for key in headers.keys()])

        return np.array(matrix)

    def rows(self, expanded=False):
        for item in self:
            yield {field.name:getattr(item, field.name) for field in self.model._meta.fields}
