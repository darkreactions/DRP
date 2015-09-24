"""A module containing abstract base classes for descriptors
Reaction and Molecular Descriptors should inherit from these"""

from django.db import models
from django.template.defaultfilters import slugify as _slugify
from django.core.validators import RegexValidator
from django.core.exceptions import ValidationError


def slugify(text):
    """returns a modified version of slug text
    so as to keep compatibility with some external programs
    """
    return _slugify(text).replace('-', '_')


class Descriptor(models.Model):
    """A class which describes a descriptor- a value which describes
    a system such as a compound or a reaction
    """

    class Meta:
        app_label = 'DRP'
        unique_together = (
            'heading',
            'calculatorSoftware',
            'calculatorSoftwareVersion'
        )

    heading = models.CharField(
        max_length=200,
        validators=[
                RegexValidator(
                    '[A-Za-z0-9][A-Za-z0-9_]+',
                    ('Please include only values which are limited to'
                     'alphanumeric characters and underscoresi, and must start'
                     'with an alphabetic character.')
                )
        ]
    )
    """A short label which is given to a description."""
    name = models.CharField('Full name', max_length=300)
    calculatorSoftware = models.CharField(max_length=100)
    calculatorSoftwareVersion = models.CharField(max_length=20)

    @property
    def csvHeader(self):
        return '{}_{}_{}'.format(
            self.heading,
            slugify(self.calculatorSoftware),
            self.calculatorSoftwareVersion
        )

    @property
    def arffHeader(self):
        """returns the base unit of an Arff Header, but this will not be
        sufficient and must be overridden by subclasses
        """
        return'@attribute {} ' .format(self.csvHeader)

    def __unicode__(self):
        return self.name


class CategoricalDescriptor(Descriptor):

    class Meta:
        app_label = 'DRP'

    @property
    def arffHeader(self):
        return (super(CategoricalDescriptor, self).arffHeader +
                '{{{}}}'.format(
                ','.join(str(v.value) for v in self.permittedValues.all())
               ))


class CategoricalDescriptorPermittedValue(models.Model):

    class Meta:
        app_label = 'DRP'
        unique_together = ('descriptor', 'value')

    value = models.CharField('Permitted Value', max_length=255)
    descriptor = models.ForeignKey(
        CategoricalDescriptor,
        related_name='permittedValues'
    )

    def __unicode__(self):
        return self.value


class OrdinalDescriptor(Descriptor):

    class Meta:
        app_label = 'DRP'

    maximum = models.IntegerField(null=True)
    minimum = models.IntegerField(null=True)

    def clean(self):
        if (
            self.maximum is not None and
            self.minimum is not None and
            self.maximum < self.minimum
           ):
            raise ValidationError(
             'The maximum value cannot be lower than the minimum value',
             'max_min_mix'
            )

    def save(self, *args, **kwargs):
        self.clean()
        super(OrdinalDescriptor, self).save(*args, **kwargs)

    @property
    def arffHeader(self):
        return (super(OrdinalDescriptor, self).arffHeader +
                '{{{}}}'.format(
                        ','.join(
                            str(i) for i in range(self.minimum, self.maximum+1)
                        )
                     )
                )


class NumericDescriptor(Descriptor):

    class Meta:
        app_label = 'DRP'

    maximum = models.FloatField(null=True)
    minimum = models.FloatField(null=True)

    def clean(self):
        if (
           self.maximum is not None and
           self.minimum is not None and
           self.maximum < self.minimum
           ):

            raise ValidationError(
                'The maximum value cannot be lower than the minimum value',
                'max_min_mix'
            )

    def save(self, *args, **kwargs):
        self.clean()
        super(NumericDescriptor, self).save(*args, **kwargs)

    @property
    def arffHeader(self):
        return super(NumericDescriptor, self).arffHeader + 'numeric'


class BooleanDescriptor(Descriptor):

    class Meta:
        app_label = 'DRP'

    @property
    def arffHeader(self):
        return super(BooleanDescriptor, self).arffHeader + '{True, False}'
