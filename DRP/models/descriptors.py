"""
A module containing abstract base classes for descriptors.

Reaction and Molecular Descriptors should inherit from these
classes.
"""

from django.db import models
from django.core.validators import RegexValidator
from django.core.exceptions import ValidationError
from django.db.models.functions import Concat
<<<<<<< HEAD
=======
from django.contrib.auth.models import User
>>>>>>> e08a9d8bcd64b253b8f31062a7cf280d17bb3a0e


class DescriptorQuerySet(models.query.QuerySet):
    """A queryset to manage a queryset or any of the subtypes thereof."""

    def __init__(self, model=None, **kwargs):
        """Initialise the queryset."""
        model = Descriptor if model is None else model
        super(DescriptorQuerySet, self).__init__(model=model, **kwargs)


class DescriptorManager(models.Manager):
    """A generic manager for all of the model subtypes, which annotates them with a heading generated from separated fields."""

    use_for_related_fields = True

    def get_queryset(self):
        """Generate and appropriately subtype the queryset."""
        return DescriptorQuerySet(model=self.model).annotate(csvHeader=Concat('heading', models.Value('_'), 'calculatorSoftware', models.Value('_'), 'calculatorSoftwareVersion'))


class Descriptor(models.Model):
    """A class which describes a descriptor.

    A descriptor is a classification of values which describe
    a system such as a compound or a reaction.
    """

    class Meta:
        app_label = 'DRP'
        unique_together = (
            'heading',
            'calculatorSoftware',
            'calculatorSoftwareVersion'
        )

    objects = DescriptorManager()

    heading = models.CharField(
        max_length=200,
        validators=[
            RegexValidator(
                '[A-Za-z0-9][A-Za-z0-9_]*',
                ('Please include only values which are limited to '
                 'alphanumeric characters and underscores, and must start '
                 'with an alphabetic character.')
            )
        ]
    )
    """A short label which is given to a description."""
    name = models.CharField('Full name', max_length=300)
    calculatorSoftware = models.CharField(
        max_length=100,
        blank=True,
        validators=[
            RegexValidator(
                '[A-Za-z0-9][A-Za-z0-9_]*',
                ('Please include only values which are limited to '
                 'alphanumeric characters and underscores, and must start '
                 'with an alphabetic character.')
            )
        ]
    )
    calculatorSoftwareVersion = models.CharField(
        max_length=20,
        blank=True,
        validators=[
            RegexValidator(
                '[A-Za-z0-9][A-Za-z0-9_]*',
                ('Please include only values which are limited to '
                 'alphanumeric characters, periods and underscores, and must start '
                 'with an alphabetic character.')
            )
        ]
    )

    @property
    def arffHeader(self):
        """Return the base unit of an Arff Header.

        This method is insufficient and must be overridden by subclasses.
        Details about the Arff file format can be found at
        http://www.cs.waikato.ac.nz/ml/weka/arff.html
        """
        return'@attribute {} ' .format(self.csvHeader)

    def __str__(self):
        """Unicode represenation of a descriptor is it's name."""
        return self.name


class CategoricalDescriptor(Descriptor):
    """A a class of descriptors which are broken up into categories."""

    class Meta:
        app_label = 'DRP'

    objects = DescriptorManager()

    @property
    def arffHeader(self):
        """Complete the Arff header for this descriptor."""
        return super(CategoricalDescriptor, self).arffHeader + '{{{}}}'.format(','.join(str(v.value) for v in self.permittedValues.all()))


class CategoricalDescriptorPermittedValue(models.Model):
    """Each instance is a value that a given descriptor may take."""

    class Meta:
        app_label = 'DRP'
        unique_together = ('descriptor', 'value')

    value = models.CharField('Permitted Value', max_length=255)
    descriptor = models.ForeignKey(
        CategoricalDescriptor,
        related_name='permittedValues'
    )

    def __str__(self):
        """Return the literal value the instance represents."""
        return self.value


class OrdinalDescriptor(Descriptor):
    """A descriptor which is ordinal in nature.

    It may be valued by discrete categories, but can be
    set in an order, such as big, medium and small.

    Values in the DRP database of this kind are
    represented by integers within a limited range.
    """

    class Meta:
        app_label = 'DRP'

    objects = DescriptorManager()

    maximum = models.IntegerField()
    """The maximal permitted value for a given descriptor instance."""
    minimum = models.IntegerField()
    """The minimal permitted value for a given descriptor instance."""

<<<<<<< HEAD
=======
    user = models.ForeignKey(User)
    # changed performedReactions to ratedReactions because of related_name issue (clashing with PerformedReactions FOreignKey)

    ratedBy = models.ForeignKey(User, related_name='ratedReactions', null=True, blank=True, default=None, verbose_name="Rated By")

>>>>>>> e08a9d8bcd64b253b8f31062a7cf280d17bb3a0e
    def clean(self):
        """Special cleaning method. Ensures max < min."""
        if self.maximum is not None and self.minimum is not None and self.maximum < self.minimum:
            raise ValidationError(
                'The maximum value cannot be lower than the minimum value',
                'max_min_mix'
            )

    def save(self, *args, **kwargs):
        """Force cleaning to be run on save."""
        self.clean()
        super(OrdinalDescriptor, self).save(*args, **kwargs)

    @property
    def arffHeader(self):
        """Complete the Arff header for this descriptor."""
        return super(OrdinalDescriptor, self).arffHeader + '{{{}}}'.format(','.join(str(i) for i in range(self.minimum, self.maximum + 1)))


class NumericDescriptor(Descriptor):
    """A descriptor which is numeric in nature.

    Numeric descriptors are stored as floating
    point numbers, and can be either positive
    or negative.
    """

    class Meta:
        app_label = 'DRP'

    objects = DescriptorManager()

    maximum = models.FloatField(null=True)
    """The maximum allowed value for a given descriptor."""
    minimum = models.FloatField(null=True)
    """The minimum allowed value for a given descriptor."""

    def clean(self):
        """Special cleaning method. Ensures max < min."""
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
        """Force cleaning to be run on save."""
        self.clean()
        super(NumericDescriptor, self).save(*args, **kwargs)

    @property
    def arffHeader(self):
        """Complete the Arff header for this descriptor."""
        return super(NumericDescriptor, self).arffHeader + 'numeric'


class BooleanDescriptor(Descriptor):
    """A descriptor which can be represented by either True or False."""

    class Meta:
        app_label = 'DRP'

    objects = DescriptorManager()

    @property
    def arffHeader(self):
        """Complete the Arff header for this descriptor."""
        return super(BooleanDescriptor, self).arffHeader + '{True, False}'


class Predictable(models.Model):
    """Abstract class to add functionality to descriptors about which predictions can be made."""

    class Meta:
        app_label = 'DRP'
        abstract = True

    def createPredictionDescriptor(self, modelContainer, modelComponent=None):
        """Create a descriptor with which to associated predicted values."""
        if not self.pk:
            raise self.DoesNotExist(
                'Cannot create a prediction descriptor of a descriptor which has not yet been saved.')
        else:
            try:
                return self.predictedDescriptorType.objects.get(modelContainer=modelContainer, statsModel=modelComponent, predictionOf=self)
            except self.predictedDescriptorType.DoesNotExist:
                pred = self.predictedDescriptorType()

                if modelComponent is None:
                    headingSuffix = '_prediction_{}_summative'.format(
                        modelContainer.pk)
                    nameSuffix = ' prediction for modelContainer {}'.format(
                        modelContainer.pk)
                else:
                    headingSuffix = '_prediction_{}_component_{}'.format(
                        modelContainer.pk, modelComponent.pk)
                    nameSuffix = ' prediction for modelcontainer {} component {}'.format(
                        modelContainer.pk, modelComponent.pk)

                pred.heading = self.heading + headingSuffix
                pred.name = self.name + nameSuffix
                pred.predictionOf = self
                pred.modelContainer = modelContainer
                pred.statsModel = modelComponent
                return pred
