"""A module containing classes which represent concrete descriptor values.

Descriptor values are not instantiated in and of themselves because they
must be related to another entity- for instance a reaction or a compound.
"""

from descriptors import CategoricalDescriptorPermittedValue
from descriptors import CategoricalDescriptor, NumericDescriptor
from descriptors import BooleanDescriptor, OrdinalDescriptor
from django.core.exceptions import ValidationError
from django.db import models


class CategoricalDescriptorValue(models.Model):

    """A concrete value for a categorical descriptor."""

    descriptor = models.ForeignKey(CategoricalDescriptor)
    """The categorical descriptor to which this pertains."""
    value = models.ForeignKey(CategoricalDescriptorPermittedValue,
                              null=True, blank=True,
                              on_delete=models.PROTECT)
    """The value of the categorical descriptor value.

    Must be an instance of a permitted descriptor value.
    """

    class Meta:
        app_label = 'DRP'
        abstract = True

    def __eq__(self, other):
        """A method for assessing equality conditions.

        In particular this method contains conditions to manage
        Permitted categorical descriptor values which are
        'standing alone', and values which are not
        categorical descriptor values but are other literals
        in python.
        """
        if isinstance(other, CategoricalDescriptorValue):
            return self.value.value == other.value.value
        elif isinstance(other, CategoricalDescriptorPermittedValue):
            return self.value.value == other.value
        else:
            return self.value.value == other

    def clean(self):
        """Ensure the chosen value is a permitted value for the descriptor."""
        if (
           self.value is not None and
           self.value not in self.descriptor.permittedValues.all()
           ):
            raise ValidationError(
                'Invalid Category Described for this Categorical Descriptor',
                'invalid_category')

    def save(self, *args, **kwargs):
        """Ensure cleaning is run on item save."""
        self.clean()
        super(CategoricalDescriptorValue, self).save(*args, **kwargs)

    def __unicode__(self):
        return self.value.value


class BooleanDescriptorValue(models.Model):

    """A concrete value for a boolean descriptor."""

    class Meta:
        app_label = 'DRP'
        abstract = True

    value = models.NullBooleanField('Value for descriptor', null=True)
    """Set to true, false or none (missing value) for instances."""
    descriptor = models.ForeignKey(BooleanDescriptor)
    """The descriptor to which this value pertains."""

    def __nonzero__(self):
        """Correctly assess instances in a boolean context."""
        return self.value

    def clean(self):
        if not (isinstance(self.value, bool) or self.value is None):
            raise ValidationError(
                'Only boolean values are allowed for numeric descriptors',
                'value_wrong_type')


class NumericDescriptorValue(models.Model):

    """A concrete value for a numeric descriptor."""

    class Meta:
        app_label = 'DRP'
        abstract = True

    value = models.FloatField(null=True, blank=True)
    """Set to the floating point value for instances."""
    descriptor = models.ForeignKey(NumericDescriptor)
    """The descriptor to which a value pertains."""

    def clean(self):
        """Ensure that the value is within the prescribed bounds."""
        if not (isinstance(self.value, float) or isinstance(self.value, int) or self.value is None):
            raise ValidationError(
                'Only float or integer values are allowed for numeric descriptors.',
                'value_wrong_type')
        if self.value is not None:
            if self.descriptor.maximum is not None and self.value > self.descriptor.maximum:
                raise ValidationError(
                    'The provided value is higher than the descriptor maximum',
                    'value_too_high')
            if self.descriptor.minimum is not None and self.value < self.descriptor.minimum:
                raise ValidationError(
                    'The provided value is lower than the descriptor minimum',
                    'value_too_low')

    def save(self, *args, **kwargs):
        """Ensure cleaning is run on save."""
        self.clean()
        super(NumericDescriptorValue, self).save(*args, **kwargs)

    def __nonzero__(self):
        """Define boolean context."""
        return bool(self.value)

    def __eq__(self, other):
        """Ensure that equivalent instances are matched as such."""
        if isinstance(other, NumericDescriptorValue):
            return self.value == other.value
        else:
            return self.value == other

    def __gt__(self, other):
        """Greater than comparison to work for instances."""
        if isinstance(other, NumericDescriptorValue):
            return self.value > other.value
        else:
            return self.value > other

    def __lt__(self, other):
        """Less than comparison to work for instances."""
        if isinstance(other, NumericDescriptorValue):
            return self.value < other.value
        else:
            return self.value < other


class OrdinalDescriptorValue(models.Model):

    """Describes a concrete value for an ordinal descriptor."""

    class Meta:
        app_label = 'DRP'
        abstract = True

    value = models.IntegerField(null=True, blank=True)
    """The integer value for the specified descriptor."""
    descriptor = models.ForeignKey(OrdinalDescriptor)

    def clean(self):
        """Ensure the value is within the prescribed bounds."""
        if not (isinstance(self.value, int) or self.value is None):
            raise ValidationError(
                'Only integer values are allowed for numeric descriptors',
                'value_wrong_type')

        if self.value is not None:
            if (
               self.descriptor.maximum is not None and
               self.value > self.descriptor.maximum
               ):
                raise ValidationError(
                    'The provided value is higher than the descriptor maximum',
                    'value_too_high')
            if (
               self.descriptor.minimum is not None and
               self.value < self.descriptor.minimum
               ):
                raise ValidationError(
                    'The provided value is lower than the descriptor minimum',
                    'value_too_low')

    def save(self, *args, **kwargs):
        """Ensure cleaning run on save."""
        self.clean()
        super(OrdinalDescriptorValue, self).save(*args, **kwargs)

    def __eq__(self, other):
        """Ensure that equivalent instances are matched as such."""
        if isinstance(other, OrdinalDescriptorValue):
            return self.value == other.value
        else:
            return self.value == other

    def __gt__(self, other):
        """Greater than comparison to work for instances."""
        if isinstance(other, OrdinalDescriptorValue):
            return self.value > other.value
        else:
            return self.value > other

    def __lt__(self, other):
        """Less than comparison to work for instances."""
        if isinstance(other, OrdinalDescriptorValue):
            return self.value < other.value
        else:
            return self.value < other
