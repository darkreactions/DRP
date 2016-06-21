"""An assortment of custom validators for model fields."""
from django.core.exceptions import ValidationError
import datetime


def notInTheFuture(value):
    """Verify that a provided date is not in the future."""
    if value is not None and value.date() > datetime.datetime.today().date():
        raise ValidationError('This date cannot be in the future')


class GreaterThanValidator(object):

    """Verify that a value is exclusively greater than a limit."""

    def __init__(self, floorValue):
        """Initialise with the floor value."""
        self.floorValue = floorValue

    def __call__(self, value):
        """Make this object callable."""
        if value <= self.floorValue:
            raise ValidationError('This value must be greater than %(value)s', params={
                                  'value': self.floorValue})

    def __eq__(self, other):
        """Test for equality; allows django some optimisation."""
        return self.floorValue == other.floorValue

    def deconstruct(self):
        """Allow for serialization into django migrations."""
        return ('DRP.models.validators.GreaterThanValidator', [self.floorValue], {})
