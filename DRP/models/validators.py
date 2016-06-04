from django.core.exceptions import ValidationError
import datetime
def notInTheFuture(value):
    if value is not None and value.date() > datetime.datetime.today().date():
        raise ValidationError('This date cannot be in the future')

class GreaterThanValidator(object):

    def __init__(self, floorValue):
        self.floorValue = floorValue

    def __call__(self, value):
        if value <= self.floorValue:
            raise ValidationError('This value must be greater than %(value)s', params={'value':self.floorValue})

    def __eq__(self, other):
        return self.floorValue == other.floorValue

    def deconstruct(self):
        return ('DRP.models.validators.GreaterThanValidator', [self.floorValue], {})
