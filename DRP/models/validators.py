from django.core.exceptions import ValidationError
import datetime

def notInTheFuture(value):
    if value is not None and value.date() > datetime.datetime.today().date():
        raise ValidationError('This date cannot be in the future')
