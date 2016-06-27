"""A module containing the reactions descriptors."""
from .descriptors import CategoricalDescriptor, OrdinalDescriptor, BooleanDescriptor
from .descriptors import CategoricalDescriptorPermittedValue, NumericDescriptor, Predictable, DescriptorManager
import DRP.models


class CatRxnDescriptor(CategoricalDescriptor, Predictable):
    """A class which describes a descriptor- a value which describes a system such as a compound or a reaction."""

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Categorical Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        """Initialise a new instance of a reaction descriptor."""
        super(CatRxnDescriptor, self).__init__(*args, **kwargs)
        # because of python's flawed dependency resolution, this is what I've
        # been reduced to.
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredCatRxnDescriptor

    def createValue(self, reaction, value):
        """Create a new reaction descriptor value object."""
        try:
            v = DRP.models.rxnDescriptorValues.CatRxnDescriptorValue.objects.get(
                descriptor=self, reaction=reaction)
        except DRP.models.rxnDescriptorValues.CatRxnDescriptorValue.doesnotExist:
            v = DRP.models.rxnDescriptorValues.CatRxnDescriptorValue(
                descriptor=self, reaction=reaction)
        v.value = CategoricalDescriptorPermittedValue.objects.get(value=value)
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Update the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.

        This allows later bulk creation.
        Return the new value object or None if no object was created (only updated).
        """
        qs = DRP.models.rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(
            descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return DRP.models.rxnDescriptorValues.CatRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)


class OrdRxnDescriptor(OrdinalDescriptor, Predictable):
    """A class which represents an ordinal descriptor."""

    class Meta:
        verbose_name = 'Ordinal Reaction Descriptor'
        app_label = 'DRP'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        """Initialise a new instance of a reaction descriptor."""
        super(OrdRxnDescriptor, self).__init__(*args, **kwargs)
        # because of python's flawed dependency resolution, this is what I've
        # been reduced to.
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredOrdRxnDescriptor

    def createValue(self, reaction, value):
        """Create a new reaction descriptor value object."""
        if not isinstance(value, int) and value is not None:
            raise TypeError(
                "You cannot create a ordinal value with non-integer type {}".format(type(value)))
        try:
            v = DRP.models.rxnDescriptorValues.OrdRxnDescriptorValue.objects.get(
                descriptor=self, reaction=reaction)
        except DRP.models.rxnDescriptorValues.OrdRxnDescriptorValue.DoesNotExist:
            v = DRP.models.rxnDescriptorValues.OrdRxnDescriptorValue(
                descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Update the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.

        This allows later bulk creation.
        Return the new value object or None if no object was created (only updated).
        """
        if not isinstance(value, int) and value is not None:
            raise TypeError(
                "You cannot create a ordinal value with non-integer type {}".format(type(value)))

        qs = DRP.models.rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(
            descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return DRP.models.rxnDescriptorValues.OrdRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)

    def createPredictionDescriptor(self, *args, **kwargs):
        """Create a new predicted descriptor which matches this one's parameters."""
        pred = super(OrdRxnDescriptor, self).createPredictionDescriptor(
            *args, **kwargs)
        pred.maximum = self.maximum
        pred.minimum = self.minimum
        return pred


class NumRxnDescriptor(NumericDescriptor, Predictable):
    """A class which represents a numerical descriptor."""

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Numerical Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        """Initialise a new instance of a reaction descriptor."""
        super(NumRxnDescriptor, self).__init__(*args, **kwargs)
        # because of python's flawed dependency resolution, this is what I've
        # been reduced to.
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredNumRxnDescriptor

    def createValue(self, reaction, value):
        """Create a new reaction descriptor value object."""
        # type thing
        if not (isinstance(value, float) or isinstance(value, int)) and value is not None:
            raise TypeError(
                "You cannot create a numerical value with non-float type {}".format(type(value)))
        try:
            v = DRP.models.rxnDescriptorValues.NumRxnDescriptorValue.objects.get(
                descriptor=self, reaction=reaction)
        except DRP.models.rxnDescriptorValues.NumRxnDescriptorValue.DoesNotExist:
            v = DRP.models.rxnDescriptorValues.NumRxnDescriptorValue(
                descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Update the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.

        This allows later bulk creation.
        Return the new value object or None if no object was created (only updated).
        """
        if not (isinstance(value, float) or isinstance(value, int)) and value is not None:
            raise TypeError(
                "You cannot create a numerical value with non-float type {}".format(type(value)))

        qs = DRP.models.rxnDescriptorValues.NumRxnDescriptorValue.objects.filter(
            descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return DRP.models.rxnDescriptorValues.NumRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)

    def createPredictionDescriptor(self, *args, **kwargs):
        """Create a new predicted descriptor which matches this one's parameters."""
        pred = super(NumRxnDescriptor, self).createPredictionDescriptor(
            *args, **kwargs)
        pred.maximum = self.maximum
        pred.minimum = self.minimum
        return pred


class BoolRxnDescriptor(BooleanDescriptor, Predictable):
    """A class which represents a boolean descriptors."""

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Boolean Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        """Initialise a new instance of a reaction descriptor."""
        super(BoolRxnDescriptor, self).__init__(*args, **kwargs)
        # because of python's flawed dependency resolution, this is what I've
        # been reduced to.
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredBoolRxnDescriptor

    def createValue(self, reaction, value):
        """Create a new reaction descriptor value object."""
        if not isinstance(value, bool) and value is not None:
            raise TypeError(
                "You cannot create a boolean value with non-boolean type {}".format(type(value)))
        try:
            v = DRP.models.rxnDescriptorValues.BoolRxnDescriptorValue.objects.get(
                descriptor=self, reaction=reaction)
        except DRP.models.rxnDescriptorValues.BoolRxnDescriptorValue.DoesNotExist:
            v = DRP.models.rxnDescriptorValues.BoolRxnDescriptorValue(
                descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Update the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.

        This allows later bulk creation.
        Return a tuple of the value object and whether it is new (needs to be saved).
        """
        if not isinstance(value, bool) and value is not None:
            raise TypeError(
                "You cannot create a boolean value with non-boolean type {}".format(type(value)))
        qs = DRP.models.rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(
            descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return DRP.models.rxnDescriptorValues.BoolRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)
