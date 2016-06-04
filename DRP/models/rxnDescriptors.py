'''A module containing the reactions descriptors'''
from descriptors import CategoricalDescriptor, OrdinalDescriptor, BooleanDescriptor
from descriptors import CategoricalDescriptorPermittedValue, NumericDescriptor, Predictable, DescriptorManager
import rxnDescriptorValues
import DRP.models


class CatRxnDescriptor(CategoricalDescriptor, Predictable):
    '''A class which describes a descriptor- a value which describes a system such as a compound or a reaction'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Categorical Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        super(CatRxnDescriptor, self).__init__(*args, **kwargs)
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredCatRxnDescriptor  # because of python's flawed dependency resolution, this is what I've been reduced to.

    def createValue(self, reaction, value):
        """Create a new reaction value object"""
        try:
            v = rxnDescriptorValues.CatRxnDescriptorValue.objects.get(descriptor=self, reaction=reaction)
        except rxnDescriptorValues.CatRxnDescriptorValue.doesnotExist:
            v = rxnDescriptorValues.CatRxnDescriptorValue(descriptor=self, reaction=reaction)
        v.value = CategoricalDescriptorPermittedValue.objects.get(value=value)
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Updates the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.
        This allows later bulk creation.
        Returns the new value object or None if no object was created (only updated).
        """
        qs = rxnDescriptorValues.CatRxnDescriptorValue.objects.filter(descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return rxnDescriptorValues.CatRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)


class OrdRxnDescriptor(OrdinalDescriptor, Predictable):
    '''A class which represents an ordinal descriptor'''

    class Meta:
        verbose_name = 'Ordinal Reaction Descriptor'
        app_label = 'DRP'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        super(OrdRxnDescriptor, self).__init__(*args, **kwargs)
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredOrdRxnDescriptor  # because of python's flawed dependency resolution, this is what I've been reduced to.

    def createValue(self, reaction, value):
        if not isinstance(value, int) and value is not None:
            raise TypeError("You cannot create a ordinal value with non-integer type {}".format(type(value)))
        try:
            v = rxnDescriptorValues.OrdRxnDescriptorValue.objects.get(descriptor=self, reaction=reaction)
        except rxnDescriptorValues.OrdRxnDescriptorValue.DoesNotExist:
            v = rxnDescriptorValues.OrdRxnDescriptorValue(descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Updates the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.
        This allows later bulk creation.
        Returns the new value object or None if no object was created (only updated).
        """
        if not isinstance(value, int) and value is not None:
            raise TypeError("You cannot create a ordinal value with non-integer type {}".format(type(value)))

        qs = rxnDescriptorValues.OrdRxnDescriptorValue.objects.filter(descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return rxnDescriptorValues.OrdRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)

    def createPredictionDescriptor(self, *args, **kwargs):
        pred = super(OrdRxnDescriptor, self).createPredictionDescriptor(*args, **kwargs)
        pred.maximum = self.maximum
        pred.minimum = self.minimum
        return pred


class NumRxnDescriptor(NumericDescriptor, Predictable):
    '''A class which represents a numerical descriptor'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Numerical Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        super(NumRxnDescriptor, self).__init__(*args, **kwargs)
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredNumRxnDescriptor  # because of python's flawed dependency resolution, this is what I've been reduced to.

    def createValue(self, reaction, value):
        # TODO These checks should be part of a more standard 'allowed value' type thing
        if not (isinstance(value, float) or isinstance(value, int)) and value is not None:
            raise TypeError("You cannot create a numerical value with non-float type {}".format(type(value)))
        try:
            v = rxnDescriptorValues.NumRxnDescriptorValue.objects.get(descriptor=self, reaction=reaction)
        except rxnDescriptorValues.NumRxnDescriptorValue.DoesNotExist:
            v = rxnDescriptorValues.NumRxnDescriptorValue(descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Updates the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.
        This allows later bulk creation.
        Returns the new value object or None if no object was created (only updated).
        """
        if not (isinstance(value, float) or isinstance(value, int)) and value is not None:
            raise TypeError("You cannot create a numerical value with non-float type {}".format(type(value)))

        qs = rxnDescriptorValues.NumRxnDescriptorValue.objects.filter(descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return rxnDescriptorValues.NumRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)

    def createPredictionDescriptor(self, *args, **kwargs):
        pred = super(NumRxnDescriptor, self).createPredictionDescriptor(*args, **kwargs)
        pred.maximum = self.maximum
        pred.minimum = self.minimum
        return pred


class BoolRxnDescriptor(BooleanDescriptor, Predictable):
    '''A class which represents a boolean descriptors'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Boolean Reaction Descriptor'

    objects = DescriptorManager()

    def __init__(self, *args, **kwargs):
        super(BoolRxnDescriptor, self).__init__(*args, **kwargs)
        self.predictedDescriptorType = DRP.models.predRxnDescriptors.PredBoolRxnDescriptor  # because of python's flawed dependency resolution, this is what I've been reduced to.

    def createValue(self, reaction, value):
        if not isinstance(value, bool) and value is not None:
            raise TypeError("You cannot create a boolean value with non-boolean type {}".format(type(value)))
        try:
            v = rxnDescriptorValues.BoolRxnDescriptorValue.objects.get(descriptor=self, reaction=reaction)
        except rxnDescriptorValues.BoolRxnDescriptorValue.DoesNotExist:
            v = rxnDescriptorValues.BoolRxnDescriptorValue(descriptor=self, reaction=reaction)
        v.value = value
        return v

    def updateOrNewValue(self, reaction, value):
        """
        Updates the value for the given reaction, saving if necessary, or creates an *unsaved* value for the reaction.
        This allows later bulk creation.
        Returns a tuple of the value object and whether it is new (needs to be saved).
        """
        if not isinstance(value, bool) and value is not None:
            raise TypeError("You cannot create a boolean value with non-boolean type {}".format(type(value)))
        qs = rxnDescriptorValues.BoolRxnDescriptorValue.objects.filter(descriptor=self, reaction=reaction)

        if qs.exists():
            qs.exclude(value=value).update(value=value)
            return None
        else:
            return rxnDescriptorValues.BoolRxnDescriptorValue(descriptor=self, reaction=reaction, value=value)
