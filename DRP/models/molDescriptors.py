'''A module containing Classes permitting the representation of molecular descriptors'''
from django.db import models
from descriptors import Descriptor, CategoricalDescriptor, OrdinalDescriptor, BooleanDescriptor
from descriptors import CategoricalDescriptorPermittedValue, NumericDescriptor, DescriptorManager


class CatMolDescriptor(CategoricalDescriptor):
    '''A class which describes a categorical molecular descriptors'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Categorical Molecular Descriptor'

    objects = DescriptorManager()


class OrdMolDescriptor(OrdinalDescriptor):
    '''A class which represents an ordinal descriptor'''

    class Meta:
        verbose_name = 'Ordinal Molecular Descriptor'
        app_label = 'DRP'

    objects = DescriptorManager()


class NumMolDescriptor(NumericDescriptor):
    '''A class which represents a numerical descriptor'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Numerical Molecular Descriptor'

    objects = DescriptorManager()


class BoolMolDescriptor(BooleanDescriptor):
    '''A class which represents a boolean descriptors'''

    class Meta:
        app_label = 'DRP'
        verbose_name = 'Boolean Molecular Descriptor'

    objects = DescriptorManager()
