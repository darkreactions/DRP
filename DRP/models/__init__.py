'''Models for the DRP app'''
#Anything imported here gets imported on `from DRP.models import *`
#Not that such behaviour is encouraged!

# Classes that should get included.

from dataSets import DataSet, DataSetRelation
from descriptors import Descriptor, CategoricalDescriptor, BooleanDescriptor, NumericDescriptor, OrdinalDescriptor, CategoricalDescriptorPermittedValue
from rxnDescriptors import CatRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptor
from predRxnDescriptors import PredCatRxnDescriptor, PredOrdRxnDescriptor, PredBoolRxnDescriptor, PredNumRxnDescriptor
from molDescriptors import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor
from molDescriptorValues import CatMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, OrdMolDescriptorValue
from rxnDescriptorValues import RxnDescriptorValue, CatRxnDescriptorValue, BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue
from LabGroup import LabGroup
from LicenseAgreement import LicenseAgreement
from License import License
from Reaction import Reaction
from PerformedReaction import PerformedReaction
from Compound import Compound
from CompoundQuantity import CompoundQuantity
from RecommendedReaction import RecommendedReaction
from StatsModel import StatsModel
from ChemicalClass import ChemicalClass
from ConfirmationCode import ConfirmationCode
from CompoundRole import CompoundRole
from ModelContainer import ModelContainer
