"""Models for the DRP app."""
# Anything imported here gets imported on `from DRP.models import *`
# Not that such behaviour is encouraged!

# Classes that should get included.

from .dataSets import DataSet, DataSetRelation
from .descriptors import Descriptor, CategoricalDescriptor, BooleanDescriptor, NumericDescriptor, OrdinalDescriptor, CategoricalDescriptorPermittedValue
from .rxnDescriptors import CatRxnDescriptor, OrdRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptor
from .predRxnDescriptors import PredCatRxnDescriptor, PredOrdRxnDescriptor, PredBoolRxnDescriptor, PredNumRxnDescriptor
from .molDescriptors import CatMolDescriptor, BoolMolDescriptor, NumMolDescriptor, OrdMolDescriptor
from .molDescriptorValues import CatMolDescriptorValue, BoolMolDescriptorValue, NumMolDescriptorValue, OrdMolDescriptorValue
from .rxnDescriptorValues import RxnDescriptorValue, CatRxnDescriptorValue, BoolRxnDescriptorValue, NumRxnDescriptorValue, OrdRxnDescriptorValue
from .labGroup import LabGroup
from .licenseAgreement import LicenseAgreement
from .license import License
from .reaction import Reaction
from .performedReaction import PerformedReaction
from .compound import Compound, CompoundGuideEntry
from .compoundQuantity import CompoundQuantity
from .recommendedReaction import RecommendedReaction
from .statsModel import StatsModel
from .chemicalClass import ChemicalClass
from .confirmationCode import ConfirmationCode
from .compoundRole import CompoundRole
from .modelContainer import ModelContainer
from .metricContainer import MetricContainer
from .featureSelectionContainer import FeatureSelectionContainer

