'''Models for the DRP app'''
#Anything imported here gets imported on `from DRP.models import *`
#Not that such behaviour is encouraged!

# Classes that should get included.
from Compound import Compound
from CompoundQuantity import CompoundQuantity
from DataSet import DataSet
from MolDescriptor import MolDescriptor
from RxnDescriptor import RxnDescriptor
from MolDescriptorValue import MolDescriptorValue
from RxnDescriptorValue import RxnDescriptorValue
from LabGroup import LabGroup
from LegacyStatsModel import LegacyStatsModel
from LicenseAgreement import LicenseAgreement
from License import License
from Reaction import Reaction
from PerformedReaction import PerformedReaction
from RecommendedReaction import RecommendedReaction
from StatsModel import StatsModel
from StatsModelTag import StatsModelTag
from ChemicalClass import ChemicalClass
from ConfirmationCode import ConfirmationCode
