'''The django forms and model forms used in the DRP
Classes.

'''
from LabGroup import LabGroupForm, LabGroupJoiningForm, LabGroupSelectionForm
from Contact import ContactForm
from FilterForm import FilterForm, FilterFormSet, filterFormSetFactory
from compound import CompoundForm, CompoundAdminForm, CompoundEditForm, CompoundDeleteForm, CompoundUploadForm, CompoundFilterForm
from compound import CompoundFilterFormSet, AdvancedCompoundFilterForm, AdvancedCompoundFilterFormSet
from authentication import ConfirmationForm, LicenseAgreementForm, UserCreationForm
from PerformedReaction import PerformedRxnAdminForm, PerformedRxnForm
from descriptor import CatRxnDescriptorForm, BoolRxnDescriptorForm, NumRxnDescriptorForm, OrdRxnDescriptorForm, CatDescPermittedValueForm
from descriptorValues import NumRxnDescValFormSet, OrdRxnDescValFormSet, BoolRxnDescValFormSet, CatRxnDescValFormSet  
from FormSet import FormSet, ModelFormSet, FormSetManagerForm 
