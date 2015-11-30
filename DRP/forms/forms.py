'''The django forms and model forms used in the DRP
Classes:

LabGroupForm: for creating Lab Groups in the django admin.
ContactForm: A very simple form for the contact page.
'''
#from django.contrib.auth.models import User
import django.forms as forms
from django.contrib.auth.models import User
from django.contrib.auth.hashers import make_password
from django.contrib.auth import authenticate
from DRP.models import LabGroup, License, LicenseAgreement
from django.core.exceptions import ValidationError
from django.contrib.auth.forms import UserCreationForm as DjangoUserCreationForm
from django.contrib.auth.forms import AuthenticationForm as DjangoAuthenticationForm
from django.utils.safestring import mark_safe
from django.utils.html import conditional_escape
from DRP.models import License

#from DRP.models import Data, Lab_Group, Lab_Member, CompoundEntry
#from DRP.settings import ACCESS_CODE_LENGTH
#from DRP.validation import validate_CG, full_validation
#from DRP.validation import TYPE_CHOICES, UNIT_CHOICES
#from DRP.validation import BOOL_CHOICES, PURITY_CHOICES, OUTCOME_CHOICES



#class CompoundGuideForm(forms.ModelForm):
#  compound = forms.CharField(widget=forms.TextInput(
#              attrs={'class':'form_text',
#              "title":"What the abbreviation stands for."}))
#  abbrev = forms.CharField(widget=forms.TextInput(
#            attrs={'class':'form_text form_text_short',
#            "title":"The abbreviation you want to type."}))
#  CAS_ID = forms.CharField(label="CAS ID", widget=forms.TextInput(
#          attrs={'class':'form_text',
#          "title":"The CAS ID of the compound if available."}),
#          required=False)
#  compound_type = forms.ChoiceField(label="Type", choices = TYPE_CHOICES,
#                  widget=forms.Select(attrs={'class':'form_text dropDownMenu',
#                  "title":"Choose the compound type: <br/> --Organic <br/> --Inorganic<br/>--pH Changing<br/>--Oxalate-like<br/>--Solute<br/>--Water"}))
#  custom = forms.CharField(widget=forms.HiddenInput(attrs={"id":"input_custom"}), required=False)
#
#  class Meta:
#    model = CompoundEntry
#    exclude = ("lab_group", "calculations", "image_url", "smiles",
#              "mw", "calculations", "calculations_failed")
#
#  def __init__(self, lab_group=None, *args, **kwargs):
#    super(CompoundGuideForm, self).__init__(*args, **kwargs)
#    self.lab_group = lab_group
#
#  def save(self, commit=True):
#    entry = super(CompoundGuideForm, self).save(commit=False)
#    entry.lab_group = self.lab_group
#
#    if commit:
#      entry.save()
#    return entry
#
#  def clean(self):
#    #Initialize the variables needed for the cleansing process.
#    dirty_data = super(CompoundGuideForm, self).clean() #Get the available raw (dirty) data
#    #Gather the clean_data and any errors found.
#    clean_data, gathered_errors = validate_CG(dirty_data, self.lab_group)
#    form_errors = {field: self.error_class([message]) for (field, message) in gathered_errors.iteritems()}
#
#    #Apply the errors to the form.
#    self._errors.update(form_errors)
#    return clean_data
#
#
#class DataEntryForm(forms.ModelForm):
#
#  reactant_attrs = {'class':'form_text autocomplete_reactant',
#   'title':'Enter the name of the reactant.'}
#  quantity_attrs = {'class':'form_text form_text_short', 'placeholder':'Amount',
#   'title':'Enter the amount of reactant.'}
#  unit_attrs = {'class':'form_text dropDownMenu',
#   'title':'\"g\": gram <br/> \"mL\": milliliter <br/> \"d\": drop'}
#
#
#  reactant_fk_1 = forms.CharField(required=True, label='Reactant 1',
#                    widget=forms.TextInput( attrs=reactant_attrs))
#  quantity_1 = forms.CharField(required=True, label='Quantity 1',
#                     widget=forms.TextInput( attrs=quantity_attrs))
#  unit_1 = forms.ChoiceField(required=True, choices = UNIT_CHOICES,
#                 widget=forms.Select( attrs=unit_attrs))
#
#  reactant_fk_2 = forms.CharField(required=True, label='Reactant 2',
#                    widget=forms.TextInput( attrs=reactant_attrs))
#  quantity_2 = forms.CharField(required=True, label='Quantity 2',
#                     widget=forms.TextInput( attrs=quantity_attrs))
#  unit_2 = forms.ChoiceField(required=True, choices = UNIT_CHOICES,
#                 widget=forms.Select( attrs=unit_attrs))
#
#  reactant_fk_3 = forms.CharField(required=True, label='Reactant 3',
#                    widget=forms.TextInput( attrs=reactant_attrs))
#  quantity_3 = forms.CharField(required=True, label='Quantity 3',
#                     widget=forms.TextInput( attrs=quantity_attrs))
#  unit_3 = forms.ChoiceField(required=True, choices = UNIT_CHOICES,
#                 widget=forms.Select( attrs=unit_attrs))
#
#  reactant_fk_4 = forms.CharField(required=True, label='Reactant 4',
#                    widget=forms.TextInput( attrs=reactant_attrs))
#  quantity_4 = forms.CharField(required=True, label='Quantity 4',
#                     widget=forms.TextInput( attrs=quantity_attrs))
#  unit_4 = forms.ChoiceField(required=True, choices = UNIT_CHOICES,
#                 widget=forms.Select( attrs=unit_attrs))
#
#  reactant_fk_5 = forms.CharField(required=True, label='Reactant 5',
#                    widget=forms.TextInput( attrs=reactant_attrs))
#  quantity_5 = forms.CharField(required=True, label='Quantity 5',
#                     widget=forms.TextInput( attrs=quantity_attrs))
#  unit_5 = forms.ChoiceField(required=True, choices = UNIT_CHOICES,
#                 widget=forms.Select( attrs=unit_attrs))
#
#
#  ref = forms.CharField(label="Ref.", widget=forms.TextInput(
#      attrs={'class':'form_text form_text_short',
#      "title":"The lab notebook and page number where the data can be found."}))
#  temp = forms.CharField(label="Temp.", widget=forms.TextInput(
#      attrs={'class':'form_text form_text_short', 'placeholder':'Celsius',
#      "title":"The temperature at which the reaction took place."}))
#  time = forms.CharField(widget=forms.TextInput(
#    attrs={'class':'form_text form_text_short', 'placeholder':'Hours',
#    "title":"How long the reaction was allowed to occur."}))
#  pH = forms.CharField(label="pH", required=False, widget=forms.TextInput(
#      attrs={'class':'form_text form_text_short', 'placeholder':'0 - 14',
#      "title":"The pH at which the reaction occurred."}))
#  slow_cool = forms.ChoiceField(label="Slow Cool", choices = BOOL_CHOICES,
#      widget=forms.Select(attrs={'class':'form_text dropDownMenu',
#      "title":"Was the reaction allowed to slow-cool?"}))
#  leak = forms.ChoiceField(choices = BOOL_CHOICES, widget=forms.Select(
#      attrs={'class':'form_text dropDownMenu',
#      "title":"Was a leak present during the reaction?"}))
#  outcome = forms.ChoiceField(choices = OUTCOME_CHOICES, widget=forms.Select(
#      attrs={'class':'form_text dropDownMenu',
#      "title":"?: Missing<br/> 1: No Solid<br/> 2: Noncrystalline/Brown<br/>3: Powder/Crystallites<br/>4: Large Single Crystals"}))
#  purity = forms.ChoiceField(choices = PURITY_CHOICES, widget=forms.Select(
#      attrs={'class':'form_text dropDownMenu',
#      "title":"?: Missing<br/> 1: Multiphase<br/> 2: Single Phase"}))
#  notes = forms.CharField(required = False, widget=forms.TextInput(
#      attrs={'class':'form_text form_text_long',
#      "title":"Additional notes about the reaction."}))
#
#  duplicate_of = forms.CharField(required = False,
#      label="Duplicate of", widget=forms.TextInput(
#      attrs={'class':'form_text form_text_short',
#      "title":"The reaction reference of which this data is a duplicate."}))
#  recommended = forms.ChoiceField(choices = BOOL_CHOICES, widget=forms.Select(
#      attrs={'class':'form_text dropDownMenu',
#      "title":"Did we recommend this reaction to you?"}))
#
#  class Meta:
#    model = Data
#    exclude = ("user","lab_group", "creation_time_dt")
#    #Set the field order.
#    fields = [
#    "reactant_fk_1", "quantity_1", "unit_1",
#    "reactant_fk_2", "quantity_2", "unit_2",
#    "reactant_fk_3", "quantity_3", "unit_3",
#    "reactant_fk_4", "quantity_4", "unit_4",
#    "reactant_fk_5", "quantity_5", "unit_5",
#    "ref", "temp", "time", "pH", "slow_cool",
#    "leak", "outcome", "purity",
#    "duplicate_of", "recommended","notes"
#    ]
#
#  def __init__(self, user=None, *args, **kwargs):
#    super(DataEntryForm, self).__init__(*args, **kwargs)
#
#    import datetime
#
#    if user:
#      self.user = user
#      self.lab_group = user.get_profile().lab_group
#    self.creation_time_dt = datetime.datetime.now()
#
#  def save(self, commit=True):
#    datum = super(DataEntryForm, self).save(commit=False)
#    datum.user = self.user
#    datum.lab_group = self.lab_group
#    datum.creation_time_dt = self.creation_time_dt
#
#    try:
#      datum.get_calculations_list()
#      datum.is_valid = True
#    except:
#      datum.is_valid = False
#
#
#    if commit:
#      datum.save()
#
#    return datum
#
#
#  #Clean the data that is input using the form.
#  def clean(self):
#    #Initialize the variables needed for the cleansing process.
#    dirty_data = super(DataEntryForm, self).clean() #Get the available raw (dirty) data
#    print dirty_data
#    #Gather the clean_data and any errors found.
#    clean_data, gathered_errors = full_validation(dirty_data, self.lab_group)
#    form_errors = {field: self.error_class([message]) for (field, message) in gathered_errors.iteritems()}
#
#    #Apply the errors to the form.
#    self._errors = form_errors
#
#    #Add the non-input information to the clean data package:
#    clean_data["lab_group"] = self.lab_group
#    clean_data["user"] = self.user
#    clean_data["creation_time_dt"] = self.creation_time_dt
#
#    return clean_data
