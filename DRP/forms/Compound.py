'''A module containing form for creating compounds'''

import django.forms as forms 
from DRP.models import Compound
from chemspipy import ChemSpider
from django.conf import settings
from django.core.exceptions import ValidationError

class CompoundAdminForm(forms.ModelForm):
  '''Form for the django admin; permits the overriding of CSID absence but forces the existence of the custom flag'''

  class Meta:
    model=Compound
    exclude=('descriptors','custom')

  def save(self, commit=True):
    compound = super(CompoundAdminForm, self).save(commit=False)
    compound.custom=True
    if commit:
      compound.save()
      self.save_m2m()
    return compound


class CompoundForm(forms.ModelForm):
  '''A form for users to add compounds to the compound guide. Forces a check against the chemspider database
  to ensure no spurious compounds make their way into the compound guide.
  '''

  CAS_ID = forms.CharField(label='CAS ID', required=False)
  '''Adding this field, not in the database, allows users to match compounds to a CAS_ID without us incuring issues for storing them'''
  CSID = forms.IntegerField(label='Chemspider ID', min_value=1, error_messages={'required':'This value must be set or selected'})
  '''If the user already knows the right value for this it allows them to skip a step'''


  class Meta:
    fields=('labGroup', 'abbrev', 'CSID', 'name', 'CAS_ID', 'chemicalClass')
    model=Compound
    help_texts = {
      'abbrev':'A local abbreviation by which the compound is known.',
      'name':'A common or IUPAC name for the compound.',
      'CAS_ID':'The CAS number for the compound. Optional.',
      'CSID':'The Chemspider ID for the compound. If this is not included, a list will be provided for you to choose from.' 
    }

  def __init__(self, user, *args, **kwargs):
    '''Overridden version of the init method allows us to place the user's lab groups as a restricted set'''
    super(CompoundForm, self).__init__(*args, **kwargs)
    self.compound = None
    self.chemSpider = ChemSpider(settings.CHEMSPIDER_TOKEN)
    self.fields['labGroup'].queryset = user.labgroup_set.all()
    if user.labgroup_set.all().count() == 1:
      self.fields['labGroup'].empty_label = None

  def clean_CSID(self):
    '''Checks that the CSID is actually a valid id from chemspider'''
    searchResults = self.chemSpider.simple_search(self.cleaned_data['CSID'])
    if(len(searchResults) < 1):
      raise ValidationError('The CSID you have provided is invalid', code='invalid_csid')
    else: 
      self.compound=searchResults[0]
    return self.cleaned_data['CSID']

  def clean(self):
    '''This method verifies that the CSID, CAS_ID (where supplied) and name are consistent'''
    self.cleaned_data = super(CompoundForm, self).clean()
    if self.cleaned_data.get('name'):
      nameResults = self.chemSpider.simple_search(self.cleaned_data['name'])
      if self.cleaned_data.get('CAS_ID') != '':
        CAS_IDResults = self.chemSpider.simple_search(self.cleaned_data['CAS_ID'])
        compoundChoices = [compound for compound in nameResults if compound in CAS_IDResults][0:10]
        #the CAS_ID always generates a more restrictive set
      else:
        compoundChoices = nameResults[0:10]
        #if the CAS_ID is not supplied, then we just create a subset based on the name search alone
  
      if self.compound is None and len(compoundChoices) > 0:
        self.fields['CSID'] = forms.ChoiceField(choices=((choice.csid, choice.common_name) for choice in compoundChoices), widget=forms.widgets.RadioSelect)
        #in essence, if a CSID was not supplied, but the chemspider search returned chemspider results, then we offer those results to the user to make a selection.
        return self.cleaned_data
      elif self.compound is None:
        raise ValidationError('Your search terms failed to validate against the Chemspider database. Please contact a local administrator.', code='no_compounds')
      else:
        if self.compound not in nameResults:
          raise ValidationError('The name provided was not valid for the CSID provided. Please change the entry, or contact your local administrator.', code='name_csid_conflict')
        elif self.cleaned_data.get('CAS_ID') and self.compound not in CAS_IDResults:
          raise ValidationError('The CAS ID provided is not valid for the CSID provided. Remove, replace, or contact your local administrator.', 'name_cas_id_conflict')
        else:
          return self.cleaned_data
    else:
      if self.compound is not None:
        #this is probably some of the most horrible code I have written, but it is the only way to get this to work - Phil.
        data = self.data.copy() #because otherwise the query dict is immutable
        data['name'] = self.compound.common_name #replace the data directly, as bad as that is...
        self._errors['name'] = self.error_class(['Please review this suggestion']) #manually input an error message which is less demanding (this is actually canonical method)
        self.data = data #override the old data
      return self.cleaned_data

  def save(self, commit=True):
    '''Creates (and if appropriate, saves) the compound instance, and adds Inchi and smiles from chemspider'''
    compound = super(CompoundForm, self).save(commit=False)
    csCompound = self.chemSpider.get_compound(compound.CSID)
    compound.INCHI = csCompound.inchi
    compound.smiles = csCompound.smiles
    if commit:
      compound.save()
      self.save_m2m()
    return compound
    
