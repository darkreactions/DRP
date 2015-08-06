'''A module containing form for creating compounds'''

import django.forms as forms 
from DRP.models import Compound
from chemspipy import ChemSpider
from django.conf import settings
from django.core.exceptions import ValidationError

class CompoundForm(forms.ModelForm):

  CAS_ID = forms.CharField(label='CAS ID', required=False)
  CSID = forms.IntegerField(min_value=1, error_messages={'required':'This value must be set or selected'})

  chemSpider = ChemSpider(settings.CHEMSPIDER_TOKEN)
  compound = None

  class Meta:
    fields=('abbrev', 'name', 'CSID', 'CAS_ID', 'chemicalClass')
    model=Compound
    help_texts = {
      'abbrev':'A local abbreviation by which the compound is known.',
      'name':'A common or IUPAC name for the compound.',
      'CAS_ID':'The CAS number for the compound. Optional.',
      'CSID':'The Chemspider ID for the compound. If this is not included, a list will be provided for you to choose from.' 
    }

  def clean_CSID(self):
    searchResults = self.chemSpider.simple_search(self.cleaned_data['CSID'])
    if(len(searchResults) == 0):
      raise ValidationError('The CSID you have provided is invalid', code='invalid_csid')
    else: 
      self.compound=searchResults[0]
    return self.cleaned_data['CSID']

  def clean(self):
    if 'name' in self.cleaned_data.keys():
      nameResults = self.chemSpider.simple_search(self.cleaned_data['name'])
      if self.cleaned_data['CAS_ID']:
        CAS_IDResults = self.chemSpider.simple_search(self.cleaned_data['CAS_ID'])
        compoundChoices = [compound for compound in nameResults if compound in CAS_IDResults][1:10]
      else:
        compoundChoices = nameResults[1:10]
        CAS_IDResults = []
  
      if self.compound is None and len(compoundChoices) > 0:
        self.fields['CSID'] = forms.ChoiceField(choices=((choice.csid, choice.common_name) for choice in compoundChoices), widget=forms.widgets.RadioSelect)
        return self.cleaned_data
      elif self.compound is None:
        raise ValidationError('Your search terms failed to validate against the Chemspider database. Please contact a local administrator.', code='no_compounds')
      else:
        if self.compound not in nameResults:
          raise ValidationError('The name provided was not valid for the CSID provided. Please change the entry, or contact your local administrator.', code='name_csid_conflict')
        elif self.cleaned_data['CAS_ID'] and self.compound not in CAS_IDResults:
          raise ValidationError('The CAS ID provided is not valid for the CSID provided. Remove, replace, or contact your local administrator.', 'name_cas_id_conflict')
        else:
          return self.cleaned_data
    else:
      return self.cleaned_data        
