from models import *

from django.contrib.auth.forms import UserCreationForm

class UserForm(ModelForm):
 username = CharField(label="Username", required=True,
  widget=TextInput(attrs={'class':'form_text'}))
 password = CharField(label="Password", required=True,
  widget=PasswordInput(attrs={'class':'form_text'}))
 first_name = CharField(label="First Name", required=True,
  widget=TextInput(attrs={'class':'form_text'}))
 last_name = CharField(label="Last Name", required=True,
  widget=TextInput(attrs={'class':'form_text'}))
 email = EmailField(label="Email", required=True,
  widget=TextInput(attrs={'class':'form_text email_text'}))

 class Meta:
  model = User
  fields = ("username", "email",
   "first_name", "last_name",
   "password")

 #Hash the user's password upon save.
 def save(self, commit=True):
  user = super(UserForm, self).save(commit=False)
  user.set_password(self.cleaned_data["password"])
  if commit:
   user.save()
  return user

class UserProfileForm(ModelForm):
 #Enumerate all of the lab titles.
 lab_group = ModelChoiceField(queryset=Lab_Group.objects.all(),
  label="Lab Group", required=True,
  widget=Select(attrs={'class':'form_text', 'title':'Choose which lab you would like to join.'}))
 access_code = CharField(label="Access Code", required=True,
  max_length=ACCESS_CODE_MAX_LENGTH,
  widget=TextInput(attrs={'class':'form_text', 'title':'The unique code given to you by a lab administrator.'}))

 class Meta:
  model = Lab_Member
  app_label = "formsite"
  fields = ["lab_group"]

class LabForm(ModelForm):
 lab_title = CharField(label="Lab Name", required=True,
  widget=TextInput(attrs={'class':'form_text', 'title':'The unique name of your lab.'}))
 lab_address = CharField(label="Address", required=True,
  widget=TextInput(attrs={'class':'form_text', 'title':'This helps us see how the Dark Reactions Project is spreading.'}))
 lab_email = EmailField(label="Contact Email", required=True,
  widget=TextInput(attrs={'class':'form_text email_text', 'title':'We will send you a verification code through this email.'}))


 class Meta:
  model = Lab_Group
  fields = ("lab_title", "lab_address", "lab_email")

 def save(self, commit=True):
  lab_group = super(LabForm, self).save(commit=False)
  lab_group.access_code = get_random_code() 
  if commit:
   lab_group.save()
  return lab_group


class CompoundGuideForm(ModelForm):
 compound = CharField(widget=TextInput(
  attrs={'class':'form_text',
  "title":"What the abbreviation stands for."}))
 abbrev = CharField(widget=TextInput(
  attrs={'class':'form_text form_text_short',
  "title":"The abbreviation you want to type."}))
 CAS_ID = CharField(label="CAS ID", widget=TextInput(
  attrs={'class':'form_text',
  "title":"The CAS ID of the compound if available."}),
  required=False)
 compound_type = ChoiceField(label="Type", choices = TYPE_CHOICES,
  widget=Select(attrs={'class':'form_text dropDownMenu',
  "title":"Choose the compound type: <br/> --Organic <br/> --Inorganic<br/>--pH Changing<br/>--Oxalate-like<br/>--Solute<br/>--Water"}))
 custom = CharField(widget=HiddenInput(attrs={"id":"input_custom"}), required=False)

 class Meta:
  model = CompoundEntry
  exclude = ("lab_group", "calculations", "image_url", "smiles", "mw", "calculations")

 def __init__(self, lab_group=None, *args, **kwargs):
  super(CompoundGuideForm, self).__init__(*args, **kwargs)
  self.lab_group = lab_group

 def save(self, commit=True):
  entry = super(CompoundGuideForm, self).save(commit=False)
  entry.lab_group = self.lab_group

  if commit:
   entry.save()
  return entry

 def clean(self):
  print "starting_cleaning"
  #Initialize the variables needed for the cleansing process.
  dirty_data = super(CompoundGuideForm, self).clean() #Get the available raw (dirty) data

  #Gather the clean_data and any errors found.
  clean_data, gathered_errors = validate_CG(dirty_data, self.lab_group)
  form_errors = {field: self.error_class([message]) for (field, message) in gathered_errors.iteritems()}

  #Apply the errors to the form.
  self._errors.update(form_errors)
  print "finished cleaning!"
  return clean_data


class DataEntryForm(ModelForm):
 #List Fields
 for i in CONFIG.reactant_range():
  if i > 2:
   required="required=False,"
  else:
   required=""
  exec("reactant_{0} = CharField({1} label='Reactant {0}', widget=TextInput(".format(i, required) +
   "attrs={'class':'form_text autocomplete_reactant',"+
   "'title':'Enter the name of the reactant.'}))")
  exec("quantity_{0} = CharField({1} label='Quantity {0}', widget=TextInput(".format(i, required) +
   "attrs={'class':'form_text form_text_short', 'placeholder':'Amount',"+
   "'title':'Enter the amount of reactant.'}))")
  exec("unit_{0} = ChoiceField({1} choices = UNIT_CHOICES, widget=Select(".format(i, required) +
   "attrs={'class':'form_text dropDownMenu',"+
   "'title':'\"g\": gram <br/> \"mL\": milliliter <br/> \"d\": drop'}))")
 ref = CharField(label="Ref.", widget=TextInput(
  attrs={'class':'form_text form_text_short',
  "title":"The lab notebook and page number where the data can be found."}))
 temp = CharField(label="Temp.", widget=TextInput(
  attrs={'class':'form_text form_text_short', 'placeholder':'Celsius',
  "title":"The temperature at which the reaction took place."}))
 time = CharField(widget=TextInput(
  attrs={'class':'form_text form_text_short', 'placeholder':'Hours',
  "title":"How long the reaction was allowed to occur."})) ###
 pH = CharField(label="pH", widget=TextInput(
  attrs={'class':'form_text form_text_short', 'placeholder':'0 - 14',
  "title":"The pH at which the reaction occurred."}))
 slow_cool = ChoiceField(label="Slow Cool", choices = BOOL_CHOICES,
  widget=Select(attrs={'class':'form_text dropDownMenu',
  "title":"Was the reaction allowed to slow-cool?"}))
 leak = ChoiceField(choices = BOOL_CHOICES, widget=Select(
  attrs={'class':'form_text dropDownMenu',
  "title":"Was a leak present during the reaction?"}))
 outcome = ChoiceField(choices = OUTCOME_CHOICES, widget=Select(
  attrs={'class':'form_text dropDownMenu',
  "title":"0: No Data Available <br/> 1: No Solid<br/> 2: Noncrystalline/Brown<br/>3: Powder/Crystallites<br/>4: Large Single Crystals"}))
 purity = ChoiceField(choices = PURITY_CHOICES, widget=Select(
  attrs={'class':'form_text dropDownMenu',
  "title":"0: No Data Available<br/> 1: Multiphase<br/> 2: Single Phase"}))
 notes = CharField(required = False, widget=TextInput(
  attrs={'class':'form_text form_text_long',
  "title":"Additional notes about the reaction."}))

 duplicate_of = CharField(required = False,
  label="Duplicate of", widget=TextInput(
  attrs={'class':'form_text form_text_short',
  "title":"The reaction reference of which this data is a duplicate."}))
 recommended = ChoiceField(choices = BOOL_CHOICES, widget=Select(
  attrs={'class':'form_text dropDownMenu',
  "title":"Did we recommend this reaction to you?"}))

 class Meta:
  model = Data
  exclude = ("user","lab_group", "creation_time")
  #Set the field order.
  fields = [
   "reactant_1", "quantity_1", "unit_1",
   "reactant_2", "quantity_2", "unit_2",
   "reactant_3", "quantity_3", "unit_3",
   "reactant_4", "quantity_4", "unit_4",
   "reactant_5", "quantity_5", "unit_5",
   "ref", "temp", "time", "pH", "slow_cool",
   "leak", "outcome", "purity",
   "duplicate_of", "recommended","notes"
  ]

 def __init__(self, user=None, *args, **kwargs):
  super(DataEntryForm, self).__init__(*args, **kwargs)

  if user:
   self.user = user
   self.lab_group = user.get_profile().lab_group
  self.creation_time = str(datetime.datetime.now())

 def save(self, commit=True):
  print "1"
  datum = super(DataEntryForm, self).save(commit=False)
  print "2"
  datum.user = self.user
  datum.lab_group = self.lab_group
  datum.creation_time = self.creation_time
  datum.is_valid = True #If validation succeeded and data is saved, then it is_valid.
  datum.atoms = "".join(get_atom_set_from_reaction(datum))

  if commit:
   datum.save()
  return datum


 #Clean the data that is input using the form.
 def clean(self):
  #Initialize the variables needed for the cleansing process.
  dirty_data = super(DataEntryForm, self).clean() #Get the available raw (dirty) data

  #Gather the clean_data and any errors found.
  clean_data, gathered_errors = full_validation(dirty_data, self.lab_group)
  form_errors = {field: self.error_class([message]) for (field, message) in gathered_errors.iteritems()}

  #Apply the errors to the form.
  self._errors.update(form_errors)

  #Add the non-input information to the clean data package:
  clean_data["lab_group"] = self.lab_group
  clean_data["user"] = self.user
  clean_data["creation_time"] = self.creation_time

  return clean_data
