import datetime
strptime = datetime.datetime.strptime

class ConfigManager(object):
 def __init__(self):
  #Database Setup Variables
  self.num_reactants = 5 #The number of reactants supported.
  self.reactants_required = 2 #The number of reactants required.
  self.fields_per_reactant = 3 #Each reactant has a name, a quantity, and a unit.
   #Note: if num_reactants changes, database must be migrated.

  #Database Page View Variables
  self.current_page_radius = 3 #The number of pages to show "around" the current page.
  self.data_per_page = 15 #The number of reactions to show per page.

  #Data Upload Variables
  self.blacklist = {"x", "-1", -1, "z", "?", "", " "} #Implies absence of data.
  self.unknown_label = "?" #The label that blacklisted values will inherit.
  self.not_required_label = "" #The label that auto-added values will inherit if empty.

  #Licensing/legal Setup
  self.raw_license_date = "2014-01-20 01:00:00.000000"
  self.current_license_date = strptime(self.raw_license_date, "%Y-%m-%d %X.%f")
  self.current_license_file = "/licenses/01_20_14.pdf" #Assumed to be in STATIC directory.

  #General Variable Setup
  self.jchem_path = "/home/drp/ChemAxon/JChem/bin"
  self.weka_path = "/home/drp/weka/weka.jar"

 def reactant_range(self):
  return xrange(1,self.num_reactants+1)

CONFIG = ConfigManager()
