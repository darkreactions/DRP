class ConfigManager(object):
 def __init__(self):
  #Database Setup Variables
  self.num_reactants = 5 #The number of reactants supported.
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
  self.current_license_date = "2013-11-30 11:33:49.612035"
  self.current_license_file = "/license/11_30_13.pdf" #Assumed to be in STATIC directory.

  #General Variable Setup
  self.jchem_path = "/home/ubuntu/ChemAxon/JChem/bin"
  self.weka_path = "/home/ubuntu/weka-3-7-9/weka.jar"

 def reactant_range(self):
  return xrange(1,self.num_reactants+1)

CONFIG = ConfigManager()
