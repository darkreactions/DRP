class ConfigManager(object):
 def __init__(self):
  self.num_reactants = 5 #The number of reactants supported.
  self.fields_per_reactant = 3 #Each reactant has a name, a quantity, and a unit.

  #Database Page Variable
  self.current_page_radius = 4 #The number of pages to show "around" the current page.
  self.data_per_page = 15 #The number of reactions to show per page.
   #Note: if num_reactants changed, South needs to migrate database.
  self.blacklist = {"x", "-1", -1, "z", "?", "", " "} #Implies absence of data.
  self.unknown_label = "?" #The label that blacklisted values will inherit.
  self.not_required_label = "" #The label that auto-added values will inherit if empty.
  self.jchem_path = "/home/ubuntu/ChemAxon/JChem/bin"
  self.weka_path = "/home/ubuntu/weka-3-7-9/weka.jar"

 def reactant_range(self):
  return xrange(1,self.num_reactants+1)

CONFIG = ConfigManager()
