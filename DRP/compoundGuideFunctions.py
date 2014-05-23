
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # # # # #  Compound Guide Helper Functions  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

from DRP.retrievalFunctions import *
from DRP.database_construction import *

#Translate any compounds to the corresponding abbrevs. 
def translate_reactants(lab_group, dataList):
  #Create a map from compounds to abbrevs.
  #(Note that duplicate compounds end up being given the "latest" abbrev.)
  compoundGuide = get_lab_CG(lab_group)
  compoundToAbbrevMap = {entry.compound:entry.abbrev for entry in compoundGuide}

  #Actually "translate" the compounds to abbrevs.
  fields = get_model_field_names(model="Recommendation")
  for (field, i) in zip(fields, range(len(dataList[2][1:]))):
    if field[:8]=="reactant":
        value = dataList[2][i+1]
        if value in compoundToAbbrevMap:
          dataList[2][i+1] = compoundToAbbrevMap[value]

  return dataList
