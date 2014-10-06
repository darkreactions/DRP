
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # # # # #  Compound Guide Helper Functions  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Translate any compounds to the corresponding abbrevs. 
def translate_reactants(lab_group, dataList):
  from DRP.models import get_lab_CG, get_model_field_names

  #Create a map from compounds to abbrevs.
  #(Note that duplicate compounds end up being given the "latest" abbrev.)
  compoundGuide = get_lab_CG(lab_group)
  compoundToAbbrevMap = {entry.compound:entry.abbrev for entry in compoundGuide}

  #Actually "translate" the compounds to abbrevs.
  fields = get_model_field_names(model="Recommendation")

  for tup in dataList:
    if len(tup)<3: 
      reactionTuple = ["","",tup]
    else:
      reactionTuple = tup
    for (field, i) in zip(fields, range(len(reactionTuple[2][1:]))):
      if field[:8]=="reactant":
        value = reactionTuple[2][i+1]
        if value in compoundToAbbrevMap:
          reactionTuple[2][i+1] = compoundToAbbrevMap[value]

  return dataList
