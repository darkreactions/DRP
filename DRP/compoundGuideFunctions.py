
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # #  Compound Guide Helper Functions  # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Translate any compounds to the corresponding abbrevs.
def translate_reactants(lab_group, dataList, single=False, onlyAbbrevs=False, direction="compound to abbrev"):
  from DRP.models import get_lab_CG, get_model_field_names

  #Create a map from compounds to abbrevs.
  #(Note that duplicate compounds end up being given the "latest" abbrev.)
  compoundGuide = get_lab_CG(lab_group)

  #TODO: Don't grab EVERYTHING if not needed...
  if direction=="compound to abbrev":
    translation_table = {entry.compound:entry.abbrev for entry in compoundGuide}
  else:
    translation_table = {entry.abbrev:entry.compound for entry in compoundGuide}

  #Actually "translate" the compounds to abbrevs.
  fields = get_model_field_names(model="Recommendation")

  if single:
    dataList = [["","", dataList]]

  for i, reactionTuple in enumerate(dataList):
    if onlyAbbrevs:
        abbrev = reactionTuple
        if abbrev in translation_table:
          dataList[i] = translation_table[abbrev]
        else:
          print "Compound '{}' not found in translation_table".format(abbrev)

    else:
      for (field, j) in zip(fields, range(len(reactionTuple[2][1:]))):
        if field[:8]=="reactant":
          value = reactionTuple[2][j+1]
          if value in translation_table:
            reactionTuple[2][j+1] = translation_table[value]

  if single:
    return dataList[0][2]

  return dataList

def getMoles(mass, compound):
  from DRP.models import CompoundEntry
  try:
    molar_mass = CompoundEntry.objects.filter(compound=compound)[0].mw
    return mass/float(molar_mass)
  except Exception as e:
    print e
    raise Exception("getMoles: No molar mass available for {}".format(compound))

def getMass(moles, compound):
  from DRP.models import CompoundEntry
  try:
    molar_mass = CompoundEntry.objects.filter(compound=compound)[0].mw
    value = moles*float(molar_mass)
    return float("{:.5f}".format(value))
  except Exception as e:
    print e
    raise Exception("getMass: No molar mass available for {}".format(compound))

