def update_all_reactions(lab_group):
  from DRP.models import get_lab_Data, update_reaction
  data = get_lab_Data(lab_group)
  for entry in data:
    try:
      update_reaction(reaction, lab_group)
    except:
      print "--could not update reaction: {}".format(entry)
  print "Finished data validation."

