from DRP.models import OrdRxnDescriptorValue, BoolRxnDescriptorValue


class ReactionSieve(object):
    def __init__(self, model_container, desired_desc_dict):
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

    def filter(self, reactions):
        print 'TYPE', type(reactions)
        print '''ENTERING FILTER'''
        a = self.model_container.predict(reactions)  # , verbose=True)
        print '''PREDICTIONS DONE'''
        plausible_reactions = []
        print '''DESIRED DESCRIPTORS:''', self.desired_desc_dict
        for desc, desc_vals in self.desired_desc_dict.items():
            for i in desc_vals:
                for a in a[i]:
                    plausible_reactions.append(a[0])
                    #         print '''PLAUSIBLE REACTIONS:''', plausible_reactions
        plausible_reactions_ids = [pr.id for pr in plausible_reactions]
        reactions = reactions.filter(pk__range=(min(plausible_reactions_ids), max(plausible_reactions_ids)))
        print '''EXITING FILTER'''
        return reactions