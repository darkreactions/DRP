class ReactionSieve(object):
    def __init__(self, model_container, desired_desc_dict):
        """Create a sieve to remove reactions that do not conform to the desired descriptor dictionary. As these will often be outcome descriptors, requires a model for predictions"""
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

    def filter(self, reactions, verbose=True):
        """Given an iterable of reactions, make model predictions and return reactions that conform to desired descriptor dictionary"""
        if verbose:
            print('''ENTERING FILTER''')

        a = self.model_container.predict(reactions)  # , verbose=True)
        if verbose:
            print('''PREDICTIONS DONE''')

        plausible_reactions = []
        if verbose:
            print('''DESIRED DESCRIPTORS:''')
            print((self.desired_desc_dict))
        for desc, desc_vals in list(self.desired_desc_dict.items()):
            for i in desc_vals:
                for a in a[i]:
                    plausible_reactions.append(a[0])
        plausible_reactions_ids = [pr.id for pr in plausible_reactions]
        reactions = reactions.filter(pk__in=plausible_reactions_ids)

        if verbose:
            print('''EXITING FILTER''')
        return reactions