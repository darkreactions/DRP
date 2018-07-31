class ReactionSieve(object):
    def __init__(self, model_container, desired_desc_dict):
        """Create a sieve to remove reactions that do not conform to the desired descriptor dictionary. As these will often be outcome descriptors, requires a model for predictions"""
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

    def filter(self, reactions, verbose=True):
        """Given an iterable of reactions, make model predictions and returns a dictionary of reaction id to outcome according
            to the desired descriptor dictionary"""
        verbose = True
        if verbose:
            print('''ENTERING FILTER''')

        a = self.model_container.predict(reactions)  # , verbose=True)
        print("--------")
        print(a)
        print("--------")
        if verbose:
            print('''PREDICTIONS DONE''')

        plausible_reactions = {}
        if verbose:
            print('''DESIRED DESCRIPTORS:''')
            print((self.desired_desc_dict))
        # reactions_from_dict = list(a.items())[0][1]
        # for tuple in reactions_from_dict:
        #     plausible_reactions.append(tuple[0])

        # # NOT SURE WHAT THIS DOES
        # # for desc, desc_vals in list(self.desired_desc_dict.items()):
        # #     for i in desc_vals:
        # #         for a in a[i]:
        # #             plausible_reactions.append(a[0])
        # plausible_reactions_ids = [pr.id for pr in plausible_reactions]
        # reactions = reactions.filter(pk__in=plausible_reactions_ids)

        for descriptor in self.desired_desc_dict:
            for prediction in a[descriptor]:
                if prediction[1].value in self.desired_desc_dict[descriptor]:
                    plausible_reactions[prediction[0].id] = prediction[1].value
        # plausible_reactions_ids = [pr.id for pr in plausible_reactions]
        # reactions = reactions.filter(pk__in=plausible_reactions_ids)

        

        if verbose:
            print('''EXITING FILTER''')
        return plausible_reactions