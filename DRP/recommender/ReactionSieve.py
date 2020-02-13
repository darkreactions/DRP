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

        predictions = self.model_container.predict(reactions)  # , verbose=True)
        print("--------")
        print(predictions)
        print("--------")
        if verbose:
            print('''PREDICTIONS DONE''')

        plausible_reactions = {}
        if verbose:
            print('''DESIRED DESCRIPTORS:''')
            print((self.desired_desc_dict))

        for descriptor in self.desired_desc_dict:
            for prediction in predictions[descriptor]:
                if prediction[1].value in self.desired_desc_dict[descriptor]:
                    plausible_reactions[prediction[0].id] = prediction[1].value
        if verbose:
            print('''EXITING FILTER''')
        return plausible_reactions