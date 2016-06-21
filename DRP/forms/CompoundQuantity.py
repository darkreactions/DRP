"""Form related things for Compound Quantity creation."""
from django import forms
from DRP.models import CompoundQuantity, Reaction


def compoundQuantityFormFactory(reactionId):
    """Return a class of compound quantity model form specified to a reaction."""
    class CompoundQuantityForm(forms.ModelForm):

        """A form for creating compound quantities."""

        class Meta:
            model = CompoundQuantity
            fields = ('reaction', 'compound', 'role', 'amount')

        def __init__(self, *args, **kwargs):
            """Restrict the form to a single reaction. Useful for formsets."""
            super(CompoundQuantityForm, self).__init__(*args, **kwargs)
            self.fields['reaction'].widget = forms.HiddenInput()
            self.fields['reaction'].initial = reactionId
            self.fields['reaction'].queryset = Reaction.objects.filter(
                id=reactionId)

    return CompoundQuantityForm
