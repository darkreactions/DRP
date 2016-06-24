"""A module containing code pertaining to forms for the reaction classes."""
from django import forms
from django.core.exceptions import ValidationError
from DRP.models import PerformedReaction, RecommendedReaction
from DRP.models import DataSetRelation
from django.contrib.auth.models import User
from django.forms.widgets import HiddenInput


class PerformedRxnAdminForm(forms.ModelForm):
    """
    A form for the performed reactions in the Django admin.

    The only things we want editing from there are whether a Reaction is invalid or if it should have the legacy recommended flag set.
    """

    class Meta:
        fields = ('valid', 'legacyRecommendedFlag')
        model = PerformedReaction


class PerformedRxnForm(forms.ModelForm):
    """A form for creating performed reaction instances in teh databases."""

    class Meta:
        fields = ('reference', 'notes', 'performedBy', 'labGroup', 'duplicateOf',
                  'performedDateTime', 'public', 'valid', 'recommendation')
        model = PerformedReaction

    def __init__(self, user, *args, **kwargs):
        """Overridden __init__ method; requires the user as the first argument so that choice of lab group etc can be validated, as well as to track who enters what."""
        super(PerformedRxnForm, self).__init__(*args, **kwargs)
        self.user = user
        labGroups = user.labgroup_set.all()
        self.fields['labGroup'].queryset = labGroups
        self.fields['recommendation'].queryset = RecommendedReaction.objects.filter(
            labGroup__in=labGroups)
        self.fields['recommendation'].widget = forms.HiddenInput()
        self.fields['duplicateOf'].queryset = PerformedReaction.objects.filter(
            labGroup__in=labGroups) | PerformedReaction.objects.filter(public=True)
        self.fields['performedBy'].queryset = User.objects.filter(
            labgroup__in=labGroups)
        if not kwargs.get('instance', False):
            self.fields['valid'].initial = False
            # a little hacky, but this is faster than making another form...
            self.fields['valid'].widget = forms.HiddenInput()
        if labGroups.exists():
            self.fields['labGroup'].empty_label = None

    def save(self, commit=True, *args, **kwargs):
        """Overriden save method automates addition of user that created this instance."""
        rxn = super(PerformedRxnForm, self).save(commit=False)
        rxn.user = self.user
        if commit:
            rxn.save()
        return rxn


class PerformedRxnDeleteForm(forms.Form):
    """A form for deleting a reaction."""

    def __init__(self, user, instance=None, *args, **kwargs):
        """Limit the form to this one reaction."""
        super(PerformedRxnDeleteForm, self).__init__(*args, **kwargs)
        self.fields['id'] = forms.ModelChoiceField(queryset=PerformedReaction.objects.filter(
            labGroup__in=user.labgroup_set.all()), initial=None if instance is None else instance.pk, widget=HiddenInput)

    def clean_id(self):
        """Protect the reaction if it is in use in the model pipeline."""
        if DataSetRelation.objects.filter(reaction__pk=self.cleaned_data['id']).exists():
            raise ValidationError(
                "This reaction is protected from deletion because it is used in one or more reactions or recommendations.")
        return self.cleaned_data['id']

    def save(self):
        """Save changes in an ironic way."""
        self.cleaned_data['id'].delete()
        return self.cleaned_data['id']


class PerformedRxnInvalidateForm(forms.Form):
    """A form for Invalidating a reaction."""

    def __init__(self, user, instance=None, *args, **kwargs):
        """Limit the form to this one reaction."""
        super(PerformedRxnInvalidateForm, self).__init__(*args, **kwargs)
        self.fields['id'] = forms.ModelChoiceField(queryset=PerformedReaction.objects.filter(
            labGroup__in=user.labgroup_set.all()), initial=None if instance is None else instance.pk, widget=HiddenInput)

    def save(self):
        """Save the change."""
        self.cleaned_data['id'].valid = False
        self.cleaned_data['id'].save()
        return self.cleaned_data['id']
