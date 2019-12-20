"""Forms pertinent to user creation and authentication."""
from django.core.exceptions import ValidationError
from django.contrib.auth.models import User
import django.forms as forms
from DRP.models import LicenseAgreement
from django.contrib.auth.forms import UserCreationForm as DjangoUserCreationForm
from django.contrib.auth.forms import AuthenticationForm as DjangoAuthenticationForm
from django.utils.safestring import mark_safe
from django.utils.html import conditional_escape
from django.contrib.auth import authenticate
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit


class UserCreationForm(DjangoUserCreationForm):
    """A form for creating users."""

    email = forms.EmailField(required=True)

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_method = 'post'
        self.helper.add_input(Submit('submit', 'Register'))

    


class ConfirmationForm(DjangoAuthenticationForm):
    """A form for confirming a user's credentials, without checking if they are 'active'."""

    def clean(self):
        """
        A very close rewrite of the DjangoAuthenticationForm method.

        Raise an error if the user is already active rather than if it is inactive.
        """
        username = self.cleaned_data.get('username')
        password = self.cleaned_data.get('password')

        if username and password:
            self.user_cache = authenticate(
                username=username, password=password)
            if self.user_cache is None:
                raise forms.ValidationError(self.error_messages['invalid_login'],
                                            code='invalid_login',
                                            params={'username': self.username_field.verbose_name})
            elif self.user_cache.is_active:
                raise forms.ValidationError('Your account has already been activated',
                                            code='active_user')
        return self.cleaned_data


class LicenseAgreementForm(DjangoAuthenticationForm):
    """A re-authentication form for the signing of site license agreements for DRP deployments."""

    licenseId = forms.IntegerField(widget=forms.widgets.HiddenInput)

    def __init__(self, user, license, *args, **kwargs):
        """Initilaiser."""
        super(LicenseAgreementForm, self).__init__(*args, **kwargs)
        self.user = user
        self.license = license
        self.fields['licenseId'].initial = license.id

    def clean(self):
        """A slightly adjusted clean method which checks that the correct license is being signed and checks that the right user is signing."""
        supercleaned = super(LicenseAgreementForm, self).clean()
        if self.user != self.user_cache:
            raise forms.ValidationError(
                'Incorrect user details entered. Please enter your own user credentials')
        if self.license.id != self.cleaned_data.get('licenseId'):
            raise forms.ValidationError(
                'Whilst you were signing the agreement, a more up-to-date agreement has been created. Please read the new agreement and sign again.')
        return supercleaned.update(self.cleaned_data)

    def as_ol(self):
        """Present this form as an ordered list."""
        text = mark_safe(
            '<pre>{0}</pre>'.format(conditional_escape(self.license.text)))
        text += mark_safe('<ol>{0}</ol>'.format(
            super(LicenseAgreementForm, self).as_ul()))
        return text

    def save(self, commit=True):
        """Save the agreement to the license."""
        agreement = LicenseAgreement(user=self.user, text=self.license)
        if commit:
            agreement.save()
        return agreement
