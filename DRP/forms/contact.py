"""A module containing forms pertinent to contacting Managers."""
import django.forms as forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit


class ContactForm(forms.Form):
    """A very simple form for the contacting of site Admins by all people viewing the DRP site."""

    email = forms.EmailField(label="Your Email Address",
                             initial="youremail@example.com")
    content = forms.CharField(label="Your Message", widget=forms.Textarea)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_method = 'post'
        self.helper.add_input(Submit('submit', 'Send'))
