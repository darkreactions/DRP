'''A module containing forms pertinent to contacting Managers'''
import django.forms as forms


class ContactForm(forms.Form):
    '''A very simple form for the contacting of site Admins by all people viewing the DRP site'''

    email = forms.EmailField(label="Your Email Address", initial="youremail@example.com")
    content = forms.CharField(label="Your Message", widget=forms.Textarea)
