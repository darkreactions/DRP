'''The django forms and model forms used in the DRP
Classes:

LabGroupForm: for creating Lab Groups in the django admin.
ContactForm: A very simple form for the contact page.
'''
from LabGroup import LabGroupForm, LabGroupJoiningForm
from Contact import ContactForm
from authentication import ConfirmationForm, LicenseAgreementForm, UserCreationForm
