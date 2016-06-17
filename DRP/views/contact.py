"""A module containing only the default contact view."""
from DRP.forms import ContactForm
from django.template import RequestContext
from DRP.Email import EmailToAdmins
from django.shortcuts import render


def contact(request):
    """The contact view.

    Validate a ContactForm against postdata, including CSRF (as part of the template).
    If a user is authenticated, prompts but does not enforce the authenticated user's email
    """
    success = False
    if request.method == "POST":
        form = ContactForm(request.POST)
        if form.is_valid():
            mail = EmailToAdmins('Contact from DRP Contact Us Page', 'Contact sent from {}\n'.format(form.cleaned_data['email']) + form.cleaned_data['content'], includeManagers=True, sender=form.cleaned_data['email'])
            if mail.send() == 1:
                success = True
            else:
                success = False
    else:
        u = request.user
        if u.is_authenticated():
            form = ContactForm(initial={'email': u.email})
        else:
            form = ContactForm()
    return render(request, 'contact.html', RequestContext(request, {'success': success, 'form': form}))
