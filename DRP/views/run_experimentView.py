from DRP.forms import ContactForm
from django.template import RequestContext
from DRP.email import EmailToAdmins
from django.shortcuts import render


def run_experiment(request):
	return render(request, 'run.html', RequestContext(request, {'success': success, 'form': form}), status=status)