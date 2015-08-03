'''contains the view which controls license signing'''
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from DRP.models import License, LicenseAgreement
from DRP.forms import LicenseAgreementForm 
from django.template.loader import get_template
from django.template import RequestContext
from django.http import HttpResponseNotFound
from django.shortcuts import render
import datetime

@login_required
def license(request):
  '''Controls requests pertaining to the signing of deployment license agreements'''
  if licenses.all().count() < 1:
    template = get_template('license_404.html')
    return HttpResponseNotFound(template.render(RequestContext(request)))
  else:
    latestLicense = License.objects.latest()
    currentSignedLicenseQ = LicenseAgreement.objects.filter(user=request.user, text=latestLicense)
    if currentSignedLicenseQ.count()==0:
      if request.method == 'POST':
        form = LicenseAgreementForm(request.user, latestLicense, request.POST)
        if form.is_valid():
          form.save()
          if 'next' in request.GET.keys():
            return redirect(request.GET['next']) 
          else:
            return render('license_signed.html')
        else:
          return render('license.html', RequestContext(request, {'form':form}))
      else:
        form = LicenseAgreementForm()
        return render('license.html', RequestContext(request, {'form':form}))
    elif currentSignedLicenseQ.count()==1:
      if 'next' in request.GET.keys():
        return redirect(request.GET['next'])
      else:
        return render('license_up_to_date.html', RequestContext(request, {'license':latestLicense}))
    else:
      raise RuntimeError('Impossible condition occured. Please contact an administrator or developer')
