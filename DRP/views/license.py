'''contains the view which controls license signing'''
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from DRP.models import License, LicenseAgreement
from DRP.forms import LicenseAgreementForm 
from django.template.loader import get_template
from django.template import RequestContext
from django.http import HttpResponseNotFound
from django.shortcuts import render, redirect
import datetime

@login_required
def license(request):
  '''Controls requests pertaining to the signing of deployment license agreements'''
  if License.objects.all().count() < 1:
    template = get_template('license_404.html')
    return HttpResponseNotFound(template.render(RequestContext(request)))
  else:
    latestLicense = License.objects.latest()
    currentSignedLicenseQ = LicenseAgreement.objects.filter(user=request.user, text=latestLicense)
    nextPage = request.GET['next'] if 'next' in request.GET.keys() else None
    if currentSignedLicenseQ.count()==0:
      if request.method == 'POST':
        form = LicenseAgreementForm(request.user, latestLicense, data=request.POST)
        if form.is_valid():
          form.save()
          if nextPage:
            return redirect(nextPage) 
          else:
            return render(request, 'license_signed.html')
        else:
          return render(request, 'license.html', RequestContext(request, {'form':form, 'next':nextPage}))
      else:
        form = LicenseAgreementForm(request.user, latestLicense)
        return render(request, 'license.html', RequestContext(request, {'form':form, 'next':nextPage}))
    elif currentSignedLicenseQ.count()==1:
      if nextPage:
        return redirect(nextPage)
      else:
        return render(request, 'license_up_to_date.html', RequestContext(request, {'license':latestLicense}))
    else:
      raise RuntimeError('Impossible condition occured. Please contact an administrator or developer')
