from django.shortcuts import render
from django.contrib.auth.decorators import login_required

def page(request, template):
  return render(request, 'global_page.html', {"template":template})


@login_required
def locked_page(request, template):
  return page(request, template)
