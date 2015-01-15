from django.shortcuts import render

def page(request, template):
  return render(request, 'global_page.html', {"template":template})
