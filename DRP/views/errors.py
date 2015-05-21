# # # # # # # # # # # # # # # # # # # 
 # # # # # # JSON Views  # # # # # #
# # # # # # # # # # # # # # # # # # # 

#Necessary Imports:
from django.shortcuts import render

######################  Error Messages  ################################
def display_404_error(request):
 response = render(request, '404_error.html')
 response.status_code = 404
 return response

def display_500_error(request):
 response = render(request, '500_error.html')
 response.status_code = 500
 return response
