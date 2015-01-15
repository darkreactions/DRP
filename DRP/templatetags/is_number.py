from django import template
register = template.Library()

@register.filter()
def is_number(value):
  try:
    int(value)
    return True
  except:
    return False
