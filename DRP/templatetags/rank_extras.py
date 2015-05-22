from django import template
from django.utils.safestring import mark_safe

#Used to create custom template tags.
register = template.Library()

@register.filter()
def prettify_reaction_list(value):
  HTML_reactants = "+".join(["""
	<div class="reactantField">
	<span class="dataField type_reactant">{0}</span>
	<span class="dataDecorator">(</span>
	<span class="dataField type_quantity">{1}</span>
	<span class="dataField type.unit">{2}</span>
	<span class="dataDecorator">)</span>
	</div>""".format(value[i], value[i+1], value[i+2]) for i in range(1, len(value[:16]), 3)])

  HTML_extras = """
	<div>
	<span class="dataField type_temp">{0}</span>
	<span class="dataField type_time">{1}</span>
	<span class="dataField type_pH">{2}</span>
	</div>
	""".format(*value[16:])

  HTML_total = HTML_reactants + HTML_extras

  return mark_safe(HTML_total)

