"Views which correspond to the data extraction api"

from django.contrib.auth.decorators import login_required
from django.http import HttpResponseForbidden, HttpResponse, Http404
from django.core.exceptions import PermissionDenied
from django.core import serializers
from django.contrib.auth.models import User
from DRP.models import License, LicenseAgreement
from DRP.models import Compound, CompoundRole, CompoundQuantity
from DRP.models import ChemicalClass, Reaction
from DRP.models import LabGroup, PerformedReaction


@login_required
def api1(request, component):
    "Version 1 of the data export api"
    if request.user.is_superuser:
        if 'limit' in request.GET:
            limit = int(request.GET['limit'])
        else:
            limit = PerformedReaction.objects.all().count()
        performedReactions = PerformedReaction.objects.all()[:limit]
        if component == 'performed_reactions':
            return HttpResponse(serializers.serialize('xml', performedReactions), content_type='Application/xml')
        else:
            if component == 'reactions':
                return HttpResponse(serializers.serialize('xml', Reaction.objects.filter(performedreaction__in=(p for p in performedReactions))), content_type='Application/xml')
            else:
                labGroups = LabGroup.objects.filter(reaction__in=[p for p in performedReactions]).distinct()
                if component == 'lab_groups':
                    return HttpResponse(serializers.serialize('xml', labGroups), content_type='Application/xml')
                else:
                    users = (User.objects.filter(labgroup__in=[l for l in labGroups]) |
                             User.objects.filter(performedReactions__in=(p for p in performedReactions)) |
                             User.objects.filter(performedreaction__in=(p for p in performedReactions)) |
                             User.objects.filter(pk=request.user.pk))
                    users = users.distinct()
                    if component == 'users':
                        return HttpResponse(serializers.serialize('xml', users, fields=('username', 'first_name', 'last_name', 'email', 'is_staff', 'is_active', 'is_superuser')), content_type='Application/xml')
                    else:
                        licenseAgreements = LicenseAgreement.objects.filter(user__in=users).distinct()
                        if component == 'license_agreements':
                            return HttpResponse(serializers.serialize('xml', licenseAgreements), content_type='Application/xml')
                        else:
                            licenses = License.objects.filter(licenseagreement__in=licenseAgreements).distinct()
                            if component == 'licenses':
                                return HttpResponse(serializers.serialize('xml', licenses), content_type='Application/xml')
                            else:
                                compoundQuantities = CompoundQuantity.objects.filter(reaction__in=[p for p in performedReactions]).distinct()
                                if component == 'compound_quantities':
                                    return HttpResponse(serializers.serialize('xml', compoundQuantities), content_type='Application/xml')
                                else:
                                    compoundRoles = CompoundRole.objects.filter(compoundquantity__in=compoundQuantities).distinct()
                                    if component == 'compound_roles':
                                        return HttpResponse(serializers.serialize('xml', compoundRoles), content_type='Application/xml')
                                    else:
                                        compounds = Compound.objects.filter(compoundquantity__in=compoundQuantities).distinct()
                                        if component == 'compounds':
                                            return HttpResponse(serializers.serialize('xml', compounds), content_type='Application/xml')
                                        else:
                                            chemicalClasses = ChemicalClass.objects.filter(compound__in=compounds).distinct()
                                            if component == 'chemical_classes':
                                                return HttpResponse(serializers.serialize('xml', chemicalClasses), content_type='Application/xml')
                                            else:
                                                raise Http404
    else:
        raise PermissionDenied
