from django.db import models
from Lab_Group import Lab_Group
from django.contrib.auth.models import User
import datetime

class Lab_Member(models.Model):
  class Meta:
    app_label = "DRP"

  user = models.OneToOneField(User, related_name="profile", unique=True)
  license_agreement_date_dt = models.DateTimeField("License Agreement Date",
                                                    null=True, blank=True)
  lab_group = models.ForeignKey(Lab_Group)

  def update_license(self):
    self.license_agreement_date_dt = datetime.datetime.now()
    self.save()

  def __unicode__(self):
    return self.user.username


def get_users(lab):
  return User.objects.filter(profile__lab_group=lab).order_by("first_name")
