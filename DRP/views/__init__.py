"""Module containing the views for the DRP project."""

from .contactView import contact
from .registerViews import register, confirm
from .licenseView import license
from .joinGroupView import joinGroup
from .selectGroupView import selectGroup
import .compound
import .decorators
import .reaction
from .leaveGroupView import leaveGroup
from api import api1
