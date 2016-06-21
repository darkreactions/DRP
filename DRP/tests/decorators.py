"""A module containing decorators which are useful in most test cases for the DRP."""

from DRP.models import Compound, LabGroup, ChemicalClass, License, LicenseAgreement, PerformedReaction, CompoundQuantity, CompoundRole
from DRP.models.rxnDescriptorValues import BoolRxnDescriptorValue, OrdRxnDescriptorValue, NumRxnDescriptorValue, CatRxnDescriptorValue
from DRP.models.rxnDescriptors import BoolRxnDescriptor, OrdRxnDescriptor, CatRxnDescriptor, NumRxnDescriptor
from DRP.models import DataSet, CompoundGuideEntry
from django.contrib.auth.models import User
from django.conf import settings
from datetime import date, timedelta
import os


def createsOrdRxnDescriptor(heading, minimum, maximum, calculatorSoftware='manual', calculatorSoftwareVersion='0'):
    """A class decorator that creates an ordinal reaction descriptor."""
    def _createsOrdRxnDescriptor(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            OrdRxnDescriptor.objects.create(name=heading, heading=heading, maximum=maximum, minimum=minimum,
                                            calculatorSoftware=calculatorSoftware, calculatorSoftwareVersion=calculatorSoftwareVersion)
            _oldSetup(self)

        def tearDown(self):
            OrdRxnDescriptor.objects.filter(name=heading, heading=heading, calculatorSoftware=calculatorSoftware,
                                            calculatorSoftwareVersion=calculatorSoftwareVersion).delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown

        return c

    return _createsOrdRxnDescriptor


def createsOrdRxnDescriptorValue(labGroupTitle, rxnRef, descHeading, value):
    """A class decorator that creates an ordinal reaction descriptor value."""
    def _createsOrdRxnDescriptorValue(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            OrdRxnDescriptorValue.objects.create(descriptor=OrdRxnDescriptor.objects.get(heading=descHeading), reaction=PerformedReaction.objects.get(
                reference=rxnRef, labGroup=LabGroup.objects.get(title=labGroupTitle)), value=value)
            _oldSetup(self)

        c.setUp = setUp

        return c

    return _createsOrdRxnDescriptorValue


def createsCompoundQuantity(rxnRef, compRef, CompRoleAbbrev, mmols):
    """Create a compound quantity for use in a test."""
    def _createsCompoundQuantity(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            CompoundQuantity.objects.create(reaction=PerformedReaction.objects.get(reference=rxnRef), compound=Compound.objects.get(
                abbrev=compRef), role=CompoundRole.objects.get(label=CompRoleAbbrev), amount=mmols)
            _oldSetup(self)

        def tearDown(self):
            if PerformedReaction.objects.filter(reference=rxnRef).exists():
                CompoundQuantity.objects.get(reaction=PerformedReaction.objects.get(reference=rxnRef), compound=Compound.objects.get(
                    abbrev=compRef), role=CompoundRole.objects.get(label=CompRoleAbbrev), amount=mmols).delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown

        return c

    return _createsCompoundQuantity


# def createsRxnDescriptor(heading, descriptorClass, options={}):
#    """A class decorator that creates a reaction descriptor."""
#    def _createsRxnDescriptor(c):
#        _oldSetup = c.setUp
#        _oldTearDown = c.tearDown
#
#        if descriptorType == "BoolRxnDescriptor":
#          descriptor = BoolRxnDescriptor()
#        elif descriptorType == "NumRxnDescriptor":
#          descriptor = NumRxnDescriptor()
#        elif descriptorType == "OrdRxnDescriptor":
#          descriptor = OrdRxnDescriptor()
#        elif descriptorType == "CatRxnDescriptor":
#          descriptor = CatRxnDescriptor()
#        else:
#            error = "Descriptor type \"{}\" unknown to descriptor.".format(descriptorType)
#            raise NotImplementedError(error)
#
#        def setUp(self):
#            descriptor.heading = heading
#            descriptor.name = heading
#
#            for key, val in options.items():
#                setattr(descriptor, key, val)
#
#            descriptor.save()
#
#            _oldSetup(self)
#
#        def tearDown(self):
#            descriptor.delete()
#            _oldTearDown(self)
#
#        c.setUp = setUp
#        c.tearDown = tearDown
#        return c
#    return _createsRxnDescriptor

def createsPerformedReaction(labTitle, username, reference, valid=True):
    """A class decorator that creates a very minimal reaction with no compounds or reactants."""
    def _createsPerformedReaction(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            labGroup = LabGroup.objects.get(title=labTitle)
            user = User.objects.get(username=username)
            reaction = PerformedReaction.objects.create(
                labGroup=labGroup, user=user, reference=reference, valid=valid)
            _oldSetup(self)

        def tearDown(self):
            labGroup = LabGroup.objects.get(title=labTitle)
            user = User.objects.get(username=username)
            PerformedReaction.objects.filter(
                labGroup=labGroup, reference=reference).delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown

        return c

    return _createsPerformedReaction


# TODO: finish replacing this
# def createsPerformedReaction(labTitle, username, reference, compoundAbbrevs=[], compoundRoles=[], compoundAmounts=[], descriptorDict={}, duplicateRef=None):
#    """A class decorator that creates a reaction using pre-existing compounds
#          with pre-existing compoundRoles."""
#    def _createsPerformedReaction(c):
#        _oldSetup = c.setUp
#        _oldTearDown = c.tearDown
#
#        reaction = PerformedReaction()
#        compoundQuantities = []
#        descriptorVals = []
#
#        def setUp(self):
#            labGroup=LabGroup.objects.get(title=labTitle)
#            reaction.labGroup = labGroup
#
#            user=User.objects.get(username=username)
#            reaction.user = user
#            reaction.reference = reference
#            reaction.public = False
#            if PerformedReaction.objects.filter(labGroup=labGroup, reference=duplicateRef).exists():
#                reaction.duplicateOf = PerformedReaction.objects.get(labGroup=labGroup, reference=duplicateRef)
#
#            reaction.save()
#
#            for abbrev, role, quantity in zip(compoundAbbrevs, compoundRoles, compoundAmounts):
#                compound = Compound.objects.get(labGroup=labGroup, abbrev=abbrev)
#                compoundRole = CompoundRole.objects.get(label=role)
#                compoundQuantity = CompoundQuantity(compound=compound, reaction=reaction,
#                                                                                        role=compoundRole, amount=quantity)
#
#                # TODO XXX bulk_create? Can't use the special save
#                compoundQuantity.save()
#
#                compoundQuantities.append(compoundQuantity)
#
#            #TODO: This is hideous and I'm not proud of it.
#            for desc_heading,val in descriptorDict.items():
#                descriptor = None
#                try:
#                    descriptor = BoolRxnDescriptor.objects.get(heading=desc_heading)
#                    descriptorVal = BoolRxnDescriptorValue()
#                except BoolRxnDescriptor.DoesNotExist:
#                  pass
#
#                try:
#                    descriptor = OrdRxnDescriptor.objects.get(heading=desc_heading)
#                    descriptorVal = OrdRxnDescriptorValue()
#                except OrdRxnDescriptor.DoesNotExist:
#                  pass
#
#                try:
#                    descriptor = CatRxnDescriptor.objects.get(heading=desc_heading)
#                    descriptorVal = CatRxnDescriptorValue()
#                except CatRxnDescriptor.DoesNotExist:
#                  pass
#
#                try:
#                    descriptor = NumRxnDescriptor.objects.get(heading=desc_heading)
#                    descriptorVal = NumRxnDescriptorValue()
#                except NumRxnDescriptor.DoesNotExist:
#                  pass
#
#                if descriptor == None:
#                    error = "Unknown descriptorValue type for '{}'".format(descriptor)
#                    raise NotImplementedError(error)
#
#                # TODO XXX: bulk_create?
#                descriptorVal.descriptor = descriptor
#                descriptorVal.value = val
#                descriptorVal.reaction = reaction
#                descriptorVal.save()
#
#                descriptorVals.append(descriptorVal)
#
#            _oldSetup(self)
#
#        def tearDown(self):
#            _oldTearDown(self)
#
#            # TODO XXX bulk deletion?
#            for cq in compoundQuantities:
#                cq.delete()
#
#            for descriptorVal in descriptorVals:
#                descriptorVal.delete()
#
#            DataSet.objects.filter(reactions__in=[reaction]).delete()
#            reaction.delete()
#
#
#        c.setUp = setUp
#        c.tearDown = tearDown
#        return c
#    return _createsPerformedReaction


def createsUser(username, password, is_superuser=False):
    """A class decorator that creates a user."""
    def _createsUser(c):

        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            user = User.objects.create_user(
                username=username, password=password)
            user.is_superuser = is_superuser
            user.save()
            _oldSetup(self)

        def tearDown(self):
            _oldTearDown(self)
            User.objects.filter(username=username).delete()

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _createsUser


def createsCompoundRole(label, description):
    """Create a compound role in a test."""
    def _createsCompoundRole(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        role = CompoundRole()

        def setUp(self):
            role.label = label
            role.description = description
            role.save()
            _oldSetup(self)

        def tearDown(self):
            _oldTearDown(self)
            role.delete()

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _createsCompoundRole


def createsCompound(abbrev, csid, classLabel, labTitle, custom=False):
    """Create a compound in a test."""
    def _createsCompound(c):

        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        compound = Compound(CSID=csid, custom=custom)

        def setUp(self):
            compound.save()
            CompoundGuideEntry.objects.create(labGroup = LabGroup.objects.get(title=labTitle), abbrev=abbrev, compound=compound)
            for c in ChemicalClass.objects.filter(label=classLabel):
                compound.chemicalClasses.add(c)
            compound.save()
            _oldSetup(self)

        def tearDown(self):
            _oldTearDown(self)
            compound.delete()

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _createsCompound


def createsChemicalClass(label, description):
    """A class decorator that creates a test chemical class for the addition of compounds into the database."""
    def _createsChemicalClass(c):

        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        chemicalClass = ChemicalClass(label=label, description=description)

        def setUp(self):
            chemicalClass.save()
            _oldSetup(self)

        def tearDown(self):
            chemicalClass.delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _createsChemicalClass


def joinsLabGroup(username, labGroupTitle):
    """
    A class decorator that creates a test lab group.

    labGroupTitle is it's title and assigns user identified by username to that lab group.
    """
    def _joinsLabGroup(c):
        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        labGroup = LabGroup(title=labGroupTitle, address='War drobe',
                            email='Aslan@example.com', access_code='new_magic')

        def setUp(self):
            user = User.objects.get(username=username)
            labGroup.save()
            user.labgroup_set.add(labGroup)
            _oldSetup(self)

        def tearDown(self):
            _oldTearDown(self)
            labGroup.delete()

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _joinsLabGroup


def signsExampleLicense(username):
    """A class decorator that creates a test license and makes the user specified by username sign it on setUp."""
    def _signsExampleLicense(c):

        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        license = License(text='This is an example license used in a test',
                          effectiveDate=date.today() - timedelta(1))

        def setUp(self):
            user = User.objects.get(username=username)
            license.save()
            agreement = LicenseAgreement(user=user, text=license)
            agreement.save()
            _oldSetup(self)

        def tearDown(self):
            LicenseAgreement.objects.filter(
                user__username=username, text=license).delete()
            if license.licenseagreement_set.count() < 1:
                license.delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _signsExampleLicense


def loadsCompoundsFromCsv(labGroupTitle, csvFileName):
    """A class decorator that creates a test set of compounds using the csvFileName, which should be stored in the tests directory resource folder."""
    def _loadsCompoundsFromCsv(c):

        _oldSetup = c.setUp
        _oldTearDown = c.tearDown

        def setUp(self):
            labGroup = LabGroup.objects.get(title=labGroupTitle)
            compounds = labGroup.compound_set.fromCsv(os.path.join(
                settings.APP_DIR, 'tests', 'resource', csvFileName))

            # TODO XXX bulk_create? Can't use our custom save then
            for compound in compounds:
                compound.csConsistencyCheck()
                compound.save()

            _oldSetup(self)

        def tearDown(self):
            CompoundQuantity.objects.all().delete()
            _oldTearDown(self)

        c.setUp = setUp
        c.tearDown = tearDown
        return c
    return _loadsCompoundsFromCsv


def createsPerformedReactionSetOrd(c):
    """a broken decorator for creating a bunch of sample reactions."""
    # Create a bunch of simple sample reactions with ordinal outcomes.
    c = createsPerformedReaction("Watchmen", "Rorschach", "R01", ["EtOH"], ["Org"], [
                                 0.13], {"outcome": 1, "testNumber": 5.04}, duplicateRef='R15')(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R02", ["EtOH"], [
                                 "Org"], [0.71], {"outcome": 1, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R03", ["dmed"], [
                                 "Org"], [0.1], {"outcome": 1, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R04", ["dmed"], [
                                 "Org"], [0.18], {"outcome": 1, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R05", ["dabco"], [
                                 "Org"], [0.71], {"outcome": 1, "testNumber": 5.04})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R06", ["dabco"], [
                                 "Org"], [0.14], {"outcome": 1, "testNumber": 5.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R07", ["Water"], [
                                 "Water"], [0.14], {"outcome": 1, "testNumber": 5.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R08", ["EtOH", "dmed", "Water"], [
                                 "Org", "Org", "Water"], [0.31, 0.3, 0.5], {"outcome": 4, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R09", ["dmed", "EtOH", "Water"], [
                                 "Org", "Org", "Water"], [0.2, 0.34, 0.3], {"outcome": 4, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R10", ["dabco", "Water"], [
                                 "Org", "Water"], [0.3, 0.32], {"outcome": 3, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R11", ["Water", "dabco"], [
                                 "Water", "Org"], [0.3, 0.34], {"outcome": 3, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R12", ["dmed", "Water"], [
                                 "Org", "Water"], [0.3, 0.34], {"outcome": 3, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R13", ["EtOH", "Water"], [
                                 "Org", "Water"], [0.3, 0.35], {"outcome": 3, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R14", ["dmed", "Water"], [
                                 "Org", "Water"], [0.4, 0.56], {"outcome": 4, "testNumber": 0.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R15", ["Water", "EtOH"], [
                                 "Water", "Org"], [0.4, 0.57], {"outcome": 4, "testNumber": 0.02})(c)

    c = createsRxnDescriptor("outcome", "OrdRxnDescriptor", options={
                             "maximum": 4, "minimum": 1})(c)
    c = createsRxnDescriptor("testNumber", "NumRxnDescriptor")(c)

    c = createsCompoundRole('Org', 'Organic')(c)
    c = createsCompoundRole('Water', 'Water')(c)
    c = createsCompound('EtOH', 682, 'Org', 'Watchmen')(c)
    c = createsCompound('dmed', 67600, 'Org', 'Watchmen')(c)
    c = createsCompound('dabco', 8882, 'Org', 'Watchmen')(c)
    c = createsCompound('Water', 937, 'Water', 'Watchmen')(c)
    c = createsChemicalClass('Org', 'Organic')(c)
    c = createsChemicalClass('Water', 'Water')(c)

    c = signsExampleLicense("Rorschach")(c)
    c = joinsLabGroup('Rorschach', 'Watchmen')(c)
    c = createsUser('Rorschach', 'whatareyouwaitingfor')(c)
    return c


def createsPerformedReactionSetBool(c):
    """Broken."""
    # Create a bunch of simple sample reactions with Boolean outcomes.
    c = createsPerformedReaction("Watchmen", "Rorschach", "R01", ["EtOH"], ["Org"], [
                                 0.13], {"outcome": False, "testNumber": 5.04}, duplicateRef='R15')(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R02", ["EtOH"], ["Org"], [
                                 0.71], {"outcome": False, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R03", ["dmed"], [
                                 "Org"], [0.1], {"outcome": False, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R04", ["dmed"], ["Org"], [
                                 0.18], {"outcome": False, "testNumber": 5.0})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R05", ["dabco"], [
                                 "Org"], [0.71], {"outcome": False, "testNumber": 5.04})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R06", ["dabco"], [
                                 "Org"], [0.14], {"outcome": False, "testNumber": 5.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R07", ["Water"], [
                                 "Water"], [0.14], {"outcome": False, "testNumber": 5.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R08", ["EtOH", "dmed", "Water"], [
                                 "Org", "Org", "Water"], [0.31, 0.3, 0.5], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R09", ["dmed", "EtOH", "Water"], [
                                 "Org", "Org", "Water"], [0.2, 0.34, 0.3], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R10", ["dabco", "Water"], [
                                 "Org", "Water"], [0.3, 0.32], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R11", ["Water", "dabco"], [
                                 "Water", "Org"], [0.3, 0.34], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R12", ["dmed", "Water"], [
                                 "Org", "Water"], [0.3, 0.34], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R13", ["EtOH", "Water"], [
                                 "Org", "Water"], [0.3, 0.35], {"outcome": True, "testNumber": 0.1})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R14", ["dmed", "Water"], [
                                 "Org", "Water"], [0.4, 0.56], {"outcome": True, "testNumber": 0.01})(c)
    c = createsPerformedReaction("Watchmen", "Rorschach", "R15", ["Water", "EtOH"], [
                                 "Water", "Org"], [0.4, 0.57], {"outcome": True, "testNumber": 0.02})(c)

    c = createsRxnDescriptor("outcome", "BoolRxnDescriptor", options={})(c)
    c = createsRxnDescriptor("testNumber", "NumRxnDescriptor")(c)

    c = createsCompoundRole('Org', 'Organic')(c)
    c = createsCompoundRole('Water', 'Water')(c)
    c = createsCompound('EtOH', 682, 'Org', 'Watchmen')(c)
    c = createsCompound('dmed', 67600, 'Org', 'Watchmen')(c)
    c = createsCompound('dabco', 8882, 'Org', 'Watchmen')(c)
    c = createsCompound('Water', 937, 'Water', 'Watchmen')(c)
    c = createsChemicalClass('Org', 'Organic')(c)
    c = createsChemicalClass('Water', 'Water')(c)

    c = signsExampleLicense("Rorschach")(c)
    c = joinsLabGroup('Rorschach', 'Watchmen')(c)
    c = createsUser('Rorschach', 'whatareyouwaitingfor')(c)
    return c
