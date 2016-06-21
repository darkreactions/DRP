"""A module containing code pertinent to manipulating Compound objects."""

from django.db import models, transaction
from descriptors import BooleanDescriptor, NumericDescriptor
from descriptors import CategoricalDescriptor, OrdinalDescriptor
from molDescriptorValues import BoolMolDescriptorValue, NumMolDescriptorValue
from molDescriptorValues import CatMolDescriptorValue, OrdMolDescriptorValue
from ChemicalClass import ChemicalClass
from LabGroup import LabGroup
import csv
from querysets import CsvQuerySet, ArffQuerySet
from chemspipy import ChemSpider
from django.conf import settings
from django.core.exceptions import ValidationError
from itertools import chain
import importlib
from collections import OrderedDict
import PerformedReaction
from django.core.validators import RegexValidator
from decimal import Decimal

descriptorPlugins = [importlib.import_module(plugin) for
                     plugin in settings.MOL_DESCRIPTOR_PLUGINS]
# This prevents a cyclic dependency problem

def elementsFormatValidator(molFormula):
    """A validator for molecular formulae."""
    elements = {}
    inBrackets = False
    currentElement = ''
    for char in molFormula:
        if inBrackets:
            if currentElement not in elements:
                elements[currentElement] = {'stoichiometry': 0}
            if char in (str(x) for x in range(0, 10)) or char == '.':
                strStoichiometry += char
            elif char == '}':
                inBrackets = False
                elements[currentElement][
                    'stoichiometry'] += float(strStoichiometry)
                currentElement = ''
                strStoichiometry = ''
            else:
                raise ValidationError(
                    'Invalid molecular formula format.', 'mol_malform')
        elif char.isalpha():
            if char.isupper():
                if currentElement != '':
                    if currentElement in elements:
                        elements[currentElement]['stoichiometry'] += 1
                    else:
                        elements[currentElement] = {'stoichiometry': 1}
                currentElement = char
            else:
                currentElement += char
        elif char == '{':
            strStoichiometry = ''
            inBrackets = True
        elif char == '_':
            pass
        else:
            raise ValidationError(
                'Invalid molecular formula format.', code='mol_malform')


class CompoundQuerySet(CsvQuerySet, ArffQuerySet):

    """A specialised queryset for outputting Compounds in specific formats."""

    def __init__(self, **kwargs):
        """Initialise teh queryset."""
        kwargs.pop('model', None)
        super(CompoundQuerySet, self).__init__(Compound, **kwargs)

    def maxChemicalClassCount(self):
        """Give a count of the maximum number of chemical classes.

        Calculated over all molecules in the given queryset.
        """
        annotated = self.annotate(
            chemicalClassCount=models.Count('chemicalClasses')
        )
        m = annotated.aggregate(max=models.Max('chemicalClassCount'))['max']
        return 0 if m is None else m

    def csvHeaders(self, whitelist=None):
        """Generate the header row information for the CSV."""
        headers = super(CompoundQuerySet, self).csvHeaders()
        m = Compound.objects.all().maxChemicalClassCount()
        headers += ['chemicalClass_{}'.format(x + 1) for x in range(0, m)]
        return headers

    def arffHeaders(self, whitelist=None):
        """Generate headers for the arff file."""
        headers = super(CompoundQuerySet, self).arffHeaders(whitelist)
        m = Compound.objects.all().maxChemicalClassCount()
        for x in range(0, m):
            label = 'chemicalClass_{0}'.format(x + 1)
            clsStrings = ('"{}"'.format(chemicalClass)
                          for chemicalClass in ChemicalClass.objects.all())
            headers[label] = '@attribute {} {{{}}}'.format(
                label, ','.join(clsStrings))
        return headers

    def expandedArffHeaders(self, whitelist=None):
        """Generate expanded headers for the arff file."""
        headers = self.arffHeaders(whitelist)
        headers.update(OrderedDict(((d.csvHeader, d.arffHeader)
                                    for d in self.descriptors)))
        return headers

    def expandedCsvHeaders(self, whitelist=None):
        """Generate the expanded header for the csv."""
        return self.csvHeaders() + [d.csvHeader for d in self.descriptors]

    @property
    def descriptors(self):
        """Return the descriptor which have relationship to the queryset."""
        return chain(
            BooleanDescriptor.objects.filter(
                boolmoldescriptorvalue__in=BoolMolDescriptorValue.objects.filter(
                    compound__in=self
                )
            ).distinct(),
            NumericDescriptor.objects.filter(
                nummoldescriptorvalue__in=NumMolDescriptorValue.objects.filter(
                    compound__in=self
                )
            ).distinct(),
            OrdinalDescriptor.objects.filter(
                ordmoldescriptorvalue__in=OrdMolDescriptorValue.objects.filter(
                    compound__in=self
                )
            ).distinct(),
            CategoricalDescriptor.objects.filter(
                catmoldescriptorvalue__in=CatMolDescriptorValue.objects.filter(
                    compound__in=self
                )
            ).distinct()
        )

    def rows(self, expanded, whitelist=None):
        """Generate 'row' (list) for each row of the file."""
        if expanded:
            compounds = self.prefetch_related(
                'boolmoldescriptorvalue_set__descriptor')
            compounds = compounds.prefetch_related(
                'catmoldescriptorvalue_set__descriptor')
            compounds = compounds.prefetch_related(
                'ordmoldescriptorvalue_set__descriptor')
            compounds = compounds.prefetch_related(
                'nummoldescriptorvalue_set__descriptor')
            compounds = compounds.prefetch_related('chemicalClasses')
            for item in compounds:
                row = {field.name: getattr(item, field.name)
                       for field in self.model._meta.fields}
                row.update(
                    {dv.descriptor.csvHeader: dv.value for dv in item.descriptorValues})
                i = 1
                for cc in item.chemicalClasses.all():
                    row['chemicalClass_{}'.format(i)] = cc.label
                    i += 1
                yield row
        else:
            for row in super(CompoundQuerySet, self).rows(expanded):
                yield row

    def calculate_descriptors(self, verbose=False, plugins=None, **kwargs):
        """Calculate descriptors for the current molecule queryset."""
        for plugin in descriptorPlugins:
            if plugins is None or plugin.__name__ in plugins:
                if verbose:
                    print "Calculating for plugin: {}".format(plugin)
                plugin.calculate_many(self, verbose=verbose, **kwargs)
                if verbose:
                    print "Done with plugin: {}\n".format(plugin)


class CompoundManager(models.Manager):

    """A custom manager for the Compound Class which permits the creation of entries to and from CSVs."""

    # NOTE:This doesn't actually work, but no-one's sure which way django is
    # going to jump on this.
    use_for_related_fields = True

    def get_queryset(self):
        """Return the default queryset."""
        return CompoundQuerySet()

    def fromCsv(self, fileName, labGroup=None):
        """Read a CSV into the creating objects, returning a list of compounds which have not yet been saved.

        This assumes that the uploaded csv will have headers which map to the names of the fields and that compound classes are
        stored as comma separated lists of the chemicalClass LABEL only.

        Each compound will perform a chemspider-based consistency check on the information it has been created with to ensure
        information is consistent- this throws an ValidationError if it is not.
        """
        if labGroup is None and hasattr(self, 'instance'):
            # we presume that if this is being called without a labgroup that's
            # because this manager belongs to a lab group
            labGroup = self.instance

        compoundsList = []
        cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
        with open(fileName) as f:
            reader = csv.DictReader(f, restkey='restKey')
            rowCount = 0
            errors = []
            for row in reader:
                try:
                    rowCount += 1
                    if 'chemicalClasses' in row:
                        classes = (c.strip()
                                   for c in row['chemicalClasses'].split(','))
                        chemicalClasses = []
                        for c in classes:
                            chemicalClass, created = ChemicalClass.objects.get_or_create(
                                label=c)
                            chemicalClasses.append(chemicalClass)
                    if row.get('CAS') not in ('', None) and row.get('CSID') in ('', None):
                        CASResults = cs.simple_search(row['CAS'])
                        if len(CASResults) < 1:
                            errors.append(ValidationError(
                                'CAS Number returned no results from ChemSpider on row %(rowCount)d of uploaded csv.', params={'rowCount': rowCount}))
                        elif len(CASResults) == 1:
                            # a little hacky, but it gets the job done
                            row['CSID'] = CASResults[0].csid
                        else:
                            errors.append(ValidationError(
                                'CAS number returns more than one ChemSpider ID on row %(rowCount)d of uploaded csv.', params={'rowCount': rowCount}))
                    elif row.get('CSID') in ('', None):
                        errors.append(ValidationError(
                            'No CSID provided on row %(rowCount)d of uploaded csv.', params={'rowCount': rowCount}))
                    kwargs = {}
                    kwargs['CSID'] = row.get('CSID')
                    kwargs['abbrev'] = row.get('abbrev')
                    kwargs['smiles'] = row.get('smiles')
                    kwargs['name'] = row.get('name')
                    kwargs['INCHI'] = row.get('INCHI')
                    compound = Compound(labGroup=labGroup, **kwargs)
                    for chemicalClass in chemicalClasses:
                        compound.lazyChemicalClasses.append(chemicalClass)
                    compoundsList.append(compound)
                except ValidationError as e:
                    for message in e.messages:
                        errors.append(ValidationError(
                            message + ' on row %(rowCount)d of uploaded csv', params={'rowCount': rowCount}))
            if len(errors) > 0:
                raise ValidationError(errors)
        return compoundsList


class Compound(models.Model):

    """
    A class for containing data about Compounds used in chemical reactions.

    The assumption is made that all chemicals used are single-species.
    """

    class Meta:
        app_label = "DRP"

    name = models.CharField('Name', max_length=400)
    """Normally the IUPAC name of the compound, however this may not be the most parsable name (which is preferable)."""
    chemicalClasses = models.ManyToManyField(
        ChemicalClass, verbose_name="Chemical Class")
    """The class of the compound- examples include Inorganic Salt."""
    CSID = models.PositiveIntegerField('Chemspider ID', null=True, unique=True)
    """The chemspider ID for the compound- preferable to the CAS_ID since it is not subject to licensing restrictions."""
    custom = models.BooleanField("Custom", default=False)
    """This flag denotes whether a compound has been added irrespective of other validation.
    This should be restricted to superusers."""
    INCHI = models.TextField('InCHI key', blank=True, default='')
    """The Inchi key for a compound- a canonical representation of a molecule which is also unique."""

    smiles = models.TextField('Smiles', blank=True, default='')
    """A non-canonical string representation of a molecule which cannot be directly used to test for identity
    but is nevertheless useful for calculating descriptors
    """

    labGroups = models.ManyToManyField(LabGroup, verbose_name="Lab Groups", through="DRP.CompoundGuideEntry")
    """Tells us whose compound guide this appears in."""

    formula = models.CharField(
        max_length=500,
        blank=True,
        help_text="A formula should be made up of element names. C_{4}H_{8} type notation should be use for subscript digits.",
        validators=[elementsFormatValidator]
    )

    objects = CompoundManager()

    def __init__(self, *args, **kwargs):
        """Instantiate an object."""
        super(Compound, self).__init__(*args, **kwargs)
        self.lazyChemicalClasses = []

    def __unicode__(self):
        """Unicode representation of a compound is it's name and abbreviation."""
        return unicode("{}".format(self.name), 'utf-8')

    def csConsistencyCheck(self):
        """Perform a consistency check of this record against chemspider. Raise a ValidationError on error."""
        if not self.custom:
            errorList = []
            cs = ChemSpider(settings.CHEMSPIDER_TOKEN)
            if self.CSID is None or self.CSID is '':
                raise ValidationError('No CSID set', 'no_csid')
            else:
                csCompound = cs.get_compound(self.CSID)
                if self.name not in ('', None):
                    nameResults = cs.simple_search(self.name)
                    if csCompound not in nameResults:
                        errorList.append(ValidationError(
                            'A compound was consistency checked and was found to have an invalid name', code='invalid_inchi'))
                else:
                    self.name = csCompound.common_name
                if self.INCHI == '':
                    self.INCHI = csCompound.stdinchi
                elif self.INCHI != csCompound.stdinchi:
                    errorList.append(ValidationError(
                        'A compound was consistency checked and was found to have an invalid InChi', code='invalid_inchi'))
                if self.smiles == '':
                    self.smiles = csCompound.smiles
                elif self.smiles != csCompound.smiles:
                    errorList.append(ValidationError(
                        'A compound was consistency checked and was found to have an invalid smiles string', code='invalid_smiles'))
                if self.formula == '':
                    self.formula = csCompound.molecular_formula
                elif self.formula != csCompound.molecular_formula:
                    errorsList.append(ValidationError(
                        'A compound was consistency checked and was found to have an invalid formula', code="invalid_formula"))
                if len(errorList) > 0:
                    raise ValidationError(errorList)

    @transaction.atomic
    def save(self, calcDescriptors=True, invalidateReactions=True, *args, **kwargs):
        """Save the compound, invalidating any consequent objects like reactions and models."""
        if self.pk is not None and invalidateReactions:
            for reaction in self.reaction_set.all():
                reaction.save()  # descriptor recalculation
                try:
                    reaction.performedreaction.save()  # invalidate models
                except PerformedReaction.PerformedReaction.DoesNotExist:
                    pass  # it doesn't matter
        super(Compound, self).save(*args, **kwargs)
        if self.pk is not None:
            # coping mechanism for compounds loaded from csv files; not to be
            # used by other means
            for lcc in self.lazyChemicalClasses:
                self.chemicalClasses.add(lcc)
            if calcDescriptors:  # not generally done, but useful for debugging
                for descriptorPlugin in descriptorPlugins:
                    descriptorPlugin.calculate(self)

    @property
    def descriptorValues(self):
        """Return an iterable of all descriptor values for this compound."""
        return chain(self.boolmoldescriptorvalue_set.all(), self.nummoldescriptorvalue_set.all(), self.ordmoldescriptorvalue_set.all(), self.catmoldescriptorvalue_set.all())

    @property
    def elements(self):
        """
        Return a dictionary of elemental symbols and their stoichiometry.

        Note that this method does not validate the data contained in the database.
        """
        elements = {}
        inBrackets = False
        currentElement = ''
        strStoichiometry = ''
        for char in self.formula:
            if inBrackets:
                if currentElement not in elements:
                    elements[currentElement] = {'stoichiometry': 0}
                if char in (str(x) for x in range(0, 10)) or char == '.':
                    strStoichiometry += char
                elif char == '}':
                    inBrackets = False
                    elements[currentElement][
                        'stoichiometry'] += Decimal(strStoichiometry)
                    currentElement = ''
                    strStoichiometry = ''
                else:
                    raise ElementException('Invalid molecular formula format.')
            elif char.isalpha():
                if char.isupper():
                    if currentElement != '':
                        if currentElement in elements:
                            elements[currentElement]['stoichiometry'] += 1
                        else:
                            elements[currentElement] = {'stoichiometry': 1}
                    currentElement = char
                else:
                    currentElement += char
            elif char == '{':
                inBrackets = True
        if currentElement != '':
            elements[currentElement] = {'stoichiometry': 1}
        return elements


class CompoundGuideEntry(models.Model):

    class Meta:
        app_label='DRP'
        unique_together=(('compound', 'labGroup'), ('abbrev', 'labGroup'))

    compound = models.ForeignKey(Compound)
    labGroup = models.ForeignKey(LabGroup)
    abbrev = models.CharField("Abbreviation", max_length=100)
    


class ElementsException(Exception):

    """An exception for element formats."""

    pass
