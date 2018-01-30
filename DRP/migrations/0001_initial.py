# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import DRP.models.rxnDescriptorValues
import django.core.validators
import DRP.models.molDescriptorValues
import DRP.models.validators
from django.conf import settings
import DRP.models.compound
import DRP.models.performedReaction
import django.db.models.deletion
import DRP.models.fileStorage


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='BoolMolDescriptorValue',
            fields=[
                ('value', models.NullBooleanField(verbose_name='Value for descriptor')),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.molDescriptorValues.molUid, max_length=36)),
            ],
            options={
                'verbose_name': 'Boolean Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='BoolRxnDescriptorValue',
            fields=[
                ('value', models.NullBooleanField(verbose_name='Value for descriptor')),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36)),
            ],
            options={
                'verbose_name': 'Boolean Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='CategoricalDescriptorPermittedValue',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('value', models.CharField(verbose_name='Permitted Value', max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='CatMolDescriptorValue',
            fields=[
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.molDescriptorValues.molUid, max_length=36)),
            ],
            options={
                'verbose_name': 'Categorical Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='CatRxnDescriptorValue',
            fields=[
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36)),
            ],
            options={
                'verbose_name': 'Categorical Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='ChemicalClass',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('label', models.CharField(unique=True, error_messages={'unique': 'A chemical class with this label already exists'}, max_length=30)),
                ('description', models.CharField(max_length=20)),
            ],
            options={
                'verbose_name_plural': 'Chemical Classes',
            },
        ),
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(verbose_name='Name', max_length=400)),
                ('CSID', models.PositiveIntegerField(unique=True, null=True, verbose_name='Chemspider ID')),
                ('custom', models.BooleanField(verbose_name='Custom', default=False)),
                ('INCHI', models.TextField(blank=True, verbose_name='InCHI key', default='')),
                ('smiles', models.TextField(blank=True, verbose_name='Smiles', default='')),
                ('dirty', models.BooleanField(default=True)),
                ('calculating', models.BooleanField(default=False)),
                ('recalculate', models.BooleanField(default=False)),
                ('formula', models.CharField(blank=True, help_text='A formula should be made up of element names. C_{4}H_{8} type notation should be use for subscript digits.', max_length=500, validators=[DRP.models.compound.elementsFormatValidator])),
                ('chemicalClasses', models.ManyToManyField(verbose_name='Chemical Class', to='DRP.ChemicalClass')),
            ],
        ),
        migrations.CreateModel(
            name='CompoundGuideEntry',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('abbrev', models.CharField(verbose_name='Abbreviation', max_length=100)),
                ('compound', models.ForeignKey(to='DRP.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='CompoundQuantity',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('amount', models.DecimalField(blank=True, decimal_places=5, help_text='(in mmoles, up to 5 decimal places)', null=True, max_digits=12, validators=[DRP.models.validators.GreaterThanValidator(0)])),
                ('compound', models.ForeignKey(to='DRP.Compound', on_delete=django.db.models.deletion.PROTECT)),
            ],
        ),
        migrations.CreateModel(
            name='CompoundRole',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('label', models.CharField(unique=True, max_length=255)),
                ('description', models.TextField()),
            ],
            options={
                'verbose_name': 'Compound Role Category',
                'verbose_name_plural': 'Compound Role Categories',
            },
        ),
        migrations.CreateModel(
            name='ConfirmationCode',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=36)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='DataSet',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=200)),
            ],
        ),
        migrations.CreateModel(
            name='DataSetRelation',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('dataSet', models.ForeignKey(to='DRP.DataSet')),
            ],
        ),
        migrations.CreateModel(
            name='Descriptor',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('heading', models.CharField(max_length=200, validators=[django.core.validators.RegexValidator('[A-Za-z0-9][A-Za-z0-9_]*', 'Please include only values which are limited to alphanumeric characters and underscores, and must start with an alphabetic character.')])),
                ('name', models.CharField(verbose_name='Full name', max_length=300)),
                ('calculatorSoftware', models.CharField(blank=True, max_length=100, validators=[django.core.validators.RegexValidator('[A-Za-z0-9][A-Za-z0-9_]*', 'Please include only values which are limited to alphanumeric characters and underscores, and must start with an alphabetic character.')])),
                ('calculatorSoftwareVersion', models.CharField(blank=True, max_length=20, validators=[django.core.validators.RegexValidator('[A-Za-z0-9][A-Za-z0-9_]*', 'Please include only values which are limited to alphanumeric characters, periods and underscores, and must start with an alphabetic character.')])),
            ],
        ),
        migrations.CreateModel(
            name='FeatureSelectionContainer',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('description', models.TextField(blank=True, default='')),
                ('featureVisitorLibrary', models.CharField(blank=True, max_length=200, default='')),
                ('featureVisitorTool', models.CharField(blank=True, max_length=200, default='')),
                ('featureVisitorOptions', models.TextField(blank=True, default='{}')),
                ('startTime', models.DateTimeField(blank=True, null=True, default=None)),
                ('endTime', models.DateTimeField(blank=True, null=True, default=None)),
                ('built', models.BooleanField(verbose_name='Has the build procedure been called with this container?', default=False, editable=False)),
                ('trainingSet', models.ForeignKey(related_name='trainingSetForFeatureSelection', null=True, to='DRP.DataSet')),
            ],
        ),
        migrations.CreateModel(
            name='LabGroup',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('title', models.CharField(unique=True, error_messages={'unique': 'This name is already taken.'}, max_length=200)),
                ('address', models.CharField(max_length=200)),
                ('email', models.CharField(max_length=254, default='')),
                ('access_code', models.CharField(max_length=128)),
                ('legacy_access_code', models.CharField(max_length=20)),
            ],
            options={
                'verbose_name': 'Lab Group',
            },
        ),
        migrations.CreateModel(
            name='License',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('text', models.TextField()),
                ('effectiveDate', models.DateField(verbose_name='Effective Date', help_text='The license will become effective on midnight of the provided date.')),
            ],
            options={
                'get_latest_by': 'effectiveDate',
            },
        ),
        migrations.CreateModel(
            name='LicenseAgreement',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('signedDateTime', models.DateTimeField(auto_now=True)),
                ('text', models.ForeignKey(to='DRP.License')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='MetricContainer',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('description', models.TextField(blank=True, default='')),
                ('metricVisitor', models.CharField(max_length=255)),
                ('startTime', models.DateTimeField(blank=True, null=True, default=None)),
                ('endTime', models.DateTimeField(blank=True, null=True, default=None)),
                ('fileName', models.FileField(blank=True, upload_to='metrics', max_length=200)),
                ('invalid', models.BooleanField(default=False)),
                ('built', models.BooleanField(verbose_name='Has the build procedure been called with this container?', default=False, editable=False)),
                ('trainingSet', models.ForeignKey(related_name='trainingSetForMetric', null=True, to='DRP.DataSet')),
            ],
        ),
        migrations.CreateModel(
            name='ModelContainer',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('description', models.TextField(blank=True)),
                ('active', models.BooleanField(verbose_name='Is this the active model?', default=False)),
                ('modelVisitorLibrary', models.CharField(max_length=200)),
                ('modelVisitorTool', models.CharField(max_length=200)),
                ('splitter', models.CharField(blank=True, max_length=200, default='')),
                ('modelVisitorOptions', models.TextField(blank=True, default='{}')),
                ('splitterOptions', models.TextField(blank=True, default='{}')),
                ('built', models.BooleanField(verbose_name='Has the build procedure been called with this container?', default=False, editable=False)),
            ],
        ),
        migrations.CreateModel(
            name='NumMolDescriptorValue',
            fields=[
                ('value', models.FloatField(blank=True, null=True)),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.molDescriptorValues.molUid, max_length=36)),
                ('compound', models.ForeignKey(to='DRP.Compound')),
            ],
            options={
                'verbose_name': 'Numeric Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='NumRxnDescriptorValue',
            fields=[
                ('value', models.FloatField(blank=True, null=True)),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36)),
            ],
            options={
                'verbose_name': 'Numeric Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='OrdMolDescriptorValue',
            fields=[
                ('value', models.IntegerField(blank=True, null=True)),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.molDescriptorValues.molUid, max_length=36)),
                ('compound', models.ForeignKey(to='DRP.Compound')),
                ('rater', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Ordinal Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='OrdRxnDescriptorValue',
            fields=[
                ('value', models.IntegerField(blank=True, null=True)),
                ('uid', models.CharField(serialize=False, primary_key=True, default=DRP.models.rxnDescriptorValues.rxnUid, max_length=36)),
                ('rater', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'verbose_name': 'Ordinal Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('notes', models.TextField(blank=True)),
                ('dirty', models.BooleanField(default=True)),
                ('calculating', models.BooleanField(default=False)),
                ('recalculate', models.BooleanField(default=False)),
            ],
        ),
        migrations.CreateModel(
            name='StatsModel',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('outputFile', models.FileField(blank=True, upload_to='models', max_length=200)),
                ('inputFile', models.FileField(blank=True, upload_to='model_inputs', max_length=255)),
                ('startTime', models.DateTimeField(null=True, default=None)),
                ('endTime', models.DateTimeField(null=True, default=None)),
                ('invalid', models.BooleanField(default=False)),
                ('container', models.ForeignKey(to='DRP.ModelContainer')),
                ('regenerationOf', models.ForeignKey(blank=True, to='DRP.StatsModel', null=True, default=None)),
                ('testSets', models.ManyToManyField(to='DRP.DataSet', related_name='testSetsFor')),
                ('trainingSet', models.ForeignKey(to='DRP.DataSet', related_name='trainingSetFor')),
            ],
        ),
        migrations.CreateModel(
            name='BooleanDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Descriptor')),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='CategoricalDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Descriptor')),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='NumericDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Descriptor')),
                ('maximum', models.FloatField(null=True)),
                ('minimum', models.FloatField(null=True)),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='OrdinalDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Descriptor')),
                ('maximum', models.IntegerField()),
                ('minimum', models.IntegerField()),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='PerformedReaction',
            fields=[
                ('reaction_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Reaction')),
                ('performedDateTime', models.DateTimeField(blank=True, default=None, null=True, help_text='Timezone assumed EST, Date in format YYYY-MM-DD', verbose_name='Date Reaction Performed', validators=[DRP.models.validators.notInTheFuture])),
                ('insertedDateTime', models.DateTimeField(verbose_name='Date Reaction Saved', auto_now_add=True)),
                ('legacyRecommendedFlag', models.NullBooleanField(default=None)),
                ('valid', models.BooleanField(default=True)),
                ('public', models.BooleanField(default=False)),
                ('legacyID', models.IntegerField(blank=True, unique=True, null=True)),
                ('legacyRef', models.CharField(blank=True, null=True, max_length=40)),
                ('convertedLegacyRef', models.CharField(blank=True, null=True, max_length=40, validators=[django.core.validators.RegexValidator('^[a-z0-9._]*[a-z][a-z0-9._]*$', 'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')])),
                ('reference', models.CharField(max_length=40, validators=[django.core.validators.RegexValidator('^[a-z0-9\\._]*[a-z][a-z0-9\\._]*$', 'Please include only values which are limited to alphanumeric characters, underscores, periods, and must include at least one alphabetic character.')])),
                ('labBookPage', models.ImageField(upload_to=DRP.models.performedReaction.UploadLabNotesTo(), null=True, verbose_name='Lab Book Image', storage=DRP.models.fileStorage.OverwriteStorage(location='/home/h205c/DRP/sec_media/lab_notes', base_url='/database/lab_notes/'))),
                ('duplicateOf', models.ForeignKey(blank=True, related_name='duplicatedBy', null=True, to='DRP.PerformedReaction', verbose_name='Duplicate Of', default=None)),
                ('performedBy', models.ForeignKey(blank=True, related_name='performedReactions', null=True, to=settings.AUTH_USER_MODEL, verbose_name='Performed By', default=None)),
            ],
            bases=('DRP.reaction',),
        ),
        migrations.CreateModel(
            name='RecommendedReaction',
            fields=[
                ('reaction_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.Reaction')),
                ('score', models.FloatField()),
                ('nonsense', models.BooleanField(default=None)),
                ('hidden', models.BooleanField(default=None)),
                ('saved', models.BooleanField(default=None)),
                ('reference', models.CharField(verbose_name='Text Reference', max_length=200)),
            ],
            bases=('DRP.reaction',),
        ),
        migrations.AddField(
            model_name='reaction',
            name='compounds',
            field=models.ManyToManyField(through='DRP.CompoundQuantity', to='DRP.Compound'),
        ),
        migrations.AddField(
            model_name='reaction',
            name='labGroup',
            field=models.ForeignKey(verbose_name='Lab Group', to='DRP.LabGroup'),
        ),
        migrations.AddField(
            model_name='ordrxndescriptorvalue',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='numrxndescriptorvalue',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='fully_trained',
            field=models.ForeignKey(blank=True, to='DRP.StatsModel', null=True),
        ),
        migrations.AddField(
            model_name='labgroup',
            name='defaultDescriptors',
            field=models.ManyToManyField(blank=True, to='DRP.Descriptor', related_name='isDefaultForLabGroups'),
        ),
        migrations.AddField(
            model_name='labgroup',
            name='users',
            field=models.ManyToManyField(blank=True, to=settings.AUTH_USER_MODEL),
        ),
        migrations.AlterUniqueTogether(
            name='descriptor',
            unique_together=set([('heading', 'calculatorSoftware', 'calculatorSoftwareVersion')]),
        ),
        migrations.AddField(
            model_name='compoundquantity',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='compoundquantity',
            name='role',
            field=models.ForeignKey(to='DRP.CompoundRole'),
        ),
        migrations.AddField(
            model_name='compoundguideentry',
            name='labGroup',
            field=models.ForeignKey(to='DRP.LabGroup'),
        ),
        migrations.AddField(
            model_name='compound',
            name='labGroups',
            field=models.ManyToManyField(verbose_name='Lab Groups', through='DRP.CompoundGuideEntry', to='DRP.LabGroup'),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='value',
            field=models.ForeignKey(blank=True, to='DRP.CategoricalDescriptorPermittedValue', on_delete=django.db.models.deletion.PROTECT, null=True),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='compound',
            field=models.ForeignKey(to='DRP.Compound'),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='value',
            field=models.ForeignKey(blank=True, to='DRP.CategoricalDescriptorPermittedValue', on_delete=django.db.models.deletion.PROTECT, null=True),
        ),
        migrations.AddField(
            model_name='boolrxndescriptorvalue',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='boolmoldescriptorvalue',
            name='compound',
            field=models.ForeignKey(to='DRP.Compound'),
        ),
        migrations.CreateModel(
            name='BoolMolDescriptor',
            fields=[
                ('booleandescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.BooleanDescriptor')),
            ],
            options={
                'verbose_name': 'Boolean Molecular Descriptor',
            },
            bases=('DRP.booleandescriptor',),
        ),
        migrations.CreateModel(
            name='BoolRxnDescriptor',
            fields=[
                ('booleandescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.BooleanDescriptor')),
            ],
            options={
                'verbose_name': 'Boolean Reaction Descriptor',
            },
            bases=('DRP.booleandescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='CatMolDescriptor',
            fields=[
                ('categoricaldescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.CategoricalDescriptor')),
            ],
            options={
                'verbose_name': 'Categorical Molecular Descriptor',
            },
            bases=('DRP.categoricaldescriptor',),
        ),
        migrations.CreateModel(
            name='CatRxnDescriptor',
            fields=[
                ('categoricaldescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.CategoricalDescriptor')),
            ],
            options={
                'verbose_name': 'Categorical Reaction Descriptor',
            },
            bases=('DRP.categoricaldescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='NumMolDescriptor',
            fields=[
                ('numericdescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.NumericDescriptor')),
            ],
            options={
                'verbose_name': 'Numerical Molecular Descriptor',
            },
            bases=('DRP.numericdescriptor',),
        ),
        migrations.CreateModel(
            name='NumRxnDescriptor',
            fields=[
                ('numericdescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.NumericDescriptor')),
            ],
            options={
                'verbose_name': 'Numerical Reaction Descriptor',
            },
            bases=('DRP.numericdescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='OrdMolDescriptor',
            fields=[
                ('ordinaldescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.OrdinalDescriptor')),
            ],
            options={
                'verbose_name': 'Ordinal Molecular Descriptor',
            },
            bases=('DRP.ordinaldescriptor',),
        ),
        migrations.CreateModel(
            name='OrdRxnDescriptor',
            fields=[
                ('ordinaldescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.OrdinalDescriptor')),
            ],
            options={
                'verbose_name': 'Ordinal Reaction Descriptor',
            },
            bases=('DRP.ordinaldescriptor', models.Model),
        ),
        migrations.AddField(
            model_name='recommendedreaction',
            name='seed',
            field=models.ForeignKey(related_name='seeded', null=True, to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='performedreaction',
            name='recommendation',
            field=models.ForeignKey(blank=True, related_name='resultantExperiment', null=True, to='DRP.RecommendedReaction', default=None),
        ),
        migrations.AddField(
            model_name='performedreaction',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='ordrxndescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.OrdinalDescriptor'),
        ),
        migrations.AddField(
            model_name='ordmoldescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.OrdinalDescriptor'),
        ),
        migrations.AddField(
            model_name='numrxndescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.NumericDescriptor'),
        ),
        migrations.AddField(
            model_name='nummoldescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.NumericDescriptor'),
        ),
        migrations.AddField(
            model_name='datasetrelation',
            name='reaction',
            field=models.ForeignKey(to='DRP.PerformedReaction', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='dataset',
            name='reactions',
            field=models.ManyToManyField(through='DRP.DataSetRelation', to='DRP.PerformedReaction'),
        ),
        migrations.AlterUniqueTogether(
            name='compoundquantity',
            unique_together=set([('reaction', 'role', 'amount')]),
        ),
        migrations.AlterUniqueTogether(
            name='compoundguideentry',
            unique_together=set([('abbrev', 'labGroup'), ('compound', 'labGroup')]),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.CategoricalDescriptor'),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.CategoricalDescriptor'),
        ),
        migrations.AddField(
            model_name='categoricaldescriptorpermittedvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.CategoricalDescriptor', related_name='permittedValues'),
        ),
        migrations.AddField(
            model_name='boolrxndescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.BooleanDescriptor'),
        ),
        migrations.AddField(
            model_name='boolmoldescriptorvalue',
            name='descriptor',
            field=models.ForeignKey(to='DRP.BooleanDescriptor'),
        ),
        migrations.CreateModel(
            name='PredBoolRxnDescriptor',
            fields=[
                ('boolrxndescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.BoolRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Boolean Rxn Descriptor',
            },
            bases=('DRP.boolrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredCatRxnDescriptor',
            fields=[
                ('catrxndescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.CatRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Categorical Rxn Descriptor',
            },
            bases=('DRP.catrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredNumRxnDescriptor',
            fields=[
                ('numrxndescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.NumRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Numeric Rxn Descriptor',
            },
            bases=('DRP.numrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredOrdRxnDescriptor',
            fields=[
                ('ordrxndescriptor_ptr', models.OneToOneField(serialize=False, primary_key=True, parent_link=True, auto_created=True, to='DRP.OrdRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Ordinal Rxn Descriptor',
            },
            bases=('DRP.ordrxndescriptor', models.Model),
        ),
        migrations.AlterUniqueTogether(
            name='ordrxndescriptorvalue',
            unique_together=set([('reaction', 'descriptor')]),
        ),
        migrations.AlterUniqueTogether(
            name='ordmoldescriptorvalue',
            unique_together=set([('descriptor', 'compound')]),
        ),
        migrations.AlterUniqueTogether(
            name='numrxndescriptorvalue',
            unique_together=set([('reaction', 'descriptor')]),
        ),
        migrations.AlterUniqueTogether(
            name='nummoldescriptorvalue',
            unique_together=set([('descriptor', 'compound')]),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='boolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='catRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='numRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='ordRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeBoolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor', related_name='outcomeForModels'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeCatRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor', related_name='outcomeForModels'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeNumRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor', related_name='outcomeForModels'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeOrdRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor', related_name='outcomeForModels'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='boolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='catRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='numRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='ordRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeBoolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor', related_name='outcomeForMetrics'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeCatRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor', related_name='outcomeForMetrics'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeNumRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor', related_name='outcomeForMetrics'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeOrdRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor', related_name='outcomeForMetrics'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='transformedRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor', related_name='transformedByMetric'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='boolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='catRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='chosenBoolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor', related_name='chosenForFeatureSelection'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='chosenCatRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor', related_name='chosenForFeatureSelection'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='chosenNumRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor', related_name='chosenForFeatureSelection'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='chosenOrdRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor', related_name='chosenForFeatureSelection'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='numRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='ordRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='outcomeBoolRxnDescriptors',
            field=models.ManyToManyField(to='DRP.BoolRxnDescriptor', related_name='outcomeForFeatureSelection'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='outcomeCatRxnDescriptors',
            field=models.ManyToManyField(to='DRP.CatRxnDescriptor', related_name='outcomeForFeatureSelections'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='outcomeNumRxnDescriptors',
            field=models.ManyToManyField(to='DRP.NumRxnDescriptor', related_name='outcomeForFeatureSelections'),
        ),
        migrations.AddField(
            model_name='featureselectioncontainer',
            name='outcomeOrdRxnDescriptors',
            field=models.ManyToManyField(to='DRP.OrdRxnDescriptor', related_name='outcomeForFeatureSelections'),
        ),
        migrations.AlterUniqueTogether(
            name='datasetrelation',
            unique_together=set([('dataSet', 'reaction')]),
        ),
        migrations.AlterUniqueTogether(
            name='catrxndescriptorvalue',
            unique_together=set([('reaction', 'descriptor')]),
        ),
        migrations.AlterUniqueTogether(
            name='catmoldescriptorvalue',
            unique_together=set([('descriptor', 'compound')]),
        ),
        migrations.AlterUniqueTogether(
            name='categoricaldescriptorpermittedvalue',
            unique_together=set([('descriptor', 'value')]),
        ),
        migrations.AlterUniqueTogether(
            name='boolrxndescriptorvalue',
            unique_together=set([('reaction', 'descriptor')]),
        ),
        migrations.AlterUniqueTogether(
            name='boolmoldescriptorvalue',
            unique_together=set([('descriptor', 'compound')]),
        ),
        migrations.AddField(
            model_name='predordrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='predordrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(to='DRP.OrdRxnDescriptor', related_name='predition_of'),
        ),
        migrations.AddField(
            model_name='predordrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(null=True, to='DRP.StatsModel'),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(to='DRP.NumRxnDescriptor', related_name='prediction_of'),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(null=True, to='DRP.StatsModel'),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(to='DRP.CatRxnDescriptor', related_name='prediction_of'),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(null=True, to='DRP.StatsModel'),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(to='DRP.BoolRxnDescriptor', related_name='prediction_of'),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(null=True, to='DRP.StatsModel'),
        ),
    ]
