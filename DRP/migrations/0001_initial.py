# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from DRP.models.Compound import elementsFormatValidator
import django.db.models.deletion
from django.conf import settings
import django.core.validators


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='BoolMolDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.NullBooleanField(
                    verbose_name=b'Value for descriptor')),
            ],
            options={
                'verbose_name': 'Boolean Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='BoolRxnDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.NullBooleanField(
                    verbose_name=b'Value for descriptor')),
            ],
            options={
                'verbose_name': 'Boolean Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='CategoricalDescriptorPermittedValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.CharField(
                    max_length=255, verbose_name=b'Permitted Value')),
            ],
        ),
        migrations.CreateModel(
            name='CatMolDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'verbose_name': 'Categorical Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='CatRxnDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'verbose_name': 'Categorical Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='ChemicalClass',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('label', models.CharField(unique=True, max_length=30, error_messages={
                 b'unique': b'A chemical class with this label already exists'})),
                ('description', models.CharField(max_length=20)),
            ],
            options={
                'verbose_name_plural': 'Chemical Classes',
            },
        ),
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('abbrev', models.CharField(
                    max_length=100, verbose_name=b'Abbreviation')),
                ('name', models.CharField(max_length=400, verbose_name=b'Name')),
                ('CSID', models.PositiveIntegerField(
                    null=True, verbose_name=b'Chemspider ID')),
                ('custom', models.BooleanField(
                    default=False, verbose_name=b'Custom')),
                ('INCHI', models.TextField(default=b'',
                                           verbose_name=b'InCHI key', blank=True)),
                ('smiles', models.TextField(default=b'',
                                            verbose_name=b'Smiles', blank=True)),
                ('formula', models.CharField(
                    blank=True, help_text=b'A formula should be made up of element names. C_{4}H_{8} type notation should be use for subscript digits.', max_length=500, validators=[elementsFormatValidator])),
                ('chemicalClasses', models.ManyToManyField(
                    to='DRP.ChemicalClass', verbose_name=b'Chemical Class')),
            ],
        ),
        migrations.CreateModel(
            name='CompoundQuantity',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('amount', models.FloatField(null=True, blank=True)),
                ('compound', models.ForeignKey(to='DRP.Compound',
                                               on_delete=django.db.models.deletion.PROTECT)),
            ],
        ),
        migrations.CreateModel(
            name='CompoundQuantityManager',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
            ],
        ),
        migrations.CreateModel(
            name='CompoundRole',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
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
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('code', models.CharField(unique=True, max_length=36)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='DataSet',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=200)),
            ],
        ),
        migrations.CreateModel(
            name='DataSetRelation',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('dataSet', models.ForeignKey(to='DRP.DataSet')),
            ],
        ),
        migrations.CreateModel(
            name='Descriptor',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('heading', models.CharField(max_length=200, validators=[django.core.validators.RegexValidator(
                    b'[A-Za-z0-9][A-Za-z0-9_]+', b'Please include only values which are limited toalphanumeric characters and underscores, and must startwith an alphabetic character.')])),
                ('name', models.CharField(max_length=300, verbose_name=b'Full name')),
                ('calculatorSoftware', models.CharField(max_length=100)),
                ('calculatorSoftwareVersion', models.CharField(max_length=20)),
            ],
        ),
        migrations.CreateModel(
            name='LabGroup',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(unique=True, max_length=200,
                                           error_messages={b'unique': b'This name is already taken.'})),
                ('address', models.CharField(max_length=200)),
                ('email', models.CharField(default=b'', max_length=254)),
                ('access_code', models.CharField(max_length=128)),
                ('legacy_access_code', models.CharField(max_length=20)),
                ('users', models.ManyToManyField(
                    to=settings.AUTH_USER_MODEL, blank=True)),
            ],
            options={
                'verbose_name': 'Lab Group',
            },
        ),
        migrations.CreateModel(
            name='License',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('text', models.TextField()),
                ('effectiveDate', models.DateField(
                    help_text=b'The license will become effective on midnight of the provided date.', verbose_name=b'Effective Date')),
            ],
            options={
                'get_latest_by': 'effectiveDate',
            },
        ),
        migrations.CreateModel(
            name='LicenseAgreement',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('signedDateTime', models.DateTimeField(auto_now=True)),
                ('text', models.ForeignKey(to='DRP.License')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='MetricContainer',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('metricVisitor', models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='ModelContainer',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('description', models.TextField()),
                ('active', models.BooleanField(default=False,
                                               verbose_name=b'Is this the active model?')),
                ('modelVisitorLibrary', models.CharField(
                    max_length=200, choices=[(b'weka', b'weka')])),
                ('modelVisitorTool', models.CharField(max_length=200, choices=[
                 ((b'weka', (b'SVM_PUK_basic', b'SVM_PUK_BCR', b'KNN', b'NaiveBayes', b'J48')), (b'weka', (b'SVM_PUK_basic', b'SVM_PUK_BCR', b'KNN', b'NaiveBayes', b'J48')))])),
                ('featureLibrary', models.CharField(default=b'',
                                                    max_length=200, choices=[(b'weka', b'weka')])),
                ('featureTool', models.CharField(default=b'', max_length=200, choices=[
                 ((b'weka', (b'SVM_PUK_basic', b'SVM_PUK_BCR', b'KNN', b'NaiveBayes', b'J48')), (b'weka', (b'SVM_PUK_basic', b'SVM_PUK_BCR', b'KNN', b'NaiveBayes', b'J48')))])),
                ('splitter', models.CharField(blank=True, max_length=200, null=True, choices=[(b'KFoldSplitter', b'KFoldSplitter'), (
                    b'MutualInfoSplitter', b'MutualInfoSplitter'), (b'NoSplitter', b'NoSplitter'), (b'RandomSplitter', b'RandomSplitter')])),
                ('built', models.BooleanField(default=False,
                                              verbose_name=b'Has the build procedure been called with this container?', editable=False)),
            ],
        ),
        migrations.CreateModel(
            name='NumMolDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.FloatField(null=True)),
                ('compound', models.ForeignKey(to='DRP.Compound')),
            ],
            options={
                'verbose_name': 'Numeric Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='NumRxnDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.FloatField(null=True)),
            ],
            options={
                'verbose_name': 'Numeric Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='OrdMolDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.IntegerField(null=True)),
                ('compound', models.ForeignKey(to='DRP.Compound')),
            ],
            options={
                'verbose_name': 'Ordinal Molecular Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='OrdRxnDescriptorValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('value', models.IntegerField(null=True)),
            ],
            options={
                'verbose_name': 'Ordinal Reaction Descriptor Value',
            },
        ),
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('notes', models.TextField(blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='StatsModel',
            fields=[
                ('id', models.AutoField(verbose_name='ID',
                                        serialize=False, auto_created=True, primary_key=True)),
                ('fileName', models.FileField(
                    max_length=200, upload_to=b'models', blank=True)),
                ('startTime', models.DateTimeField(default=None, null=True)),
                ('endTime', models.DateTimeField(default=None, null=True)),
                ('invalid', models.BooleanField(default=False)),
                ('container', models.ForeignKey(to='DRP.ModelContainer')),
                ('regenerationOf', models.ForeignKey(default=None,
                                                     blank=True, to='DRP.StatsModel', null=True)),
                ('testSets', models.ManyToManyField(
                    related_name='testSetsFor', to='DRP.DataSet')),
                ('trainingSet', models.ForeignKey(
                    related_name='trainingSetFor', to='DRP.DataSet')),
            ],
        ),
        migrations.CreateModel(
            name='BooleanDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                        primary_key=True, serialize=False, to='DRP.Descriptor')),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='CategoricalDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                        primary_key=True, serialize=False, to='DRP.Descriptor')),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='NumericDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                        primary_key=True, serialize=False, to='DRP.Descriptor')),
                ('maximum', models.FloatField(null=True)),
                ('minimum', models.FloatField(null=True)),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='OrdinalDescriptor',
            fields=[
                ('descriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                        primary_key=True, serialize=False, to='DRP.Descriptor')),
                ('maximum', models.IntegerField()),
                ('minimum', models.IntegerField()),
            ],
            bases=('DRP.descriptor',),
        ),
        migrations.CreateModel(
            name='PerformedReaction',
            fields=[
                ('reaction_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                      primary_key=True, serialize=False, to='DRP.Reaction')),
                ('reference', models.CharField(unique=True, max_length=40)),
                ('performedDateTime', models.DateTimeField(
                    default=None, help_text=b'Date in format YYYY-MM-DD', null=True, verbose_name=b'Date Reaction Performed')),
                ('insertedDateTime', models.DateTimeField(
                    auto_now_add=True, verbose_name=b'Date Reaction Saved')),
                ('legacyRecommendedFlag', models.NullBooleanField(default=None)),
                ('valid', models.BooleanField(default=True)),
                ('public', models.BooleanField(default=False)),
                ('duplicateOf', models.ForeignKey(related_name='duplicatedBy',
                                                  default=None, blank=True, to='DRP.PerformedReaction', null=True)),
                ('performedBy', models.ForeignKey(related_name='performedReactions',
                                                  default=None, to=settings.AUTH_USER_MODEL, null=True)),
            ],
            bases=('DRP.reaction',),
        ),
        migrations.CreateModel(
            name='RecommendedReaction',
            fields=[
                ('reaction_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                      primary_key=True, serialize=False, to='DRP.Reaction')),
                ('score', models.FloatField()),
                ('nonsense', models.BooleanField(default=None)),
                ('hidden', models.BooleanField(default=None)),
                ('saved', models.BooleanField(default=None)),
                ('reference', models.CharField(
                    max_length=200, verbose_name=b'Text Reference')),
            ],
            bases=('DRP.reaction',),
        ),
        migrations.AddField(
            model_name='reaction',
            name='compounds',
            field=models.ManyToManyField(
                to='DRP.Compound', through='DRP.CompoundQuantity'),
        ),
        migrations.AddField(
            model_name='reaction',
            name='labGroup',
            field=models.ForeignKey(to='DRP.LabGroup'),
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
            field=models.ForeignKey(to='DRP.StatsModel', null=True),
        ),
        migrations.AlterUniqueTogether(
            name='descriptor',
            unique_together=set(
                [('heading', 'calculatorSoftware', 'calculatorSoftwareVersion')]),
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
            model_name='compound',
            name='labGroup',
            field=models.ForeignKey(
                verbose_name=b'Lab Group', to='DRP.LabGroup'),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='reaction',
            field=models.ForeignKey(to='DRP.Reaction'),
        ),
        migrations.AddField(
            model_name='catrxndescriptorvalue',
            name='value',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT,
                                    to='DRP.CategoricalDescriptorPermittedValue', null=True),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='compound',
            field=models.ForeignKey(to='DRP.Compound'),
        ),
        migrations.AddField(
            model_name='catmoldescriptorvalue',
            name='value',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT,
                                    to='DRP.CategoricalDescriptorPermittedValue', null=True),
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
                ('booleandescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.BooleanDescriptor')),
            ],
            options={
                'verbose_name': 'Boolean Molecular Descriptor',
            },
            bases=('DRP.booleandescriptor',),
        ),
        migrations.CreateModel(
            name='BoolRxnDescriptor',
            fields=[
                ('booleandescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.BooleanDescriptor')),
            ],
            options={
                'verbose_name': 'Boolean Reaction Descriptor',
            },
            bases=('DRP.booleandescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='CatMolDescriptor',
            fields=[
                ('categoricaldescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                                   primary_key=True, serialize=False, to='DRP.CategoricalDescriptor')),
            ],
            options={
                'verbose_name': 'Categorical Molecular Descriptor',
            },
            bases=('DRP.categoricaldescriptor',),
        ),
        migrations.CreateModel(
            name='CatRxnDescriptor',
            fields=[
                ('categoricaldescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                                   primary_key=True, serialize=False, to='DRP.CategoricalDescriptor')),
            ],
            options={
                'verbose_name': 'Categorical Reaction Descriptor',
            },
            bases=('DRP.categoricaldescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='NumMolDescriptor',
            fields=[
                ('numericdescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.NumericDescriptor')),
            ],
            options={
                'verbose_name': 'Numerical Molecular Descriptor',
            },
            bases=('DRP.numericdescriptor',),
        ),
        migrations.CreateModel(
            name='NumRxnDescriptor',
            fields=[
                ('numericdescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.NumericDescriptor')),
            ],
            options={
                'verbose_name': 'Numerical Reaction Descriptor',
            },
            bases=('DRP.numericdescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='OrdMolDescriptor',
            fields=[
                ('ordinaldescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.OrdinalDescriptor')),
            ],
            options={
                'verbose_name': 'Ordinal Molecular Descriptor',
            },
            bases=('DRP.ordinaldescriptor',),
        ),
        migrations.CreateModel(
            name='OrdRxnDescriptor',
            fields=[
                ('ordinaldescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.OrdinalDescriptor')),
            ],
            options={
                'verbose_name': 'Ordinal Reaction Descriptor',
            },
            bases=('DRP.ordinaldescriptor', models.Model),
        ),
        migrations.AddField(
            model_name='recommendedreaction',
            name='seed',
            field=models.ForeignKey(
                related_name='seeded', to='DRP.Reaction', null=True),
        ),
        migrations.AddField(
            model_name='performedreaction',
            name='recommendation',
            field=models.ForeignKey(related_name='resultantExperiment', default=None,
                                    blank=True, to='DRP.RecommendedReaction', null=True),
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
            field=models.ForeignKey(
                to='DRP.PerformedReaction', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='dataset',
            name='reactions',
            field=models.ManyToManyField(
                to='DRP.PerformedReaction', through='DRP.DataSetRelation'),
        ),
        migrations.AlterUniqueTogether(
            name='compoundquantity',
            unique_together=set([('reaction', 'role', 'amount')]),
        ),
        migrations.AlterUniqueTogether(
            name='compound',
            unique_together=set(
                [('CSID', 'labGroup'), ('abbrev', 'labGroup')]),
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
            field=models.ForeignKey(
                related_name='permittedValues', to='DRP.CategoricalDescriptor'),
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
                ('boolrxndescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                               primary_key=True, serialize=False, to='DRP.BoolRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Boolean Rxn Descriptor',
            },
            bases=('DRP.boolrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredCatRxnDescriptor',
            fields=[
                ('catrxndescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                              primary_key=True, serialize=False, to='DRP.CatRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Categorical Rxn Descriptor',
            },
            bases=('DRP.catrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredNumRxnDescriptor',
            fields=[
                ('numrxndescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                              primary_key=True, serialize=False, to='DRP.NumRxnDescriptor')),
            ],
            options={
                'verbose_name': 'Predicted Numeric Rxn Descriptor',
            },
            bases=('DRP.numrxndescriptor', models.Model),
        ),
        migrations.CreateModel(
            name='PredOrdRxnDescriptor',
            fields=[
                ('ordrxndescriptor_ptr', models.OneToOneField(parent_link=True, auto_created=True,
                                                              primary_key=True, serialize=False, to='DRP.OrdRxnDescriptor')),
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
            field=models.ManyToManyField(
                related_name='outcomeForModels', to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeCatRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForModels', to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeNumRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForModels', to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='modelcontainer',
            name='outcomeOrdRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForModels', to='DRP.OrdRxnDescriptor'),
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
            field=models.ManyToManyField(
                related_name='outcomeForMetrics', to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeCatRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForMetrics', to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeNumRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForMetrics', to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='outcomeOrdRxnDescriptors',
            field=models.ManyToManyField(
                related_name='outcomeForMetrics', to='DRP.OrdRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='metriccontainer',
            name='transformedRxnDescriptors',
            field=models.ManyToManyField(
                related_name='transformedByMetric', to='DRP.NumRxnDescriptor'),
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
            field=models.ForeignKey(
                related_name='predition_of', to='DRP.OrdRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='predordrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(to='DRP.StatsModel', null=True),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(
                related_name='prediction_of', to='DRP.NumRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='prednumrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(to='DRP.StatsModel', null=True),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(
                related_name='prediction_of', to='DRP.CatRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='predcatrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(to='DRP.StatsModel', null=True),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='modelContainer',
            field=models.ForeignKey(to='DRP.ModelContainer'),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='predictionOf',
            field=models.ForeignKey(
                related_name='prediction_of', to='DRP.BoolRxnDescriptor'),
        ),
        migrations.AddField(
            model_name='predboolrxndescriptor',
            name='statsModel',
            field=models.ForeignKey(to='DRP.StatsModel', null=True),
        ),
    ]
