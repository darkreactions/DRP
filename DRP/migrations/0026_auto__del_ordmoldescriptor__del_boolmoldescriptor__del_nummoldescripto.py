# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Removing unique constraint on 'MolDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.delete_unique(u'DRP_moldescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Removing unique constraint on 'CatMolDescriptorPermitted', fields ['descriptor', 'value']
        db.delete_unique(u'DRP_catmoldescriptorpermitted', ['descriptor_id', 'value'])

        # Deleting model 'OrdMolDescriptor'
        db.delete_table(u'DRP_ordmoldescriptor')

        # Deleting model 'BoolMolDescriptor'
        db.delete_table(u'DRP_boolmoldescriptor')

        # Deleting model 'NumMolDescriptor'
        db.delete_table(u'DRP_nummoldescriptor')

        # Deleting model 'CatMolDescriptorPermitted'
        db.delete_table(u'DRP_catmoldescriptorpermitted')

        # Deleting model 'CatMolDescriptor'
        db.delete_table(u'DRP_catmoldescriptor')

        # Deleting model 'MolDescriptor'
        db.delete_table(u'DRP_moldescriptor')

        # Adding model 'CategoricalDescriptorPermittedValue'
        db.create_table(u'DRP_categoricaldescriptorpermittedvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=255)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(related_name='permittedValues', to=orm['DRP.CategoricalDescriptor'])),
        ))
        db.send_create_signal('DRP', ['CategoricalDescriptorPermittedValue'])

        # Adding unique constraint on 'CategoricalDescriptorPermittedValue', fields ['descriptor', 'value']
        db.create_unique(u'DRP_categoricaldescriptorpermittedvalue', ['descriptor_id', 'value'])

        # Adding model 'OrdinalDescriptor'
        db.create_table(u'DRP_ordinaldescriptor', (
            (u'descriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Descriptor'], unique=True, primary_key=True)),
            ('maximum', self.gf('django.db.models.fields.IntegerField')(null=True)),
            ('minimum', self.gf('django.db.models.fields.IntegerField')(null=True)),
        ))
        db.send_create_signal('DRP', ['OrdinalDescriptor'])

        # Adding model 'CategoricalDescriptor'
        db.create_table(u'DRP_categoricaldescriptor', (
            (u'descriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Descriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['CategoricalDescriptor'])

        # Adding model 'Descriptor'
        db.create_table(u'DRP_descriptor', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('heading', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            ('calculatorSoftware', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('calculatorSoftwareVersion', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['Descriptor'])

        # Adding unique constraint on 'Descriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.create_unique(u'DRP_descriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Adding model 'NumericDescriptor'
        db.create_table(u'DRP_numericdescriptor', (
            (u'descriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Descriptor'], unique=True, primary_key=True)),
            ('maximum', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('minimum', self.gf('django.db.models.fields.FloatField')(null=True)),
        ))
        db.send_create_signal('DRP', ['NumericDescriptor'])

        # Adding model 'BooleanDescriptor'
        db.create_table(u'DRP_booleandescriptor', (
            (u'descriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Descriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['BooleanDescriptor'])


        # Changing field 'BoolMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_boolmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.BooleanDescriptor']))

        # Changing field 'NumMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_nummoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.NumericDescriptor']))

        # Changing field 'OrdMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_ordmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.OrdinalDescriptor']))

        # Changing field 'CatMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_catmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CategoricalDescriptor']))

        # Changing field 'CatMolDescriptorValue.value'
        db.alter_column(u'DRP_catmoldescriptorvalue', 'value_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CategoricalDescriptorPermittedValue'], null=True))

    def backwards(self, orm):
        # Removing unique constraint on 'Descriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.delete_unique(u'DRP_descriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Removing unique constraint on 'CategoricalDescriptorPermittedValue', fields ['descriptor', 'value']
        db.delete_unique(u'DRP_categoricaldescriptorpermittedvalue', ['descriptor_id', 'value'])

        # Adding model 'OrdMolDescriptor'
        db.create_table(u'DRP_ordmoldescriptor', (
            (u'moldescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.MolDescriptor'], unique=True, primary_key=True)),
            ('minimum', self.gf('django.db.models.fields.IntegerField')(null=True)),
            ('maximum', self.gf('django.db.models.fields.IntegerField')(null=True)),
        ))
        db.send_create_signal('DRP', ['OrdMolDescriptor'])

        # Adding model 'BoolMolDescriptor'
        db.create_table(u'DRP_boolmoldescriptor', (
            (u'moldescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.MolDescriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['BoolMolDescriptor'])

        # Adding model 'NumMolDescriptor'
        db.create_table(u'DRP_nummoldescriptor', (
            (u'moldescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.MolDescriptor'], unique=True, primary_key=True)),
            ('minimum', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('maximum', self.gf('django.db.models.fields.FloatField')(null=True)),
        ))
        db.send_create_signal('DRP', ['NumMolDescriptor'])

        # Adding model 'CatMolDescriptorPermitted'
        db.create_table(u'DRP_catmoldescriptorpermitted', (
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(related_name='permittedValues', to=orm['DRP.CatMolDescriptor'])),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal('DRP', ['CatMolDescriptorPermitted'])

        # Adding unique constraint on 'CatMolDescriptorPermitted', fields ['descriptor', 'value']
        db.create_unique(u'DRP_catmoldescriptorpermitted', ['descriptor_id', 'value'])

        # Adding model 'CatMolDescriptor'
        db.create_table(u'DRP_catmoldescriptor', (
            (u'moldescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.MolDescriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['CatMolDescriptor'])

        # Adding model 'MolDescriptor'
        db.create_table(u'DRP_moldescriptor', (
            ('calculatorSoftwareVersion', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('calculatorSoftware', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('heading', self.gf('django.db.models.fields.CharField')(max_length=200)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
        ))
        db.send_create_signal('DRP', ['MolDescriptor'])

        # Adding unique constraint on 'MolDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.create_unique(u'DRP_moldescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Deleting model 'CategoricalDescriptorPermittedValue'
        db.delete_table(u'DRP_categoricaldescriptorpermittedvalue')

        # Deleting model 'OrdinalDescriptor'
        db.delete_table(u'DRP_ordinaldescriptor')

        # Deleting model 'CategoricalDescriptor'
        db.delete_table(u'DRP_categoricaldescriptor')

        # Deleting model 'Descriptor'
        db.delete_table(u'DRP_descriptor')

        # Deleting model 'NumericDescriptor'
        db.delete_table(u'DRP_numericdescriptor')

        # Deleting model 'BooleanDescriptor'
        db.delete_table(u'DRP_booleandescriptor')


        # Changing field 'BoolMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_boolmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.BoolMolDescriptor']))

        # Changing field 'NumMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_nummoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.NumMolDescriptor']))

        # Changing field 'OrdMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_ordmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.OrdMolDescriptor']))

        # Changing field 'CatMolDescriptorValue.descriptor'
        db.alter_column(u'DRP_catmoldescriptorvalue', 'descriptor_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CatMolDescriptor']))

        # Changing field 'CatMolDescriptorValue.value'
        db.alter_column(u'DRP_catmoldescriptorvalue', 'value_id', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CatMolDescriptorPermitted'], null=True))

    models = {
        'DRP.booleandescriptor': {
            'Meta': {'object_name': 'BooleanDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.boolmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'BoolMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.BooleanDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'})
        },
        'DRP.categoricaldescriptor': {
            'Meta': {'object_name': 'CategoricalDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.categoricaldescriptorpermittedvalue': {
            'Meta': {'unique_together': "(('descriptor', 'value'),)", 'object_name': 'CategoricalDescriptorPermittedValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'permittedValues'", 'to': "orm['DRP.CategoricalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '255'})
        },
        'DRP.catmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'CatMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptorPermittedValue']", 'null': 'True'})
        },
        'DRP.chemicalclass': {
            'Meta': {'object_name': 'ChemicalClass'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        'DRP.compound': {
            'CSID': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True'}),
            'INCHI': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'}),
            'Meta': {'unique_together': "(('abbrev', 'labGroup'), ('CSID', 'labGroup'))", 'object_name': 'Compound'},
            'abbrev': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'chemicalClasses': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.ChemicalClass']", 'symmetrical': 'False'}),
            'custom': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'labGroup': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.LabGroup']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'}),
            'smiles': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'})
        },
        'DRP.compoundquantity': {
            'Meta': {'object_name': 'CompoundQuantity'},
            'amount': ('django.db.models.fields.FloatField', [], {}),
            'amountUnit': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']", 'on_delete': 'models.PROTECT'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"})
        },
        'DRP.confirmationcode': {
            'Meta': {'object_name': 'ConfirmationCode'},
            'code': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '36'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        },
        'DRP.dataset': {
            'Meta': {'object_name': 'DataSet'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isTestSet': ('django.db.models.fields.BooleanField', [], {}),
            'isTrainingSet': ('django.db.models.fields.BooleanField', [], {}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']"}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.PerformedReaction']"})
        },
        'DRP.descriptor': {
            'Meta': {'unique_together': "(('heading', 'calculatorSoftware', 'calculatorSoftwareVersion'),)", 'object_name': 'Descriptor'},
            'calculatorSoftware': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'calculatorSoftwareVersion': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'heading': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'})
        },
        'DRP.labgroup': {
            'Meta': {'object_name': 'LabGroup'},
            'access_code': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'address': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'email': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '254'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'legacy_access_code': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'title': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.User']", 'symmetrical': 'False', 'blank': 'True'})
        },
        'DRP.legacystatsmodel': {
            'Meta': {'object_name': 'LegacyStatsModel'},
            'active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'confusion_table': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            'correct_vals': ('django.db.models.fields.CharField', [], {'default': '\'["3","4"]\'', 'max_length': '100'}),
            'description': ('django.db.models.fields.TextField', [], {'default': "''"}),
            'end_time': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'filename': ('django.db.models.fields.CharField', [], {'default': "'/home/padler1/programming/drp/DRP/modelsuntitled.model'", 'max_length': '128'}),
            'headers': ('django.db.models.fields.TextField', [], {'default': "'[]'"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'iterations': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'library': ('django.db.models.fields.CharField', [], {'default': "'weka'", 'max_length': '128'}),
            'response': ('django.db.models.fields.CharField', [], {'default': "'outcomoe'", 'max_length': '128'}),
            'start_time': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'tags': ('django.db.models.fields.TextField', [], {'default': "''"}),
            'title': ('django.db.models.fields.CharField', [], {'default': "'untitled'", 'max_length': '100'}),
            'tmp_confusion_table': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            'tool': ('django.db.models.fields.CharField', [], {'default': "'svc'", 'max_length': '128'}),
            'train_confusion_table': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            'usable': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'DRP.license': {
            'Meta': {'object_name': 'License'},
            'effectiveDate': ('django.db.models.fields.DateField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'text': ('django.db.models.fields.TextField', [], {})
        },
        'DRP.licenseagreement': {
            'Meta': {'object_name': 'LicenseAgreement'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'signedDateTime': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'text': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.License']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        'DRP.numericdescriptor': {
            'Meta': {'object_name': 'NumericDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'maximum': ('django.db.models.fields.FloatField', [], {'null': 'True'}),
            'minimum': ('django.db.models.fields.FloatField', [], {'null': 'True'})
        },
        'DRP.nummoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'NumMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.NumericDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {'null': 'True'})
        },
        'DRP.ordinaldescriptor': {
            'Meta': {'object_name': 'OrdinalDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'maximum': ('django.db.models.fields.IntegerField', [], {'null': 'True'}),
            'minimum': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.ordmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'OrdMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.OrdinalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.performedreaction': {
            'Meta': {'object_name': 'PerformedReaction', '_ormbases': ['DRP.Reaction']},
            'duplicateOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'duplicatedBy'", 'to': "orm['DRP.Reaction']"}),
            'inTestSetFor': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'testSet'", 'symmetrical': 'False', 'to': "orm['DRP.StatsModel']"}),
            'inTrainingSetFor': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'trainingSet'", 'symmetrical': 'False', 'to': "orm['DRP.StatsModel']"}),
            'legacyRecommendedFlag': ('django.db.models.fields.NullBooleanField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'performedDateTime': ('django.db.models.fields.DateTimeField', [], {}),
            'public': ('django.db.models.fields.BooleanField', [], {}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'recommendation': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'resultantExperiment'", 'null': 'True', 'to': "orm['DRP.RecommendedReaction']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"}),
            'valid': ('django.db.models.fields.BooleanField', [], {})
        },
        'DRP.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'compounds': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.Compound']", 'through': "orm['DRP.CompoundQuantity']", 'symmetrical': 'False'}),
            'descriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.RxnDescriptor']", 'through': "orm['DRP.RxnDescriptorValue']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'labGroup': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.LabGroup']"}),
            'leak': ('django.db.models.fields.BooleanField', [], {}),
            'notes': ('django.db.models.fields.TextField', [], {}),
            'purity': ('django.db.models.fields.IntegerField', [], {}),
            'slowCool': ('django.db.models.fields.BooleanField', [], {}),
            'temp': ('django.db.models.fields.IntegerField', [], {}),
            'time': ('django.db.models.fields.IntegerField', [], {})
        },
        'DRP.recommendedreaction': {
            'Meta': {'object_name': 'RecommendedReaction', '_ormbases': ['DRP.Reaction']},
            'hidden': ('django.db.models.fields.BooleanField', [], {}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'nonsense': ('django.db.models.fields.BooleanField', [], {}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'reference': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'saved': ('django.db.models.fields.BooleanField', [], {}),
            'score': ('django.db.models.fields.FloatField', [], {}),
            'seed': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'seeded'", 'null': 'True', 'to': "orm['DRP.Reaction']"})
        },
        'DRP.rxndescriptor': {
            'Meta': {'unique_together': "(('heading', 'calculatorSoftware', 'calculatorSoftwareVersion'),)", 'object_name': 'RxnDescriptor'},
            'calculatorSoftware': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'calculatorSoftwareVersion': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'heading': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kind': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'})
        },
        'DRP.rxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'RxnDescriptorValue'},
            'booleanValue': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'catValue': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True'}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.RxnDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isPredicted': ('django.db.models.fields.BooleanField', [], {}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'numValue': ('django.db.models.fields.FloatField', [], {'null': 'True'}),
            'ordValue': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'})
        },
        'DRP.statsmodel': {
            'Meta': {'object_name': 'StatsModel'},
            'active': ('django.db.models.fields.BooleanField', [], {}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'end_time': ('django.db.models.fields.DateTimeField', [], {'default': 'None', 'null': 'True'}),
            'fileName': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'iterations': ('django.db.models.fields.IntegerField', [], {}),
            'library': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'start_time': ('django.db.models.fields.DateTimeField', [], {}),
            'tags': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.StatsModelTag']", 'symmetrical': 'False'}),
            'tool': ('django.db.models.fields.CharField', [], {'max_length': '200'})
        },
        'DRP.statsmodeltag': {
            'Meta': {'object_name': 'StatsModelTag'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'text': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'})
        },
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['DRP']