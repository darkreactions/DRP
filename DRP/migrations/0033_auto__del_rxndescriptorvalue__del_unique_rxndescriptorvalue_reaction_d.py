# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Removing unique constraint on 'RxnDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.delete_unique(u'DRP_rxndescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Removing unique constraint on 'RxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_rxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Deleting model 'RxnDescriptorValue'
        db.delete_table(u'DRP_rxndescriptorvalue')

        # Deleting model 'RxnDescriptor'
        db.delete_table(u'DRP_rxndescriptor')

        # Deleting model 'DataSet'
        db.delete_table(u'DRP_dataset')

        # Adding model 'CatRxnDescriptor'
        db.create_table(u'DRP_catrxndescriptor', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['CatRxnDescriptor'])

        # Adding model 'BoolRxnDescriptorValue'
        db.create_table(u'DRP_boolrxndescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.BooleanDescriptor'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['BoolRxnDescriptorValue'])

        # Adding unique constraint on 'BoolRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_boolrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'OrdRxnDescriptorValue'
        db.create_table(u'DRP_ordrxndescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.IntegerField')(null=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.OrdinalDescriptor'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['OrdRxnDescriptorValue'])

        # Adding unique constraint on 'OrdRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_ordrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'BoolRxnDescriptor'
        db.create_table(u'DRP_boolrxndescriptor', (
            (u'booleandescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.BooleanDescriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['BoolRxnDescriptor'])

        # Adding model 'NumRxnDescriptorValue'
        db.create_table(u'DRP_numrxndescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.NumericDescriptor'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['NumRxnDescriptorValue'])

        # Adding unique constraint on 'NumRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_numrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'CatRxnDescriptorValue'
        db.create_table(u'DRP_catrxndescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CategoricalDescriptor'])),
            ('value', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.CategoricalDescriptorPermittedValue'], null=True)),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['CatRxnDescriptorValue'])

        # Adding unique constraint on 'CatRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_catrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'OrdRxnDescriptor'
        db.create_table(u'DRP_ordrxndescriptor', (
            (u'ordinaldescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.OrdinalDescriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['OrdRxnDescriptor'])

        # Adding model 'NumRxnDescriptor'
        db.create_table(u'DRP_numrxndescriptor', (
            (u'numericdescriptor_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.NumericDescriptor'], unique=True, primary_key=True)),
        ))
        db.send_create_signal('DRP', ['NumRxnDescriptor'])


    def backwards(self, orm):
        # Removing unique constraint on 'CatRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_catrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Removing unique constraint on 'NumRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_numrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Removing unique constraint on 'OrdRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_ordrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Removing unique constraint on 'BoolRxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_boolrxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'RxnDescriptorValue'
        db.create_table(u'DRP_rxndescriptorvalue', (
            ('ordValue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.RxnDescriptor'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('catValue', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('numValue', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
            ('booleanValue', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('isPredicted', self.gf('django.db.models.fields.BooleanField')()),
        ))
        db.send_create_signal('DRP', ['RxnDescriptorValue'])

        # Adding unique constraint on 'RxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_rxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'RxnDescriptor'
        db.create_table(u'DRP_rxndescriptor', (
            ('kind', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('calculatorSoftwareVersion', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('calculatorSoftware', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('heading', self.gf('django.db.models.fields.CharField')(max_length=200, unique=True)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
        ))
        db.send_create_signal('DRP', ['RxnDescriptor'])

        # Adding unique constraint on 'RxnDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.create_unique(u'DRP_rxndescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Adding model 'DataSet'
        db.create_table(u'DRP_dataset', (
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.PerformedReaction'])),
            ('isTrainingSet', self.gf('django.db.models.fields.BooleanField')()),
            ('isTestSet', self.gf('django.db.models.fields.BooleanField')()),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'])),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['DataSet'])

        # Deleting model 'CatRxnDescriptor'
        db.delete_table(u'DRP_catrxndescriptor')

        # Deleting model 'BoolRxnDescriptorValue'
        db.delete_table(u'DRP_boolrxndescriptorvalue')

        # Deleting model 'OrdRxnDescriptorValue'
        db.delete_table(u'DRP_ordrxndescriptorvalue')

        # Deleting model 'BoolRxnDescriptor'
        db.delete_table(u'DRP_boolrxndescriptor')

        # Deleting model 'NumRxnDescriptorValue'
        db.delete_table(u'DRP_numrxndescriptorvalue')

        # Deleting model 'CatRxnDescriptorValue'
        db.delete_table(u'DRP_catrxndescriptorvalue')

        # Deleting model 'OrdRxnDescriptor'
        db.delete_table(u'DRP_ordrxndescriptor')

        # Deleting model 'NumRxnDescriptor'
        db.delete_table(u'DRP_numrxndescriptor')


    models = {
        'DRP.booleandescriptor': {
            'Meta': {'object_name': 'BooleanDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.boolmoldescriptor': {
            'Meta': {'object_name': 'BoolMolDescriptor', '_ormbases': ['DRP.BooleanDescriptor']},
            u'booleandescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.BooleanDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.boolmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'BoolMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.BooleanDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'})
        },
        'DRP.boolrxndescriptor': {
            'Meta': {'object_name': 'BoolRxnDescriptor', '_ormbases': ['DRP.BooleanDescriptor']},
            u'booleandescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.BooleanDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.boolrxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'BoolRxnDescriptorValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.BooleanDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'}),
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
        'DRP.catmoldescriptor': {
            'Meta': {'object_name': 'CatMolDescriptor', '_ormbases': ['DRP.CategoricalDescriptor']},
            u'categoricaldescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.CategoricalDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.catmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'CatMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptorPermittedValue']", 'null': 'True'})
        },
        'DRP.catrxndescriptor': {
            'Meta': {'object_name': 'CatRxnDescriptor'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'DRP.catrxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'CatRxnDescriptorValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'}),
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
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']", 'on_delete': 'models.PROTECT'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"}),
            'role': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CompoundRole']"})
        },
        'DRP.compoundrole': {
            'Meta': {'object_name': 'CompoundRole'},
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
        'DRP.confirmationcode': {
            'Meta': {'object_name': 'ConfirmationCode'},
            'code': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '36'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
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
        'DRP.nummoldescriptor': {
            'Meta': {'object_name': 'NumMolDescriptor', '_ormbases': ['DRP.NumericDescriptor']},
            u'numericdescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.NumericDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.nummoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'NumMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.NumericDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {'null': 'True'})
        },
        'DRP.numrxndescriptor': {
            'Meta': {'object_name': 'NumRxnDescriptor', '_ormbases': ['DRP.NumericDescriptor']},
            u'numericdescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.NumericDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.numrxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'NumRxnDescriptorValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.NumericDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'}),
            'value': ('django.db.models.fields.FloatField', [], {'null': 'True'})
        },
        'DRP.ordinaldescriptor': {
            'Meta': {'object_name': 'OrdinalDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'maximum': ('django.db.models.fields.IntegerField', [], {'null': 'True'}),
            'minimum': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.ordmoldescriptor': {
            'Meta': {'object_name': 'OrdMolDescriptor', '_ormbases': ['DRP.OrdinalDescriptor']},
            u'ordinaldescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.OrdinalDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.ordmoldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'OrdMolDescriptorValue'},
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.OrdinalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.ordrxndescriptor': {
            'Meta': {'object_name': 'OrdRxnDescriptor', '_ormbases': ['DRP.OrdinalDescriptor']},
            u'ordinaldescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.OrdinalDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.ordrxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'OrdRxnDescriptorValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.OrdinalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'}),
            'value': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.performedreaction': {
            'Meta': {'object_name': 'PerformedReaction', '_ormbases': ['DRP.Reaction']},
            'duplicateOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'duplicatedBy'", 'to': "orm['DRP.PerformedReaction']"}),
            'inTestSetFor': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'testSet'", 'symmetrical': 'False', 'to': "orm['DRP.StatsModel']"}),
            'inTrainingSetFor': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'trainingSet'", 'symmetrical': 'False', 'to': "orm['DRP.StatsModel']"}),
            'legacyRecommendedFlag': ('django.db.models.fields.NullBooleanField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'performedDateTime': ('django.db.models.fields.DateTimeField', [], {}),
            'public': ('django.db.models.fields.BooleanField', [], {}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'recommendation': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'resultantExperiment'", 'null': 'True', 'to': "orm['DRP.RecommendedReaction']"}),
            'reference': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '40'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"}),
            'valid': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'DRP.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'compounds': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.Compound']", 'through': "orm['DRP.CompoundQuantity']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'labGroup': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.LabGroup']"}),
            'notes': ('django.db.models.fields.TextField', [], {})
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