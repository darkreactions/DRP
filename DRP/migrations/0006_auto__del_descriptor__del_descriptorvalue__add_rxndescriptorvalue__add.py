# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting model 'Descriptor'
        db.delete_table(u'DRP_descriptor')

        # Deleting model 'DescriptorValue'
        db.delete_table(u'DRP_descriptorvalue')

        # Adding model 'RxnDescriptorValue'
        db.create_table(u'DRP_rxndescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.RxnDescriptor'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('booleanValue', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('ordValue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True)),
            ('catValue', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('numValue', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('isPredicted', self.gf('django.db.models.fields.BooleanField')()),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['RxnDescriptorValue'])

        # Adding unique constraint on 'RxnDescriptorValue', fields ['reaction', 'descriptor']
        db.create_unique(u'DRP_rxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'MolDescriptorValue'
        db.create_table(u'DRP_moldescriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.MolDescriptor'])),
            ('compound', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Compound'], null=True)),
            ('booleanValue', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('ordValue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True)),
            ('catValue', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('numValue', self.gf('django.db.models.fields.FloatField')(null=True)),
        ))
        db.send_create_signal('DRP', ['MolDescriptorValue'])

        # Adding unique constraint on 'MolDescriptorValue', fields ['descriptor', 'compound']
        db.create_unique(u'DRP_moldescriptorvalue', ['descriptor_id', 'compound_id'])

        # Adding model 'RxnDescriptor'
        db.create_table(u'DRP_rxndescriptor', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('heading', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            ('kind', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('calculatorSoftware', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('calculatorSoftwareVersion', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['RxnDescriptor'])

        # Adding unique constraint on 'RxnDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.create_unique(u'DRP_rxndescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Adding model 'MolDescriptor'
        db.create_table(u'DRP_moldescriptor', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('heading', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            ('kind', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('calculatorSoftware', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('calculatorSoftwareVersion', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['MolDescriptor'])

        # Adding unique constraint on 'MolDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.create_unique(u'DRP_moldescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])


    def backwards(self, orm):
        # Removing unique constraint on 'MolDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.delete_unique(u'DRP_moldescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Removing unique constraint on 'RxnDescriptor', fields ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion']
        db.delete_unique(u'DRP_rxndescriptor', ['heading', 'calculatorSoftware', 'calculatorSoftwareVersion'])

        # Removing unique constraint on 'MolDescriptorValue', fields ['descriptor', 'compound']
        db.delete_unique(u'DRP_moldescriptorvalue', ['descriptor_id', 'compound_id'])

        # Removing unique constraint on 'RxnDescriptorValue', fields ['reaction', 'descriptor']
        db.delete_unique(u'DRP_rxndescriptorvalue', ['reaction_id', 'descriptor_id'])

        # Adding model 'Descriptor'
        db.create_table(u'DRP_descriptor', (
            ('kind', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('heading', self.gf('django.db.models.fields.CharField')(max_length=200, unique=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['Descriptor'])

        # Adding model 'DescriptorValue'
        db.create_table(u'DRP_descriptorvalue', (
            ('Reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('catValue', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('booleanValue', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('isPredicted', self.gf('django.db.models.fields.BooleanField')()),
            ('ordValue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Descriptor'])),
            ('Compound', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Compound'], null=True)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
            ('numValue', self.gf('django.db.models.fields.FloatField')(null=True)),
        ))
        db.send_create_signal('DRP', ['DescriptorValue'])

        # Deleting model 'RxnDescriptorValue'
        db.delete_table(u'DRP_rxndescriptorvalue')

        # Deleting model 'MolDescriptorValue'
        db.delete_table(u'DRP_moldescriptorvalue')

        # Deleting model 'RxnDescriptor'
        db.delete_table(u'DRP_rxndescriptor')

        # Deleting model 'MolDescriptor'
        db.delete_table(u'DRP_moldescriptor')


    models = {
        'DRP.chemicalclass': {
            'Meta': {'object_name': 'ChemicalClass'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        'DRP.compound': {
            'CAS_ID': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '13', 'blank': 'True'}),
            'CHEBI_ID': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'CSID': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'ChemicalClass': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.ChemicalClass']", 'symmetrical': 'False'}),
            'INCHI': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'}),
            'Meta': {'object_name': 'Compound'},
            'abbrev': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'custom': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'descriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.MolDescriptor']", 'through': "orm['DRP.MolDescriptorValue']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'}),
            'smiles': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'})
        },
        'DRP.compoundquantity': {
            'Meta': {'object_name': 'CompoundQuantity'},
            'amount': ('django.db.models.fields.FloatField', [], {}),
            'amountUnit': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Compound']"}),
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
            'filename': ('django.db.models.fields.CharField', [], {'default': "'/home/padler1/programming/drp/DRP/models/untitled.model'", 'max_length': '128'}),
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
        'DRP.moldescriptor': {
            'Meta': {'unique_together': "(('heading', 'calculatorSoftware', 'calculatorSoftwareVersion'),)", 'object_name': 'MolDescriptor'},
            'calculatorSoftware': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'calculatorSoftwareVersion': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'heading': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kind': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'})
        },
        'DRP.moldescriptorvalue': {
            'Meta': {'unique_together': "(('descriptor', 'compound'),)", 'object_name': 'MolDescriptorValue'},
            'booleanValue': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'catValue': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True'}),
            'compound': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Compound']", 'null': 'True'}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.MolDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'numValue': ('django.db.models.fields.FloatField', [], {'null': 'True'}),
            'ordValue': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True'})
        },
        'DRP.performedreaction': {
            'Meta': {'object_name': 'PerformedReaction', '_ormbases': ['DRP.Reaction']},
            'duplicateOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'duplicatedBy'", 'to': "orm['DRP.Reaction']"}),
            'legacyRecommendedFlag': ('django.db.models.fields.NullBooleanField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'performedDateTime': ('django.db.models.fields.DateTimeField', [], {}),
            'public': ('django.db.models.fields.BooleanField', [], {}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'recommendation': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'resultantExperiment'", 'null': 'True', 'to': "orm['DRP.RecommendedReaction']"}),
            'usedForModel': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.StatsModel']", 'through': "orm['DRP.DataSet']", 'symmetrical': 'False'}),
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