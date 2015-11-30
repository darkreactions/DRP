# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Descriptor'
        db.create_table(u'DRP_descriptor', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('heading', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            ('kind', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['Descriptor'])

        # Adding model 'ChemicalClass'
        db.create_table(u'DRP_chemicalclass', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('label', self.gf('django.db.models.fields.CharField')(unique=True, max_length=30)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['ChemicalClass'])

        # Adding model 'LabGroup'
        db.create_table(u'DRP_labgroup', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('title', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('address', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('email', self.gf('django.db.models.fields.CharField')(default='', max_length=254)),
            ('access_code', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('legacy_access_code', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal('DRP', ['LabGroup'])

        # Adding M2M table for field users on 'LabGroup'
        db.create_table(u'DRP_labgroup_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('labgroup', models.ForeignKey(orm['DRP.labgroup'], null=False)),
            ('user', models.ForeignKey(orm[u'auth.user'], null=False))
        ))
        db.create_unique(u'DRP_labgroup_users', ['labgroup_id', 'user_id'])

        # Adding model 'CG_calculations'
        db.create_table(u'DRP_cg_calculations', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('json_data', self.gf('django.db.models.fields.TextField')()),
            ('compound', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('smiles', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('json', self.gf('django.db.models.fields.TextField')(default='{}', null=True)),
        ))
        db.send_create_signal('DRP', ['CG_calculations'])

        # Adding model 'Compound'
        db.create_table(u'DRP_compound', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('abbrev', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=300)),
            ('CAS_ID', self.gf('django.db.models.fields.CharField')(default='', max_length=13, blank=True)),
            ('CHEBI_ID', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('CSID', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('custom', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('INCHI', self.gf('django.db.models.fields.TextField')(default='', blank=True)),
            ('smiles', self.gf('django.db.models.fields.TextField')(default='', blank=True)),
        ))
        db.send_create_signal('DRP', ['Compound'])

        # Adding M2M table for field ChemicalClass on 'Compound'
        db.create_table(u'DRP_compound_ChemicalClass', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('compound', models.ForeignKey(orm['DRP.compound'], null=False)),
            ('chemicalclass', models.ForeignKey(orm['DRP.chemicalclass'], null=False))
        ))
        db.create_unique(u'DRP_compound_ChemicalClass', ['compound_id', 'chemicalclass_id'])

        # Adding model 'Reaction'
        db.create_table(u'DRP_reaction', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('temp', self.gf('django.db.models.fields.IntegerField')()),
            ('slowCool', self.gf('django.db.models.fields.BooleanField')()),
            ('time', self.gf('django.db.models.fields.IntegerField')()),
            ('leak', self.gf('django.db.models.fields.BooleanField')()),
            ('purity', self.gf('django.db.models.fields.IntegerField')()),
            ('notes', self.gf('django.db.models.fields.TextField')()),
            ('labGroup', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.LabGroup'])),
        ))
        db.send_create_signal('DRP', ['Reaction'])

        # Adding model 'CompoundQuantity'
        db.create_table(u'DRP_compoundquantity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('compound', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Compound'])),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Reaction'])),
            ('amount', self.gf('django.db.models.fields.FloatField')()),
            ('amountUnit', self.gf('django.db.models.fields.CharField')(max_length=10)),
        ))
        db.send_create_signal('DRP', ['CompoundQuantity'])

        # Adding model 'StatsModelTag'
        db.create_table(u'DRP_statsmodeltag', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('text', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
        ))
        db.send_create_signal('DRP', ['StatsModelTag'])

        # Adding model 'StatsModel'
        db.create_table(u'DRP_statsmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('fileName', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('active', self.gf('django.db.models.fields.BooleanField')()),
            ('start_time', self.gf('django.db.models.fields.DateTimeField')()),
            ('end_time', self.gf('django.db.models.fields.DateTimeField')(default=None, null=True)),
            ('iterations', self.gf('django.db.models.fields.IntegerField')()),
            ('library', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('tool', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('DRP', ['StatsModel'])

        # Adding M2M table for field tags on 'StatsModel'
        db.create_table(u'DRP_statsmodel_tags', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('statsmodeltag', models.ForeignKey(orm['DRP.statsmodeltag'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_tags', ['statsmodel_id', 'statsmodeltag_id'])

        # Adding model 'RecommendedReaction'
        db.create_table(u'DRP_recommendedreaction', (
            (u'reaction_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Reaction'], unique=True, primary_key=True)),
            ('score', self.gf('django.db.models.fields.FloatField')()),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True)),
            ('seed', self.gf('django.db.models.fields.related.ForeignKey')(related_name='seeded', null=True, to=orm['DRP.Reaction'])),
            ('nonsense', self.gf('django.db.models.fields.BooleanField')()),
            ('hidden', self.gf('django.db.models.fields.BooleanField')()),
            ('saved', self.gf('django.db.models.fields.BooleanField')()),
            ('reference', self.gf('django.db.models.fields.CharField')(max_length=200)),
        ))
        db.send_create_signal('DRP', ['RecommendedReaction'])

        # Adding model 'PerformedReaction'
        db.create_table(u'DRP_performedreaction', (
            (u'reaction_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['DRP.Reaction'], unique=True, primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('performedDateTime', self.gf('django.db.models.fields.DateTimeField')()),
            ('recommendation', self.gf('django.db.models.fields.related.ForeignKey')(default=None, related_name='resultantExperiment', null=True, to=orm['DRP.RecommendedReaction'])),
            ('legacyRecommendedFlag', self.gf('django.db.models.fields.NullBooleanField')(default=None, null=True, blank=True)),
            ('valid', self.gf('django.db.models.fields.BooleanField')()),
            ('public', self.gf('django.db.models.fields.BooleanField')()),
            ('duplicateOf', self.gf('django.db.models.fields.related.ForeignKey')(related_name='duplicatedBy', to=orm['DRP.Reaction'])),
        ))
        db.send_create_signal('DRP', ['PerformedReaction'])

        # Adding model 'DataSet'
        db.create_table(u'DRP_dataset', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.PerformedReaction'])),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'])),
            ('isTestSet', self.gf('django.db.models.fields.BooleanField')()),
            ('isTrainingSet', self.gf('django.db.models.fields.BooleanField')()),
        ))
        db.send_create_signal('DRP', ['DataSet'])

        # Adding model 'DescriptorValue'
        db.create_table(u'DRP_descriptorvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('descriptor', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Descriptor'])),
            ('Reaction', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Reaction'], null=True)),
            ('Compound', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Compound'], null=True)),
            ('booleanValue', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('ordValue', self.gf('django.db.models.fields.PositiveIntegerField')(null=True)),
            ('catValue', self.gf('django.db.models.fields.CharField')(max_length=200, null=True)),
            ('numValue', self.gf('django.db.models.fields.FloatField')(null=True)),
            ('isPredicted', self.gf('django.db.models.fields.BooleanField')()),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True)),
        ))
        db.send_create_signal('DRP', ['DescriptorValue'])

        # Adding model 'LegacyStatsModel'
        db.create_table(u'DRP_legacystatsmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            ('train_confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            ('tmp_confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            ('headers', self.gf('django.db.models.fields.TextField')(default='[]')),
            ('correct_vals', self.gf('django.db.models.fields.CharField')(default='["3","4"]', max_length=100)),
            ('title', self.gf('django.db.models.fields.CharField')(default='untitled', max_length=100)),
            ('description', self.gf('django.db.models.fields.TextField')(default='')),
            ('tags', self.gf('django.db.models.fields.TextField')(default='')),
            ('iterations', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('filename', self.gf('django.db.models.fields.CharField')(default='/home/padler1/programming/drp/DRP/models/untitled.model', max_length=128)),
            ('active', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('start_time', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('end_time', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('usable', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('library', self.gf('django.db.models.fields.CharField')(default='weka', max_length=128)),
            ('tool', self.gf('django.db.models.fields.CharField')(default='svc', max_length=128)),
            ('response', self.gf('django.db.models.fields.CharField')(default='outcomoe', max_length=128)),
        ))
        db.send_create_signal('DRP', ['LegacyStatsModel'])

        # Adding model 'License'
        db.create_table(u'DRP_license', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('text', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal('DRP', ['License'])

        # Adding model 'LicenseAgreement'
        db.create_table(u'DRP_licenseagreement', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('text', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.License'])),
            ('signedDateTime', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal('DRP', ['LicenseAgreement'])

        # Adding M2M table for field users on 'LicenseAgreement'
        db.create_table(u'DRP_licenseagreement_users', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('licenseagreement', models.ForeignKey(orm['DRP.licenseagreement'], null=False)),
            ('user', models.ForeignKey(orm[u'auth.user'], null=False))
        ))
        db.create_unique(u'DRP_licenseagreement_users', ['licenseagreement_id', 'user_id'])


    def backwards(self, orm):
        # Deleting model 'Descriptor'
        db.delete_table(u'DRP_descriptor')

        # Deleting model 'ChemicalClass'
        db.delete_table(u'DRP_chemicalclass')

        # Deleting model 'LabGroup'
        db.delete_table(u'DRP_labgroup')

        # Removing M2M table for field users on 'LabGroup'
        db.delete_table('DRP_labgroup_users')

        # Deleting model 'CG_calculations'
        db.delete_table(u'DRP_cg_calculations')

        # Deleting model 'Compound'
        db.delete_table(u'DRP_compound')

        # Removing M2M table for field ChemicalClass on 'Compound'
        db.delete_table('DRP_compound_ChemicalClass')

        # Deleting model 'Reaction'
        db.delete_table(u'DRP_reaction')

        # Deleting model 'CompoundQuantity'
        db.delete_table(u'DRP_compoundquantity')

        # Deleting model 'StatsModelTag'
        db.delete_table(u'DRP_statsmodeltag')

        # Deleting model 'StatsModel'
        db.delete_table(u'DRP_statsmodel')

        # Removing M2M table for field tags on 'StatsModel'
        db.delete_table('DRP_statsmodel_tags')

        # Deleting model 'RecommendedReaction'
        db.delete_table(u'DRP_recommendedreaction')

        # Deleting model 'PerformedReaction'
        db.delete_table(u'DRP_performedreaction')

        # Deleting model 'DataSet'
        db.delete_table(u'DRP_dataset')

        # Deleting model 'DescriptorValue'
        db.delete_table(u'DRP_descriptorvalue')

        # Deleting model 'LegacyStatsModel'
        db.delete_table(u'DRP_legacystatsmodel')

        # Deleting model 'License'
        db.delete_table(u'DRP_license')

        # Deleting model 'LicenseAgreement'
        db.delete_table(u'DRP_licenseagreement')

        # Removing M2M table for field users on 'LicenseAgreement'
        db.delete_table('DRP_licenseagreement_users')


    models = {
        'DRP.cg_calculations': {
            'Meta': {'object_name': 'CG_calculations'},
            'compound': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json': ('django.db.models.fields.TextField', [], {'default': "'{}'", 'null': 'True'}),
            'json_data': ('django.db.models.fields.TextField', [], {}),
            'smiles': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
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
            'descriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.Descriptor']", 'through': "orm['DRP.DescriptorValue']", 'symmetrical': 'False'}),
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
        'DRP.dataset': {
            'Meta': {'object_name': 'DataSet'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isTestSet': ('django.db.models.fields.BooleanField', [], {}),
            'isTrainingSet': ('django.db.models.fields.BooleanField', [], {}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']"}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.PerformedReaction']"})
        },
        'DRP.descriptor': {
            'Meta': {'object_name': 'Descriptor'},
            'heading': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kind': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '300'})
        },
        'DRP.descriptorvalue': {
            'Compound': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Compound']", 'null': 'True'}),
            'Meta': {'object_name': 'DescriptorValue'},
            'Reaction': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Reaction']", 'null': 'True'}),
            'booleanValue': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'catValue': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True'}),
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Descriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isPredicted': ('django.db.models.fields.BooleanField', [], {}),
            'model': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            'numValue': ('django.db.models.fields.FloatField', [], {'null': 'True'}),
            'ordValue': ('django.db.models.fields.PositiveIntegerField', [], {'null': 'True'})
        },
        'DRP.labgroup': {
            'Meta': {'object_name': 'LabGroup'},
            'access_code': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'address': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'email': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '254'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'legacy_access_code': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'title': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.User']", 'symmetrical': 'False'})
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
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'text': ('django.db.models.fields.TextField', [], {})
        },
        'DRP.licenseagreement': {
            'Meta': {'object_name': 'LicenseAgreement'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'signedDateTime': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'text': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.License']"}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.User']", 'symmetrical': 'False'})
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
            'descriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.Descriptor']", 'through': "orm['DRP.DescriptorValue']", 'symmetrical': 'False'}),
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