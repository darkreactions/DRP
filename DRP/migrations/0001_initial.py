# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Lab_Group'
        db.create_table(u'DRP_lab_group', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('lab_title', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
            ('lab_address', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('lab_email', self.gf('django.db.models.fields.CharField')(max_length=254)),
            ('access_code', self.gf('django.db.models.fields.CharField')(default='IcHcOYhg9P2l5zDlK6Yl', max_length=20)),
        ))
        db.send_create_signal('DRP', ['Lab_Group'])

        # Adding model 'Lab_Member'
        db.create_table(u'DRP_lab_member', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(related_name='profile', unique=True, to=orm['auth.User'])),
            ('license_agreement_date_dt', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('lab_group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Lab_Group'])),
        ))
        db.send_create_signal('DRP', ['Lab_Member'])

        # Adding model 'ModelStats'
        db.create_table(u'DRP_modelstats', (
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
            ('response', self.gf('django.db.models.fields.CharField')(default='outcome', max_length=128)),
        ))
        db.send_create_signal('DRP', ['ModelStats'])

        # Adding model 'DataCalc'
        db.create_table(u'DRP_datacalc', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('contents', self.gf('django.db.models.fields.TextField')(default='{}')),
        ))
        db.send_create_signal('DRP', ['DataCalc'])

        # Adding model 'CG_calculations'
        db.create_table(u'DRP_cg_calculations', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('json_data', self.gf('django.db.models.fields.TextField')()),
            ('compound', self.gf('django.db.models.fields.CharField')(max_length=200)),
            ('smiles', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('json', self.gf('django.db.models.fields.TextField')(default='{}', null=True)),
        ))
        db.send_create_signal('DRP', ['CG_calculations'])

        # Adding model 'CompoundEntry'
        db.create_table(u'DRP_compoundentry', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('abbrev', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('compound', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('CAS_ID', self.gf('django.db.models.fields.CharField')(default='', max_length=13, blank=True)),
            ('compound_type', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('image_url', self.gf('django.db.models.fields.CharField')(default='', max_length=100, blank=True)),
            ('smiles', self.gf('django.db.models.fields.CharField')(default='', max_length=255, blank=True)),
            ('mw', self.gf('django.db.models.fields.CharField')(default='', max_length=20)),
            ('custom', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('lab_group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Lab_Group'])),
            ('calculations', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.CG_calculations'], null=True)),
            ('calculations_failed', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal('DRP', ['CompoundEntry'])

        # Adding model 'Data'
        db.create_table(u'DRP_data', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ref', self.gf('django.db.models.fields.CharField')(unique=True, max_length=30)),
            ('reactant_fk_1', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant_key_1', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_1', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_1', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_2', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant_key_2', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_2', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_2', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_3', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant_key_3', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_3', self.gf('django.db.models.fields.CharField')(max_length=10, blank=True)),
            ('unit_3', self.gf('django.db.models.fields.CharField')(max_length=4, blank=True)),
            ('reactant_fk_4', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant_key_4', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_4', self.gf('django.db.models.fields.CharField')(max_length=10, blank=True)),
            ('unit_4', self.gf('django.db.models.fields.CharField')(max_length=4, blank=True)),
            ('reactant_fk_5', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant_key_5', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_5', self.gf('django.db.models.fields.CharField')(max_length=10, blank=True)),
            ('unit_5', self.gf('django.db.models.fields.CharField')(max_length=4, blank=True)),
            ('temp', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('time', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('pH', self.gf('django.db.models.fields.CharField')(max_length=5)),
            ('slow_cool', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('leak', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('outcome', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('purity', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('notes', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('calculations', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.DataCalc'], null=True, on_delete=models.SET_NULL, blank=True)),
            ('calculated_pH', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('calculated_temp', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('calculated_time', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('atoms', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('lab_group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Lab_Group'])),
            ('creation_time_dt', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('is_valid', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('public', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('duplicate_of', self.gf('django.db.models.fields.CharField')(max_length=12, null=True, blank=True)),
            ('recommended', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('persistent_homologies', self.gf('django.db.models.fields.TextField')(default='[]', blank=True)),
        ))
        db.send_create_signal('DRP', ['Data'])

        # Adding model 'Recommendation'
        db.create_table(u'DRP_recommendation', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('ref', self.gf('django.db.models.fields.CharField')(default='', max_length=30)),
            ('reactant_fk_1', self.gf('django.db.models.fields.related.ForeignKey')(related_name='rec_reactant_key_1', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_1', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_1', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_2', self.gf('django.db.models.fields.related.ForeignKey')(related_name='rec_reactant_key_2', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_2', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_2', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_3', self.gf('django.db.models.fields.related.ForeignKey')(related_name='rec_reactant_key_3', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_3', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_3', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_4', self.gf('django.db.models.fields.related.ForeignKey')(related_name='rec_reactant_key_4', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_4', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_4', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('reactant_fk_5', self.gf('django.db.models.fields.related.ForeignKey')(related_name='rec_reactant_key_5', default=None, to=orm['DRP.CompoundEntry'], max_length=30, blank=True, null=True)),
            ('quantity_5', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('unit_5', self.gf('django.db.models.fields.CharField')(max_length=4)),
            ('score', self.gf('django.db.models.fields.FloatField')()),
            ('temp', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('time', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('pH', self.gf('django.db.models.fields.CharField')(max_length=5)),
            ('slow_cool', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('leak', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('outcome', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('purity', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('atoms', self.gf('django.db.models.fields.CharField')(max_length=30, blank=True)),
            ('lab_group', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.Lab_Group'])),
            ('model_version', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.ModelStats'])),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(default=None, related_name='last_user', null=True, blank=True, to=orm['auth.User'])),
            ('assigned_user', self.gf('django.db.models.fields.related.ForeignKey')(default=None, related_name='assigned_user', null=True, blank=True, to=orm['auth.User'])),
            ('seed', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.Data'], null=True, blank=True)),
            ('date_dt', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('complete', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('seeded', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('saved', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('nonsense', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('hidden', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('notes', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('calculations', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.DataCalc'], null=True, on_delete=models.SET_NULL, blank=True)),
        ))
        db.send_create_signal('DRP', ['Recommendation'])

        # Adding model 'RankedReactionList'
        db.create_table(u'DRP_rankedreactionlist', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('original_list', self.gf('django.db.models.fields.TextField')()),
            ('seed', self.gf('django.db.models.fields.TextField')()),
            ('ranked_list', self.gf('django.db.models.fields.TextField')()),
            ('ranker', self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['auth.User'], null=True)),
        ))
        db.send_create_signal('DRP', ['RankedReactionList'])


    def backwards(self, orm):
        # Deleting model 'Lab_Group'
        db.delete_table(u'DRP_lab_group')

        # Deleting model 'Lab_Member'
        db.delete_table(u'DRP_lab_member')

        # Deleting model 'ModelStats'
        db.delete_table(u'DRP_modelstats')

        # Deleting model 'DataCalc'
        db.delete_table(u'DRP_datacalc')

        # Deleting model 'CG_calculations'
        db.delete_table(u'DRP_cg_calculations')

        # Deleting model 'CompoundEntry'
        db.delete_table(u'DRP_compoundentry')

        # Deleting model 'Data'
        db.delete_table(u'DRP_data')

        # Deleting model 'Recommendation'
        db.delete_table(u'DRP_recommendation')

        # Deleting model 'RankedReactionList'
        db.delete_table(u'DRP_rankedreactionlist')


    models = {
        'DRP.cg_calculations': {
            'Meta': {'object_name': 'CG_calculations'},
            'compound': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json': ('django.db.models.fields.TextField', [], {'default': "'{}'", 'null': 'True'}),
            'json_data': ('django.db.models.fields.TextField', [], {}),
            'smiles': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'})
        },
        'DRP.compoundentry': {
            'CAS_ID': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '13', 'blank': 'True'}),
            'Meta': {'object_name': 'CompoundEntry'},
            'abbrev': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'calculations': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.CG_calculations']", 'null': 'True'}),
            'calculations_failed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'compound': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'compound_type': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'custom': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'image_url': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '100', 'blank': 'True'}),
            'lab_group': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Lab_Group']"}),
            'mw': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '20'}),
            'smiles': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '255', 'blank': 'True'})
        },
        'DRP.data': {
            'Meta': {'object_name': 'Data'},
            'atoms': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'calculated_pH': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'calculated_temp': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'calculated_time': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'calculations': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.DataCalc']", 'null': 'True', 'on_delete': 'models.SET_NULL', 'blank': 'True'}),
            'creation_time_dt': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'duplicate_of': ('django.db.models.fields.CharField', [], {'max_length': '12', 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_valid': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'lab_group': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Lab_Group']"}),
            'leak': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'notes': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'outcome': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'pH': ('django.db.models.fields.CharField', [], {'max_length': '5'}),
            'persistent_homologies': ('django.db.models.fields.TextField', [], {'default': "'[]'", 'blank': 'True'}),
            'public': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'purity': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'quantity_1': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_2': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_3': ('django.db.models.fields.CharField', [], {'max_length': '10', 'blank': 'True'}),
            'quantity_4': ('django.db.models.fields.CharField', [], {'max_length': '10', 'blank': 'True'}),
            'quantity_5': ('django.db.models.fields.CharField', [], {'max_length': '10', 'blank': 'True'}),
            'reactant_fk_1': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant_key_1'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_2': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant_key_2'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_3': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant_key_3'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_4': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant_key_4'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_5': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant_key_5'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'recommended': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'ref': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'}),
            'slow_cool': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'temp': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'time': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'unit_1': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_2': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_3': ('django.db.models.fields.CharField', [], {'max_length': '4', 'blank': 'True'}),
            'unit_4': ('django.db.models.fields.CharField', [], {'max_length': '4', 'blank': 'True'}),
            'unit_5': ('django.db.models.fields.CharField', [], {'max_length': '4', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        'DRP.datacalc': {
            'Meta': {'object_name': 'DataCalc'},
            'contents': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        'DRP.lab_group': {
            'Meta': {'object_name': 'Lab_Group'},
            'access_code': ('django.db.models.fields.CharField', [], {'default': "'wlLMosRQ0L6QFqeblZwB'", 'max_length': '20'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lab_address': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'lab_email': ('django.db.models.fields.CharField', [], {'max_length': '254'}),
            'lab_title': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'})
        },
        'DRP.lab_member': {
            'Meta': {'object_name': 'Lab_Member'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lab_group': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Lab_Group']"}),
            'license_agreement_date_dt': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'related_name': "'profile'", 'unique': 'True', 'to': u"orm['auth.User']"})
        },
        'DRP.modelstats': {
            'Meta': {'object_name': 'ModelStats'},
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
            'response': ('django.db.models.fields.CharField', [], {'default': "'outcome'", 'max_length': '128'}),
            'start_time': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'tags': ('django.db.models.fields.TextField', [], {'default': "''"}),
            'title': ('django.db.models.fields.CharField', [], {'default': "'untitled'", 'max_length': '100'}),
            'tmp_confusion_table': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            'tool': ('django.db.models.fields.CharField', [], {'default': "'svc'", 'max_length': '128'}),
            'train_confusion_table': ('django.db.models.fields.TextField', [], {'default': "'{}'"}),
            'usable': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'DRP.rankedreactionlist': {
            'Meta': {'object_name': 'RankedReactionList'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'original_list': ('django.db.models.fields.TextField', [], {}),
            'ranked_list': ('django.db.models.fields.TextField', [], {}),
            'ranker': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': u"orm['auth.User']", 'null': 'True'}),
            'seed': ('django.db.models.fields.TextField', [], {})
        },
        'DRP.recommendation': {
            'Meta': {'object_name': 'Recommendation'},
            'assigned_user': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'assigned_user'", 'null': 'True', 'blank': 'True', 'to': u"orm['auth.User']"}),
            'atoms': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'calculations': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.DataCalc']", 'null': 'True', 'on_delete': 'models.SET_NULL', 'blank': 'True'}),
            'complete': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'date_dt': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'hidden': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lab_group': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Lab_Group']"}),
            'leak': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'model_version': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelStats']"}),
            'nonsense': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'notes': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'outcome': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'pH': ('django.db.models.fields.CharField', [], {'max_length': '5'}),
            'purity': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'quantity_1': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_2': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_3': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_4': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'quantity_5': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'reactant_fk_1': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'rec_reactant_key_1'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_2': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'rec_reactant_key_2'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_3': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'rec_reactant_key_3'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_4': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'rec_reactant_key_4'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'reactant_fk_5': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'rec_reactant_key_5'", 'default': 'None', 'to': "orm['DRP.CompoundEntry']", 'max_length': '30', 'blank': 'True', 'null': 'True'}),
            'ref': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '30'}),
            'saved': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'score': ('django.db.models.fields.FloatField', [], {}),
            'seed': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.Data']", 'null': 'True', 'blank': 'True'}),
            'seeded': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'slow_cool': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'temp': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'time': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'unit_1': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_2': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_3': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_4': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'unit_5': ('django.db.models.fields.CharField', [], {'max_length': '4'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'last_user'", 'null': 'True', 'blank': 'True', 'to': u"orm['auth.User']"})
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