# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting model 'StatsModelTag'
        db.delete_table(u'DRP_statsmodeltag')

        # Deleting model 'LegacyStatsModel'
        db.delete_table(u'DRP_legacystatsmodel')

        # Adding model 'MetricContainer'
        db.create_table(u'DRP_metriccontainer', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('metricVisitor', self.gf('django.db.models.fields.CharField')(max_length=255)),
        ))
        db.send_create_signal('DRP', ['MetricContainer'])

        # Adding M2M table for field boolRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_boolRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('boolrxndescriptor', models.ForeignKey(orm['DRP.boolrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_boolRxnDescriptors', ['metriccontainer_id', 'boolrxndescriptor_id'])

        # Adding M2M table for field ordRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_ordRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('ordrxndescriptor', models.ForeignKey(orm['DRP.ordrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_ordRxnDescriptors', ['metriccontainer_id', 'ordrxndescriptor_id'])

        # Adding M2M table for field catRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_catRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('catrxndescriptor', models.ForeignKey(orm['DRP.catrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_catRxnDescriptors', ['metriccontainer_id', 'catrxndescriptor_id'])

        # Adding M2M table for field numRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_numRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('numrxndescriptor', models.ForeignKey(orm['DRP.numrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_numRxnDescriptors', ['metriccontainer_id', 'numrxndescriptor_id'])

        # Adding M2M table for field outcomeBoolRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_outcomeBoolRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('boolrxndescriptor', models.ForeignKey(orm['DRP.boolrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_outcomeBoolRxnDescriptors', ['metriccontainer_id', 'boolrxndescriptor_id'])

        # Adding M2M table for field outcomeOrdRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_outcomeOrdRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('ordrxndescriptor', models.ForeignKey(orm['DRP.ordrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_outcomeOrdRxnDescriptors', ['metriccontainer_id', 'ordrxndescriptor_id'])

        # Adding M2M table for field outcomeCatRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_outcomeCatRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('catrxndescriptor', models.ForeignKey(orm['DRP.catrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_outcomeCatRxnDescriptors', ['metriccontainer_id', 'catrxndescriptor_id'])

        # Adding M2M table for field outcomeNumRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_outcomeNumRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('numrxndescriptor', models.ForeignKey(orm['DRP.numrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_outcomeNumRxnDescriptors', ['metriccontainer_id', 'numrxndescriptor_id'])

        # Adding M2M table for field transformedRxnDescriptors on 'MetricContainer'
        db.create_table(u'DRP_metriccontainer_transformedRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('metriccontainer', models.ForeignKey(orm['DRP.metriccontainer'], null=False)),
            ('numrxndescriptor', models.ForeignKey(orm['DRP.numrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_metriccontainer_transformedRxnDescriptors', ['metriccontainer_id', 'numrxndescriptor_id'])


    def backwards(self, orm):
        # Adding model 'StatsModelTag'
        db.create_table(u'DRP_statsmodeltag', (
            ('text', self.gf('django.db.models.fields.CharField')(max_length=200, unique=True)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['StatsModelTag'])

        # Adding model 'LegacyStatsModel'
        db.create_table(u'DRP_legacystatsmodel', (
            ('usable', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('description', self.gf('django.db.models.fields.TextField')(default='')),
            ('tags', self.gf('django.db.models.fields.TextField')(default='')),
            ('tmp_confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            ('start_time', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('title', self.gf('django.db.models.fields.CharField')(default='untitled', max_length=100)),
            ('tool', self.gf('django.db.models.fields.CharField')(default='svc', max_length=128)),
            ('correct_vals', self.gf('django.db.models.fields.CharField')(default='["3","4"]', max_length=100)),
            ('library', self.gf('django.db.models.fields.CharField')(default='weka', max_length=128)),
            ('filename', self.gf('django.db.models.fields.CharField')(default='/home/padler1/programming/drp/DRP/modelsuntitled.model', max_length=128)),
            ('headers', self.gf('django.db.models.fields.TextField')(default='[]')),
            ('confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            ('end_time', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('iterations', self.gf('django.db.models.fields.IntegerField')(default=1)),
            ('active', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('response', self.gf('django.db.models.fields.CharField')(default='outcomoe', max_length=128)),
            ('train_confusion_table', self.gf('django.db.models.fields.TextField')(default='{}')),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['LegacyStatsModel'])

        # Deleting model 'MetricContainer'
        db.delete_table(u'DRP_metriccontainer')

        # Removing M2M table for field boolRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_boolRxnDescriptors')

        # Removing M2M table for field ordRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_ordRxnDescriptors')

        # Removing M2M table for field catRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_catRxnDescriptors')

        # Removing M2M table for field numRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_numRxnDescriptors')

        # Removing M2M table for field outcomeBoolRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_outcomeBoolRxnDescriptors')

        # Removing M2M table for field outcomeOrdRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_outcomeOrdRxnDescriptors')

        # Removing M2M table for field outcomeCatRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_outcomeCatRxnDescriptors')

        # Removing M2M table for field outcomeNumRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_outcomeNumRxnDescriptors')

        # Removing M2M table for field transformedRxnDescriptors on 'MetricContainer'
        db.delete_table('DRP_metriccontainer_transformedRxnDescriptors')


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
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"}),
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
            'value': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptorPermittedValue']", 'null': 'True', 'on_delete': 'models.PROTECT'})
        },
        'DRP.catrxndescriptor': {
            'Meta': {'object_name': 'CatRxnDescriptor', '_ormbases': ['DRP.CategoricalDescriptor']},
            u'categoricaldescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.CategoricalDescriptor']", 'unique': 'True', 'primary_key': 'True'})
        },
        'DRP.catrxndescriptorvalue': {
            'Meta': {'unique_together': "(('reaction', 'descriptor'),)", 'object_name': 'CatRxnDescriptorValue'},
            'descriptor': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptor']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"}),
            'value': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.CategoricalDescriptorPermittedValue']", 'null': 'True', 'on_delete': 'models.PROTECT'})
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
            'formula': ('django.db.models.fields.CharField', [], {'max_length': '500', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'labGroup': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.LabGroup']"}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '400'}),
            'smiles': ('django.db.models.fields.TextField', [], {'default': "''", 'blank': 'True'})
        },
        'DRP.compoundquantity': {
            'Meta': {'unique_together': "(('reaction', 'role', 'amount'),)", 'object_name': 'CompoundQuantity'},
            'amount': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
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
        'DRP.dataset': {
            'Meta': {'object_name': 'DataSet'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '200'}),
            'reactions': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.PerformedReaction']", 'through': "orm['DRP.DataSetRelation']", 'symmetrical': 'False'})
        },
        'DRP.datasetrelation': {
            'Meta': {'unique_together': "(('dataSet', 'reaction'),)", 'object_name': 'DataSetRelation'},
            'dataSet': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.DataSet']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.PerformedReaction']", 'on_delete': 'models.PROTECT'})
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
        'DRP.metriccontainer': {
            'Meta': {'object_name': 'MetricContainer'},
            'boolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.BoolRxnDescriptor']", 'symmetrical': 'False'}),
            'catRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.CatRxnDescriptor']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'metricVisitor': ('django.db.models.fields.CharField', [], {'max_length': '255'}),
            'numRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.NumRxnDescriptor']", 'symmetrical': 'False'}),
            'ordRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.OrdRxnDescriptor']", 'symmetrical': 'False'}),
            'outcomeBoolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForMetrics'", 'symmetrical': 'False', 'to': "orm['DRP.BoolRxnDescriptor']"}),
            'outcomeCatRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForMetrics'", 'symmetrical': 'False', 'to': "orm['DRP.CatRxnDescriptor']"}),
            'outcomeNumRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForMetrics'", 'symmetrical': 'False', 'to': "orm['DRP.NumRxnDescriptor']"}),
            'outcomeOrdRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForMetrics'", 'symmetrical': 'False', 'to': "orm['DRP.OrdRxnDescriptor']"}),
            'transformedRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'transformedByMetric'", 'symmetrical': 'False', 'to': "orm['DRP.NumRxnDescriptor']"})
        },
        'DRP.modelcontainer': {
            'Meta': {'object_name': 'ModelContainer'},
            'active': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'boolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.BoolRxnDescriptor']", 'symmetrical': 'False'}),
            'built': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'catRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.CatRxnDescriptor']", 'symmetrical': 'False'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'featureLibrary': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '200'}),
            'featureTool': ('django.db.models.fields.CharField', [], {'default': "''", 'max_length': '200'}),
            'fully_trained': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'modelVisitorLibrary': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'modelVisitorTool': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'numRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.NumRxnDescriptor']", 'symmetrical': 'False'}),
            'ordRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.OrdRxnDescriptor']", 'symmetrical': 'False'}),
            'outcomeBoolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.BoolRxnDescriptor']"}),
            'outcomeCatRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.CatRxnDescriptor']"}),
            'outcomeNumRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.NumRxnDescriptor']"}),
            'outcomeOrdRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.OrdRxnDescriptor']"}),
            'splitter': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'})
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
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"}),
            'value': ('django.db.models.fields.FloatField', [], {'null': 'True'})
        },
        'DRP.ordinaldescriptor': {
            'Meta': {'object_name': 'OrdinalDescriptor', '_ormbases': ['DRP.Descriptor']},
            u'descriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Descriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'maximum': ('django.db.models.fields.IntegerField', [], {}),
            'minimum': ('django.db.models.fields.IntegerField', [], {})
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
            'reaction': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.Reaction']"}),
            'value': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        'DRP.performedreaction': {
            'Meta': {'object_name': 'PerformedReaction', '_ormbases': ['DRP.Reaction']},
            'duplicateOf': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'duplicatedBy'", 'null': 'True', 'blank': 'True', 'to': "orm['DRP.PerformedReaction']"}),
            'insertedDateTime': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'legacyRecommendedFlag': ('django.db.models.fields.NullBooleanField', [], {'default': 'None', 'null': 'True', 'blank': 'True'}),
            'performedBy': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'performedReactions'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'performedDateTime': ('django.db.models.fields.DateTimeField', [], {'default': 'None', 'null': 'True'}),
            'public': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'recommendation': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'related_name': "'resultantExperiment'", 'null': 'True', 'blank': 'True', 'to': "orm['DRP.RecommendedReaction']"}),
            'reference': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '40'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"}),
            'valid': ('django.db.models.fields.BooleanField', [], {'default': 'True'})
        },
        'DRP.predboolrxndescriptor': {
            'Meta': {'object_name': 'PredBoolRxnDescriptor', '_ormbases': ['DRP.BoolRxnDescriptor']},
            u'boolrxndescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.BoolRxnDescriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'modelContainer': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelContainer']"}),
            'predictionOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'prediction_of'", 'to': "orm['DRP.BoolRxnDescriptor']"}),
            'statsModel': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'})
        },
        'DRP.predcatrxndescriptor': {
            'Meta': {'object_name': 'PredCatRxnDescriptor', '_ormbases': ['DRP.CatRxnDescriptor']},
            u'catrxndescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.CatRxnDescriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'modelContainer': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelContainer']"}),
            'predictionOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'prediction_of'", 'to': "orm['DRP.CatRxnDescriptor']"}),
            'statsModel': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'})
        },
        'DRP.prednumrxndescriptor': {
            'Meta': {'object_name': 'PredNumRxnDescriptor', '_ormbases': ['DRP.NumRxnDescriptor']},
            'modelContainer': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelContainer']"}),
            u'numrxndescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.NumRxnDescriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'predictionOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'prediction_of'", 'to': "orm['DRP.NumRxnDescriptor']"}),
            'statsModel': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'})
        },
        'DRP.predordrxndescriptor': {
            'Meta': {'object_name': 'PredOrdRxnDescriptor', '_ormbases': ['DRP.OrdRxnDescriptor']},
            'modelContainer': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelContainer']"}),
            u'ordrxndescriptor_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.OrdRxnDescriptor']", 'unique': 'True', 'primary_key': 'True'}),
            'predictionOf': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'predition_of'", 'to': "orm['DRP.OrdRxnDescriptor']"}),
            'statsModel': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'})
        },
        'DRP.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'compounds': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.Compound']", 'through': "orm['DRP.CompoundQuantity']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'labGroup': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.LabGroup']"}),
            'notes': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        'DRP.recommendedreaction': {
            'Meta': {'object_name': 'RecommendedReaction', '_ormbases': ['DRP.Reaction']},
            'hidden': ('django.db.models.fields.BooleanField', [], {}),
            'nonsense': ('django.db.models.fields.BooleanField', [], {}),
            u'reaction_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['DRP.Reaction']", 'unique': 'True', 'primary_key': 'True'}),
            'reference': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'saved': ('django.db.models.fields.BooleanField', [], {}),
            'score': ('django.db.models.fields.FloatField', [], {}),
            'seed': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'seeded'", 'null': 'True', 'to': "orm['DRP.Reaction']"})
        },
        'DRP.statsmodel': {
            'Meta': {'object_name': 'StatsModel'},
            'container': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.ModelContainer']"}),
            'endTime': ('django.db.models.fields.DateTimeField', [], {'default': 'None', 'null': 'True'}),
            'fileName': ('django.db.models.fields.files.FileField', [], {'max_length': '200', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'invalid': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'regenerationOf': ('django.db.models.fields.related.ForeignKey', [], {'default': 'None', 'to': "orm['DRP.StatsModel']", 'null': 'True', 'blank': 'True'}),
            'startTime': ('django.db.models.fields.DateTimeField', [], {'default': 'None', 'null': 'True'}),
            'testSets': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'testSetsFor'", 'symmetrical': 'False', 'to': "orm['DRP.DataSet']"}),
            'trainingSet': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'trainingSetFor'", 'to': "orm['DRP.DataSet']"})
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