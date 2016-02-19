# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Removing unique constraint on 'TestSetRelation', fields ['test_set', 'reaction']
        db.delete_unique(u'DRP_testsetrelation', ['test_set_id', 'reaction_id'])

        # Deleting model 'TestSetRelation'
        db.delete_table(u'DRP_testsetrelation')

        # Deleting model 'TrainingSet'
        db.delete_table(u'DRP_trainingset')

        # Deleting model 'TestSet'
        db.delete_table(u'DRP_testset')

        # Adding model 'DataSet'
        db.create_table(u'DRP_dataset', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=200)),
        ))
        db.send_create_signal('DRP', ['DataSet'])

        # Adding model 'DataSetRelation'
        db.create_table(u'DRP_datasetrelation', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.PerformedReaction'], on_delete=models.PROTECT)),
            ('dataSet', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.DataSet'])),
        ))
        db.send_create_signal('DRP', ['DataSetRelation'])

        # Adding unique constraint on 'DataSetRelation', fields ['dataSet', 'reaction']
        db.create_unique(u'DRP_datasetrelation', ['dataSet_id', 'reaction_id'])

        # Deleting field 'PredOrdRxnDescriptor.prediction_of'
        db.delete_column(u'DRP_predordrxndescriptor', 'prediction_of_id')

        # Deleting field 'PredOrdRxnDescriptor.stats_model'
        db.delete_column(u'DRP_predordrxndescriptor', 'stats_model_id')

        # Adding field 'PredOrdRxnDescriptor.modelContainer'
        db.add_column(u'DRP_predordrxndescriptor', 'modelContainer',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['DRP.ModelContainer']),
                      keep_default=False)

        # Adding field 'PredOrdRxnDescriptor.statsModel'
        db.add_column(u'DRP_predordrxndescriptor', 'statsModel',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'PredOrdRxnDescriptor.predictionOf'
        db.add_column(u'DRP_predordrxndescriptor', 'predictionOf',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, related_name='predition_of', to=orm['DRP.OrdRxnDescriptor']),
                      keep_default=False)

        # Deleting field 'PredNumRxnDescriptor.prediction_of'
        db.delete_column(u'DRP_prednumrxndescriptor', 'prediction_of_id')

        # Deleting field 'PredNumRxnDescriptor.stats_model'
        db.delete_column(u'DRP_prednumrxndescriptor', 'stats_model_id')

        # Adding field 'PredNumRxnDescriptor.modelContainer'
        db.add_column(u'DRP_prednumrxndescriptor', 'modelContainer',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['DRP.ModelContainer']),
                      keep_default=False)

        # Adding field 'PredNumRxnDescriptor.statsModel'
        db.add_column(u'DRP_prednumrxndescriptor', 'statsModel',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'PredNumRxnDescriptor.predictionOf'
        db.add_column(u'DRP_prednumrxndescriptor', 'predictionOf',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, related_name='prediction_of', to=orm['DRP.NumRxnDescriptor']),
                      keep_default=False)

        # Deleting field 'PredCatRxnDescriptor.prediction_of'
        db.delete_column(u'DRP_predcatrxndescriptor', 'prediction_of_id')

        # Deleting field 'PredCatRxnDescriptor.stats_model'
        db.delete_column(u'DRP_predcatrxndescriptor', 'stats_model_id')

        # Adding field 'PredCatRxnDescriptor.modelContainer'
        db.add_column(u'DRP_predcatrxndescriptor', 'modelContainer',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['DRP.ModelContainer']),
                      keep_default=False)

        # Adding field 'PredCatRxnDescriptor.statsModel'
        db.add_column(u'DRP_predcatrxndescriptor', 'statsModel',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'PredCatRxnDescriptor.predictionOf'
        db.add_column(u'DRP_predcatrxndescriptor', 'predictionOf',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, related_name='prediction_of', to=orm['DRP.CatRxnDescriptor']),
                      keep_default=False)

        # Adding field 'ModelContainer.description'
        db.add_column(u'DRP_modelcontainer', 'description',
                      self.gf('django.db.models.fields.TextField')(default=1),
                      keep_default=False)

        # Adding field 'ModelContainer.active'
        db.add_column(u'DRP_modelcontainer', 'active',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)

        # Adding field 'ModelContainer.built'
        db.add_column(u'DRP_modelcontainer', 'built',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)

        # Adding M2M table for field boolRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_boolRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('boolrxndescriptor', models.ForeignKey(orm['DRP.boolrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_boolRxnDescriptors', ['modelcontainer_id', 'boolrxndescriptor_id'])

        # Adding M2M table for field ordRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_ordRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('ordrxndescriptor', models.ForeignKey(orm['DRP.ordrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_ordRxnDescriptors', ['modelcontainer_id', 'ordrxndescriptor_id'])

        # Adding M2M table for field catRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_catRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('catrxndescriptor', models.ForeignKey(orm['DRP.catrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_catRxnDescriptors', ['modelcontainer_id', 'catrxndescriptor_id'])

        # Adding M2M table for field numRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_numRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('numrxndescriptor', models.ForeignKey(orm['DRP.numrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_numRxnDescriptors', ['modelcontainer_id', 'numrxndescriptor_id'])

        # Adding M2M table for field outcomeBoolRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_outcomeBoolRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('boolrxndescriptor', models.ForeignKey(orm['DRP.boolrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_outcomeBoolRxnDescriptors', ['modelcontainer_id', 'boolrxndescriptor_id'])

        # Adding M2M table for field outcomeOrdRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_outcomeOrdRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('ordrxndescriptor', models.ForeignKey(orm['DRP.ordrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_outcomeOrdRxnDescriptors', ['modelcontainer_id', 'ordrxndescriptor_id'])

        # Adding M2M table for field outcomeCatRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_outcomeCatRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('catrxndescriptor', models.ForeignKey(orm['DRP.catrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_outcomeCatRxnDescriptors', ['modelcontainer_id', 'catrxndescriptor_id'])

        # Adding M2M table for field outcomeNumRxnDescriptors on 'ModelContainer'
        db.create_table(u'DRP_modelcontainer_outcomeNumRxnDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('modelcontainer', models.ForeignKey(orm['DRP.modelcontainer'], null=False)),
            ('numrxndescriptor', models.ForeignKey(orm['DRP.numrxndescriptor'], null=False))
        ))
        db.create_unique(u'DRP_modelcontainer_outcomeNumRxnDescriptors', ['modelcontainer_id', 'numrxndescriptor_id'])


        # Changing field 'ModelContainer.splitter'
        db.alter_column(u'DRP_modelcontainer', 'splitter', self.gf('django.db.models.fields.CharField')(max_length=200, null=True))
        # Deleting field 'NumRxnDescriptorValue.model'
        db.delete_column(u'DRP_numrxndescriptorvalue', 'model_id')

        # Deleting field 'BoolRxnDescriptorValue.model'
        db.delete_column(u'DRP_boolrxndescriptorvalue', 'model_id')

        # Deleting field 'PredBoolRxnDescriptor.prediction_of'
        db.delete_column(u'DRP_predboolrxndescriptor', 'prediction_of_id')

        # Deleting field 'PredBoolRxnDescriptor.stats_model'
        db.delete_column(u'DRP_predboolrxndescriptor', 'stats_model_id')

        # Adding field 'PredBoolRxnDescriptor.modelContainer'
        db.add_column(u'DRP_predboolrxndescriptor', 'modelContainer',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, to=orm['DRP.ModelContainer']),
                      keep_default=False)

        # Adding field 'PredBoolRxnDescriptor.statsModel'
        db.add_column(u'DRP_predboolrxndescriptor', 'statsModel',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'PredBoolRxnDescriptor.predictionOf'
        db.add_column(u'DRP_predboolrxndescriptor', 'predictionOf',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, related_name='prediction_of', to=orm['DRP.BoolRxnDescriptor']),
                      keep_default=False)

        # Deleting field 'OrdRxnDescriptorValue.model'
        db.delete_column(u'DRP_ordrxndescriptorvalue', 'model_id')

        # Deleting field 'CatRxnDescriptorValue.model'
        db.delete_column(u'DRP_catrxndescriptorvalue', 'model_id')

        # Deleting field 'StatsModel.description'
        db.delete_column(u'DRP_statsmodel', 'description')

        # Deleting field 'StatsModel.start_time'
        db.delete_column(u'DRP_statsmodel', 'start_time')

        # Deleting field 'StatsModel.end_time'
        db.delete_column(u'DRP_statsmodel', 'end_time')

        # Deleting field 'StatsModel.iterations'
        db.delete_column(u'DRP_statsmodel', 'iterations')

        # Deleting field 'StatsModel.active'
        db.delete_column(u'DRP_statsmodel', 'active')

        # Deleting field 'StatsModel.snapShot'
        db.delete_column(u'DRP_statsmodel', 'snapShot')

        # Adding field 'StatsModel.startTime'
        db.add_column(u'DRP_statsmodel', 'startTime',
                      self.gf('django.db.models.fields.DateTimeField')(default=None, null=True),
                      keep_default=False)

        # Adding field 'StatsModel.endTime'
        db.add_column(u'DRP_statsmodel', 'endTime',
                      self.gf('django.db.models.fields.DateTimeField')(default=None, null=True),
                      keep_default=False)

        # Adding field 'StatsModel.trainingSet'
        db.add_column(u'DRP_statsmodel', 'trainingSet',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=1, related_name='trainingSetFor', to=orm['DRP.DataSet']),
                      keep_default=False)

        # Removing M2M table for field tags on 'StatsModel'
        db.delete_table('DRP_statsmodel_tags')

        # Removing M2M table for field descriptors on 'StatsModel'
        db.delete_table('DRP_statsmodel_descriptors')

        # Removing M2M table for field outcomeDescriptors on 'StatsModel'
        db.delete_table('DRP_statsmodel_outcomeDescriptors')

        # Removing M2M table for field predictsDescriptors on 'StatsModel'
        db.delete_table('DRP_statsmodel_predictsDescriptors')

        # Adding M2M table for field testSets on 'StatsModel'
        db.create_table(u'DRP_statsmodel_testSets', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('dataset', models.ForeignKey(orm['DRP.dataset'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_testSets', ['statsmodel_id', 'dataset_id'])

        # Deleting field 'RecommendedReaction.model'
        db.delete_column(u'DRP_recommendedreaction', 'model_id')


    def backwards(self, orm):
        # Removing unique constraint on 'DataSetRelation', fields ['dataSet', 'reaction']
        db.delete_unique(u'DRP_datasetrelation', ['dataSet_id', 'reaction_id'])

        # Adding model 'TestSetRelation'
        db.create_table(u'DRP_testsetrelation', (
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.PerformedReaction'], on_delete=models.PROTECT)),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('test_set', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.TestSet'])),
        ))
        db.send_create_signal('DRP', ['TestSetRelation'])

        # Adding unique constraint on 'TestSetRelation', fields ['test_set', 'reaction']
        db.create_unique(u'DRP_testsetrelation', ['test_set_id', 'reaction_id'])

        # Adding model 'TrainingSet'
        db.create_table(u'DRP_trainingset', (
            ('reaction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.PerformedReaction'], on_delete=models.PROTECT)),
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'])),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal('DRP', ['TrainingSet'])

        # Adding model 'TestSet'
        db.create_table(u'DRP_testset', (
            ('model', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'])),
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, unique=True)),
        ))
        db.send_create_signal('DRP', ['TestSet'])

        # Deleting model 'DataSet'
        db.delete_table(u'DRP_dataset')

        # Deleting model 'DataSetRelation'
        db.delete_table(u'DRP_datasetrelation')


        # User chose to not deal with backwards NULL issues for 'PredOrdRxnDescriptor.prediction_of'
        raise RuntimeError("Cannot reverse this migration. 'PredOrdRxnDescriptor.prediction_of' and its values cannot be restored.")

        # User chose to not deal with backwards NULL issues for 'PredOrdRxnDescriptor.stats_model'
        raise RuntimeError("Cannot reverse this migration. 'PredOrdRxnDescriptor.stats_model' and its values cannot be restored.")
        # Deleting field 'PredOrdRxnDescriptor.modelContainer'
        db.delete_column(u'DRP_predordrxndescriptor', 'modelContainer_id')

        # Deleting field 'PredOrdRxnDescriptor.statsModel'
        db.delete_column(u'DRP_predordrxndescriptor', 'statsModel_id')

        # Deleting field 'PredOrdRxnDescriptor.predictionOf'
        db.delete_column(u'DRP_predordrxndescriptor', 'predictionOf_id')


        # User chose to not deal with backwards NULL issues for 'PredNumRxnDescriptor.prediction_of'
        raise RuntimeError("Cannot reverse this migration. 'PredNumRxnDescriptor.prediction_of' and its values cannot be restored.")

        # User chose to not deal with backwards NULL issues for 'PredNumRxnDescriptor.stats_model'
        raise RuntimeError("Cannot reverse this migration. 'PredNumRxnDescriptor.stats_model' and its values cannot be restored.")
        # Deleting field 'PredNumRxnDescriptor.modelContainer'
        db.delete_column(u'DRP_prednumrxndescriptor', 'modelContainer_id')

        # Deleting field 'PredNumRxnDescriptor.statsModel'
        db.delete_column(u'DRP_prednumrxndescriptor', 'statsModel_id')

        # Deleting field 'PredNumRxnDescriptor.predictionOf'
        db.delete_column(u'DRP_prednumrxndescriptor', 'predictionOf_id')


        # User chose to not deal with backwards NULL issues for 'PredCatRxnDescriptor.prediction_of'
        raise RuntimeError("Cannot reverse this migration. 'PredCatRxnDescriptor.prediction_of' and its values cannot be restored.")

        # User chose to not deal with backwards NULL issues for 'PredCatRxnDescriptor.stats_model'
        raise RuntimeError("Cannot reverse this migration. 'PredCatRxnDescriptor.stats_model' and its values cannot be restored.")
        # Deleting field 'PredCatRxnDescriptor.modelContainer'
        db.delete_column(u'DRP_predcatrxndescriptor', 'modelContainer_id')

        # Deleting field 'PredCatRxnDescriptor.statsModel'
        db.delete_column(u'DRP_predcatrxndescriptor', 'statsModel_id')

        # Deleting field 'PredCatRxnDescriptor.predictionOf'
        db.delete_column(u'DRP_predcatrxndescriptor', 'predictionOf_id')

        # Deleting field 'ModelContainer.description'
        db.delete_column(u'DRP_modelcontainer', 'description')

        # Deleting field 'ModelContainer.active'
        db.delete_column(u'DRP_modelcontainer', 'active')

        # Deleting field 'ModelContainer.built'
        db.delete_column(u'DRP_modelcontainer', 'built')

        # Removing M2M table for field boolRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_boolRxnDescriptors')

        # Removing M2M table for field ordRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_ordRxnDescriptors')

        # Removing M2M table for field catRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_catRxnDescriptors')

        # Removing M2M table for field numRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_numRxnDescriptors')

        # Removing M2M table for field outcomeBoolRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_outcomeBoolRxnDescriptors')

        # Removing M2M table for field outcomeOrdRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_outcomeOrdRxnDescriptors')

        # Removing M2M table for field outcomeCatRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_outcomeCatRxnDescriptors')

        # Removing M2M table for field outcomeNumRxnDescriptors on 'ModelContainer'
        db.delete_table('DRP_modelcontainer_outcomeNumRxnDescriptors')


        # User chose to not deal with backwards NULL issues for 'ModelContainer.splitter'
        raise RuntimeError("Cannot reverse this migration. 'ModelContainer.splitter' and its values cannot be restored.")
        # Adding field 'NumRxnDescriptorValue.model'
        db.add_column(u'DRP_numrxndescriptorvalue', 'model',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'BoolRxnDescriptorValue.model'
        db.add_column(u'DRP_boolrxndescriptorvalue', 'model',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)


        # User chose to not deal with backwards NULL issues for 'PredBoolRxnDescriptor.prediction_of'
        raise RuntimeError("Cannot reverse this migration. 'PredBoolRxnDescriptor.prediction_of' and its values cannot be restored.")

        # User chose to not deal with backwards NULL issues for 'PredBoolRxnDescriptor.stats_model'
        raise RuntimeError("Cannot reverse this migration. 'PredBoolRxnDescriptor.stats_model' and its values cannot be restored.")
        # Deleting field 'PredBoolRxnDescriptor.modelContainer'
        db.delete_column(u'DRP_predboolrxndescriptor', 'modelContainer_id')

        # Deleting field 'PredBoolRxnDescriptor.statsModel'
        db.delete_column(u'DRP_predboolrxndescriptor', 'statsModel_id')

        # Deleting field 'PredBoolRxnDescriptor.predictionOf'
        db.delete_column(u'DRP_predboolrxndescriptor', 'predictionOf_id')

        # Adding field 'OrdRxnDescriptorValue.model'
        db.add_column(u'DRP_ordrxndescriptorvalue', 'model',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)

        # Adding field 'CatRxnDescriptorValue.model'
        db.add_column(u'DRP_catrxndescriptorvalue', 'model',
                      self.gf('django.db.models.fields.related.ForeignKey')(default=None, to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)


        # User chose to not deal with backwards NULL issues for 'StatsModel.description'
        raise RuntimeError("Cannot reverse this migration. 'StatsModel.description' and its values cannot be restored.")
        # Adding field 'StatsModel.start_time'
        db.add_column(u'DRP_statsmodel', 'start_time',
                      self.gf('django.db.models.fields.DateTimeField')(default=None, null=True),
                      keep_default=False)

        # Adding field 'StatsModel.end_time'
        db.add_column(u'DRP_statsmodel', 'end_time',
                      self.gf('django.db.models.fields.DateTimeField')(default=None, null=True),
                      keep_default=False)


        # User chose to not deal with backwards NULL issues for 'StatsModel.iterations'
        raise RuntimeError("Cannot reverse this migration. 'StatsModel.iterations' and its values cannot be restored.")
        # Adding field 'StatsModel.active'
        db.add_column(u'DRP_statsmodel', 'active',
                      self.gf('django.db.models.fields.BooleanField')(default=False),
                      keep_default=False)

        # Adding field 'StatsModel.snapShot'
        db.add_column(u'DRP_statsmodel', 'snapShot',
                      self.gf('django.db.models.fields.files.FileField')(default=None, max_length=200, null=True),
                      keep_default=False)

        # Deleting field 'StatsModel.startTime'
        db.delete_column(u'DRP_statsmodel', 'startTime')

        # Deleting field 'StatsModel.endTime'
        db.delete_column(u'DRP_statsmodel', 'endTime')

        # Deleting field 'StatsModel.trainingSet'
        db.delete_column(u'DRP_statsmodel', 'trainingSet_id')

        # Adding M2M table for field tags on 'StatsModel'
        db.create_table(u'DRP_statsmodel_tags', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('statsmodeltag', models.ForeignKey(orm['DRP.statsmodeltag'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_tags', ['statsmodel_id', 'statsmodeltag_id'])

        # Adding M2M table for field descriptors on 'StatsModel'
        db.create_table(u'DRP_statsmodel_descriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('descriptor', models.ForeignKey(orm['DRP.descriptor'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_descriptors', ['statsmodel_id', 'descriptor_id'])

        # Adding M2M table for field outcomeDescriptors on 'StatsModel'
        db.create_table(u'DRP_statsmodel_outcomeDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('descriptor', models.ForeignKey(orm['DRP.descriptor'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_outcomeDescriptors', ['statsmodel_id', 'descriptor_id'])

        # Adding M2M table for field predictsDescriptors on 'StatsModel'
        db.create_table(u'DRP_statsmodel_predictsDescriptors', (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('statsmodel', models.ForeignKey(orm['DRP.statsmodel'], null=False)),
            ('descriptor', models.ForeignKey(orm['DRP.descriptor'], null=False))
        ))
        db.create_unique(u'DRP_statsmodel_predictsDescriptors', ['statsmodel_id', 'descriptor_id'])

        # Removing M2M table for field testSets on 'StatsModel'
        db.delete_table('DRP_statsmodel_testSets')

        # Adding field 'RecommendedReaction.model'
        db.add_column(u'DRP_recommendedreaction', 'model',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['DRP.StatsModel'], null=True),
                      keep_default=False)


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
        'DRP.modelcontainer': {
            'Meta': {'object_name': 'ModelContainer'},
            'active': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'boolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.BoolRxnDescriptor']", 'symmetrical': 'False'}),
            'built': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'catRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.CatRxnDescriptor']", 'symmetrical': 'False'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'fully_trained': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['DRP.StatsModel']", 'null': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'library': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'numRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.NumRxnDescriptor']", 'symmetrical': 'False'}),
            'ordRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['DRP.OrdRxnDescriptor']", 'symmetrical': 'False'}),
            'outcomeBoolRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.BoolRxnDescriptor']"}),
            'outcomeCatRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.CatRxnDescriptor']"}),
            'outcomeNumRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.NumRxnDescriptor']"}),
            'outcomeOrdRxnDescriptors': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'outcomeForModels'", 'symmetrical': 'False', 'to': "orm['DRP.OrdRxnDescriptor']"}),
            'splitter': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True', 'blank': 'True'}),
            'tool': ('django.db.models.fields.CharField', [], {'max_length': '200'})
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