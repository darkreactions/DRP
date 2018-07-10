from copy import deepcopy

from sklearn.metrics.cluster import adjusted_mutual_info_score

from DRP.plugins.rxndescriptors.rxnhash import calculate
from DRP.recommender.AbstractRecommender import AbstractRecommender
from DRP.recommender.ReactionGenerator import ReactionGenerator
from DRP.recommender.ReactionSieve import ReactionSieve
from DRP.models import Descriptor, PerformedReaction, BoolRxnDescriptorValue, CatRxnDescriptorValue, CompoundQuantity, Reaction, Compound, RecommendedReaction, LabGroup

headers_to_use = ["C_mols","Cr_mols","Ga_mols","H_mols","I_mols","Inorg_amount_molarity","Inorg_boolean_group_10_DRP_1_5_False_molarity","Inorg_boolean_group_11_DRP_1_5_False_molarity","Inorg_boolean_group_12_DRP_1_5_any","Inorg_boolean_group_12_DRP_1_5_False_molarity","Inorg_boolean_group_13_DRP_1_5_any","Inorg_boolean_group_13_DRP_1_5_False_count","Inorg_boolean_group_13_DRP_1_5_False_molarity","Inorg_boolean_group_13_DRP_1_5_True_count","Inorg_boolean_group_13_DRP_1_5_True_molarity","Inorg_boolean_group_14_DRP_1_5_any","Inorg_boolean_group_14_DRP_1_5_False_molarity","Inorg_boolean_group_14_DRP_1_5_True_count","Inorg_boolean_group_15_DRP_1_5_any","Inorg_boolean_group_15_DRP_1_5_False_molarity","Inorg_boolean_group_16_DRP_1_5_True_molarity","Inorg_boolean_group_17_DRP_1_5_any","Inorg_boolean_group_17_DRP_1_5_False_molarity","Inorg_boolean_group_18_DRP_1_5_False_molarity","Inorg_boolean_group_1_DRP_1_5_any","Inorg_boolean_group_1_DRP_1_5_False_count","Inorg_boolean_group_1_DRP_1_5_False_molarity","Inorg_boolean_group_1_DRP_1_5_True_count","Inorg_boolean_group_1_DRP_1_5_True_molarity","Inorg_boolean_group_2_DRP_1_5_False_molarity","Inorg_boolean_group_3_DRP_1_5_False_molarity","Inorg_boolean_group_4_DRP_1_5_False_molarity","Inorg_boolean_group_5_DRP_1_5_any","Inorg_boolean_group_5_DRP_1_5_False_count","Inorg_boolean_group_5_DRP_1_5_False_molarity","Inorg_boolean_group_5_DRP_1_5_True_molarity","Inorg_boolean_group_6_DRP_1_5_any","Inorg_boolean_group_6_DRP_1_5_False_count","Inorg_boolean_group_6_DRP_1_5_False_molarity","Inorg_boolean_group_7_DRP_1_5_False_molarity","Inorg_boolean_group_8_DRP_1_5_False_molarity","Inorg_boolean_group_9_DRP_1_5_any","Inorg_boolean_group_9_DRP_1_5_False_molarity","Inorg_boolean_period_1_DRP_1_5_any","Inorg_boolean_period_1_DRP_1_5_False_molarity","Inorg_boolean_period_2_DRP_1_5_True_molarity","Inorg_boolean_period_3_DRP_1_5_any","Inorg_boolean_period_3_DRP_1_5_False_count","Inorg_boolean_period_3_DRP_1_5_False_molarity","Inorg_boolean_period_3_DRP_1_5_True_count","Inorg_boolean_period_3_DRP_1_5_True_molarity","Inorg_boolean_period_4_DRP_1_5_any","Inorg_boolean_period_4_DRP_1_5_False_count","Inorg_boolean_period_4_DRP_1_5_False_molarity","Inorg_boolean_period_4_DRP_1_5_True_count","Inorg_boolean_period_4_DRP_1_5_True_molarity","Inorg_boolean_period_5_DRP_1_5_any","Inorg_boolean_period_5_DRP_1_5_False_count","Inorg_boolean_period_5_DRP_1_5_False_molarity","Inorg_boolean_period_5_DRP_1_5_True_count","Inorg_boolean_period_5_DRP_1_5_True_molarity","Inorg_boolean_period_6_DRP_1_5_False_molarity","Inorg_boolean_period_7_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_12_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_True_count","Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_group_9_DRP_1_5_any","Inorg_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_any","Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_True_count","Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_any","Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_True_count","Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_2_DRP_1_5_any","Inorg_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_any","Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_True_count","Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_valence_4_DRP_1_5_any","Inorg_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_any","Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_False_count","Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_any","Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity","Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_True_count","Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_True_molarity","Inorg_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity","Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_gmean_molarity","Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_gmean_molarity","Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_gmean_molarity","Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_Max","Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_Range","Inorg_drpInorgAtomElectronAffinity_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomElectronAffinity_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomElectronAffinity_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_Max","Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_Range","Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomHardness_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomHardness_max_DRP_0_02_Max","Inorg_drpInorgAtomHardness_max_DRP_0_02_Range","Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_Max","Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_Range","Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_Max","Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_Range","Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_gmean_count","Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_Max","Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_Range","Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_gmean_count","Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_Max","Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_Range","Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_gmean_count","Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_Max","Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_Range","Inorg_mw_DRP_rdkit_0_02_gmean_count","Inorg_mw_DRP_rdkit_0_02_gmean_molarity","Inorg_mw_DRP_rdkit_0_02_Max","Inorg_mw_DRP_rdkit_0_02_Range","K_mols","leak","Mo_mols","N_mols","Na_mols","O_mols","Org_amount_molarity","Org_boolean_group_10_DRP_1_5_False_molarity","Org_boolean_group_11_DRP_1_5_False_molarity","Org_boolean_group_12_DRP_1_5_False_molarity","Org_boolean_group_13_DRP_1_5_False_molarity","Org_boolean_group_14_DRP_1_5_True_molarity","Org_boolean_group_15_DRP_1_5_any","Org_boolean_group_15_DRP_1_5_False_count","Org_boolean_group_15_DRP_1_5_True_molarity","Org_boolean_group_16_DRP_1_5_any","Org_boolean_group_16_DRP_1_5_False_molarity","Org_boolean_group_17_DRP_1_5_any","Org_boolean_group_17_DRP_1_5_False_molarity","Org_boolean_group_18_DRP_1_5_False_molarity","Org_boolean_group_1_DRP_1_5_True_molarity","Org_boolean_group_2_DRP_1_5_False_molarity","Org_boolean_group_3_DRP_1_5_False_molarity","Org_boolean_group_4_DRP_1_5_False_molarity","Org_boolean_group_5_DRP_1_5_False_molarity","Org_boolean_group_6_DRP_1_5_False_molarity","Org_boolean_group_7_DRP_1_5_False_molarity","Org_boolean_group_8_DRP_1_5_False_molarity","Org_boolean_group_9_DRP_1_5_False_molarity","Org_boolean_period_1_DRP_1_5_True_molarity","Org_boolean_period_2_DRP_1_5_True_molarity","Org_boolean_period_3_DRP_1_5_any","Org_boolean_period_3_DRP_1_5_False_molarity","Org_boolean_period_4_DRP_1_5_any","Org_boolean_period_4_DRP_1_5_False_molarity","Org_boolean_period_5_DRP_1_5_False_molarity","Org_boolean_period_6_DRP_1_5_False_molarity","Org_boolean_period_7_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity","Org_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity","Org_mw_DRP_rdkit_0_02_gmean_count","Org_mw_DRP_rdkit_0_02_gmean_molarity","Org_mw_DRP_rdkit_0_02_Max","Ox_amount_count","Ox_amount_molarity","Ox_boolean_group_10_DRP_1_5_False_molarity","Ox_boolean_group_11_DRP_1_5_False_molarity","Ox_boolean_group_12_DRP_1_5_False_molarity","Ox_boolean_group_13_DRP_1_5_False_molarity","Ox_boolean_group_14_DRP_1_5_True_molarity","Ox_boolean_group_15_DRP_1_5_False_molarity","Ox_boolean_group_16_DRP_1_5_True_molarity","Ox_boolean_group_17_DRP_1_5_False_molarity","Ox_boolean_group_18_DRP_1_5_False_molarity","Ox_boolean_group_1_DRP_1_5_True_molarity","Ox_boolean_group_2_DRP_1_5_False_molarity","Ox_boolean_group_3_DRP_1_5_False_molarity","Ox_boolean_group_4_DRP_1_5_False_molarity","Ox_boolean_group_5_DRP_1_5_False_molarity","Ox_boolean_group_6_DRP_1_5_False_molarity","Ox_boolean_group_7_DRP_1_5_False_molarity","Ox_boolean_group_8_DRP_1_5_False_molarity","Ox_boolean_group_9_DRP_1_5_False_molarity","Ox_boolean_period_1_DRP_1_5_False_molarity","Ox_boolean_period_2_DRP_1_5_True_molarity","Ox_boolean_period_3_DRP_1_5_True_molarity","Ox_boolean_period_4_DRP_1_5_False_molarity","Ox_boolean_period_5_DRP_1_5_False_molarity","Ox_boolean_period_6_DRP_1_5_False_molarity","Ox_boolean_period_7_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity","Ox_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity","reaction_pH","reaction_temperature","reaction_time","Se_mols","slow_cool","Solv_amount_molarity","Solv_boolean_group_10_DRP_1_5_False_molarity","Solv_boolean_group_11_DRP_1_5_False_molarity","Solv_boolean_group_12_DRP_1_5_False_molarity","Solv_boolean_group_13_DRP_1_5_False_molarity","Solv_boolean_group_14_DRP_1_5_False_molarity","Solv_boolean_group_15_DRP_1_5_False_molarity","Solv_boolean_group_16_DRP_1_5_True_molarity","Solv_boolean_group_17_DRP_1_5_False_molarity","Solv_boolean_group_18_DRP_1_5_False_molarity","Solv_boolean_group_1_DRP_1_5_True_molarity","Solv_boolean_group_2_DRP_1_5_False_molarity","Solv_boolean_group_3_DRP_1_5_False_molarity","Solv_boolean_group_4_DRP_1_5_False_molarity","Solv_boolean_group_5_DRP_1_5_False_molarity","Solv_boolean_group_6_DRP_1_5_False_molarity","Solv_boolean_group_7_DRP_1_5_False_molarity","Solv_boolean_group_8_DRP_1_5_False_molarity","Solv_boolean_group_9_DRP_1_5_False_molarity","Solv_boolean_period_1_DRP_1_5_True_molarity","Solv_boolean_period_2_DRP_1_5_True_molarity","Solv_boolean_period_3_DRP_1_5_False_molarity","Solv_boolean_period_4_DRP_1_5_False_molarity","Solv_boolean_period_5_DRP_1_5_False_molarity","Solv_boolean_period_6_DRP_1_5_False_molarity","Solv_boolean_period_7_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity","Solv_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity","Te_mols","V_mols", "rxnSpaceHash1"]

# headers_to_use = ['Pu_mols_DRP_0.02', 'Tb_mols_DRP_0.02', 'Lu_mols_DRP_0.02', 'Ru_mols_DRP_0.02', 'Ga_mols_DRP_0.02', 'Au_mols_DRP_0.02', 'B_mols_DRP_0.02', 'Li_mols_DRP_0.02', 'Ca_mols_DRP_0.02', 'Hf_mols_DRP_0.02', 'Pd_mols_DRP_0.02', 'In_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Range_DRP_DRP', 'Pr_mols_DRP_0.02', 'Ho_mols_DRP_0.02', 'Na_mols_DRP_0.02', 'Cu_mols_DRP_0.02', 'reaction_pH_manual_0', 'Po_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Zn_mols_DRP_0.02', 'Tm_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Range_DRP_DRP', 'H_mols_DRP_0.02', 'reaction_temperature_manual_0', 'At_mols_DRP_0.02', 'Te_mols_DRP_0.02', 'Rh_mols_DRP_0.02', 'Yb_mols_DRP_0.02', 'Inorg_amount_count_DRP_0.02', 'Pt_mols_DRP_0.02', 'Pm_mols_DRP_0.02', 'Sb_mols_DRP_0.02', 'Al_mols_DRP_0.02', 'F_mols_DRP_0.02', 'S_mols_DRP_0.02', 'Mn_mols_DRP_0.02', 'Cr_mols_DRP_0.02', 'Nd_mols_DRP_0.02', 'Co_mols_DRP_0.02', 'Fe_mols_DRP_0.02', 'C_mols_DRP_0.02', 'Os_mols_DRP_0.02', 'Ti_mols_DRP_0.02', 'Ir_mols_DRP_0.02', 'Si_mols_DRP_0.02', 'Re_mols_DRP_0.02', 'Ar_mols_DRP_0.02', 'Solv_amount_count_DRP_0.02', 'W_mols_DRP_0.02', 'Dy_mols_DRP_0.02', 'Y_mols_DRP_0.02', 'Cd_mols_DRP_0.02', 'Ce_mols_DRP_0.02', 'As_mols_DRP_0.02', 'Ag_mols_DRP_0.02', 'La_mols_DRP_0.02', 'U_mols_DRP_0.02', 'Np_mols_DRP_0.02', 'Eu_mols_DRP_0.02', 'Er_mols_DRP_0.02', 'K_mols_DRP_0.02', 'Pa_mols_DRP_0.02', 'Ni_mols_DRP_0.02', 'pH_amount_molarity_DRP_0.02', 'Ox_amount_molarity_DRP_0.02', 'N_mols_DRP_0.02', 'Am_mols_DRP_0.02', 'Tl_mols_DRP_0.02', 'Kr_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Max_DRP_DRP', 'Ra_mols_DRP_0.02', 'Org_amount_count_DRP_0.02', 'Se_mols_DRP_0.02', 'Be_mols_DRP_0.02', 'Th_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Mg_mols_DRP_0.02', 'Rn_mols_DRP_0.02', 'Tc_mols_DRP_0.02', 'Cl_mols_DRP_0.02', 'Sn_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Sc_mols_DRP_0.02', 'Ac_mols_DRP_0.02', 'Gd_mols_DRP_0.02', 'Xe_mols_DRP_0.02', 'Hg_mols_DRP_0.02', 'Ge_mols_DRP_0.02', 'Mo_mols_DRP_0.02', 'Cs_mols_DRP_0.02', 'Sr_mols_DRP_0.02', 'O_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Nb_mols_DRP_0.02', 'I_mols_DRP_0.02', 'Ta_mols_DRP_0.02', 'Fr_mols_DRP_0.02', 'Rb_mols_DRP_0.02', 'Ox_amount_count_DRP_0.02', 'pH_amount_count_DRP_0.02', 'Zr_mols_DRP_0.02', 'He_mols_DRP_0.02', 'Ne_mols_DRP_0.02', 'reaction_time_manual_0', 'Pb_mols_DRP_0.02', 'Br_mols_DRP_0.02', 'Inorg_amount_molarity_DRP_0.02', 'V_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Max_DRP_DRP', 'Bi_mols_DRP_0.02', 'Sm_mols_DRP_0.02', 'Solv_amount_molarity_DRP_0.02', 'Ba_mols_DRP_0.02', 'Org_amount_molarity_DRP_0.02', 'P_mols_DRP_0.02']
# headers_to_use = ['Pu_mols_DRP_0.02', 'Tb_mols_DRP_0.02', 'Lu_mols_DRP_0.02', 'Ru_mols_DRP_0.02', 'Ga_mols_DRP_0.02', 'Au_mols_DRP_0.02', 'B_mols_DRP_0.02', 'Li_mols_DRP_0.02', 'Ca_mols_DRP_0.02', 'Hf_mols_DRP_0.02', 'Pd_mols_DRP_0.02', 'In_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Range_DRP_DRP', 'Pr_mols_DRP_0.02', 'Ho_mols_DRP_0.02', 'Na_mols_DRP_0.02', 'Cu_mols_DRP_0.02', 'reaction_pH_manual_0', 'Po_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Zn_mols_DRP_0.02', 'Tm_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Range_DRP_DRP', 'H_mols_DRP_0.02', 'reaction_temperature_manual_0', 'At_mols_DRP_0.02', 'Te_mols_DRP_0.02', 'Rh_mols_DRP_0.02', 'Yb_mols_DRP_0.02', 'Inorg_amount_count_DRP_0.02', 'Pt_mols_DRP_0.02', 'Pm_mols_DRP_0.02', 'Sb_mols_DRP_0.02', 'Al_mols_DRP_0.02', 'F_mols_DRP_0.02', 'S_mols_DRP_0.02', 'Mn_mols_DRP_0.02', 'Cr_mols_DRP_0.02', 'Nd_mols_DRP_0.02', 'Co_mols_DRP_0.02', 'Fe_mols_DRP_0.02', 'C_mols_DRP_0.02', 'Os_mols_DRP_0.02', 'Ti_mols_DRP_0.02', 'Ir_mols_DRP_0.02', 'Si_mols_DRP_0.02', 'Re_mols_DRP_0.02', 'Ar_mols_DRP_0.02', 'Solv_amount_count_DRP_0.02', 'W_mols_DRP_0.02', 'Dy_mols_DRP_0.02', 'Y_mols_DRP_0.02', 'Cd_mols_DRP_0.02', 'Ce_mols_DRP_0.02', 'As_mols_DRP_0.02', 'Ag_mols_DRP_0.02', 'La_mols_DRP_0.02', 'U_mols_DRP_0.02', 'Np_mols_DRP_0.02', 'Eu_mols_DRP_0.02', 'Er_mols_DRP_0.02', 'K_mols_DRP_0.02', 'Pa_mols_DRP_0.02', 'Ni_mols_DRP_0.02', 'pH_amount_molarity_DRP_0.02', 'Ox_amount_molarity_DRP_0.02', 'N_mols_DRP_0.02', 'Am_mols_DRP_0.02', 'Tl_mols_DRP_0.02', 'Kr_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_Max_DRP_DRP', 'Ra_mols_DRP_0.02', 'Org_amount_count_DRP_0.02', 'Se_mols_DRP_0.02', 'Be_mols_DRP_0.02', 'Th_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Mg_mols_DRP_0.02', 'Rn_mols_DRP_0.02', 'Tc_mols_DRP_0.02', 'Cl_mols_DRP_0.02', 'Sn_mols_DRP_0.02', 'Inorg_mw_drp/rdkit_0_gmean_molarity_DRP_DRP', 'Sc_mols_DRP_0.02', 'Ac_mols_DRP_0.02', 'Gd_mols_DRP_0.02', 'Xe_mols_DRP_0.02', 'Hg_mols_DRP_0.02', 'Ge_mols_DRP_0.02', 'Mo_mols_DRP_0.02', 'Cs_mols_DRP_0.02', 'Sr_mols_DRP_0.02', 'O_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_gmean_count_DRP_DRP', 'Nb_mols_DRP_0.02', 'I_mols_DRP_0.02', 'Ta_mols_DRP_0.02', 'Fr_mols_DRP_0.02', 'Rb_mols_DRP_0.02', 'Ox_amount_count_DRP_0.02', 'pH_amount_count_DRP_0.02', 'Zr_mols_DRP_0.02', 'He_mols_DRP_0.02', 'Ne_mols_DRP_0.02', 'reaction_time_manual_0', 'Pb_mols_DRP_0.02', 'Br_mols_DRP_0.02', 'Inorg_amount_molarity_DRP_0.02', 'V_mols_DRP_0.02', 'Org_mw_drp/rdkit_0_Max_DRP_DRP', 'Bi_mols_DRP_0.02', 'Sm_mols_DRP_0.02', 'Solv_amount_molarity_DRP_0.02', 'Ba_mols_DRP_0.02', 'Org_amount_molarity_DRP_0.02', 'P_mols_DRP_0.02']
# headers_to_use = ['reaction_temperature_manual_0', 'reaction_time_manual_0']


class Recommender(AbstractRecommender):
    def __init__(self, model_container, grid_params, desired_desc_dict, seed = None): # Reaction.objects.get(notes = 'Placeholder reaction')):
        """
        model_container : trained statistical model from database
        grid_params : dictionary of decscriptor objects and corresponding values to input to the generator
        desired_desc_dict : dictionary of descriptor objects and values one is looking for in the final recommendation

        High level class for making recommendation. Takes a specific model, parameters given from a gridparameter object, a dictionary of desired descriptors and an optional seed.
        All Descriptor Dicts use (Descriptor, Value List) pairs, where the
        Value List is a list of any values this Descriptor should have.
        """

        # print(grid_params)
        self.generator = ReactionGenerator(grid_params)

        # self.plausible_reactions]
        self.model_container = model_container
        self.desired_desc_dict = desired_desc_dict

        # TODO remove
        self.desc_dict = grid_params
        self.hash_MI = {}

        self.lshash = None
        self.headers = None

        # Might remove?
        self.seed = seed


    def build_hash_outcome_dict(self):
        """Build a dictionary of the form 'hash' : [outcome1, outcome2, outcome3].
           The hash is a reaction descriptor in the database, outcomes correspond to our predictors.
        """
        outcome_dict = {}
        count = 0
        for rxn in PerformedReaction.objects.all()[:10]: # <<<< this is messed around with for debugging sake
            count +=1
            try:
                # TODO TEST TO MAKE SURE HASH CATEGORICAL DESCRIPTOR PK DOesN'T CHANGE
                print(Descriptor.objects.filter(heading="rxnSpaceHash1")[2].calculatorSoftwareVersion)
                hash_descriptor_id = Descriptor.objects.filter(heading="rxnSpaceHash1")[2].id
                rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=hash_descriptor_id,reaction=rxn).value
            except CatRxnDescriptorValue.DoesNotExist:
                calculate(rxn)
                try:
                    rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=hash_descriptor_id,reaction=rxn).value
                except CatRxnDescriptorValue.DoesNotExist: # missing or incorrect reactant set up
                    continue

            try:
                # TODO: change the magic number for descriptor__pk in the next line (2 is the manual boolean_crystallisation_outcome descriptor)
                reaction_outcome = BoolRxnDescriptorValue.objects.get(descriptor__pk=2,reaction=rxn)
            except BoolRxnDescriptorValue.DoesNotExist:
                continue

            # TODO: Check why reaction_outcome is solely determined if the reaction_outcome object exists rather than the value
            #       contained in reaction_outcome
            if reaction_outcome:
                reaction_outcome = 1
            else:
                reaction_outcome = 0

            if rxn_hash in outcome_dict:
                outcome_dict[rxn_hash].append(reaction_outcome)
            else:
                outcome_dict[rxn_hash] = [reaction_outcome]

        to_remove = []

        for i in outcome_dict:
            if 1 not in outcome_dict[i] or 0 not in outcome_dict[i]:
                to_remove.append(i)
            else:
                continue

        for i in to_remove:
            del outcome_dict[i]

        return outcome_dict

    def get_mutual_info(self, outcome_dict):
        """Gets mutual information given hash strings and outcomes using sklearn function."""
        hashes, outcomes = list(zip(*list(outcome_dict.items())))

        outcomes = [element for outcomeval in outcomes for element in outcomeval]
        print("Outcomes: ", outcomes)
        # TODO: check why we call this function with hashes and outcomes when it expects predicted outcomes vs. actual outcomes
        return adjusted_mutual_info_score(hashes, outcomes) # <<<< problem here, outcomes isn't of correct shape

    def hash_perform_generator(self, compound_dict, anash):
        """Determine the minimum change in information given a compound_dict and a hash to add."""
        deepdictlow = deepcopy(compound_dict)
        deepdicthigh = deepcopy(compound_dict)
        if anash in deepdictlow:
            deepdictlow[anash].append(0)
            deepdicthigh[anash].append(1)
        else:
            deepdictlow[anash] = [0]
            deepdicthigh[anash] = [1]
        return (anash, min(self.get_mutual_info(deepdictlow), self.get_mutual_info(deepdicthigh)))

    def get_best_delta_MI(self, k=100):
        """Get the compound triples that lead to the k hashes with the highest guaranteed mutual information change"""
        base_dict = self.build_hash_outcome_dict()
        triples_dict = self.generator.get_all_triple_hashes()
        hash_min_delta_MI = []
        for rxn_hash in triples_dict:
            result = self.hash_perform_generator(base_dict, rxn_hash)
            hash_min_delta_MI.append(result)
            self.hash_MI[rxn_hash] = result[1]

        hash_min_delta_MI.sort(key=(lambda x: x[1]), reverse=True)

        greatest_delta_MI_compounds = []
        for a_hash in hash_min_delta_MI[:k]:
            greatest_delta_MI_compounds.append(triples_dict[a_hash[0]])

        # greatest_delta_MI_compounds is a list of queryset of compound triples (ReactionGenerator line 65)
        return greatest_delta_MI_compounds

    def get_lshash(self):
        """Index all existing reactions based on specified headers into an lshash."""
        from lshash.lshash import LSHash
        headers = headers_to_use
        count = 0

        lsh = LSHash(1, len(headers))
        for i in PerformedReaction.objects.all().rows(True):
            to__index = []
            for header in headers:
                try:
                    to__index.append(i[header])
                except KeyError:
                    continue
            if len(to__index) == len(headers):
                lsh.index(to__index)
                count += 1
            else:
                pass

        print('count', count)
        self.lshash = lsh
        self.headers = headers

    def normalize_dist_dict(self, dist_dict):
        """Normalize all distances on a 0 to 1 scale"""
        min_knn_dist = min(dist_dict.values())
        max_knn_dist = max(dist_dict.values())
        split = max_knn_dist - min_knn_dist
        if split != 0:
            for key in dist_dict:
                dist_dict[key] = (dist_dict[key] - min_knn_dist) / split
        else:
            for key in dist_dict:
                dist_dict[key] = 1

    def get_distance_dict(self, number_of_neighbors=20, normalize=True):
        """Get the average distance of each plausible reaction to the number_of_neighbors nearest neighbors."""
        dist_dict = {}

        if self.lshash is None:
            self.get_lshash()
        # UNSURE OF WHAT CODE BELOW DOES
        # for rxn, rxn_row in zip(self.plausible_reactions, self.plausible_reactions.rows(True)):
        #     point = []
        #     for header in self.headers:
        #        point.append(rxn_row[header])
        #     # Point : Average distance from 20 nearest neighbors
        #     dist_dict[rxn] = sum([query_result[1] for query_result in self.lshash.query(point, num_results=number_of_neighbors)])/number_of_neighbors
        for rxn in self.plausible_reactions:
            dist_dict[rxn] = 1

        if normalize:
            self.normalize_dist_dict(dist_dict)

        return dist_dict


    def score_candidates(self):
        """Determine the score for all plausible reactions."""

        score_dict = {}
        dist_dict = self.get_distance_dict()
        hash_descriptor_id = Descriptor.objects.filter(heading="rxnSpaceHash1")[2].id
        for rxn in self.plausible_reactions:
            rxn_hash = CatRxnDescriptorValue.objects.get(descriptor__pk=hash_descriptor_id,reaction=rxn).value.value
            # TODO: hash_MI has xxhash as keys double check that the hashes generated by rxnhash.py for the generated reactions is the same!
            score_dict[rxn] = -1 #self.hash_MI[rxn_hash] * dist_dict[rxn] #*confidence

        return score_dict


    def create_recommendations(self):
        """Take scored reactions and create RecommendedReaction objects."""
        score_dict = self.score_candidates()

        for rxn, score in score_dict.items():
            print(rxn, score)
            rec_rxn = RecommendedReaction()
            rec_rxn.reaction_ptr = rxn
            rec_rxn.seed = self.seed
            rec_rxn.nonsense = True # TODO Change in production
            rec_rxn.hidden = True
            rec_rxn.reference = 'Holdingforlaterreplacethis' # What is this supposed to be?
            rec_rxn.labGroup = LabGroup.objects.get(pk=1) # Should the labgroup be 1 (Norquist's lab?)
            rec_rxn.score = score
            rec_rxn.saved = True
            rec_rxn.clean()
            rec_rxn.save()

    def recommend(self, k=1):
        """Recommend k highest scoring plausible reactions."""
        import time
        start_time = time.time()

        compound_objects = Compound.objects.filter(pk__in=(CompoundQuantity.objects.exclude(amount=None).values_list('compound__pk', flat=True)))

        self.generator.get_reasonable_compound_amounts(compound_objects)

        print((" top MI gain --- %s seconds ---" % (time.time() - start_time)))

        top_MI_gain_reactant_triples = self.get_best_delta_MI(k=k)

        print(("--- %s seconds ---" % (time.time() - start_time)))

        self.plausible_reactions = self.generator.generate(top_MI_gain_reactant_triples)

        # TODO: check next line
        self.a = top_MI_gain_reactant_triples

        print(("--- %s seconds ---" % (time.time() - start_time)))

        plausible_reaction_ids = [rxn.pk for rxn in self.plausible_reactions][:k]

        print(("--- %s seconds ---" % (time.time() - start_time)))

        self.plausible_reactions = Reaction.objects.filter(pk__in=plausible_reaction_ids)

        self.sieve = ReactionSieve(self.model_container, self.desired_desc_dict)

        self.plausible_reactions = self.sieve.filter(self.plausible_reactions)

        self.create_recommendations()

        print(str(k) + " reactions recommended")
        print(("--- %s seconds ---" % (time.time() - start_time)))
