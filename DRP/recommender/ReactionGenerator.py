from itertools import product

import numpy as np
import xxhash

import DRP
import DRP.plugins.rxndescriptors.drp as drp_rxn_plugin
import DRP.plugins.rxndescriptors.rxnhash as rxnhash_rxn_plugin
from DRP.models import Compound, CompoundRole, CompoundQuantity, Reaction, NumRxnDescriptor, BoolRxnDescriptor, NumRxnDescriptorValue, BoolRxnDescriptorValue

common_desc = ['reaction_temperature_manual_0', 'reaction_time_manual_0', 'reaction_pH_manual_0', 'Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_Range_DRP_0_02', 'Org_boolean_period_1_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_group_13_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_Range_DRP_0_02', 'Org_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_15_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_period_5_DRP_1_5_False_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_True_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'H_mols_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_Max_DRP_0_02', 'Inorg_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'I_mols_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Org_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_14_DRP_1_5_True_molarity_DRP_0_02', 'Org_boolean_group_15_DRP_1_5_False_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_Max_DRP_0_02', 'Inorg_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomHardness_max_DRP_0_02_Range_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_period_2_DRP_1_5_True_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_Max_DRP_0_02', 'Ox_amount_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Ga_mols_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Org_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_14_DRP_1_5_True_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_Range_DRP_0_02', 'Solv_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_Range_DRP_0_02', 'Inorg_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_Range_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_False_count_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_gmean_count_DRP_0_02', 'Org_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'K_mols_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_True_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_Max_DRP_0_02', 'Inorg_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_Max_DRP_0_02', 'Solv_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_period_4_DRP_1_5_True_molarity_DRP_0_02', 'Solv_boolean_group_1_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_period_4_DRP_1_5_False_count_DRP_0_02', 'Inorg_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Org_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'Na_mols_DRP_0_02', 'Org_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_2_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_boolean_group_15_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_Range_DRP_0_02', 'Org_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_True_count_DRP_0_02', 'Org_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_3_DRP_1_5_True_molarity_DRP_0_02', 'Solv_boolean_period_2_DRP_1_5_True_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity_DRP_0_02', 'C_mols_DRP_0_02', 'Org_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_True_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_group_1_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_gmean_molarity_DRP_0_02', 'Inorg_boolean_period_5_DRP_1_5_True_count_DRP_0_02', 'Org_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_16_DRP_1_5_True_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Org_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'V_mols_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_period_1_DRP_1_5_True_molarity_DRP_0_02', 'Solv_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_15_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_2_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_True_count_DRP_0_02', 'Inorg_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_Range_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_6_DRP_1_5_False_count_DRP_0_02', 'Inorg_boolean_group_14_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_False_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_Max_DRP_0_02', 'Solv_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_False_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_17_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_13_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_gmean_molarity_DRP_0_02', 'Org_boolean_group_16_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_14_DRP_1_5_False_molarity_DRP_0_02', 'Ox_amount_molarity_DRP_0_02', 'Inorg_boolean_period_5_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_period_3_DRP_1_5_False_count_DRP_0_02', 'Org_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_max_DRP_0_02_Max_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity_DRP_0_02', 'Mo_mols_DRP_0_02', 'Inorg_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_max_DRP_0_02_gmean_count_DRP_0_02', 'Org_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_gmean_count_DRP_0_02', 'Te_mols_DRP_0_02', 'Org_drpInorgAtom_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_1_DRP_1_5_False_molarity_DRP_0_02', 'Se_mols_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_gmean_molarity_DRP_0_02', 'Inorg_boolean_group_1_DRP_1_5_False_count_DRP_0_02', 'Cr_mols_DRP_0_02', 'Org_drpInorgAtom_boolean_group_15_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'O_mols_DRP_0_02', 'Org_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_False_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_1_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_max_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_False_count_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_max_DRP_0_02_Max_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_Range_DRP_0_02', 'N_mols_DRP_0_02', 'Org_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_Range_DRP_0_02', 'Inorg_amount_molarity_DRP_0_02', 'Solv_boolean_group_16_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomHardness_max_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_7_DRP_1_5_False_molarity_DRP_0_02', 'Solv_amount_molarity_DRP_0_02', 'Org_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Inorg_boolean_group_1_DRP_1_5_True_molarity_DRP_0_02', 'Solv_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Inorg_boolean_group_5_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_boolean_period_4_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_stoich_DRP_0_02_Range_DRP_0_02', 'Org_drpInorgAtom_boolean_group_7_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomIonizationEnergy_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_period_2_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_10_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_geom_unw_DRP_0_02_Max_DRP_0_02', 'Org_drpInorgAtom_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_9_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_valence_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_True_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_period_3_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_Max_DRP_0_02', 'Solv_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_1_DRP_1_5_True_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_8_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_14_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_13_DRP_1_5_False_count_DRP_0_02', 'Org_drpInorgAtom_boolean_period_7_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_6_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_valence_0_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_False_count_DRP_0_02', 'Inorg_boolean_group_5_DRP_1_5_False_count_DRP_0_02', 'Inorg_boolean_group_13_DRP_1_5_False_count_DRP_0_02', 'Org_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_True_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_1_DRP_1_5_False_molarity_DRP_0_02', 'Org_drpInorgAtom_boolean_group_12_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_6_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomHardness_geom_unw_DRP_0_02_Max_DRP_0_02', 'Solv_boolean_group_5_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomElectronAffinity_geom_stoich_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_max_DRP_0_02_gmean_count_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_False_count_DRP_0_02', 'Inorg_boolean_period_3_DRP_1_5_True_count_DRP_0_02', 'Inorg_drpInorgAtomAtomicRadius_geom_stoich_DRP_0_02_Range_DRP_0_02', 'Inorg_drpInorgAtom_boolean_group_11_DRP_1_5_False_molarity_DRP_0_02', 'Solv_boolean_group_18_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPearsonElectronegativity_geom_unw_DRP_0_02_gmean_count_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_16_DRP_1_5_False_molarity_DRP_0_02', 'Org_amount_molarity_DRP_0_02', 'Inorg_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_group_1_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_drpInorgAtomPaulingElectronegativity_geom_stoich_DRP_0_02_Max_DRP_0_02', 'Inorg_drpInorgAtomHardness_max_DRP_0_02_Max_DRP_0_02', 'Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_False_molarity_DRP_0_02', 'Inorg_boolean_period_1_DRP_1_5_False_molarity_DRP_0_02', 'Org_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'Solv_drpInorgAtom_boolean_group_2_DRP_1_5_False_molarity_DRP_0_02', 'leak_manual_0', 'slow_cool_manual_0', 'Inorg_boolean_group_12_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_17_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_group_16_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_period_5_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_5_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_group_6_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_13_DRP_1_5_any_DRP_1_5', 'Org_boolean_group_15_DRP_1_5_any_DRP_1_5', 'Org_boolean_period_3_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_period_5_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_group_5_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_valence_5_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_1_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_period_3_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_period_4_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_period_4_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_valence_4_DRP_1_5_any_DRP_1_5', 'Org_boolean_group_16_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_14_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_15_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_valence_6_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_6_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_valence_2_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_group_9_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_group_9_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_valence_3_DRP_1_5_any_DRP_1_5', 'Inorg_boolean_period_1_DRP_1_5_any_DRP_1_5', 'Inorg_drpInorgAtom_boolean_group_12_DRP_1_5_any_DRP_1_5', 'Org_boolean_period_4_DRP_1_5_any_DRP_1_5', 'Org_boolean_group_17_DRP_1_5_any_DRP_1_5']

def get_compound_set(role):
    compound_set = set()
    count = 0
    for comp in CompoundQuantity.objects.filter(role=role):
        print(comp)
        compound_set.add(comp.compound.id)
        count += 1
        if count > 5:
            break
    return list(compound_set)


class ReactionGenerator(object):
    """Generate the reaction grid space."""

    def __init__(self, desc_dict,  compound_dict=None, amounts_dict=None, compound_roles=None):
        """The desc dict lists descriptors and different values that they may have to construct the grid from"""
        self.desc_dict = desc_dict # A dictionary of (Descriptor, Value List) pairs.
        self.compound_dict = compound_dict
        self.compound_amounts_dict = amounts_dict
        self.compound_roles = compound_roles
        self.triples_dict = {}

    def get_all_compound_triples(self):
        org = get_compound_set(6)
        print("org: {}".format(org))
        inorg = get_compound_set(7)
        print("inorg: {}".format(inorg))
        # return looks like a tuple of tuple triples
        return product(inorg, inorg, org)

    def get_reasonable_compound_amounts(self, compounds):
        """Take a list of compound objects and return a dictionary of compounds to reasonable amounts"""
        compound_dict = {}
        for compound in compounds:
            quantities = CompoundQuantity.objects.filter(compound=compound).exclude(amount=None)
            quantities = [quantity.amount for quantity in quantities]
            quantities = np.array(quantities)
            mean_quantity = np.mean(quantities)
            quantity_sd = np.std(quantities)
            print("USING ONLY MEAN FOR COMPOUND AMOUNTS CHANGE IN PRODUCTION!")
            sample_quantities = [mean_quantity]
            #sample_quantities = [mean_quantity - 2 * quantity_sd, mean_quantity - quantity_sd, mean_quantity,
            #                     mean_quantity + quantity_sd, mean_quantity + 2 * quantity_sd]

            sample_quantities = [sq for sq in sample_quantities if sq > 0]
            compound_dict[compound] = sample_quantities

        # self.compound_amounts_dict is a dictionary with compound objects as keys and a list of values as the value
        self.compound_amounts_dict = compound_dict

    def get_all_triple_hashes(self):
        compound_triples = self.get_all_compound_triples()
        triples_dict = {}
        for triple in compound_triples:
            if len(triple) == len(set(triple)):
                h = xxhash.xxh64()  # generates a hash
                compounds = Compound.objects.filter(pk__in=triple).order_by('name')
                print (compounds)
                for reactant in compounds:
                    h.update(reactant.name) # problem here!!!, can't resolve keyword 'abbrev' into field, trying just order the reactions by a name
                triples_dict[h.hexdigest()] = compounds
        print("Triples_dict: ", triples_dict)
        # triples_dict has the form of {xxhash: compound triple queryset, ... }
        return triples_dict

    def get_all_triples(self):
        compound_triples = self.get_all_compound_triples()
        triples = []
        for triple in compound_triples:
            if len(triple) == len(set(triple)):
                compounds = Compound.objects.filter(pk__in=triple).order_by('name')
                triples.append(compounds)
        print("Triples: ", triples)
        return triples


    def get_triples_and_amounts(self, compound_triples):
        triples_and_amounts = {}
        # print 'gta', compound_triples
        # print 'cad', self.compound_amounts_dict
        for c1, c2, c3 in compound_triples:
            triples_and_amounts[(c1,c2,c3)] = list(product(self.compound_amounts_dict[c1],self.compound_amounts_dict[c2], self.compound_amounts_dict[c3]))
        
        # triples_and_amounts is a dictionary with compound triple tuple as key and tuples with three lists of compound values as the value
        # {(c1,c2,c3): ([values for c1],[values for c2],[values for c3]), ...}
        return triples_and_amounts

    
    def generate(self, compound_triples):
        """Populate the database with every recommended reaction within the grid"""
        descs = list(self.desc_dict.keys())
        desc_vals = list(self.desc_dict.values())

        triples_and_amounts = self.get_triples_and_amounts(compound_triples)

        # print (desc_vals)
        # print (descs)
        desc_combos = product(*desc_vals) 
        # print (desc_combos)
        reaction_set = []
        counter = 0
        for desc_values_instance in desc_combos:
            # print (desc_values_instance)
            for triple in triples_and_amounts:
                for compound_amts in triples_and_amounts[triple]:
                    new_rxn = Reaction()
                    # Currently assigning all of these to Alex's group 
                    new_rxn.labGroup_id = 1 
                    new_rxn.notes = 'Part of a grid search'
                    new_rxn.save()


                    # Generate (numeric) descriptor values for this reaction
                    for desc, dvalue in zip(descs, desc_values_instance):
                        #TODO Make this work for general descriptors
                        if isinstance(desc, NumRxnDescriptor):
                            new_desc = NumRxnDescriptorValue()
                        if isinstance(desc, BoolRxnDescriptor):
                            new_desc = BoolRxnDescriptorValue()
                        new_desc.id = None
                        new_desc.reaction = new_rxn
                        new_desc.descriptor = desc 
                        new_desc.value = dvalue
                        new_desc.save()

                    # Populate the compound quantity fields for the reaction
                    for compound, amount, role in zip(triple, compound_amts, [7,7,6]):
                        compound_quantity = CompoundQuantity()
                        compound_quantity.compound = compound
                        compound_quantity.reaction = new_rxn

                        #TODO Figure out compound role
                        compound_quantity.role = CompoundRole.objects.get(pk=role)
                        #END TODO 
                        compound_quantity.amount = amount

                        compound_quantity.save()

                    new_rxn.save()
                    counter += 1
                    print('Reaction #{} created.'.format(counter))
                    reaction_set.append(new_rxn)

                    # Save us from calculating the obscene number of reactions in the grid until necessary
                    yield new_rxn
        drp_rxn_plugin.calculate_many([reaction_set[0]])
        rxnhash_rxn_plugin.calculate_many([reaction_set[0]])
        # calculate reaction hashes here or manually set them from the ones generated in Hash_MI in ReactionRecommender?
        #rxnhash_rxn_plugin.calculate_many(reaction_set)
