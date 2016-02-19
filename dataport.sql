# This file enables the porting of data from the old database structure (pre 0.02) to databases of v0.02.
# This file will not be updated post version 0.02- so if you need to upgrade to higher databases than 0.02
# you should update to 0.02 first.
# It exports the data as a series of .tsv files. These .tsv files are then digested by the manage.py commands
# /path/to/manage.py port_database /path/to/tsv/files/
# hence, it is recommended that you direct this file to export the tsv files to their own isolated directory.

#users
SELECT 'password', 'last_login', 'is_superuser', 'username', 'first_name', 'last_name', 'email', 'is_staff'
UNION ALL
SELECT password, last_login, is_superuser, username, first_name, last_name, email, is_staff FROM auth_user
    INTO OUTFILE 'User.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#lab groups
SELECT 'title', 'legacy_access_code', 'address', 'email'
UNION ALL
SELECT lab_title, access_code, lab_address, lab_email FROM DRP_lab_group
    INTO OUTFILE 'labGroup.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#link those
SELECT 'username', 'title'
UNION ALL
SELECT auth_user.username, DRP_lab_group.lab_title FROM
    DRP_lab_member JOIN
    auth_user ON DRP_lab_member.user_id = auth_user.id JOIN
    DRP_lab_group ON DRP_lab_member.lab_group_id = DRP_lab_group.id
        INTO OUTFILE 'labgroup_users.tsv'
        FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#tie compounds to labs
SELECT 'abbrev', 'name', 'CAS_ID', 'labGroup.title'
UNION ALL
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_1_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_2_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_2_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_3_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_4_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_compoundentry.compound, DRP_compoundentry.CAS_ID, DRP_lab_group.lab_title FROM
    DRP_compoundentry JOIN 
        DRP_data ON DRP_compoundentry.id = DRP_data.reactant_fk_5_id JOIN
        DRP_lab_group ON DRP_data.lab_group_id = DRP_lab_group.id
    INTO OUTFILE 'compound_labs.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#fetch classes and join to compounds
SELECT 'chemicalClass.label', 'compound.abbrev'
UNION ALL
SELECT DRP_compoundentry.compound_type, DRP_compoundentry.abbrev FROM DRP_compoundentry
    INTO OUTFILE 'compound_chemicalClasses.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#fetch reactions and innate descriptor values
SELECT 'reference', 'labGroup.title', 'notes', 'user.username', 'valid', 'duplicateOf.reference', 'legacyRecommendedFlag', 'insertedDateTime', 'public', 'temp', 'time', 'leak', 'slow_cool', 'pH','outcome', 'purity'
UNION ALL
SELECT ref, DRP_lab_group.lab_title, notes, auth_user.username, is_valid, duplicate_of, recommended, creation_time_dt, public, temp, time, leak, slow_cool, pH, outcome, purity FROM DRP_data
    JOIN DRP_lab_group ON DRP_lab_group.id = DRP_data.lab_group_id
    JOIN auth_user ON auth_user.id = DRP_data.user_id
    INTO OUTFILE 'performedReactions.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

#fetch quantities
SELECT 'compound.abbrev', 'reaction.reference', 'compoundrole.name', 'amount', 'unit'
UNION ALL
SELECT DRP_compoundentry.abbrev, DRP_data.ref, DRP_compoundentry.compound_type, DRP_data.quantity_1, DRP_data.unit_1 FROM DRP_data
    JOIN DRP_compoundentry on DRP_data.reactant_fk_1_id = DRP_compoundentry.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_data.ref, DRP_compoundentry.compound_type, DRP_data.quantity_2, DRP_data.unit_2 FROM DRP_data
    JOIN DRP_compoundentry on DRP_data.reactant_fk_2_id = DRP_compoundentry.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_data.ref, DRP_compoundentry.compound_type, DRP_data.quantity_3, DRP_data.unit_3 FROM DRP_data
    JOIN DRP_compoundentry on DRP_data.reactant_fk_3_id = DRP_compoundentry.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_data.ref, DRP_compoundentry.compound_type, DRP_data.quantity_4, DRP_data.unit_4 FROM DRP_data
    JOIN DRP_compoundentry on DRP_data.reactant_fk_4_id = DRP_compoundentry.id
UNION
SELECT DRP_compoundentry.abbrev, DRP_data.ref, DRP_compoundentry.compound_type, DRP_data.quantity_5, DRP_data.unit_5 FROM DRP_data
    JOIN DRP_compoundentry on DRP_data.reactant_fk_5_id = DRP_compoundentry.id
    INTO OUTFILE 'compoundquantities.tsv'
    FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';
