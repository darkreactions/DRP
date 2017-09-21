CREATE OR REPLACE VIEW reactions_view AS
    SELECT DRP_data.*, cp_1.abbrev as cp1_ab, cp_2.abbrev as cp2_ab, cp_3.abbrev as cp3_ab, cp_4.abbrev as cp4_ab, cp_5.abbrev as cp5_ab FROM DRP_data
        LEFT JOIN DRP_compoundentry AS cp_1 ON DRP_data.reactant_fk_1_id = cp_1.id
        LEFT JOIN DRP_compoundentry AS cp_2 ON DRP_data.reactant_fk_2_id = cp_2.id
        LEFT JOIN DRP_compoundentry AS cp_3 ON DRP_data.reactant_fk_3_id = cp_3.id
        LEFT JOIN DRP_compoundentry AS cp_4 ON DRP_data.reactant_fk_4_id = cp_4.id
        LEFT JOIN DRP_compoundentry AS cp_5 ON DRP_data.reactant_fk_5_id = cp_5.id;


SELECT GROUP_CONCAT(CONCAT('"', COLUMN_NAME, '"')) FROM INFORMATION_SCHEMA.COLUMNS
    WHERE TABLE_NAME='reactions_view' AND
        TABLE_SCHEMA='old_drp_again'
    ORDER BY ORDINAL_POSITION
        INTO OUTFILE '/tmp/reactions_headers.tsv'
            FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';

SELECT * FROM reactions_view
    INTO OUTFILE '/tmp/reactions_out.tsv'
        FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n';
