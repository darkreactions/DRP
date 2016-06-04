"""Script for parsing descriptor headings into their underlying calculation"""
from sys import argv
from DRP import chemical_data
import csv

def strip_headers(fn, legacy=True):
    with open(fn) as f:
        headers = [l.strip() for l in f]

    if legacy:
        strip_method = legacy_strip
    else:
        strip_method = new_strip

    headers = [(header, compress(strip_method(header)).lower()) for header in headers]

    return headers
    
def strip_header(header, banned_substrings=[], banned_prefixes=[], banned_suffixes=[]):
    stripped_header = header

    for pref in banned_prefixes:
        if stripped_header.startswith(pref):
            stripped_header = stripped_header[len(pref):]
            
    for suf in banned_suffixes:
        if stripped_header.endswith(suf):
            stripped_header = stripped_header[:-len(suf)]

    for sub in banned_substrings:
        stripped_header = stripped_header.replace(sub, '_')

    stripped_header = stripped_header.replace('_', '')
    #stripped_header = stripped_header.lower()

    return stripped_header

def compress(header):
    return compress_elements(compress_groups(compress_periods(compress_valence(header))))

def compress_elements(header):
    return 'element' if header in chemical_data.elements.keys() else header

def compress_groups(header):
    return 'group' if header in ['G{}'.format(i) for i in range(1,19)] else header
    
def compress_periods(header):
    return 'period' if header in ['P{}'.format(i) for i in range(1,8)] else header
    
def compress_valence(header):
    return 'valence' if header in ['V{}'.format(i) for i in range(1,8)] else header

def legacy_strip(header):
    banned_substrings = ['Min',
                        'Max',
                        'Mean',
                        'Avg',
                        'Arith',
                        'Geom',
                        '_pHdependent',
                        'Weighted',
                        ]

    banned_prefixes = ['oxlike',
                        'org',
                        'inorg',
                        ]
                    
    banned_suffixes = ['_legacy',
                        ]

    return strip_header(header, banned_substrings, banned_prefixes, banned_suffixes)

def new_strip(header):
    banned_substrings = ['chemaxoncxcalc_15.6',
                         'drprdkit_0',
                         'drp_0.02',
                         'gmean',
                         'Max',
                         'Range',
                         'geom',
                         'stoich',
                         'drpInorgAtom',
                         'unw',
                         'max_',
                         'range',
                        ]
    banned_prefixes = ['Inorg',
                       'Org',
                       'Ox',
                       'Solv',
                       'reaction',
                      ]
    banned_suffixes = ['count',
                       'molarity',
                       'mols',
                       ]

    pH_strings = ['_pH{}_'.format(i) for i in range(1,15)] + ['_pHreaction']
    banned_substrings += pH_strings

    return strip_header(header, banned_substrings, banned_prefixes, banned_suffixes)

if __name__ == '__main__':
    fn = argv[1]
    csv_out = argv[2]

    header_tuples = strip_headers(fn, legacy=False)

    #print '\n'.join(['\t'.join(header_tuple) for header_tuple in header_tuples])

    header_dict = {}

    for header, stripped_header in header_tuples:
        try:
            header_dict[stripped_header].append(header)
        except KeyError:
            header_dict[stripped_header] = [header]

    with open(csv_out, 'w') as f:
        writer = csv.writer(f)
        for stripped_header, headers in sorted(header_dict.items()):
            writer.writerow([stripped_header] + headers)
