import math

Rv = 3.1
c = 2.9979e10
a = 1.0

prefix_in_out_file_3 = """read serr 2
cpd /xw
scr white
la x Time (MJD)
la y ZTF nufnu flux
la file
time off
ma 17 on 2
"""

list_of_filters = ['U  ', 'B  ', 'V  ', 'R  ', 'I  ', 'J  ', 'H  ', 'K  ', 'u  ', 'g  ', 'r  ', 'i  ', 'z  ', 'psg',
                   'psr', 'psi', 'psz', 'psy', 'fuv', 'nuv', 'su ', 'sv ', 'sb ', 'sw1', 'sm2', 'sw2', 'xu ', 'xv ',
                   'xb ', 'xw1', 'xm2', 'xw2', 'ww1', 'ww2', 'ww3', 'ww4', 'bbG']

filter_data = {
    'U  ': (3600, math.log10(1810) - 23),
    'B  ': (4400, math.log10(4260) - 23),
    'V  ': (5500, math.log10(3640) - 23),
    'R  ': (6400, math.log10(3080) - 23),
    'I  ': (7900, math.log10(2550) - 23),
    'J  ': (12350, math.log10(1594) - 23),
    'H  ': (16620, math.log10(1024) - 23),
    'K  ': (21590, math.log10(666.7) - 23),
    'u  ': (3568, math.log10(3631) - 23),
    'g  ': (4653, math.log10(3631) - 23),
    'r  ': (6148, math.log10(3631) - 23),
    'i  ': (7468, math.log10(3631) - 23),
    'z  ': (8863, math.log10(3631) - 23),
    'psg': (4810, math.log10(3631) - 23),
    'psr': (6170, math.log10(3631) - 23),
    'psi': (7520, math.log10(3631) - 23),
    'psz': (8660, math.log10(3631) - 23),
    'psy': (9620, math.log10(3631) - 23),
    'bbG': (6730, math.log10(2918) - 23),
    'fuv': (1538.6, math.log10(3631) - 23),
    'nuv': (2315.7, math.log10(3631) - 23),
    'su ': (3501, math.log10(3631) - 23),
    'sb ': (4329, math.log10(3631) - 23),
    'sv ': (5402, math.log10(3631) - 23),
    'sm2': (2231, math.log10(3631) - 23),
    'sw2': (2030, math.log10(3631) - 23),
    'xu ': (3440, math.log10(3631) - 23),
    'xb ': (4500, math.log10(3631) - 23),
    'xv ': (5430, math.log10(3631) - 23),
    'xw1': (2910, math.log10(3631) - 23),
    'xm2': (2310, math.log10(3631) - 23),
    'xw2': (2120, math.log10(3631) - 23),
    'ww1': (34000, math.log10(309.540) - 23),
    'ww2': (46000, math.log10(171.787) - 23),
    'ww3': (120000, math.log10(31.674) - 23),
    'ww4': (220000, math.log10(8.363) - 23),
}