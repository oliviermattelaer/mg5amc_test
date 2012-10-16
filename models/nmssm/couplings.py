# This file was automatically created by FeynRules 1.7.53
# Mathematica version: 8.0 for Linux x86 (64-bit) (February 23, 2011)
# Date: Wed 1 Aug 2012 00:44:44


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec



GC_1 = Coupling(name = 'GC_1',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '-G',
                order = {'QCD':1})

GC_7 = Coupling(name = 'GC_7',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_8 = Coupling(name = 'GC_8',
                value = 'G',
                order = {'QCD':1})

GC_9 = Coupling(name = 'GC_9',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_10 = Coupling(name = 'GC_10',
                 value = '(2*ee**2*complex(0,1)*I12x11)/9.',
                 order = {'QED':2})

GC_11 = Coupling(name = 'GC_11',
                 value = '(-2*ee*complex(0,1)*G*I12x11)/3.',
                 order = {'QCD':1,'QED':1})

GC_12 = Coupling(name = 'GC_12',
                 value = 'complex(0,1)*G**2*I12x11',
                 order = {'QCD':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '(2*ee**2*complex(0,1)*I12x22)/9.',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '(-2*ee*complex(0,1)*G*I12x22)/3.',
                 order = {'QCD':1,'QED':1})

GC_15 = Coupling(name = 'GC_15',
                 value = 'complex(0,1)*G**2*I12x22',
                 order = {'QCD':2})

GC_16 = Coupling(name = 'GC_16',
                 value = '(2*ee**2*complex(0,1)*I12x33)/9. + (2*ee**2*complex(0,1)*I13x33)/9.',
                 order = {'QED':2})

GC_17 = Coupling(name = 'GC_17',
                 value = '(-2*ee*complex(0,1)*G*I12x33)/3. - (2*ee*complex(0,1)*G*I13x33)/3.',
                 order = {'QCD':1,'QED':1})

GC_18 = Coupling(name = 'GC_18',
                 value = 'complex(0,1)*G**2*I12x33 + complex(0,1)*G**2*I13x33',
                 order = {'QCD':2})

GC_19 = Coupling(name = 'GC_19',
                 value = '(2*ee**2*complex(0,1)*I12x36)/9. + (2*ee**2*complex(0,1)*I13x36)/9.',
                 order = {'QED':2})

GC_20 = Coupling(name = 'GC_20',
                 value = '(-2*ee*complex(0,1)*G*I12x36)/3. - (2*ee*complex(0,1)*G*I13x36)/3.',
                 order = {'QCD':1,'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = 'complex(0,1)*G**2*I12x36 + complex(0,1)*G**2*I13x36',
                 order = {'QCD':2})

GC_22 = Coupling(name = 'GC_22',
                 value = '(2*ee**2*complex(0,1)*I13x44)/9.',
                 order = {'QED':2})

GC_23 = Coupling(name = 'GC_23',
                 value = '(-2*ee*complex(0,1)*G*I13x44)/3.',
                 order = {'QCD':1,'QED':1})

GC_24 = Coupling(name = 'GC_24',
                 value = 'complex(0,1)*G**2*I13x44',
                 order = {'QCD':2})

GC_25 = Coupling(name = 'GC_25',
                 value = '(2*ee**2*complex(0,1)*I13x55)/9.',
                 order = {'QED':2})

GC_26 = Coupling(name = 'GC_26',
                 value = '(-2*ee*complex(0,1)*G*I13x55)/3.',
                 order = {'QCD':1,'QED':1})

GC_27 = Coupling(name = 'GC_27',
                 value = 'complex(0,1)*G**2*I13x55',
                 order = {'QCD':2})

GC_28 = Coupling(name = 'GC_28',
                 value = '(2*ee**2*complex(0,1)*I12x63)/9. + (2*ee**2*complex(0,1)*I13x63)/9.',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '(-2*ee*complex(0,1)*G*I12x63)/3. - (2*ee*complex(0,1)*G*I13x63)/3.',
                 order = {'QCD':1,'QED':1})

GC_30 = Coupling(name = 'GC_30',
                 value = 'complex(0,1)*G**2*I12x63 + complex(0,1)*G**2*I13x63',
                 order = {'QCD':2})

GC_31 = Coupling(name = 'GC_31',
                 value = '(2*ee**2*complex(0,1)*I12x66)/9. + (2*ee**2*complex(0,1)*I13x66)/9.',
                 order = {'QED':2})

GC_32 = Coupling(name = 'GC_32',
                 value = '(-2*ee*complex(0,1)*G*I12x66)/3. - (2*ee*complex(0,1)*G*I13x66)/3.',
                 order = {'QCD':1,'QED':1})

GC_33 = Coupling(name = 'GC_33',
                 value = 'complex(0,1)*G**2*I12x66 + complex(0,1)*G**2*I13x66',
                 order = {'QCD':2})

GC_34 = Coupling(name = 'GC_34',
                 value = 'ee*complex(0,1)*I19x11',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = 'ee*complex(0,1)*I19x22',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = 'ee*complex(0,1)*I19x33 + ee*complex(0,1)*I20x33',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = 'ee*complex(0,1)*I19x36 + ee*complex(0,1)*I20x36',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = 'ee*complex(0,1)*I20x44',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = 'ee*complex(0,1)*I20x55',
                 order = {'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = '-(ee*complex(0,1)*I19x63) - ee*complex(0,1)*I20x63',
                 order = {'QED':1})

GC_41 = Coupling(name = 'GC_41',
                 value = 'ee*complex(0,1)*I19x66 + ee*complex(0,1)*I20x66',
                 order = {'QED':1})

GC_42 = Coupling(name = 'GC_42',
                 value = '2*ee**2*complex(0,1)*I25x11',
                 order = {'QED':2})

GC_43 = Coupling(name = 'GC_43',
                 value = '2*ee**2*complex(0,1)*I25x22',
                 order = {'QED':2})

GC_44 = Coupling(name = 'GC_44',
                 value = '2*ee**2*complex(0,1)*I25x33 + 2*ee**2*complex(0,1)*I26x33',
                 order = {'QED':2})

GC_45 = Coupling(name = 'GC_45',
                 value = '2*ee**2*complex(0,1)*I25x36 + 2*ee**2*complex(0,1)*I26x36',
                 order = {'QED':2})

GC_46 = Coupling(name = 'GC_46',
                 value = '2*ee**2*complex(0,1)*I26x44',
                 order = {'QED':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '2*ee**2*complex(0,1)*I26x55',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '2*ee**2*complex(0,1)*I25x63 + 2*ee**2*complex(0,1)*I26x63',
                 order = {'QED':2})

GC_49 = Coupling(name = 'GC_49',
                 value = '2*ee**2*complex(0,1)*I25x66 + 2*ee**2*complex(0,1)*I26x66',
                 order = {'QED':2})

GC_50 = Coupling(name = 'GC_50',
                 value = '(-2*ee*complex(0,1)*I39x11)/3.',
                 order = {'QED':1})

GC_51 = Coupling(name = 'GC_51',
                 value = '-(complex(0,1)*G*I39x11)',
                 order = {'QCD':1})

GC_52 = Coupling(name = 'GC_52',
                 value = '(-2*ee*complex(0,1)*I39x22)/3.',
                 order = {'QED':1})

GC_53 = Coupling(name = 'GC_53',
                 value = '-(complex(0,1)*G*I39x22)',
                 order = {'QCD':1})

GC_54 = Coupling(name = 'GC_54',
                 value = '(-2*ee*complex(0,1)*I39x33)/3. - (2*ee*complex(0,1)*I40x33)/3.',
                 order = {'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = '-(complex(0,1)*G*I39x33) - complex(0,1)*G*I40x33',
                 order = {'QCD':1})

GC_56 = Coupling(name = 'GC_56',
                 value = '(-2*ee*complex(0,1)*I39x36)/3. - (2*ee*complex(0,1)*I40x36)/3.',
                 order = {'QED':1})

GC_57 = Coupling(name = 'GC_57',
                 value = '-(complex(0,1)*G*I39x36) - complex(0,1)*G*I40x36',
                 order = {'QCD':1})

GC_58 = Coupling(name = 'GC_58',
                 value = '(-2*ee*complex(0,1)*I40x44)/3.',
                 order = {'QED':1})

GC_59 = Coupling(name = 'GC_59',
                 value = '-(complex(0,1)*G*I40x44)',
                 order = {'QCD':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '(-2*ee*complex(0,1)*I40x55)/3.',
                 order = {'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = '-(complex(0,1)*G*I40x55)',
                 order = {'QCD':1})

GC_62 = Coupling(name = 'GC_62',
                 value = '(2*ee*complex(0,1)*I39x63)/3. + (2*ee*complex(0,1)*I40x63)/3.',
                 order = {'QED':1})

GC_63 = Coupling(name = 'GC_63',
                 value = 'complex(0,1)*G*I39x63 + complex(0,1)*G*I40x63',
                 order = {'QCD':1})

GC_64 = Coupling(name = 'GC_64',
                 value = '(-2*ee*complex(0,1)*I39x66)/3. - (2*ee*complex(0,1)*I40x66)/3.',
                 order = {'QED':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '-(complex(0,1)*G*I39x66) - complex(0,1)*G*I40x66',
                 order = {'QCD':1})

GC_66 = Coupling(name = 'GC_66',
                 value = '(ee*complex(0,1)*I5x11)/3.',
                 order = {'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '-(complex(0,1)*G*I5x11)',
                 order = {'QCD':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '(ee*complex(0,1)*I5x22)/3.',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '-(complex(0,1)*G*I5x22)',
                 order = {'QCD':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '(8*ee**2*complex(0,1)*I62x11)/9.',
                 order = {'QED':2})

GC_71 = Coupling(name = 'GC_71',
                 value = '(4*ee*complex(0,1)*G*I62x11)/3.',
                 order = {'QCD':1,'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = 'complex(0,1)*G**2*I62x11',
                 order = {'QCD':2})

GC_73 = Coupling(name = 'GC_73',
                 value = '(8*ee**2*complex(0,1)*I62x22)/9.',
                 order = {'QED':2})

GC_74 = Coupling(name = 'GC_74',
                 value = '(4*ee*complex(0,1)*G*I62x22)/3.',
                 order = {'QCD':1,'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = 'complex(0,1)*G**2*I62x22',
                 order = {'QCD':2})

GC_76 = Coupling(name = 'GC_76',
                 value = '(8*ee**2*complex(0,1)*I62x33)/9. + (8*ee**2*complex(0,1)*I63x33)/9.',
                 order = {'QED':2})

GC_77 = Coupling(name = 'GC_77',
                 value = '(4*ee*complex(0,1)*G*I62x33)/3. + (4*ee*complex(0,1)*G*I63x33)/3.',
                 order = {'QCD':1,'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = 'complex(0,1)*G**2*I62x33 + complex(0,1)*G**2*I63x33',
                 order = {'QCD':2})

GC_79 = Coupling(name = 'GC_79',
                 value = '(8*ee**2*complex(0,1)*I62x36)/9. + (8*ee**2*complex(0,1)*I63x36)/9.',
                 order = {'QED':2})

GC_80 = Coupling(name = 'GC_80',
                 value = '(4*ee*complex(0,1)*G*I62x36)/3. + (4*ee*complex(0,1)*G*I63x36)/3.',
                 order = {'QCD':1,'QED':1})

GC_81 = Coupling(name = 'GC_81',
                 value = 'complex(0,1)*G**2*I62x36 + complex(0,1)*G**2*I63x36',
                 order = {'QCD':2})

GC_82 = Coupling(name = 'GC_82',
                 value = '(8*ee**2*complex(0,1)*I63x44)/9.',
                 order = {'QED':2})

GC_83 = Coupling(name = 'GC_83',
                 value = '(4*ee*complex(0,1)*G*I63x44)/3.',
                 order = {'QCD':1,'QED':1})

GC_84 = Coupling(name = 'GC_84',
                 value = 'complex(0,1)*G**2*I63x44',
                 order = {'QCD':2})

GC_85 = Coupling(name = 'GC_85',
                 value = '(8*ee**2*complex(0,1)*I63x55)/9.',
                 order = {'QED':2})

GC_86 = Coupling(name = 'GC_86',
                 value = '(4*ee*complex(0,1)*G*I63x55)/3.',
                 order = {'QCD':1,'QED':1})

GC_87 = Coupling(name = 'GC_87',
                 value = 'complex(0,1)*G**2*I63x55',
                 order = {'QCD':2})

GC_88 = Coupling(name = 'GC_88',
                 value = '(8*ee**2*complex(0,1)*I62x63)/9. + (8*ee**2*complex(0,1)*I63x63)/9.',
                 order = {'QED':2})

GC_89 = Coupling(name = 'GC_89',
                 value = '(4*ee*complex(0,1)*G*I62x63)/3. + (4*ee*complex(0,1)*G*I63x63)/3.',
                 order = {'QCD':1,'QED':1})

GC_90 = Coupling(name = 'GC_90',
                 value = 'complex(0,1)*G**2*I62x63 + complex(0,1)*G**2*I63x63',
                 order = {'QCD':2})

GC_91 = Coupling(name = 'GC_91',
                 value = '(8*ee**2*complex(0,1)*I62x66)/9. + (8*ee**2*complex(0,1)*I63x66)/9.',
                 order = {'QED':2})

GC_92 = Coupling(name = 'GC_92',
                 value = '(4*ee*complex(0,1)*G*I62x66)/3. + (4*ee*complex(0,1)*G*I63x66)/3.',
                 order = {'QCD':1,'QED':1})

GC_93 = Coupling(name = 'GC_93',
                 value = 'complex(0,1)*G**2*I62x66 + complex(0,1)*G**2*I63x66',
                 order = {'QCD':2})

GC_94 = Coupling(name = 'GC_94',
                 value = '(ee*complex(0,1)*I5x33)/3. + (ee*complex(0,1)*I6x33)/3.',
                 order = {'QED':1})

GC_95 = Coupling(name = 'GC_95',
                 value = '-(complex(0,1)*G*I5x33) - complex(0,1)*G*I6x33',
                 order = {'QCD':1})

GC_96 = Coupling(name = 'GC_96',
                 value = '(ee*complex(0,1)*I5x36)/3. + (ee*complex(0,1)*I6x36)/3.',
                 order = {'QED':1})

GC_97 = Coupling(name = 'GC_97',
                 value = '-(complex(0,1)*G*I5x36) - complex(0,1)*G*I6x36',
                 order = {'QCD':1})

GC_98 = Coupling(name = 'GC_98',
                 value = '(ee*complex(0,1)*I6x44)/3.',
                 order = {'QED':1})

GC_99 = Coupling(name = 'GC_99',
                 value = '-(complex(0,1)*G*I6x44)',
                 order = {'QCD':1})

GC_100 = Coupling(name = 'GC_100',
                  value = '(ee*complex(0,1)*I6x55)/3.',
                  order = {'QED':1})

GC_101 = Coupling(name = 'GC_101',
                  value = '-(complex(0,1)*G*I6x55)',
                  order = {'QCD':1})

GC_102 = Coupling(name = 'GC_102',
                  value = '-(ee*complex(0,1)*I5x63)/3. - (ee*complex(0,1)*I6x63)/3.',
                  order = {'QED':1})

GC_103 = Coupling(name = 'GC_103',
                  value = 'complex(0,1)*G*I5x63 + complex(0,1)*G*I6x63',
                  order = {'QCD':1})

GC_104 = Coupling(name = 'GC_104',
                  value = '(ee*complex(0,1)*I5x66)/3. + (ee*complex(0,1)*I6x66)/3.',
                  order = {'QED':1})

GC_105 = Coupling(name = 'GC_105',
                  value = '-(complex(0,1)*G*I5x66) - complex(0,1)*G*I6x66',
                  order = {'QCD':1})

GC_106 = Coupling(name = 'GC_106',
                  value = '-(complex(0,1)*G*Rd1x1*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_107 = Coupling(name = 'GC_107',
                  value = '-(complex(0,1)*G*Rd2x2*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_108 = Coupling(name = 'GC_108',
                  value = '-(complex(0,1)*G*Rd3x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_109 = Coupling(name = 'GC_109',
                  value = 'complex(0,1)*G*Rd3x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_110 = Coupling(name = 'GC_110',
                  value = 'complex(0,1)*G*Rd4x4*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_111 = Coupling(name = 'GC_111',
                  value = 'complex(0,1)*G*Rd5x5*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_112 = Coupling(name = 'GC_112',
                  value = '-(complex(0,1)*G*Rd6x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_113 = Coupling(name = 'GC_113',
                  value = 'complex(0,1)*G*Rd6x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_114 = Coupling(name = 'GC_114',
                  value = '-(complex(0,1)*G*Ru1x1*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_115 = Coupling(name = 'GC_115',
                  value = '-(complex(0,1)*G*Ru2x2*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_116 = Coupling(name = 'GC_116',
                  value = '-(complex(0,1)*G*Ru3x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_117 = Coupling(name = 'GC_117',
                  value = 'complex(0,1)*G*Ru3x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_118 = Coupling(name = 'GC_118',
                  value = 'complex(0,1)*G*Ru4x4*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_119 = Coupling(name = 'GC_119',
                  value = 'complex(0,1)*G*Ru5x5*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_120 = Coupling(name = 'GC_120',
                  value = '-(complex(0,1)*G*Ru6x3*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_121 = Coupling(name = 'GC_121',
                  value = 'complex(0,1)*G*Ru6x6*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_122 = Coupling(name = 'GC_122',
                  value = '-(ee**2*complex(0,1)) + (ee**2*complex(0,1))/sw**2',
                  order = {'QED':2})

GC_123 = Coupling(name = 'GC_123',
                  value = '(ee**2*complex(0,1))/(2.*sw**2)',
                  order = {'QED':2})

GC_124 = Coupling(name = 'GC_124',
                  value = '-((ee**2*complex(0,1))/sw**2)',
                  order = {'QED':2})

GC_125 = Coupling(name = 'GC_125',
                  value = '(ee**2*complex(0,1)*I96x11)/(2.*sw**2)',
                  order = {'QED':2})

GC_126 = Coupling(name = 'GC_126',
                  value = '(ee**2*complex(0,1)*I96x22)/(2.*sw**2)',
                  order = {'QED':2})

GC_127 = Coupling(name = 'GC_127',
                  value = '(ee**2*complex(0,1)*I96x33)/(2.*sw**2)',
                  order = {'QED':2})

GC_128 = Coupling(name = 'GC_128',
                  value = '(ee**2*complex(0,1)*I96x36)/(2.*sw**2)',
                  order = {'QED':2})

GC_129 = Coupling(name = 'GC_129',
                  value = '(ee**2*complex(0,1)*I96x63)/(2.*sw**2)',
                  order = {'QED':2})

GC_130 = Coupling(name = 'GC_130',
                  value = '(ee**2*complex(0,1)*I96x66)/(2.*sw**2)',
                  order = {'QED':2})

GC_131 = Coupling(name = 'GC_131',
                  value = '(ee**2*complex(0,1)*I97x11)/(2.*sw**2)',
                  order = {'QED':2})

GC_132 = Coupling(name = 'GC_132',
                  value = '(ee**2*complex(0,1)*I97x22)/(2.*sw**2)',
                  order = {'QED':2})

GC_133 = Coupling(name = 'GC_133',
                  value = '(ee**2*complex(0,1)*I97x33)/(2.*sw**2)',
                  order = {'QED':2})

GC_134 = Coupling(name = 'GC_134',
                  value = '(ee**2*complex(0,1)*I97x36)/(2.*sw**2)',
                  order = {'QED':2})

GC_135 = Coupling(name = 'GC_135',
                  value = '(ee**2*complex(0,1)*I97x63)/(2.*sw**2)',
                  order = {'QED':2})

GC_136 = Coupling(name = 'GC_136',
                  value = '(ee**2*complex(0,1)*I97x66)/(2.*sw**2)',
                  order = {'QED':2})

GC_137 = Coupling(name = 'GC_137',
                  value = '(ee**2*complex(0,1)*I98x11)/(2.*sw**2)',
                  order = {'QED':2})

GC_138 = Coupling(name = 'GC_138',
                  value = '(ee**2*complex(0,1)*I98x22)/(2.*sw**2)',
                  order = {'QED':2})

GC_139 = Coupling(name = 'GC_139',
                  value = '(ee**2*complex(0,1)*I98x33)/(2.*sw**2)',
                  order = {'QED':2})

GC_140 = Coupling(name = 'GC_140',
                  value = '(ee**2*complex(0,1)*I98x36)/(2.*sw**2)',
                  order = {'QED':2})

GC_141 = Coupling(name = 'GC_141',
                  value = '(ee**2*complex(0,1)*I98x63)/(2.*sw**2)',
                  order = {'QED':2})

GC_142 = Coupling(name = 'GC_142',
                  value = '(ee**2*complex(0,1)*I98x66)/(2.*sw**2)',
                  order = {'QED':2})

GC_143 = Coupling(name = 'GC_143',
                  value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_144 = Coupling(name = 'GC_144',
                  value = '(cw*ee*complex(0,1))/sw',
                  order = {'QED':1})

GC_145 = Coupling(name = 'GC_145',
                  value = '(-2*cw*ee**2*complex(0,1))/sw',
                  order = {'QED':2})

GC_146 = Coupling(name = 'GC_146',
                  value = '-((ee**2*complex(0,1)*I27x11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_147 = Coupling(name = 'GC_147',
                  value = '-((ee**2*complex(0,1)*I27x22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_148 = Coupling(name = 'GC_148',
                  value = '-((ee**2*complex(0,1)*I27x33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_149 = Coupling(name = 'GC_149',
                  value = '-((ee**2*complex(0,1)*I27x36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_150 = Coupling(name = 'GC_150',
                  value = '-((ee**2*complex(0,1)*I33x11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_151 = Coupling(name = 'GC_151',
                  value = '-((ee**2*complex(0,1)*I33x22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_152 = Coupling(name = 'GC_152',
                  value = '-((ee**2*complex(0,1)*I33x33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_153 = Coupling(name = 'GC_153',
                  value = '-((ee**2*complex(0,1)*I33x36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_154 = Coupling(name = 'GC_154',
                  value = '(ee**2*complex(0,1)*I41x11)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_155 = Coupling(name = 'GC_155',
                  value = '(ee*complex(0,1)*G*I41x11*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_156 = Coupling(name = 'GC_156',
                  value = '(ee**2*complex(0,1)*I41x22)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_157 = Coupling(name = 'GC_157',
                  value = '(ee*complex(0,1)*G*I41x22*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_158 = Coupling(name = 'GC_158',
                  value = '(ee**2*complex(0,1)*I41x33)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_159 = Coupling(name = 'GC_159',
                  value = '(ee*complex(0,1)*G*I41x33*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_160 = Coupling(name = 'GC_160',
                  value = '(ee**2*complex(0,1)*I41x36)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_161 = Coupling(name = 'GC_161',
                  value = '(ee*complex(0,1)*G*I41x36*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_162 = Coupling(name = 'GC_162',
                  value = '(ee**2*complex(0,1)*I41x63)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_163 = Coupling(name = 'GC_163',
                  value = '(ee*complex(0,1)*G*I41x63*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_164 = Coupling(name = 'GC_164',
                  value = '(ee**2*complex(0,1)*I41x66)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_165 = Coupling(name = 'GC_165',
                  value = '(ee*complex(0,1)*G*I41x66*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_166 = Coupling(name = 'GC_166',
                  value = '(ee**2*complex(0,1)*I54x11)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_167 = Coupling(name = 'GC_167',
                  value = '(ee*complex(0,1)*G*I54x11*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_168 = Coupling(name = 'GC_168',
                  value = '(ee**2*complex(0,1)*I54x22)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_169 = Coupling(name = 'GC_169',
                  value = '(ee*complex(0,1)*G*I54x22*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_170 = Coupling(name = 'GC_170',
                  value = '(ee**2*complex(0,1)*I54x33)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_171 = Coupling(name = 'GC_171',
                  value = '(ee*complex(0,1)*G*I54x33*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_172 = Coupling(name = 'GC_172',
                  value = '(ee**2*complex(0,1)*I54x36)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_173 = Coupling(name = 'GC_173',
                  value = '(ee*complex(0,1)*G*I54x36*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_174 = Coupling(name = 'GC_174',
                  value = '(ee**2*complex(0,1)*I54x63)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_175 = Coupling(name = 'GC_175',
                  value = '(ee*complex(0,1)*G*I54x63*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_176 = Coupling(name = 'GC_176',
                  value = '(ee**2*complex(0,1)*I54x66)/(3.*sw*cmath.sqrt(2))',
                  order = {'QED':2})

GC_177 = Coupling(name = 'GC_177',
                  value = '(ee*complex(0,1)*G*I54x66*cmath.sqrt(2))/sw',
                  order = {'QCD':1,'QED':1})

GC_178 = Coupling(name = 'GC_178',
                  value = '-((ee*complex(0,1)*I92x11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_179 = Coupling(name = 'GC_179',
                  value = '-((ee*complex(0,1)*I92x22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_180 = Coupling(name = 'GC_180',
                  value = '-((ee*complex(0,1)*I92x33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_181 = Coupling(name = 'GC_181',
                  value = '-((ee*complex(0,1)*I92x36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_182 = Coupling(name = 'GC_182',
                  value = '-((ee*complex(0,1)*I92x63)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_183 = Coupling(name = 'GC_183',
                  value = '-((ee*complex(0,1)*I92x66)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_184 = Coupling(name = 'GC_184',
                  value = '-((ee*complex(0,1)*I93x11)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_185 = Coupling(name = 'GC_185',
                  value = '-((ee*complex(0,1)*I93x22)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_186 = Coupling(name = 'GC_186',
                  value = '-((ee*complex(0,1)*I93x33)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_187 = Coupling(name = 'GC_187',
                  value = '-((ee*complex(0,1)*I93x36)/(sw*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_188 = Coupling(name = 'GC_188',
                  value = '(ee*complex(0,1)*I94x11)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_189 = Coupling(name = 'GC_189',
                  value = '(ee*complex(0,1)*I94x22)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_190 = Coupling(name = 'GC_190',
                  value = '(ee*complex(0,1)*I94x33)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_191 = Coupling(name = 'GC_191',
                  value = '(ee*complex(0,1)*I94x36)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_192 = Coupling(name = 'GC_192',
                  value = '(ee*complex(0,1)*I94x63)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_193 = Coupling(name = 'GC_193',
                  value = '(ee*complex(0,1)*I94x66)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_194 = Coupling(name = 'GC_194',
                  value = '(ee*complex(0,1)*I95x11)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_195 = Coupling(name = 'GC_195',
                  value = '(ee*complex(0,1)*I95x22)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_196 = Coupling(name = 'GC_196',
                  value = '(ee*complex(0,1)*I95x33)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_197 = Coupling(name = 'GC_197',
                  value = '(ee*complex(0,1)*I95x36)/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_198 = Coupling(name = 'GC_198',
                  value = '(cw*ee**2*complex(0,1)*I92x11)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_199 = Coupling(name = 'GC_199',
                  value = '(cw*ee**2*complex(0,1)*I92x22)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_200 = Coupling(name = 'GC_200',
                  value = '(cw*ee**2*complex(0,1)*I92x33)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_201 = Coupling(name = 'GC_201',
                  value = '(cw*ee**2*complex(0,1)*I92x36)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_202 = Coupling(name = 'GC_202',
                  value = '(cw*ee**2*complex(0,1)*I92x63)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_203 = Coupling(name = 'GC_203',
                  value = '(cw*ee**2*complex(0,1)*I92x66)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_204 = Coupling(name = 'GC_204',
                  value = '-((cw*ee**2*complex(0,1)*I93x11)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_205 = Coupling(name = 'GC_205',
                  value = '-((cw*ee**2*complex(0,1)*I93x22)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_206 = Coupling(name = 'GC_206',
                  value = '-((cw*ee**2*complex(0,1)*I93x33)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_207 = Coupling(name = 'GC_207',
                  value = '-((cw*ee**2*complex(0,1)*I93x36)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_208 = Coupling(name = 'GC_208',
                  value = '(cw*ee**2*complex(0,1)*I94x11)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_209 = Coupling(name = 'GC_209',
                  value = '(cw*ee**2*complex(0,1)*I94x22)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_210 = Coupling(name = 'GC_210',
                  value = '(cw*ee**2*complex(0,1)*I94x33)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_211 = Coupling(name = 'GC_211',
                  value = '(cw*ee**2*complex(0,1)*I94x36)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_212 = Coupling(name = 'GC_212',
                  value = '(cw*ee**2*complex(0,1)*I94x63)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_213 = Coupling(name = 'GC_213',
                  value = '(cw*ee**2*complex(0,1)*I94x66)/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':2})

GC_214 = Coupling(name = 'GC_214',
                  value = '-((cw*ee**2*complex(0,1)*I95x11)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_215 = Coupling(name = 'GC_215',
                  value = '-((cw*ee**2*complex(0,1)*I95x22)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_216 = Coupling(name = 'GC_216',
                  value = '-((cw*ee**2*complex(0,1)*I95x33)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_217 = Coupling(name = 'GC_217',
                  value = '-((cw*ee**2*complex(0,1)*I95x36)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)))',
                  order = {'QED':2})

GC_218 = Coupling(name = 'GC_218',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_219 = Coupling(name = 'GC_219',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_220 = Coupling(name = 'GC_220',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_221 = Coupling(name = 'GC_221',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_222 = Coupling(name = 'GC_222',
                  value = '(cw*ee*complex(0,1)*NN5x1*Rd4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_223 = Coupling(name = 'GC_223',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_224 = Coupling(name = 'GC_224',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_225 = Coupling(name = 'GC_225',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_226 = Coupling(name = 'GC_226',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_227 = Coupling(name = 'GC_227',
                  value = '(cw*ee*complex(0,1)*NN5x1*Rd5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_228 = Coupling(name = 'GC_228',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_229 = Coupling(name = 'GC_229',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_230 = Coupling(name = 'GC_230',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_231 = Coupling(name = 'GC_231',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_232 = Coupling(name = 'GC_232',
                  value = '(cw*ee*complex(0,1)*NN5x1*Rl4x4*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_233 = Coupling(name = 'GC_233',
                  value = '(cw*ee*complex(0,1)*NN1x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_234 = Coupling(name = 'GC_234',
                  value = '(cw*ee*complex(0,1)*NN2x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_235 = Coupling(name = 'GC_235',
                  value = '(cw*ee*complex(0,1)*NN3x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_236 = Coupling(name = 'GC_236',
                  value = '(cw*ee*complex(0,1)*NN4x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_237 = Coupling(name = 'GC_237',
                  value = '(cw*ee*complex(0,1)*NN5x1*Rl5x5*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_238 = Coupling(name = 'GC_238',
                  value = '(-2*cw*ee*complex(0,1)*NN1x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_239 = Coupling(name = 'GC_239',
                  value = '(-2*cw*ee*complex(0,1)*NN2x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_240 = Coupling(name = 'GC_240',
                  value = '(-2*cw*ee*complex(0,1)*NN3x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_241 = Coupling(name = 'GC_241',
                  value = '(-2*cw*ee*complex(0,1)*NN4x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_242 = Coupling(name = 'GC_242',
                  value = '(-2*cw*ee*complex(0,1)*NN5x1*Ru4x4*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_243 = Coupling(name = 'GC_243',
                  value = '(-2*cw*ee*complex(0,1)*NN1x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_244 = Coupling(name = 'GC_244',
                  value = '(-2*cw*ee*complex(0,1)*NN2x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_245 = Coupling(name = 'GC_245',
                  value = '(-2*cw*ee*complex(0,1)*NN3x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_246 = Coupling(name = 'GC_246',
                  value = '(-2*cw*ee*complex(0,1)*NN4x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_247 = Coupling(name = 'GC_247',
                  value = '(-2*cw*ee*complex(0,1)*NN5x1*Ru5x5*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_248 = Coupling(name = 'GC_248',
                  value = '-(ee**2*complex(0,1))/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_249 = Coupling(name = 'GC_249',
                  value = '-(cw*ee*complex(0,1))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_250 = Coupling(name = 'GC_250',
                  value = '(cw*ee*complex(0,1))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_251 = Coupling(name = 'GC_251',
                  value = '-(cw*ee*complex(0,1)*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_252 = Coupling(name = 'GC_252',
                  value = '(2*cw*ee*complex(0,1)*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_253 = Coupling(name = 'GC_253',
                  value = '-((cw*ee*complex(0,1)*sw)/((-1 + sw)*(1 + sw)))',
                  order = {'QED':1})

GC_254 = Coupling(name = 'GC_254',
                  value = '(cw*ee*complex(0,1)*I100x44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_255 = Coupling(name = 'GC_255',
                  value = '(cw*ee*complex(0,1)*I100x55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_256 = Coupling(name = 'GC_256',
                  value = '(cw*ee*complex(0,1)*I101x44*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_257 = Coupling(name = 'GC_257',
                  value = '(cw*ee*complex(0,1)*I101x55*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_258 = Coupling(name = 'GC_258',
                  value = '(-2*cw*ee*complex(0,1)*I102x44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_259 = Coupling(name = 'GC_259',
                  value = '(-2*cw*ee*complex(0,1)*I102x55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_260 = Coupling(name = 'GC_260',
                  value = '(2*cw*ee**2*complex(0,1)*I20x44*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_261 = Coupling(name = 'GC_261',
                  value = '(2*cw*ee**2*complex(0,1)*I20x55*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_262 = Coupling(name = 'GC_262',
                  value = '(8*cw*ee**2*complex(0,1)*I40x44*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_263 = Coupling(name = 'GC_263',
                  value = '(4*cw*ee*complex(0,1)*G*I40x44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_264 = Coupling(name = 'GC_264',
                  value = '(8*cw*ee**2*complex(0,1)*I40x55*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_265 = Coupling(name = 'GC_265',
                  value = '(4*cw*ee*complex(0,1)*G*I40x55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_266 = Coupling(name = 'GC_266',
                  value = '(2*cw*ee**2*complex(0,1)*I6x44*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_267 = Coupling(name = 'GC_267',
                  value = '(-2*cw*ee*complex(0,1)*G*I6x44*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_268 = Coupling(name = 'GC_268',
                  value = '(2*cw*ee**2*complex(0,1)*I6x55*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_269 = Coupling(name = 'GC_269',
                  value = '(-2*cw*ee*complex(0,1)*G*I6x55*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_270 = Coupling(name = 'GC_270',
                  value = '(-2*ee**2*complex(0,1)*I100x44*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_271 = Coupling(name = 'GC_271',
                  value = '(-2*ee**2*complex(0,1)*I100x55*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_272 = Coupling(name = 'GC_272',
                  value = '(-2*ee**2*complex(0,1)*I101x44*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_273 = Coupling(name = 'GC_273',
                  value = '(-2*ee**2*complex(0,1)*I101x55*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_274 = Coupling(name = 'GC_274',
                  value = '(-8*ee**2*complex(0,1)*I102x44*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_275 = Coupling(name = 'GC_275',
                  value = '(-8*ee**2*complex(0,1)*I102x55*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_276 = Coupling(name = 'GC_276',
                  value = '-((cw*ee**2*complex(0,1)*I19x11)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x11*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_277 = Coupling(name = 'GC_277',
                  value = '-((cw*ee**2*complex(0,1)*I19x22)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x22*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_278 = Coupling(name = 'GC_278',
                  value = '-((cw*ee**2*complex(0,1)*I19x33)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x33*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I20x33*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_279 = Coupling(name = 'GC_279',
                  value = '-((cw*ee**2*complex(0,1)*I19x36)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x36*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I20x36*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_280 = Coupling(name = 'GC_280',
                  value = '-((cw*ee**2*complex(0,1)*I19x63)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x63*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I20x63*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_281 = Coupling(name = 'GC_281',
                  value = '-((cw*ee**2*complex(0,1)*I19x66)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*I19x66*sw)/((-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I20x66*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_282 = Coupling(name = 'GC_282',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x11)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x11*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_283 = Coupling(name = 'GC_283',
                  value = '-((cw*ee*complex(0,1)*G*I39x11)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_284 = Coupling(name = 'GC_284',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x22)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x22*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_285 = Coupling(name = 'GC_285',
                  value = '-((cw*ee*complex(0,1)*G*I39x22)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_286 = Coupling(name = 'GC_286',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x33)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x33*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I40x33*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_287 = Coupling(name = 'GC_287',
                  value = '-((cw*ee*complex(0,1)*G*I39x33)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x33*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I40x33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_288 = Coupling(name = 'GC_288',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x36)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x36*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I40x36*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_289 = Coupling(name = 'GC_289',
                  value = '-((cw*ee*complex(0,1)*G*I39x36)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x36*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I40x36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_290 = Coupling(name = 'GC_290',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x63)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x63*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I40x63*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_291 = Coupling(name = 'GC_291',
                  value = '-((cw*ee*complex(0,1)*G*I39x63)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x63*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I40x63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_292 = Coupling(name = 'GC_292',
                  value = '(-2*cw*ee**2*complex(0,1)*I39x66)/(3.*(-1 + sw)*sw*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I39x66*sw)/(9.*(-1 + sw)*(1 + sw)) + (8*cw*ee**2*complex(0,1)*I40x66*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_293 = Coupling(name = 'GC_293',
                  value = '-((cw*ee*complex(0,1)*G*I39x66)/((-1 + sw)*sw*(1 + sw))) + (4*cw*ee*complex(0,1)*G*I39x66*sw)/(3.*(-1 + sw)*(1 + sw)) + (4*cw*ee*complex(0,1)*G*I40x66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_294 = Coupling(name = 'GC_294',
                  value = '-(cw*ee**2*complex(0,1)*I5x11)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x11*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_295 = Coupling(name = 'GC_295',
                  value = '(cw*ee*complex(0,1)*G*I5x11)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_296 = Coupling(name = 'GC_296',
                  value = '-(cw*ee**2*complex(0,1)*I5x22)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x22*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_297 = Coupling(name = 'GC_297',
                  value = '(cw*ee*complex(0,1)*G*I5x22)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_298 = Coupling(name = 'GC_298',
                  value = '-(cw*ee**2*complex(0,1)*I5x33)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x33*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I6x33*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_299 = Coupling(name = 'GC_299',
                  value = '(cw*ee*complex(0,1)*G*I5x33)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x33*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I6x33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_300 = Coupling(name = 'GC_300',
                  value = '-(cw*ee**2*complex(0,1)*I5x36)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x36*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I6x36*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_301 = Coupling(name = 'GC_301',
                  value = '(cw*ee*complex(0,1)*G*I5x36)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x36*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I6x36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_302 = Coupling(name = 'GC_302',
                  value = '-(cw*ee**2*complex(0,1)*I5x63)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x63*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I6x63*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_303 = Coupling(name = 'GC_303',
                  value = '(cw*ee*complex(0,1)*G*I5x63)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x63*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I6x63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_304 = Coupling(name = 'GC_304',
                  value = '-(cw*ee**2*complex(0,1)*I5x66)/(3.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I5x66*sw)/(9.*(-1 + sw)*(1 + sw)) + (2*cw*ee**2*complex(0,1)*I6x66*sw)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_305 = Coupling(name = 'GC_305',
                  value = '(cw*ee*complex(0,1)*G*I5x66)/((-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I5x66*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*G*I6x66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QCD':1,'QED':1})

GC_306 = Coupling(name = 'GC_306',
                  value = '-(cw*ee*complex(0,1)*I96x11)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I96x11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_307 = Coupling(name = 'GC_307',
                  value = '-(cw*ee*complex(0,1)*I96x22)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I96x22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_308 = Coupling(name = 'GC_308',
                  value = '-(cw*ee*complex(0,1)*I96x33)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I100x33*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I96x33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_309 = Coupling(name = 'GC_309',
                  value = '-(cw*ee*complex(0,1)*I96x36)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I100x36*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I96x36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_310 = Coupling(name = 'GC_310',
                  value = '(cw*ee*complex(0,1)*I96x63)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*I100x63*sw)/(3.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*I96x63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_311 = Coupling(name = 'GC_311',
                  value = '-(cw*ee*complex(0,1)*I96x66)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I100x66*sw)/(3.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I96x66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_312 = Coupling(name = 'GC_312',
                  value = '-(cw*ee*complex(0,1)*I97x11)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97x11*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_313 = Coupling(name = 'GC_313',
                  value = '-(cw*ee*complex(0,1)*I97x22)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I97x22*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_314 = Coupling(name = 'GC_314',
                  value = '-(cw*ee*complex(0,1)*I97x33)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I101x33*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I97x33*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_315 = Coupling(name = 'GC_315',
                  value = '-(cw*ee*complex(0,1)*I97x36)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I101x36*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I97x36*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_316 = Coupling(name = 'GC_316',
                  value = '(cw*ee*complex(0,1)*I97x63)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*I101x63*sw)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*I97x63*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_317 = Coupling(name = 'GC_317',
                  value = '-(cw*ee*complex(0,1)*I97x66)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*I101x66*sw)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*I97x66*sw)/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_318 = Coupling(name = 'GC_318',
                  value = '(cw*ee*complex(0,1)*I98x11)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I98x11*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_319 = Coupling(name = 'GC_319',
                  value = '(cw*ee*complex(0,1)*I98x22)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I98x22*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_320 = Coupling(name = 'GC_320',
                  value = '(cw*ee*complex(0,1)*I98x33)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I102x33*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I98x33*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_321 = Coupling(name = 'GC_321',
                  value = '(cw*ee*complex(0,1)*I98x36)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I102x36*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I98x36*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_322 = Coupling(name = 'GC_322',
                  value = '-(cw*ee*complex(0,1)*I98x63)/(2.*(-1 + sw)*sw*(1 + sw)) + (2*cw*ee*complex(0,1)*I102x63*sw)/(3.*(-1 + sw)*(1 + sw)) + (2*cw*ee*complex(0,1)*I98x63*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_323 = Coupling(name = 'GC_323',
                  value = '(cw*ee*complex(0,1)*I98x66)/(2.*(-1 + sw)*sw*(1 + sw)) - (2*cw*ee*complex(0,1)*I102x66*sw)/(3.*(-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*I98x66*sw)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_324 = Coupling(name = 'GC_324',
                  value = '(2*ee**2*complex(0,1)*I96x11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I96x11*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_325 = Coupling(name = 'GC_325',
                  value = '(2*ee**2*complex(0,1)*I96x22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I96x22*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_326 = Coupling(name = 'GC_326',
                  value = '(2*ee**2*complex(0,1)*I96x33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I100x33*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I96x33*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_327 = Coupling(name = 'GC_327',
                  value = '(2*ee**2*complex(0,1)*I96x36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I100x36*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I96x36*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_328 = Coupling(name = 'GC_328',
                  value = '(2*ee**2*complex(0,1)*I96x63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I100x63*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I96x63*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_329 = Coupling(name = 'GC_329',
                  value = '(2*ee**2*complex(0,1)*I96x66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I96x66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I100x66*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I96x66*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_330 = Coupling(name = 'GC_330',
                  value = '(2*ee**2*complex(0,1)*I97x11)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97x11*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_331 = Coupling(name = 'GC_331',
                  value = '(2*ee**2*complex(0,1)*I97x22)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I97x22*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_332 = Coupling(name = 'GC_332',
                  value = '(2*ee**2*complex(0,1)*I97x33)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I101x33*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I97x33*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_333 = Coupling(name = 'GC_333',
                  value = '(2*ee**2*complex(0,1)*I97x36)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I101x36*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I97x36*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_334 = Coupling(name = 'GC_334',
                  value = '(2*ee**2*complex(0,1)*I97x63)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I101x63*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I97x63*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_335 = Coupling(name = 'GC_335',
                  value = '(2*ee**2*complex(0,1)*I97x66)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I97x66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*I101x66*sw**2)/((-1 + sw)*(1 + sw)) - (2*ee**2*complex(0,1)*I97x66*sw**2)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_336 = Coupling(name = 'GC_336',
                  value = '(4*ee**2*complex(0,1)*I98x11)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x11)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I98x11*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_337 = Coupling(name = 'GC_337',
                  value = '(4*ee**2*complex(0,1)*I98x22)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x22)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I98x22*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_338 = Coupling(name = 'GC_338',
                  value = '(4*ee**2*complex(0,1)*I98x33)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x33)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I102x33*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I98x33*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_339 = Coupling(name = 'GC_339',
                  value = '(4*ee**2*complex(0,1)*I98x36)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x36)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I102x36*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I98x36*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_340 = Coupling(name = 'GC_340',
                  value = '(4*ee**2*complex(0,1)*I98x63)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x63)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I102x63*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I98x63*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_341 = Coupling(name = 'GC_341',
                  value = '(4*ee**2*complex(0,1)*I98x66)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I98x66)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (8*ee**2*complex(0,1)*I102x66*sw**2)/(9.*(-1 + sw)*(1 + sw)) - (8*ee**2*complex(0,1)*I98x66*sw**2)/(9.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_342 = Coupling(name = 'GC_342',
                  value = '(complex(0,1)*I10x33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_343 = Coupling(name = 'GC_343',
                  value = '(complex(0,1)*I10x36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_344 = Coupling(name = 'GC_344',
                  value = '(complex(0,1)*I21x33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_345 = Coupling(name = 'GC_345',
                  value = '(complex(0,1)*I21x36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_346 = Coupling(name = 'GC_346',
                  value = '(complex(0,1)*I49x33*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x33*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN1x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_347 = Coupling(name = 'GC_347',
                  value = '(complex(0,1)*I49x36*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x36*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN1x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_348 = Coupling(name = 'GC_348',
                  value = '(complex(0,1)*I10x33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_349 = Coupling(name = 'GC_349',
                  value = '(complex(0,1)*I10x36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_350 = Coupling(name = 'GC_350',
                  value = '(complex(0,1)*I21x33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_351 = Coupling(name = 'GC_351',
                  value = '(complex(0,1)*I21x36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_352 = Coupling(name = 'GC_352',
                  value = '(complex(0,1)*I49x33*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x33*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN2x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_353 = Coupling(name = 'GC_353',
                  value = '(complex(0,1)*I49x36*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x36*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN2x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_354 = Coupling(name = 'GC_354',
                  value = '(complex(0,1)*I10x33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_355 = Coupling(name = 'GC_355',
                  value = '(complex(0,1)*I10x36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_356 = Coupling(name = 'GC_356',
                  value = '(complex(0,1)*I21x33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_357 = Coupling(name = 'GC_357',
                  value = '(complex(0,1)*I21x36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_358 = Coupling(name = 'GC_358',
                  value = '(complex(0,1)*I49x33*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x33*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN3x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_359 = Coupling(name = 'GC_359',
                  value = '(complex(0,1)*I49x36*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x36*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN3x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_360 = Coupling(name = 'GC_360',
                  value = '(complex(0,1)*I10x33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_361 = Coupling(name = 'GC_361',
                  value = '(complex(0,1)*I10x36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_362 = Coupling(name = 'GC_362',
                  value = '(complex(0,1)*I21x33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_363 = Coupling(name = 'GC_363',
                  value = '(complex(0,1)*I21x36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_364 = Coupling(name = 'GC_364',
                  value = '(complex(0,1)*I49x33*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x33*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN4x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_365 = Coupling(name = 'GC_365',
                  value = '(complex(0,1)*I49x36*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x36*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN4x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_366 = Coupling(name = 'GC_366',
                  value = '(complex(0,1)*I10x33*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x33*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*Rd3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_367 = Coupling(name = 'GC_367',
                  value = '(complex(0,1)*I10x36*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I10x36*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*Rd6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_368 = Coupling(name = 'GC_368',
                  value = '(complex(0,1)*I21x33*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x33*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*Rl3x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_369 = Coupling(name = 'GC_369',
                  value = '(complex(0,1)*I21x36*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I21x36*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*Rl6x6*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_370 = Coupling(name = 'GC_370',
                  value = '(complex(0,1)*I49x33*NN5x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x33*NN5x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN5x1*Ru3x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_371 = Coupling(name = 'GC_371',
                  value = '(complex(0,1)*I49x36*NN5x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I49x36*NN5x4*sw**2)/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*NN5x1*Ru6x6*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_372 = Coupling(name = 'GC_372',
                  value = '(ee**2*complex(0,1)*UP1x1**2)/(2.*sw**2) + (ee**2*complex(0,1)*UP1x2**2)/(2.*sw**2)',
                  order = {'QED':2})

GC_373 = Coupling(name = 'GC_373',
                  value = '-(ee**2*complex(0,1)*UP1x1**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*UP1x2**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_374 = Coupling(name = 'GC_374',
                  value = '(cw*ee*NN1x1*NN1x3*UP1x1)/((-1 + sw)*(1 + sw)) - (ee*NN1x2*NN1x3*UP1x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN1x3*sw*UP1x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN1x4*UP1x2)/((-1 + sw)*(1 + sw)) + (ee*NN1x2*NN1x4*UP1x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN1x4*sw*UP1x2)/((-1 + sw)*(1 + sw)) - (NMl*NN1x4*NN1x5*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x4*NN1x5*sw**2*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN1x3*NN1x5*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x3*NN1x5*sw**2*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN1x3*NN1x4*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN1x5**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x3*NN1x4*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5**2*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_375 = Coupling(name = 'GC_375',
                  value = '(cw*ee*NN1x3*NN2x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN2x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN2x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN2x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN2x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN2x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN2x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN2x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN2x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN2x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN2x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN2x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN2x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN2x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN2x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN2x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN2x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN2x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN2x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN2x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN2x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN2x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN2x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN2x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN2x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN2x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_376 = Coupling(name = 'GC_376',
                  value = '(cw*ee*NN2x1*NN2x3*UP1x1)/((-1 + sw)*(1 + sw)) - (ee*NN2x2*NN2x3*UP1x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN2x3*sw*UP1x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN2x4*UP1x2)/((-1 + sw)*(1 + sw)) + (ee*NN2x2*NN2x4*UP1x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN2x4*sw*UP1x2)/((-1 + sw)*(1 + sw)) - (NMl*NN2x4*NN2x5*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x4*NN2x5*sw**2*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN2x3*NN2x5*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x3*NN2x5*sw**2*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN2x3*NN2x4*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN2x5**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x3*NN2x4*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5**2*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_377 = Coupling(name = 'GC_377',
                  value = '(cw*ee*NN1x3*NN3x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN3x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN3x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN3x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN3x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN3x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN3x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN3x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN3x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN3x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN3x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN3x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN3x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN3x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN3x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN3x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN3x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN3x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN3x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN3x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN3x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN3x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN3x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN3x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN3x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN3x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_378 = Coupling(name = 'GC_378',
                  value = '(cw*ee*NN2x3*NN3x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN3x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN3x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN3x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN3x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN3x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN3x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN3x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN3x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN3x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN3x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN3x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN3x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN3x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN3x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN3x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN3x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN3x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN3x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN3x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN3x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN3x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN3x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN3x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN3x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN3x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_379 = Coupling(name = 'GC_379',
                  value = '(cw*ee*NN3x1*NN3x3*UP1x1)/((-1 + sw)*(1 + sw)) - (ee*NN3x2*NN3x3*UP1x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN3x3*sw*UP1x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN3x4*UP1x2)/((-1 + sw)*(1 + sw)) + (ee*NN3x2*NN3x4*UP1x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN3x4*sw*UP1x2)/((-1 + sw)*(1 + sw)) - (NMl*NN3x4*NN3x5*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x4*NN3x5*sw**2*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN3x3*NN3x5*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x3*NN3x5*sw**2*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN3x3*NN3x4*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN3x5**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x3*NN3x4*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5**2*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_380 = Coupling(name = 'GC_380',
                  value = '(cw*ee*NN1x3*NN4x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN4x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN4x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN4x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN4x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN4x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN4x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN4x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN4x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN4x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN4x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN4x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN4x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN4x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN4x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN4x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN4x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN4x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN4x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN4x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN4x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN4x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN4x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN4x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN4x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN4x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_381 = Coupling(name = 'GC_381',
                  value = '(cw*ee*NN2x3*NN4x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN4x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN4x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN4x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN4x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN4x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN4x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN4x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN4x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN4x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN4x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN4x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN4x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN4x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN4x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN4x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN4x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN4x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN4x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN4x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN4x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN4x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN4x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN4x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN4x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN4x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_382 = Coupling(name = 'GC_382',
                  value = '(cw*ee*NN3x3*NN4x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN4x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x3*NN4x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN4x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x3*NN4x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x2*NN4x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x4*NN4x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN4x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN4x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN4x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN4x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN4x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN3x5*NN4x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN4x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN4x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN4x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x5*NN4x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN4x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN4x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN4x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN4x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN4x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN4x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN4x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN3x5*NN4x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5*NN4x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_383 = Coupling(name = 'GC_383',
                  value = '(cw*ee*NN4x1*NN4x3*UP1x1)/((-1 + sw)*(1 + sw)) - (ee*NN4x2*NN4x3*UP1x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN4x2*NN4x3*sw*UP1x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN4x1*NN4x4*UP1x2)/((-1 + sw)*(1 + sw)) + (ee*NN4x2*NN4x4*UP1x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN4x4*sw*UP1x2)/((-1 + sw)*(1 + sw)) - (NMl*NN4x4*NN4x5*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x4*NN4x5*sw**2*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN4x3*NN4x5*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x3*NN4x5*sw**2*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN4x3*NN4x4*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN4x5**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x3*NN4x4*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN4x5**2*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_384 = Coupling(name = 'GC_384',
                  value = '(cw*ee*NN1x3*NN5x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN5x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN5x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN5x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN5x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN5x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN5x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN5x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN5x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN5x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN5x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN5x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN5x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN5x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN5x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN5x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN5x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN5x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN5x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN5x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN5x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN5x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN5x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN5x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN5x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN5x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_385 = Coupling(name = 'GC_385',
                  value = '(cw*ee*NN2x3*NN5x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN5x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN5x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN5x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN5x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN5x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN5x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN5x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN5x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN5x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN5x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN5x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN5x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN5x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN5x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN5x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN5x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN5x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN5x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN5x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN5x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN5x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN5x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN5x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN5x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN5x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_386 = Coupling(name = 'GC_386',
                  value = '(cw*ee*NN3x3*NN5x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN5x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x3*NN5x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN5x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x3*NN5x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x2*NN5x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x4*NN5x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN5x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN5x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN5x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN5x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN5x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN3x5*NN5x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN5x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN5x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN5x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x5*NN5x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN5x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN5x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN5x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN5x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN5x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN5x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN5x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN3x5*NN5x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5*NN5x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_387 = Coupling(name = 'GC_387',
                  value = '(cw*ee*NN4x3*NN5x1*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN4x1*NN5x3*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN4x3*NN5x2*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN5x3*UP1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN4x3*NN5x2*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN4x2*NN5x3*sw*UP1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN4x4*NN5x1*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN4x1*NN5x4*UP1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN4x4*NN5x2*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN4x2*NN5x4*UP1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN4x4*NN5x2*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN4x2*NN5x4*sw*UP1x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN4x5*NN5x4*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x4*NN5x5*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x5*NN5x4*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x4*NN5x5*sw**2*UP1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x5*NN5x3*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x3*NN5x5*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x5*NN5x3*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x3*NN5x5*sw**2*UP1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x4*NN5x3*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x3*NN5x4*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x4*NN5x3*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x3*NN5x4*sw**2*UP1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN4x5*NN5x5*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN4x5*NN5x5*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_388 = Coupling(name = 'GC_388',
                  value = '(cw*ee*NN5x1*NN5x3*UP1x1)/((-1 + sw)*(1 + sw)) - (ee*NN5x2*NN5x3*UP1x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN5x2*NN5x3*sw*UP1x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN5x1*NN5x4*UP1x2)/((-1 + sw)*(1 + sw)) + (ee*NN5x2*NN5x4*UP1x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN5x2*NN5x4*sw*UP1x2)/((-1 + sw)*(1 + sw)) - (NMl*NN5x4*NN5x5*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x4*NN5x5*sw**2*UP1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN5x3*NN5x5*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x3*NN5x5*sw**2*UP1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN5x3*NN5x4*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN5x5**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x3*NN5x4*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN5x5**2*sw**2*UP1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_389 = Coupling(name = 'GC_389',
                  value = '(ee**2*complex(0,1)*UP1x1*UP2x1)/(2.*sw**2) + (ee**2*complex(0,1)*UP1x2*UP2x2)/(2.*sw**2)',
                  order = {'QED':2})

GC_390 = Coupling(name = 'GC_390',
                  value = '-(ee**2*complex(0,1)*UP1x1*UP2x1)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*UP1x2*UP2x2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_391 = Coupling(name = 'GC_391',
                  value = '(ee**2*complex(0,1)*UP2x1**2)/(2.*sw**2) + (ee**2*complex(0,1)*UP2x2**2)/(2.*sw**2)',
                  order = {'QED':2})

GC_392 = Coupling(name = 'GC_392',
                  value = '-(ee**2*complex(0,1)*UP2x1**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*UP2x2**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_393 = Coupling(name = 'GC_393',
                  value = '(cw*ee*NN1x1*NN1x3*UP2x1)/((-1 + sw)*(1 + sw)) - (ee*NN1x2*NN1x3*UP2x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN1x3*sw*UP2x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN1x4*UP2x2)/((-1 + sw)*(1 + sw)) + (ee*NN1x2*NN1x4*UP2x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN1x4*sw*UP2x2)/((-1 + sw)*(1 + sw)) - (NMl*NN1x4*NN1x5*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x4*NN1x5*sw**2*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN1x3*NN1x5*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x3*NN1x5*sw**2*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN1x3*NN1x4*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN1x5**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN1x3*NN1x4*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5**2*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_394 = Coupling(name = 'GC_394',
                  value = '(cw*ee*NN1x3*NN2x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN2x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN2x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN2x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN2x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN2x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN2x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN2x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN2x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN2x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN2x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN2x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN2x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN2x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN2x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN2x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN2x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN2x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN2x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN2x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN2x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN2x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN2x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN2x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN2x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN2x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_395 = Coupling(name = 'GC_395',
                  value = '(cw*ee*NN2x1*NN2x3*UP2x1)/((-1 + sw)*(1 + sw)) - (ee*NN2x2*NN2x3*UP2x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN2x3*sw*UP2x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN2x4*UP2x2)/((-1 + sw)*(1 + sw)) + (ee*NN2x2*NN2x4*UP2x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN2x4*sw*UP2x2)/((-1 + sw)*(1 + sw)) - (NMl*NN2x4*NN2x5*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x4*NN2x5*sw**2*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN2x3*NN2x5*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x3*NN2x5*sw**2*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN2x3*NN2x4*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN2x5**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN2x3*NN2x4*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5**2*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_396 = Coupling(name = 'GC_396',
                  value = '(cw*ee*NN1x3*NN3x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN3x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN3x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN3x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN3x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN3x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN3x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN3x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN3x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN3x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN3x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN3x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN3x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN3x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN3x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN3x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN3x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN3x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN3x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN3x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN3x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN3x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN3x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN3x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN3x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN3x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_397 = Coupling(name = 'GC_397',
                  value = '(cw*ee*NN2x3*NN3x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN3x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN3x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN3x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN3x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN3x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN3x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN3x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN3x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN3x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN3x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN3x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN3x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN3x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN3x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN3x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN3x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN3x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN3x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN3x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN3x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN3x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN3x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN3x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN3x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN3x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_398 = Coupling(name = 'GC_398',
                  value = '(cw*ee*NN3x1*NN3x3*UP2x1)/((-1 + sw)*(1 + sw)) - (ee*NN3x2*NN3x3*UP2x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN3x3*sw*UP2x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN3x4*UP2x2)/((-1 + sw)*(1 + sw)) + (ee*NN3x2*NN3x4*UP2x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN3x4*sw*UP2x2)/((-1 + sw)*(1 + sw)) - (NMl*NN3x4*NN3x5*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x4*NN3x5*sw**2*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN3x3*NN3x5*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x3*NN3x5*sw**2*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN3x3*NN3x4*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN3x5**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN3x3*NN3x4*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5**2*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_399 = Coupling(name = 'GC_399',
                  value = '(cw*ee*NN1x3*NN4x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN4x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN4x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN4x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN4x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN4x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN4x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN4x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN4x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN4x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN4x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN4x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN4x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN4x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN4x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN4x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN4x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN4x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN4x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN4x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN4x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN4x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN4x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN4x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN4x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN4x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_400 = Coupling(name = 'GC_400',
                  value = '(cw*ee*NN2x3*NN4x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN4x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN4x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN4x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN4x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN4x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN4x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN4x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN4x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN4x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN4x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN4x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN4x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN4x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN4x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN4x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN4x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN4x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN4x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN4x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN4x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN4x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN4x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN4x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN4x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN4x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_401 = Coupling(name = 'GC_401',
                  value = '(cw*ee*NN3x3*NN4x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN4x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x3*NN4x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN4x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x3*NN4x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x2*NN4x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x4*NN4x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN4x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN4x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN4x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN4x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN4x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN3x5*NN4x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN4x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN4x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN4x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x5*NN4x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN4x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN4x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN4x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN4x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN4x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN4x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN4x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN3x5*NN4x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5*NN4x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_402 = Coupling(name = 'GC_402',
                  value = '(cw*ee*NN4x1*NN4x3*UP2x1)/((-1 + sw)*(1 + sw)) - (ee*NN4x2*NN4x3*UP2x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN4x2*NN4x3*sw*UP2x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN4x1*NN4x4*UP2x2)/((-1 + sw)*(1 + sw)) + (ee*NN4x2*NN4x4*UP2x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN4x4*sw*UP2x2)/((-1 + sw)*(1 + sw)) - (NMl*NN4x4*NN4x5*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x4*NN4x5*sw**2*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN4x3*NN4x5*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x3*NN4x5*sw**2*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN4x3*NN4x4*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN4x5**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN4x3*NN4x4*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN4x5**2*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_403 = Coupling(name = 'GC_403',
                  value = '(cw*ee*NN1x3*NN5x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN1x1*NN5x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x3*NN5x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x2*NN5x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x3*NN5x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x2*NN5x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x4*NN5x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN1x1*NN5x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN1x4*NN5x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN1x2*NN5x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN1x4*NN5x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN1x2*NN5x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN1x5*NN5x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN5x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN5x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN5x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x5*NN5x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN5x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x5*NN5x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN5x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x4*NN5x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN1x3*NN5x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x4*NN5x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN1x3*NN5x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN1x5*NN5x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN1x5*NN5x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_404 = Coupling(name = 'GC_404',
                  value = '(cw*ee*NN2x3*NN5x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN2x1*NN5x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x3*NN5x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x2*NN5x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x3*NN5x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x2*NN5x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x4*NN5x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN2x1*NN5x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN2x4*NN5x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN2x2*NN5x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN2x4*NN5x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN2x2*NN5x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN2x5*NN5x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN5x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN5x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN5x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x5*NN5x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN5x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x5*NN5x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN5x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x4*NN5x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN2x3*NN5x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x4*NN5x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN2x3*NN5x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN2x5*NN5x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN2x5*NN5x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_405 = Coupling(name = 'GC_405',
                  value = '(cw*ee*NN3x3*NN5x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN3x1*NN5x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x3*NN5x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x2*NN5x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x3*NN5x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x2*NN5x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x4*NN5x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN3x1*NN5x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN3x4*NN5x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN3x2*NN5x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN3x4*NN5x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN3x2*NN5x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN3x5*NN5x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN5x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN5x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN5x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x5*NN5x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN5x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x5*NN5x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN5x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x4*NN5x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN3x3*NN5x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x4*NN5x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN3x3*NN5x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN3x5*NN5x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN3x5*NN5x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_406 = Coupling(name = 'GC_406',
                  value = '(cw*ee*NN4x3*NN5x1*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*NN4x1*NN5x3*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN4x3*NN5x2*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN4x2*NN5x3*UP2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN4x3*NN5x2*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN4x2*NN5x3*sw*UP2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN4x4*NN5x1*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*NN4x1*NN5x4*UP2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*NN4x4*NN5x2*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*NN4x2*NN5x4*UP2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*NN4x4*NN5x2*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*NN4x2*NN5x4*sw*UP2x2)/(2.*(-1 + sw)*(1 + sw)) - (NMl*NN4x5*NN5x4*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x4*NN5x5*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x5*NN5x4*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x4*NN5x5*sw**2*UP2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x5*NN5x3*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x3*NN5x5*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x5*NN5x3*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x3*NN5x5*sw**2*UP2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x4*NN5x3*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*NN4x3*NN5x4*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x4*NN5x3*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*NN4x3*NN5x4*sw**2*UP2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMk*NN4x5*NN5x5*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN4x5*NN5x5*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_407 = Coupling(name = 'GC_407',
                  value = '(cw*ee*NN5x1*NN5x3*UP2x1)/((-1 + sw)*(1 + sw)) - (ee*NN5x2*NN5x3*UP2x1)/((-1 + sw)*sw*(1 + sw)) + (ee*NN5x2*NN5x3*sw*UP2x1)/((-1 + sw)*(1 + sw)) - (cw*ee*NN5x1*NN5x4*UP2x2)/((-1 + sw)*(1 + sw)) + (ee*NN5x2*NN5x4*UP2x2)/((-1 + sw)*sw*(1 + sw)) - (ee*NN5x2*NN5x4*sw*UP2x2)/((-1 + sw)*(1 + sw)) - (NMl*NN5x4*NN5x5*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x4*NN5x5*sw**2*UP2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN5x3*NN5x5*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x3*NN5x5*sw**2*UP2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*NN5x3*NN5x4*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*NN5x5**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*NN5x3*NN5x4*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*NN5x5**2*sw**2*UP2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_408 = Coupling(name = 'GC_408',
                  value = '-(cw*ee*UP1x1*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP1x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_409 = Coupling(name = 'GC_409',
                  value = '-(cw*ee*UP2x1*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP2x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_410 = Coupling(name = 'GC_410',
                  value = '(ee**2*complex(0,1)*US1x1**2)/(2.*sw**2) + (ee**2*complex(0,1)*US1x2**2)/(2.*sw**2)',
                  order = {'QED':2})

GC_411 = Coupling(name = 'GC_411',
                  value = '-(ee**2*complex(0,1)*US1x1**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_412 = Coupling(name = 'GC_412',
                  value = '-((cw*ee*complex(0,1)*NN1x1*NN1x3*US1x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN1x2*NN1x3*US1x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x3*sw*US1x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN1x4*US1x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x4*US1x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x4*sw*US1x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x4*NN1x5*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x4*NN1x5*sw**2*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x5*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x5*sw**2*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x4*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN1x5**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x4*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5**2*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_413 = Coupling(name = 'GC_413',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN2x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN2x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN2x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN2x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN2x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN2x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN2x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN2x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN2x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN2x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN2x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_414 = Coupling(name = 'GC_414',
                  value = '-((cw*ee*complex(0,1)*NN2x1*NN2x3*US1x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN2x2*NN2x3*US1x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x3*sw*US1x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN2x4*US1x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x4*US1x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x4*sw*US1x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x4*NN2x5*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x4*NN2x5*sw**2*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x5*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x5*sw**2*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x4*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN2x5**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x4*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5**2*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_415 = Coupling(name = 'GC_415',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN3x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN3x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN3x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN3x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN3x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN3x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN3x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN3x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN3x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN3x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN3x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_416 = Coupling(name = 'GC_416',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN3x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN3x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN3x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN3x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN3x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN3x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN3x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN3x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN3x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN3x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN3x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_417 = Coupling(name = 'GC_417',
                  value = '-((cw*ee*complex(0,1)*NN3x1*NN3x3*US1x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN3x2*NN3x3*US1x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x3*sw*US1x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN3x4*US1x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x4*US1x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x4*sw*US1x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x4*NN3x5*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x4*NN3x5*sw**2*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x5*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x5*sw**2*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x4*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN3x5**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x4*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5**2*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_418 = Coupling(name = 'GC_418',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN4x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN4x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN4x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN4x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN4x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN4x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN4x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN4x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN4x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN4x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN4x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_419 = Coupling(name = 'GC_419',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN4x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN4x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN4x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN4x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN4x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN4x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN4x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN4x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN4x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN4x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN4x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_420 = Coupling(name = 'GC_420',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN4x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN4x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN4x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN4x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN4x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN4x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN4x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN4x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN4x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN4x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN4x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_421 = Coupling(name = 'GC_421',
                  value = '-((cw*ee*complex(0,1)*NN4x1*NN4x3*US1x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN4x2*NN4x3*US1x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x3*sw*US1x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN4x4*US1x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x4*US1x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x4*sw*US1x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x4*NN4x5*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x4*NN4x5*sw**2*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x5*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x5*sw**2*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x4*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN4x5**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x4*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5**2*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_422 = Coupling(name = 'GC_422',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN5x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN5x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN5x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN5x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN5x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN5x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN5x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN5x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN5x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN5x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN5x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN5x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_423 = Coupling(name = 'GC_423',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN5x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN5x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN5x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN5x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN5x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN5x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN5x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN5x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN5x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN5x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN5x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN5x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_424 = Coupling(name = 'GC_424',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN5x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN5x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN5x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN5x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN5x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN5x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN5x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN5x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN5x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN5x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN5x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN5x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_425 = Coupling(name = 'GC_425',
                  value = '-(cw*ee*complex(0,1)*NN4x3*NN5x1*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*NN5x3*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x3*NN5x2*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x3*US1x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*NN5x2*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x3*sw*US1x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*NN5x1*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN5x4*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x4*NN5x2*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x4*US1x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x4*NN5x2*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x4*sw*US1x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*NN5x4*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x5*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x4*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x5*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x5*NN5x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x5*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x3*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x5*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x3*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x4*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x3*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x4*sw**2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN4x5*NN5x5*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5*NN5x5*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_426 = Coupling(name = 'GC_426',
                  value = '-((cw*ee*complex(0,1)*NN5x1*NN5x3*US1x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN5x2*NN5x3*US1x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x3*sw*US1x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*NN5x4*US1x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x4*US1x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN5x2*NN5x4*sw*US1x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x4*NN5x5*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x4*NN5x5*sw**2*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x5*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x5*sw**2*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x4*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN5x5**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x4*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN5x5**2*sw**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_427 = Coupling(name = 'GC_427',
                  value = '-(cw*ee*UP1x1*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP1x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_428 = Coupling(name = 'GC_428',
                  value = '-(cw*ee*UP2x1*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP2x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_429 = Coupling(name = 'GC_429',
                  value = '(ee**2*complex(0,1)*US1x1*US2x1)/(2.*sw**2) + (ee**2*complex(0,1)*US1x2*US2x2)/(2.*sw**2)',
                  order = {'QED':2})

GC_430 = Coupling(name = 'GC_430',
                  value = '-(ee**2*complex(0,1)*US1x1*US2x1)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_431 = Coupling(name = 'GC_431',
                  value = '(ee**2*complex(0,1)*US2x1**2)/(2.*sw**2) + (ee**2*complex(0,1)*US2x2**2)/(2.*sw**2)',
                  order = {'QED':2})

GC_432 = Coupling(name = 'GC_432',
                  value = '-(ee**2*complex(0,1)*US2x1**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_433 = Coupling(name = 'GC_433',
                  value = '-((cw*ee*complex(0,1)*NN1x1*NN1x3*US2x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN1x2*NN1x3*US2x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x3*sw*US2x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN1x4*US2x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x4*US2x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x4*sw*US2x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x4*NN1x5*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x4*NN1x5*sw**2*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x5*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x5*sw**2*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x4*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN1x5**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x4*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5**2*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_434 = Coupling(name = 'GC_434',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN2x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN2x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN2x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN2x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN2x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN2x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN2x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN2x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN2x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN2x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN2x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_435 = Coupling(name = 'GC_435',
                  value = '-((cw*ee*complex(0,1)*NN2x1*NN2x3*US2x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN2x2*NN2x3*US2x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x3*sw*US2x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN2x4*US2x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x4*US2x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x4*sw*US2x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x4*NN2x5*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x4*NN2x5*sw**2*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x5*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x5*sw**2*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x4*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN2x5**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x4*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5**2*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_436 = Coupling(name = 'GC_436',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN3x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN3x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN3x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN3x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN3x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN3x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN3x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN3x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN3x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN3x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN3x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_437 = Coupling(name = 'GC_437',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN3x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN3x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN3x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN3x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN3x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN3x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN3x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN3x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN3x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN3x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN3x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_438 = Coupling(name = 'GC_438',
                  value = '-((cw*ee*complex(0,1)*NN3x1*NN3x3*US2x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN3x2*NN3x3*US2x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x3*sw*US2x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN3x4*US2x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x4*US2x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x4*sw*US2x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x4*NN3x5*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x4*NN3x5*sw**2*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x5*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x5*sw**2*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x4*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN3x5**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x4*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5**2*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_439 = Coupling(name = 'GC_439',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN4x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN4x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN4x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN4x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN4x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN4x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN4x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN4x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN4x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN4x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN4x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_440 = Coupling(name = 'GC_440',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN4x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN4x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN4x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN4x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN4x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN4x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN4x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN4x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN4x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN4x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN4x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_441 = Coupling(name = 'GC_441',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN4x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN4x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN4x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN4x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN4x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN4x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN4x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN4x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN4x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN4x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN4x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_442 = Coupling(name = 'GC_442',
                  value = '-((cw*ee*complex(0,1)*NN4x1*NN4x3*US2x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN4x2*NN4x3*US2x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x3*sw*US2x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN4x4*US2x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x4*US2x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x4*sw*US2x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x4*NN4x5*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x4*NN4x5*sw**2*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x5*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x5*sw**2*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x4*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN4x5**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x4*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5**2*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_443 = Coupling(name = 'GC_443',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN5x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN5x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN5x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN5x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN5x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN5x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN5x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN5x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN5x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN5x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN5x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN5x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_444 = Coupling(name = 'GC_444',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN5x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN5x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN5x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN5x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN5x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN5x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN5x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN5x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN5x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN5x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN5x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN5x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_445 = Coupling(name = 'GC_445',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN5x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN5x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN5x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN5x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN5x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN5x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN5x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN5x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN5x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN5x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN5x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN5x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_446 = Coupling(name = 'GC_446',
                  value = '-(cw*ee*complex(0,1)*NN4x3*NN5x1*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*NN5x3*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x3*NN5x2*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x3*US2x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*NN5x2*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x3*sw*US2x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*NN5x1*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN5x4*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x4*NN5x2*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x4*US2x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x4*NN5x2*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x4*sw*US2x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*NN5x4*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x5*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x4*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x5*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x5*NN5x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x5*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x3*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x5*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x3*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x4*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x3*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x4*sw**2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN4x5*NN5x5*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5*NN5x5*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_447 = Coupling(name = 'GC_447',
                  value = '-((cw*ee*complex(0,1)*NN5x1*NN5x3*US2x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN5x2*NN5x3*US2x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x3*sw*US2x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*NN5x4*US2x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x4*US2x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN5x2*NN5x4*sw*US2x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x4*NN5x5*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x4*NN5x5*sw**2*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x5*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x5*sw**2*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x4*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN5x5**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x4*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN5x5**2*sw**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_448 = Coupling(name = 'GC_448',
                  value = '-(cw*ee*UP1x1*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP1x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_449 = Coupling(name = 'GC_449',
                  value = '-(cw*ee*UP2x1*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*UP2x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_450 = Coupling(name = 'GC_450',
                  value = '(ee**2*complex(0,1)*US1x1*US3x1)/(2.*sw**2) + (ee**2*complex(0,1)*US1x2*US3x2)/(2.*sw**2)',
                  order = {'QED':2})

GC_451 = Coupling(name = 'GC_451',
                  value = '-(ee**2*complex(0,1)*US1x1*US3x1)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US3x2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_452 = Coupling(name = 'GC_452',
                  value = '(ee**2*complex(0,1)*US2x1*US3x1)/(2.*sw**2) + (ee**2*complex(0,1)*US2x2*US3x2)/(2.*sw**2)',
                  order = {'QED':2})

GC_453 = Coupling(name = 'GC_453',
                  value = '-(ee**2*complex(0,1)*US2x1*US3x1)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2*US3x2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_454 = Coupling(name = 'GC_454',
                  value = '(ee**2*complex(0,1)*US3x1**2)/(2.*sw**2) + (ee**2*complex(0,1)*US3x2**2)/(2.*sw**2)',
                  order = {'QED':2})

GC_455 = Coupling(name = 'GC_455',
                  value = '-(ee**2*complex(0,1)*US3x1**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US3x2**2)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':2})

GC_456 = Coupling(name = 'GC_456',
                  value = '-((cw*ee*complex(0,1)*NN1x1*NN1x3*US3x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN1x2*NN1x3*US3x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x3*sw*US3x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN1x4*US3x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN1x4*US3x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN1x4*sw*US3x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x4*NN1x5*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x4*NN1x5*sw**2*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x5*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x5*sw**2*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x3*NN1x4*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN1x5**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN1x3*NN1x4*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5**2*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_457 = Coupling(name = 'GC_457',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN2x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN2x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN2x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN2x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN2x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN2x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN2x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN2x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN2x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN2x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN2x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN2x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN2x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN2x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN2x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN2x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN2x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN2x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN2x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_458 = Coupling(name = 'GC_458',
                  value = '-((cw*ee*complex(0,1)*NN2x1*NN2x3*US3x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN2x2*NN2x3*US3x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x3*sw*US3x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN2x4*US3x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN2x4*US3x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN2x4*sw*US3x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x4*NN2x5*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x4*NN2x5*sw**2*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x5*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x5*sw**2*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x3*NN2x4*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN2x5**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN2x3*NN2x4*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5**2*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_459 = Coupling(name = 'GC_459',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN3x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN3x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN3x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN3x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN3x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN3x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN3x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN3x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN3x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN3x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN3x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN3x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN3x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN3x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN3x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN3x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN3x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN3x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN3x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_460 = Coupling(name = 'GC_460',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN3x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN3x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN3x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN3x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN3x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN3x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN3x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN3x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN3x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN3x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN3x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN3x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN3x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN3x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN3x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN3x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN3x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN3x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN3x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_461 = Coupling(name = 'GC_461',
                  value = '-((cw*ee*complex(0,1)*NN3x1*NN3x3*US3x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN3x2*NN3x3*US3x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x3*sw*US3x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN3x4*US3x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN3x4*US3x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN3x4*sw*US3x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x4*NN3x5*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x4*NN3x5*sw**2*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x5*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x5*sw**2*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x3*NN3x4*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN3x5**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN3x3*NN3x4*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5**2*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_462 = Coupling(name = 'GC_462',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN4x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN4x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN4x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN4x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN4x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN4x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN4x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN4x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN4x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN4x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN4x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN4x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN4x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN4x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN4x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN4x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN4x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN4x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN4x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_463 = Coupling(name = 'GC_463',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN4x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN4x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN4x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN4x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN4x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN4x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN4x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN4x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN4x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN4x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN4x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN4x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN4x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN4x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN4x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN4x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN4x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN4x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN4x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_464 = Coupling(name = 'GC_464',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN4x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN4x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN4x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN4x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN4x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN4x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN4x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN4x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN4x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN4x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN4x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN4x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN4x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN4x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN4x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN4x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN4x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN4x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN4x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_465 = Coupling(name = 'GC_465',
                  value = '-((cw*ee*complex(0,1)*NN4x1*NN4x3*US3x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN4x2*NN4x3*US3x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x3*sw*US3x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN4x4*US3x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN4x4*US3x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN4x4*sw*US3x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x4*NN4x5*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x4*NN4x5*sw**2*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x5*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x5*sw**2*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x3*NN4x4*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN4x5**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN4x3*NN4x4*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5**2*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_466 = Coupling(name = 'GC_466',
                  value = '-(cw*ee*complex(0,1)*NN1x3*NN5x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*NN5x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*NN5x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*NN5x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*NN5x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*NN5x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN1x4*NN5x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x2*NN5x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN1x4*NN5x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x2*NN5x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*NN5x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x5*NN5x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*NN5x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x4*NN5x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN1x3*NN5x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x4*NN5x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x3*NN5x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN1x5*NN5x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN1x5*NN5x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_467 = Coupling(name = 'GC_467',
                  value = '-(cw*ee*complex(0,1)*NN2x3*NN5x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*NN5x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*NN5x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*NN5x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*NN5x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*NN5x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN2x4*NN5x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x2*NN5x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN2x4*NN5x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x2*NN5x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*NN5x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x5*NN5x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*NN5x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x4*NN5x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN2x3*NN5x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x4*NN5x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x3*NN5x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN2x5*NN5x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN2x5*NN5x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_468 = Coupling(name = 'GC_468',
                  value = '-(cw*ee*complex(0,1)*NN3x3*NN5x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*NN5x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*NN5x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*NN5x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*NN5x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*NN5x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN3x4*NN5x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x2*NN5x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN3x4*NN5x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x2*NN5x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*NN5x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x5*NN5x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*NN5x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x4*NN5x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN3x3*NN5x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x4*NN5x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x3*NN5x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN3x5*NN5x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN3x5*NN5x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_469 = Coupling(name = 'GC_469',
                  value = '-(cw*ee*complex(0,1)*NN4x3*NN5x1*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*NN5x3*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x3*NN5x2*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x3*US3x1)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*NN5x2*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x3*sw*US3x1)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*NN5x1*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*NN5x4*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN4x4*NN5x2*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x2*NN5x4*US3x2)/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN4x4*NN5x2*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x2*NN5x4*sw*US3x2)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*NN5x4*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x5*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x4*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x5*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x5*NN5x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x5*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*NN5x3*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x5*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x4*NN5x3*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*NN4x3*NN5x4*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x4*NN5x3*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x3*NN5x4*sw**2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*NN4x5*NN5x5*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN4x5*NN5x5*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_470 = Coupling(name = 'GC_470',
                  value = '-((cw*ee*complex(0,1)*NN5x1*NN5x3*US3x1)/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*NN5x2*NN5x3*US3x1)/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x3*sw*US3x1)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*NN5x4*US3x2)/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*NN5x2*NN5x4*US3x2)/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*NN5x2*NN5x4*sw*US3x2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x4*NN5x5*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x4*NN5x5*sw**2*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x5*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x5*sw**2*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x3*NN5x4*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NN5x5**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*NN5x3*NN5x4*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NN5x5**2*sw**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_471 = Coupling(name = 'GC_471',
                  value = '-((ee*complex(0,1)*I82x11*UU1x1)/sw)',
                  order = {'QED':1})

GC_472 = Coupling(name = 'GC_472',
                  value = '-((ee*complex(0,1)*I82x22*UU1x1)/sw)',
                  order = {'QED':1})

GC_473 = Coupling(name = 'GC_473',
                  value = '-((ee*complex(0,1)*I85x11*UU1x1)/sw)',
                  order = {'QED':1})

GC_474 = Coupling(name = 'GC_474',
                  value = '-((ee*complex(0,1)*I85x22*UU1x1)/sw)',
                  order = {'QED':1})

GC_475 = Coupling(name = 'GC_475',
                  value = 'complex(0,1)*I32x33*UU1x2',
                  order = {'QED':1})

GC_476 = Coupling(name = 'GC_476',
                  value = 'complex(0,1)*I52x33*UU1x2',
                  order = {'QED':1})

GC_477 = Coupling(name = 'GC_477',
                  value = 'complex(0,1)*I52x36*UU1x2',
                  order = {'QED':1})

GC_478 = Coupling(name = 'GC_478',
                  value = '-((ee*complex(0,1)*I82x33*UU1x1)/sw) + complex(0,1)*I83x33*UU1x2',
                  order = {'QED':1})

GC_479 = Coupling(name = 'GC_479',
                  value = '-((ee*complex(0,1)*I82x36*UU1x1)/sw) + complex(0,1)*I83x36*UU1x2',
                  order = {'QED':1})

GC_480 = Coupling(name = 'GC_480',
                  value = '-((ee*complex(0,1)*I85x33*UU1x1)/sw) + complex(0,1)*I86x33*UU1x2',
                  order = {'QED':1})

GC_481 = Coupling(name = 'GC_481',
                  value = '-((ee*complex(0,1)*I85x36*UU1x1)/sw) + complex(0,1)*I86x36*UU1x2',
                  order = {'QED':1})

GC_482 = Coupling(name = 'GC_482',
                  value = '-((ee*complex(0,1)*I82x11*UU2x1)/sw)',
                  order = {'QED':1})

GC_483 = Coupling(name = 'GC_483',
                  value = '-((ee*complex(0,1)*I82x22*UU2x1)/sw)',
                  order = {'QED':1})

GC_484 = Coupling(name = 'GC_484',
                  value = '-((ee*complex(0,1)*I85x11*UU2x1)/sw)',
                  order = {'QED':1})

GC_485 = Coupling(name = 'GC_485',
                  value = '-((ee*complex(0,1)*I85x22*UU2x1)/sw)',
                  order = {'QED':1})

GC_486 = Coupling(name = 'GC_486',
                  value = 'complex(0,1)*I32x33*UU2x2',
                  order = {'QED':1})

GC_487 = Coupling(name = 'GC_487',
                  value = 'complex(0,1)*I52x33*UU2x2',
                  order = {'QED':1})

GC_488 = Coupling(name = 'GC_488',
                  value = 'complex(0,1)*I52x36*UU2x2',
                  order = {'QED':1})

GC_489 = Coupling(name = 'GC_489',
                  value = '-((ee*complex(0,1)*I82x33*UU2x1)/sw) + complex(0,1)*I83x33*UU2x2',
                  order = {'QED':1})

GC_490 = Coupling(name = 'GC_490',
                  value = '-((ee*complex(0,1)*I82x36*UU2x1)/sw) + complex(0,1)*I83x36*UU2x2',
                  order = {'QED':1})

GC_491 = Coupling(name = 'GC_491',
                  value = '-((ee*complex(0,1)*I85x33*UU2x1)/sw) + complex(0,1)*I86x33*UU2x2',
                  order = {'QED':1})

GC_492 = Coupling(name = 'GC_492',
                  value = '-((ee*complex(0,1)*I85x36*UU2x1)/sw) + complex(0,1)*I86x36*UU2x2',
                  order = {'QED':1})

GC_493 = Coupling(name = 'GC_493',
                  value = '-(I72x33*NMl*UP1x3*vd)/2. + (I75x33*NMl*UP1x3*vd)/2. - (I73x33*UP1x2)/cmath.sqrt(2) + (I74x33*UP1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_494 = Coupling(name = 'GC_494',
                  value = '-(I72x36*NMl*UP1x3*vd)/2. + (I75x36*NMl*UP1x3*vd)/2. - (I73x36*UP1x2)/cmath.sqrt(2) + (I74x36*UP1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_495 = Coupling(name = 'GC_495',
                  value = '-(I72x63*NMl*UP1x3*vd)/2. + (I75x63*NMl*UP1x3*vd)/2. - (I73x63*UP1x2)/cmath.sqrt(2) + (I74x63*UP1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_496 = Coupling(name = 'GC_496',
                  value = '-(I72x66*NMl*UP1x3*vd)/2. + (I75x66*NMl*UP1x3*vd)/2. - (I73x66*UP1x2)/cmath.sqrt(2) + (I74x66*UP1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_497 = Coupling(name = 'GC_497',
                  value = '-(I72x33*NMl*UP2x3*vd)/2. + (I75x33*NMl*UP2x3*vd)/2. - (I73x33*UP2x2)/cmath.sqrt(2) + (I74x33*UP2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_498 = Coupling(name = 'GC_498',
                  value = '-(I72x36*NMl*UP2x3*vd)/2. + (I75x36*NMl*UP2x3*vd)/2. - (I73x36*UP2x2)/cmath.sqrt(2) + (I74x36*UP2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_499 = Coupling(name = 'GC_499',
                  value = '-(I72x63*NMl*UP2x3*vd)/2. + (I75x63*NMl*UP2x3*vd)/2. - (I73x63*UP2x2)/cmath.sqrt(2) + (I74x63*UP2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_500 = Coupling(name = 'GC_500',
                  value = '-(I72x66*NMl*UP2x3*vd)/2. + (I75x66*NMl*UP2x3*vd)/2. - (I73x66*UP2x2)/cmath.sqrt(2) + (I74x66*UP2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_501 = Coupling(name = 'GC_501',
                  value = '-(I72x33*NMl*UP1x1*vs)/2. + (I75x33*NMl*UP1x1*vs)/2.',
                  order = {'QED':2})

GC_502 = Coupling(name = 'GC_502',
                  value = '-(I72x36*NMl*UP1x1*vs)/2. + (I75x36*NMl*UP1x1*vs)/2.',
                  order = {'QED':2})

GC_503 = Coupling(name = 'GC_503',
                  value = '-(I72x63*NMl*UP1x1*vs)/2. + (I75x63*NMl*UP1x1*vs)/2.',
                  order = {'QED':2})

GC_504 = Coupling(name = 'GC_504',
                  value = '-(I72x66*NMl*UP1x1*vs)/2. + (I75x66*NMl*UP1x1*vs)/2.',
                  order = {'QED':2})

GC_505 = Coupling(name = 'GC_505',
                  value = '-(I66x33*NMl*UP1x2*vs)/2. + (I67x33*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_506 = Coupling(name = 'GC_506',
                  value = '-(I66x36*NMl*UP1x2*vs)/2. + (I67x36*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_507 = Coupling(name = 'GC_507',
                  value = '-(I66x63*NMl*UP1x2*vs)/2. + (I67x63*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_508 = Coupling(name = 'GC_508',
                  value = '-(I66x66*NMl*UP1x2*vs)/2. + (I67x66*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_509 = Coupling(name = 'GC_509',
                  value = '-(I70x33*NMl*UP1x2*vs)/2. + (I71x33*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_510 = Coupling(name = 'GC_510',
                  value = '-(I70x36*NMl*UP1x2*vs)/2. + (I71x36*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_511 = Coupling(name = 'GC_511',
                  value = '-(I70x63*NMl*UP1x2*vs)/2. + (I71x63*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_512 = Coupling(name = 'GC_512',
                  value = '-(I70x66*NMl*UP1x2*vs)/2. + (I71x66*NMl*UP1x2*vs)/2.',
                  order = {'QED':2})

GC_513 = Coupling(name = 'GC_513',
                  value = '-(I72x33*NMl*UP2x1*vs)/2. + (I75x33*NMl*UP2x1*vs)/2.',
                  order = {'QED':2})

GC_514 = Coupling(name = 'GC_514',
                  value = '-(I72x36*NMl*UP2x1*vs)/2. + (I75x36*NMl*UP2x1*vs)/2.',
                  order = {'QED':2})

GC_515 = Coupling(name = 'GC_515',
                  value = '-(I72x63*NMl*UP2x1*vs)/2. + (I75x63*NMl*UP2x1*vs)/2.',
                  order = {'QED':2})

GC_516 = Coupling(name = 'GC_516',
                  value = '-(I72x66*NMl*UP2x1*vs)/2. + (I75x66*NMl*UP2x1*vs)/2.',
                  order = {'QED':2})

GC_517 = Coupling(name = 'GC_517',
                  value = '-(I66x33*NMl*UP2x2*vs)/2. + (I67x33*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_518 = Coupling(name = 'GC_518',
                  value = '-(I66x36*NMl*UP2x2*vs)/2. + (I67x36*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_519 = Coupling(name = 'GC_519',
                  value = '-(I66x63*NMl*UP2x2*vs)/2. + (I67x63*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_520 = Coupling(name = 'GC_520',
                  value = '-(I66x66*NMl*UP2x2*vs)/2. + (I67x66*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_521 = Coupling(name = 'GC_521',
                  value = '-(I70x33*NMl*UP2x2*vs)/2. + (I71x33*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_522 = Coupling(name = 'GC_522',
                  value = '-(I70x36*NMl*UP2x2*vs)/2. + (I71x36*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_523 = Coupling(name = 'GC_523',
                  value = '-(I70x63*NMl*UP2x2*vs)/2. + (I71x63*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_524 = Coupling(name = 'GC_524',
                  value = '-(I70x66*NMl*UP2x2*vs)/2. + (I71x66*NMl*UP2x2*vs)/2.',
                  order = {'QED':2})

GC_525 = Coupling(name = 'GC_525',
                  value = '-(complex(0,1)*I72x33*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_526 = Coupling(name = 'GC_526',
                  value = '-(complex(0,1)*I72x36*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_527 = Coupling(name = 'GC_527',
                  value = '-(complex(0,1)*I72x63*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_528 = Coupling(name = 'GC_528',
                  value = '-(complex(0,1)*I72x66*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US1x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_529 = Coupling(name = 'GC_529',
                  value = '-(complex(0,1)*I66x33*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_530 = Coupling(name = 'GC_530',
                  value = '-(complex(0,1)*I66x36*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_531 = Coupling(name = 'GC_531',
                  value = '-(complex(0,1)*I66x63*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_532 = Coupling(name = 'GC_532',
                  value = '-(complex(0,1)*I66x66*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_533 = Coupling(name = 'GC_533',
                  value = '-(complex(0,1)*I70x33*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_534 = Coupling(name = 'GC_534',
                  value = '-(complex(0,1)*I70x36*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_535 = Coupling(name = 'GC_535',
                  value = '-(complex(0,1)*I70x63*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_536 = Coupling(name = 'GC_536',
                  value = '-(complex(0,1)*I70x66*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US1x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_537 = Coupling(name = 'GC_537',
                  value = '(-2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US1x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US1x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP1x1*UP1x2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3**2*US1x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_538 = Coupling(name = 'GC_538',
                  value = '-((complex(0,1)*NMk*NMl*UP1x3*UP2x2*US1x1*vs)/((-1 + sw)*(1 + sw))) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US1x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US1x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US1x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US1x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US1x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US1x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x2*UP2x1*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x1*US1x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x1*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3*UP2x3*US1x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3*UP2x3*US1x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_539 = Coupling(name = 'GC_539',
                  value = '(-2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US1x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US1x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US1x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP2x1*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP2x3**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP2x3**2*US1x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_540 = Coupling(name = 'GC_540',
                  value = '(3*complex(0,1)*NMl**2*US1x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x1**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk*NMl*US1x1*US1x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US1x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x2**2*US1x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3**3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3**3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_541 = Coupling(name = 'GC_541',
                  value = '-(complex(0,1)*I72x33*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_542 = Coupling(name = 'GC_542',
                  value = '-(complex(0,1)*I72x36*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_543 = Coupling(name = 'GC_543',
                  value = '-(complex(0,1)*I72x63*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_544 = Coupling(name = 'GC_544',
                  value = '-(complex(0,1)*I72x66*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US2x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_545 = Coupling(name = 'GC_545',
                  value = '-(complex(0,1)*I66x33*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_546 = Coupling(name = 'GC_546',
                  value = '-(complex(0,1)*I66x36*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_547 = Coupling(name = 'GC_547',
                  value = '-(complex(0,1)*I66x63*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_548 = Coupling(name = 'GC_548',
                  value = '-(complex(0,1)*I66x66*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_549 = Coupling(name = 'GC_549',
                  value = '-(complex(0,1)*I70x33*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_550 = Coupling(name = 'GC_550',
                  value = '-(complex(0,1)*I70x36*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_551 = Coupling(name = 'GC_551',
                  value = '-(complex(0,1)*I70x63*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_552 = Coupling(name = 'GC_552',
                  value = '-(complex(0,1)*I70x66*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US2x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_553 = Coupling(name = 'GC_553',
                  value = '(-2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP1x1*UP1x2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3**2*US2x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_554 = Coupling(name = 'GC_554',
                  value = '-((complex(0,1)*NMk*NMl*UP1x3*UP2x2*US2x1*vs)/((-1 + sw)*(1 + sw))) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US2x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US2x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US2x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x2*UP2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x1*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3*UP2x3*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3*UP2x3*US2x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_555 = Coupling(name = 'GC_555',
                  value = '(-2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US2x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP2x1*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP2x3**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP2x3**2*US2x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_556 = Coupling(name = 'GC_556',
                  value = '(2*complex(0,1)*NMl**2*US1x1*US1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US1x3*US2x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US1x3*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3**2*US2x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_557 = Coupling(name = 'GC_557',
                  value = '(complex(0,1)*NMl**2*US1x3*US2x1**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x1**2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US2x1*US2x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US2x1*US2x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x2**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x2**2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US2x1*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3*US2x3**2*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3*US2x3**2*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_558 = Coupling(name = 'GC_558',
                  value = '(3*complex(0,1)*NMl**2*US2x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x1**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk*NMl*US2x1*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk*NMl*sw**2*US2x1*US2x2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US2x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x2**2*US2x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US2x3**3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US2x3**3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_559 = Coupling(name = 'GC_559',
                  value = '-(complex(0,1)*I72x33*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_560 = Coupling(name = 'GC_560',
                  value = '-(complex(0,1)*I72x36*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_561 = Coupling(name = 'GC_561',
                  value = '-(complex(0,1)*I72x63*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_562 = Coupling(name = 'GC_562',
                  value = '-(complex(0,1)*I72x66*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US3x1*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_563 = Coupling(name = 'GC_563',
                  value = '-(complex(0,1)*I66x33*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_564 = Coupling(name = 'GC_564',
                  value = '-(complex(0,1)*I66x36*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_565 = Coupling(name = 'GC_565',
                  value = '-(complex(0,1)*I66x63*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_566 = Coupling(name = 'GC_566',
                  value = '-(complex(0,1)*I66x66*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_567 = Coupling(name = 'GC_567',
                  value = '-(complex(0,1)*I70x33*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_568 = Coupling(name = 'GC_568',
                  value = '-(complex(0,1)*I70x36*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_569 = Coupling(name = 'GC_569',
                  value = '-(complex(0,1)*I70x63*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_570 = Coupling(name = 'GC_570',
                  value = '-(complex(0,1)*I70x66*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US3x2*vs)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_571 = Coupling(name = 'GC_571',
                  value = '(-2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP1x1*UP1x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3**2*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_572 = Coupling(name = 'GC_572',
                  value = '-((complex(0,1)*NMk*NMl*UP1x3*UP2x2*US3x1*vs)/((-1 + sw)*(1 + sw))) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US3x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x2*UP2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x1*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP1x3*UP2x3*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP1x3*UP2x3*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_573 = Coupling(name = 'GC_573',
                  value = '(-2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*UP2x1*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk**2*UP2x3**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk**2*sw**2*UP2x3**2*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_574 = Coupling(name = 'GC_574',
                  value = '(2*complex(0,1)*NMl**2*US1x1*US1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US1x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US1x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3**2*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_575 = Coupling(name = 'GC_575',
                  value = '(complex(0,1)*NMl**2*US1x3*US2x1*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x1*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x2*US3x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x2*US3x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x2*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x2*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x2*US3x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x2*US3x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x1*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x1*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x2*US2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x2*US2x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x1*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x1*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3*US2x3*US3x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3*US2x3*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_576 = Coupling(name = 'GC_576',
                  value = '(2*complex(0,1)*NMl**2*US2x1*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x2*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x2*US2x3*US3x1*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x1*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x1*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x2*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x2*US2x3*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x1*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x1*US2x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US2x3**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US2x3**2*US3x3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_577 = Coupling(name = 'GC_577',
                  value = '(complex(0,1)*NMl**2*US1x3*US3x1**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US3x1**2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US3x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US3x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US3x2**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US3x2**2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US1x3*US3x3**2*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US1x3*US3x3**2*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_578 = Coupling(name = 'GC_578',
                  value = '(complex(0,1)*NMl**2*US2x3*US3x1**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3*US3x1**2*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x3*US3x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x3*US3x1*US3x2*vs)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x3*US3x2**2*vs)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3*US3x2**2*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x1*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x2*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x2*US3x1*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x2*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x2*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US2x3*US3x3**2*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US2x3*US3x3**2*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_579 = Coupling(name = 'GC_579',
                  value = '(3*complex(0,1)*NMl**2*US3x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x1**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk*NMl*US3x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk*NMl*sw**2*US3x1*US3x2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US3x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x2**2*US3x3*vs)/((-1 + sw)*(1 + sw)) + (6*complex(0,1)*NMk**2*US3x3**3*vs)/((-1 + sw)*(1 + sw)) - (6*complex(0,1)*NMk**2*sw**2*US3x3**3*vs)/((-1 + sw)*(1 + sw))',
                  order = {'QED':2})

GC_580 = Coupling(name = 'GC_580',
                  value = '-(I66x33*NMl*UP1x3*vu)/2. + (I67x33*NMl*UP1x3*vu)/2. - (I64x33*UP1x1)/cmath.sqrt(2) + (I65x33*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_581 = Coupling(name = 'GC_581',
                  value = '-(I66x36*NMl*UP1x3*vu)/2. + (I67x36*NMl*UP1x3*vu)/2. - (I64x36*UP1x1)/cmath.sqrt(2) + (I65x36*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_582 = Coupling(name = 'GC_582',
                  value = '-(I66x63*NMl*UP1x3*vu)/2. + (I67x63*NMl*UP1x3*vu)/2. - (I64x63*UP1x1)/cmath.sqrt(2) + (I65x63*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_583 = Coupling(name = 'GC_583',
                  value = '-(I66x66*NMl*UP1x3*vu)/2. + (I67x66*NMl*UP1x3*vu)/2. - (I64x66*UP1x1)/cmath.sqrt(2) + (I65x66*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_584 = Coupling(name = 'GC_584',
                  value = '-(I70x33*NMl*UP1x3*vu)/2. + (I71x33*NMl*UP1x3*vu)/2. - (I68x33*UP1x1)/cmath.sqrt(2) + (I69x33*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_585 = Coupling(name = 'GC_585',
                  value = '-(I70x36*NMl*UP1x3*vu)/2. + (I71x36*NMl*UP1x3*vu)/2. - (I68x36*UP1x1)/cmath.sqrt(2) + (I69x36*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_586 = Coupling(name = 'GC_586',
                  value = '-(I70x63*NMl*UP1x3*vu)/2. + (I71x63*NMl*UP1x3*vu)/2. - (I68x63*UP1x1)/cmath.sqrt(2) + (I69x63*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_587 = Coupling(name = 'GC_587',
                  value = '-(I70x66*NMl*UP1x3*vu)/2. + (I71x66*NMl*UP1x3*vu)/2. - (I68x66*UP1x1)/cmath.sqrt(2) + (I69x66*UP1x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_588 = Coupling(name = 'GC_588',
                  value = '-(I66x33*NMl*UP2x3*vu)/2. + (I67x33*NMl*UP2x3*vu)/2. - (I64x33*UP2x1)/cmath.sqrt(2) + (I65x33*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_589 = Coupling(name = 'GC_589',
                  value = '-(I66x36*NMl*UP2x3*vu)/2. + (I67x36*NMl*UP2x3*vu)/2. - (I64x36*UP2x1)/cmath.sqrt(2) + (I65x36*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_590 = Coupling(name = 'GC_590',
                  value = '-(I66x63*NMl*UP2x3*vu)/2. + (I67x63*NMl*UP2x3*vu)/2. - (I64x63*UP2x1)/cmath.sqrt(2) + (I65x63*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_591 = Coupling(name = 'GC_591',
                  value = '-(I66x66*NMl*UP2x3*vu)/2. + (I67x66*NMl*UP2x3*vu)/2. - (I64x66*UP2x1)/cmath.sqrt(2) + (I65x66*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_592 = Coupling(name = 'GC_592',
                  value = '-(I70x33*NMl*UP2x3*vu)/2. + (I71x33*NMl*UP2x3*vu)/2. - (I68x33*UP2x1)/cmath.sqrt(2) + (I69x33*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_593 = Coupling(name = 'GC_593',
                  value = '-(I70x36*NMl*UP2x3*vu)/2. + (I71x36*NMl*UP2x3*vu)/2. - (I68x36*UP2x1)/cmath.sqrt(2) + (I69x36*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_594 = Coupling(name = 'GC_594',
                  value = '-(I70x63*NMl*UP2x3*vu)/2. + (I71x63*NMl*UP2x3*vu)/2. - (I68x63*UP2x1)/cmath.sqrt(2) + (I69x63*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_595 = Coupling(name = 'GC_595',
                  value = '-(I70x66*NMl*UP2x3*vu)/2. + (I71x66*NMl*UP2x3*vu)/2. - (I68x66*UP2x1)/cmath.sqrt(2) + (I69x66*UP2x1)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_596 = Coupling(name = 'GC_596',
                  value = '(ee**2*complex(0,1)*US1x1*vd)/(2.*sw**2) + (ee**2*complex(0,1)*US1x2*vu)/(2.*sw**2)',
                  order = {'QED':1})

GC_597 = Coupling(name = 'GC_597',
                  value = '-(ee**2*complex(0,1)*I20x44*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x44*US1x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_598 = Coupling(name = 'GC_598',
                  value = '-(ee**2*complex(0,1)*I20x55*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x55*US1x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_599 = Coupling(name = 'GC_599',
                  value = '(ee**2*complex(0,1)*I40x44*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x44*US1x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_600 = Coupling(name = 'GC_600',
                  value = '(ee**2*complex(0,1)*I40x55*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x55*US1x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_601 = Coupling(name = 'GC_601',
                  value = '-(ee**2*complex(0,1)*I6x44*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x44*US1x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_602 = Coupling(name = 'GC_602',
                  value = '-(ee**2*complex(0,1)*I6x55*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x55*US1x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_603 = Coupling(name = 'GC_603',
                  value = '(ee**2*complex(0,1)*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_604 = Coupling(name = 'GC_604',
                  value = '-(ee**2*complex(0,1)*US1x1*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_605 = Coupling(name = 'GC_605',
                  value = '(ee**2*complex(0,1)*I19x11*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x11*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_606 = Coupling(name = 'GC_606',
                  value = '(ee**2*complex(0,1)*I19x22*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x22*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_607 = Coupling(name = 'GC_607',
                  value = '-(ee**2*complex(0,1)*I39x11*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x11*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_608 = Coupling(name = 'GC_608',
                  value = '-(ee**2*complex(0,1)*I39x22*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x22*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_609 = Coupling(name = 'GC_609',
                  value = '(ee**2*complex(0,1)*I5x11*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x11*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_610 = Coupling(name = 'GC_610',
                  value = '(ee**2*complex(0,1)*I5x22*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x22*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_611 = Coupling(name = 'GC_611',
                  value = '-(ee**2*complex(0,1)*I39x33*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x33*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x33*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x33*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x33*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x33*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x33*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x33*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x33*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x33*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x33*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x33*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x33*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_612 = Coupling(name = 'GC_612',
                  value = '-(ee**2*complex(0,1)*I39x36*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x36*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x36*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x36*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x36*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x36*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x36*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x36*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x36*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x36*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x36*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x36*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x36*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_613 = Coupling(name = 'GC_613',
                  value = '-(ee**2*complex(0,1)*I39x63*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x63*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x63*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x63*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x63*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x63*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x63*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x63*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x63*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x63*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x63*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x63*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x63*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_614 = Coupling(name = 'GC_614',
                  value = '-(ee**2*complex(0,1)*I39x66*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x66*US1x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x66*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US1x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x66*US1x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x66*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x66*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x66*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x66*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x66*sw**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x66*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x66*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x66*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x66*sw**2*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_615 = Coupling(name = 'GC_615',
                  value = '(ee**2*complex(0,1)*I5x33*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x33*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x33*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x33*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x33*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x33*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x33*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x33*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x33*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x33*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x33*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x33*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x33*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_616 = Coupling(name = 'GC_616',
                  value = '(ee**2*complex(0,1)*I5x36*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x36*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x36*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x36*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x36*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x36*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x36*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x36*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x36*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x36*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x36*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x36*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x36*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_617 = Coupling(name = 'GC_617',
                  value = '(ee**2*complex(0,1)*I5x63*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x63*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x63*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x63*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x63*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x63*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x63*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x63*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x63*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x63*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x63*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x63*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x63*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_618 = Coupling(name = 'GC_618',
                  value = '(ee**2*complex(0,1)*I5x66*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x66*US1x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x66*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x66*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x66*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x66*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x66*US1x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x66*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x66*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x66*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x66*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x66*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x66*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_619 = Coupling(name = 'GC_619',
                  value = '(ee**2*complex(0,1)*I19x33*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x33*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x33*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x33*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x33*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x33*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x33*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x33*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x33*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x33*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x33*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x33*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x33*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_620 = Coupling(name = 'GC_620',
                  value = '(ee**2*complex(0,1)*I19x36*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x36*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x36*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x36*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x36*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x36*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x36*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x36*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x36*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x36*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x36*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x36*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x36*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_621 = Coupling(name = 'GC_621',
                  value = '(ee**2*complex(0,1)*I19x63*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x63*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x63*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x63*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x63*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x63*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x63*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x63*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x63*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x63*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x63*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x63*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x63*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_622 = Coupling(name = 'GC_622',
                  value = '(ee**2*complex(0,1)*I19x66*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x66*US1x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x66*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x66*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x66*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x66*sw**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x66*US1x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x66*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x66*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US1x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x66*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x66*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x66*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x66*sw**2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_623 = Coupling(name = 'GC_623',
                  value = '(ee**2*complex(0,1)*UP1x1**2*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2**2*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US1x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US1x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US1x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US1x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1**2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2**2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x2*UP1x3*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP1x3*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x3*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x3*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP1x3**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_624 = Coupling(name = 'GC_624',
                  value = '(ee**2*complex(0,1)*UP1x1*UP2x1*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2*UP2x2*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US1x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US1x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x2*US1x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US1x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US1x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US1x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1*UP2x1*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US1x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2*UP2x2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US1x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US1x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x2*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x3*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x3*US1x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x1*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x1*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x3*US1x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x1*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x1*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x2*US1x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAk*NMk*UP1x3*UP2x3*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3*UP2x3*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_625 = Coupling(name = 'GC_625',
                  value = '(ee**2*complex(0,1)*UP2x1**2*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x2**2*US1x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US1x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US1x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US1x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US1x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US1x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US1x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US1x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x1**2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP2x2**2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US1x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x2*UP2x3*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x2*UP2x3*US1x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x3*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x3*US1x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP2x3**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP2x3**2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_626 = Coupling(name = 'GC_626',
                  value = '(3*ee**2*complex(0,1)*US1x1**3*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*complex(0,1)*NMl**2*US1x1*US1x2**2*vd)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US1x1*US1x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x1*US1x2**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US1x1*US1x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x1*US1x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US1x2*US1x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US1x2*US1x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US1x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US1x1**2*US1x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x1**2*US1x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2**3*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US1x1*US1x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US1x2*US1x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US1x2*US1x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMAl*NMl*US1x1*US1x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMAl*NMl*sw**2*US1x1*US1x2*US1x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US1x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_627 = Coupling(name = 'GC_627',
                  value = '(ee**2*complex(0,1)*US2x1*vd)/(2.*sw**2) + (ee**2*complex(0,1)*US2x2*vu)/(2.*sw**2)',
                  order = {'QED':1})

GC_628 = Coupling(name = 'GC_628',
                  value = '-(ee**2*complex(0,1)*I20x44*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x44*US2x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_629 = Coupling(name = 'GC_629',
                  value = '-(ee**2*complex(0,1)*I20x55*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x55*US2x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_630 = Coupling(name = 'GC_630',
                  value = '(ee**2*complex(0,1)*I40x44*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x44*US2x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_631 = Coupling(name = 'GC_631',
                  value = '(ee**2*complex(0,1)*I40x55*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x55*US2x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_632 = Coupling(name = 'GC_632',
                  value = '-(ee**2*complex(0,1)*I6x44*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x44*US2x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_633 = Coupling(name = 'GC_633',
                  value = '-(ee**2*complex(0,1)*I6x55*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x55*US2x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_634 = Coupling(name = 'GC_634',
                  value = '(ee**2*complex(0,1)*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_635 = Coupling(name = 'GC_635',
                  value = '-(ee**2*complex(0,1)*US2x1*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_636 = Coupling(name = 'GC_636',
                  value = '(ee**2*complex(0,1)*I19x11*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x11*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_637 = Coupling(name = 'GC_637',
                  value = '(ee**2*complex(0,1)*I19x22*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x22*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_638 = Coupling(name = 'GC_638',
                  value = '-(ee**2*complex(0,1)*I39x11*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x11*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_639 = Coupling(name = 'GC_639',
                  value = '-(ee**2*complex(0,1)*I39x22*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x22*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_640 = Coupling(name = 'GC_640',
                  value = '(ee**2*complex(0,1)*I5x11*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x11*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_641 = Coupling(name = 'GC_641',
                  value = '(ee**2*complex(0,1)*I5x22*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x22*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_642 = Coupling(name = 'GC_642',
                  value = '-(ee**2*complex(0,1)*I39x33*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x33*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x33*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x33*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x33*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x33*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x33*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x33*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x33*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x33*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x33*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x33*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x33*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_643 = Coupling(name = 'GC_643',
                  value = '-(ee**2*complex(0,1)*I39x36*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x36*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x36*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x36*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x36*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x36*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x36*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x36*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x36*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x36*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x36*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x36*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x36*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_644 = Coupling(name = 'GC_644',
                  value = '-(ee**2*complex(0,1)*I39x63*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x63*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x63*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x63*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x63*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x63*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x63*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x63*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x63*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x63*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x63*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x63*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x63*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_645 = Coupling(name = 'GC_645',
                  value = '-(ee**2*complex(0,1)*I39x66*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x66*US2x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x66*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US2x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x66*US2x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x66*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x66*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x66*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x66*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x66*sw**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x66*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x66*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x66*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x66*sw**2*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_646 = Coupling(name = 'GC_646',
                  value = '(ee**2*complex(0,1)*I5x33*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x33*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x33*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x33*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x33*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x33*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x33*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x33*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x33*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x33*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x33*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x33*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x33*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_647 = Coupling(name = 'GC_647',
                  value = '(ee**2*complex(0,1)*I5x36*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x36*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x36*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x36*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x36*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x36*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x36*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x36*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x36*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x36*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x36*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x36*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x36*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_648 = Coupling(name = 'GC_648',
                  value = '(ee**2*complex(0,1)*I5x63*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x63*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x63*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x63*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x63*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x63*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x63*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x63*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x63*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x63*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x63*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x63*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x63*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_649 = Coupling(name = 'GC_649',
                  value = '(ee**2*complex(0,1)*I5x66*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x66*US2x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x66*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x66*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x66*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x66*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x66*US2x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x66*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x66*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x66*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x66*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x66*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x66*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_650 = Coupling(name = 'GC_650',
                  value = '(ee**2*complex(0,1)*I19x33*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x33*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x33*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x33*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x33*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x33*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x33*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x33*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x33*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x33*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x33*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x33*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x33*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_651 = Coupling(name = 'GC_651',
                  value = '(ee**2*complex(0,1)*I19x36*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x36*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x36*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x36*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x36*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x36*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x36*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x36*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x36*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x36*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x36*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x36*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x36*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_652 = Coupling(name = 'GC_652',
                  value = '(ee**2*complex(0,1)*I19x63*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x63*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x63*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x63*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x63*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x63*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x63*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x63*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x63*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x63*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x63*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x63*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x63*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_653 = Coupling(name = 'GC_653',
                  value = '(ee**2*complex(0,1)*I19x66*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x66*US2x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x66*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x66*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x66*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x66*sw**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x66*US2x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x66*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x66*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US2x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x66*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x66*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x66*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x66*sw**2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_654 = Coupling(name = 'GC_654',
                  value = '(ee**2*complex(0,1)*UP1x1**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x2*UP1x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP1x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP1x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_655 = Coupling(name = 'GC_655',
                  value = '(ee**2*complex(0,1)*UP1x1*UP2x1*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2*UP2x2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US2x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US2x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x2*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US2x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US2x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US2x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1*UP2x1*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US2x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2*UP2x2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US2x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US2x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x2*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x3*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x3*US2x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x1*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x1*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x3*US2x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x1*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x1*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x2*US2x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAk*NMk*UP1x3*UP2x3*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3*UP2x3*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_656 = Coupling(name = 'GC_656',
                  value = '(ee**2*complex(0,1)*UP2x1**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x2**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x1**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP2x2**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x2*UP2x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x2*UP2x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP2x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP2x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_657 = Coupling(name = 'GC_657',
                  value = '(3*ee**2*complex(0,1)*US1x1**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2**2*US2x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3**2*US2x1*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US1x2*US2x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x2*US2x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3**2*US2x2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US1x3*US2x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x2*US2x1*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US1x2*US2x1*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x2*US2x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3**2*US2x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US1x3*US2x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x2*US1x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US1x3*US2x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US1x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US1x3*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US1x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US1x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US1x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3**2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_658 = Coupling(name = 'GC_658',
                  value = '(3*ee**2*complex(0,1)*US1x1*US2x1**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US2x1*US2x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x1*US2x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US2x1*US2x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x2**2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US2x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x2**2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x3*US2x1*US2x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x3*US2x1*US2x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US2x2*US2x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US2x2*US2x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x2*US2x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x2*US2x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x1**2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x1**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x1**2*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US2x1*US2x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US2x1*US2x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US2x1*US2x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2*US2x2**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US2x1*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US2x1*US2x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x3*US2x2*US2x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x3*US2x2*US2x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x1*US2x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x1*US2x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x3*US2x1*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x3*US2x1*US2x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x2*US2x1*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US2x1*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US1x3*US2x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3*US2x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_659 = Coupling(name = 'GC_659',
                  value = '(3*ee**2*complex(0,1)*US2x1**3*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*complex(0,1)*NMl**2*US2x1*US2x2**2*vd)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US2x1*US2x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x1*US2x2**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US2x1*US2x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x1*US2x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US2x2*US2x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US2x2*US2x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US2x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US2x1**2*US2x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x1**2*US2x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US2x2**3*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US2x1*US2x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US2x1*US2x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US2x2*US2x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US2x2*US2x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMAl*NMl*US2x1*US2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMAl*NMl*sw**2*US2x1*US2x2*US2x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US2x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US2x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_660 = Coupling(name = 'GC_660',
                  value = '(ee**2*complex(0,1)*US3x1*vd)/(2.*sw**2) + (ee**2*complex(0,1)*US3x2*vu)/(2.*sw**2)',
                  order = {'QED':1})

GC_661 = Coupling(name = 'GC_661',
                  value = '-(ee**2*complex(0,1)*I20x44*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x44*US3x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_662 = Coupling(name = 'GC_662',
                  value = '-(ee**2*complex(0,1)*I20x55*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x55*US3x2*vu)/(2.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_663 = Coupling(name = 'GC_663',
                  value = '(ee**2*complex(0,1)*I40x44*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x44*US3x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_664 = Coupling(name = 'GC_664',
                  value = '(ee**2*complex(0,1)*I40x55*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x55*US3x2*vu)/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_665 = Coupling(name = 'GC_665',
                  value = '-(ee**2*complex(0,1)*I6x44*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x44*US3x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_666 = Coupling(name = 'GC_666',
                  value = '-(ee**2*complex(0,1)*I6x55*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x55*US3x2*vu)/(6.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_667 = Coupling(name = 'GC_667',
                  value = '(ee**2*complex(0,1)*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_668 = Coupling(name = 'GC_668',
                  value = '-(ee**2*complex(0,1)*US3x1*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US3x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_669 = Coupling(name = 'GC_669',
                  value = '(ee**2*complex(0,1)*I19x11*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x11*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x11*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_670 = Coupling(name = 'GC_670',
                  value = '(ee**2*complex(0,1)*I19x22*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I19x22*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x22*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_671 = Coupling(name = 'GC_671',
                  value = '-(ee**2*complex(0,1)*I39x11*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x11*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x11*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_672 = Coupling(name = 'GC_672',
                  value = '-(ee**2*complex(0,1)*I39x22*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*I39x22*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x22*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_673 = Coupling(name = 'GC_673',
                  value = '(ee**2*complex(0,1)*I5x11*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x11*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x11*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_674 = Coupling(name = 'GC_674',
                  value = '(ee**2*complex(0,1)*I5x22*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*I5x22*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x22*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                  order = {'QED':1})

GC_675 = Coupling(name = 'GC_675',
                  value = '-(ee**2*complex(0,1)*I39x33*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x33*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x33*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x33*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x33*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x33*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x33*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x33*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x33*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x33*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x33*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x33*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x33*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x33*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x33*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x33*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x33*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_676 = Coupling(name = 'GC_676',
                  value = '-(ee**2*complex(0,1)*I39x36*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x36*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x36*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x36*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x36*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x36*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x36*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x36*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x36*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x36*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x36*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x36*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x36*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x36*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x36*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x36*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x36*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_677 = Coupling(name = 'GC_677',
                  value = '-(ee**2*complex(0,1)*I39x63*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x63*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x63*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x63*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x63*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x63*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x63*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x63*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x63*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x63*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x63*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x63*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x63*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x63*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x63*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x63*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x63*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_678 = Coupling(name = 'GC_678',
                  value = '-(ee**2*complex(0,1)*I39x66*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I40x66*US3x1*vd)/(3.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I72x66*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I75x66*NMl*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I72x66*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I75x66*NMl*sw**2*US3x3*vd)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I39x66*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I40x66*US3x2*vu)/(3.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I80x66*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I81x66*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I39x66*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I80x66*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I81x66*sw**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I73x66*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I74x66*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I73x66*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I74x66*sw**2*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_679 = Coupling(name = 'GC_679',
                  value = '(ee**2*complex(0,1)*I5x33*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x33*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x33*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x33*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x33*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x33*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x33*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x33*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x33*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x33*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x33*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x33*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x33*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x33*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x33*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x33*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x33*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_680 = Coupling(name = 'GC_680',
                  value = '(ee**2*complex(0,1)*I5x36*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x36*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x36*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x36*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x36*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x36*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x36*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x36*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x36*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x36*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x36*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x36*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x36*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x36*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x36*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x36*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x36*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_681 = Coupling(name = 'GC_681',
                  value = '(ee**2*complex(0,1)*I5x63*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x63*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x63*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x63*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x63*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x63*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x63*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x63*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x63*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x63*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x63*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x63*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x63*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x63*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x63*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x63*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x63*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_682 = Coupling(name = 'GC_682',
                  value = '(ee**2*complex(0,1)*I5x66*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I6x66*US3x1*vd)/(6.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I76x66*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I77x66*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I76x66*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I77x66*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I5x66*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I6x66*US3x2*vu)/(6.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I5x66*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I66x66*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I67x66*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I66x66*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I67x66*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I64x66*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I65x66*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I64x66*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I65x66*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_683 = Coupling(name = 'GC_683',
                  value = '(ee**2*complex(0,1)*I19x33*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x33*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x33*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x33*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x33*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x33*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x33*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x33*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x33*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x33*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x33*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x33*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x33*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x33*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x33*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x33*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x33*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_684 = Coupling(name = 'GC_684',
                  value = '(ee**2*complex(0,1)*I19x36*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x36*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x36*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x36*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x36*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x36*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x36*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x36*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x36*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x36*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x36*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x36*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x36*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x36*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x36*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x36*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x36*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_685 = Coupling(name = 'GC_685',
                  value = '(ee**2*complex(0,1)*I19x63*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x63*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x63*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x63*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x63*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x63*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x63*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x63*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x63*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x63*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x63*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x63*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x63*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x63*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x63*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x63*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x63*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_686 = Coupling(name = 'GC_686',
                  value = '(ee**2*complex(0,1)*I19x66*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I20x66*US3x1*vd)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I78x66*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*I79x66*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I78x66*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I79x66*sw**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*I19x66*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I20x66*US3x2*vu)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*I19x66*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*I70x66*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*I71x66*NMl*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I70x66*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I71x66*NMl*sw**2*US3x3*vu)/(2.*(-1 + sw)*(1 + sw)) + (complex(0,1)*I68x66*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*I69x66*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I68x66*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*I69x66*sw**2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_687 = Coupling(name = 'GC_687',
                  value = '(ee**2*complex(0,1)*UP1x1**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x2*UP1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x2*UP1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP1x1*UP1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP1x1*UP1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x2*UP1x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP1x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x1*UP1x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP1x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP1x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_688 = Coupling(name = 'GC_688',
                  value = '(ee**2*complex(0,1)*UP1x1*UP2x1*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x2*UP2x2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x2*UP2x2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x2*UP2x2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x2*US3x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x2*UP2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x2*UP2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP1x3*UP2x3*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x3*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP1x1*UP2x1*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP1x1*UP2x1*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x1*UP2x1*US3x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP1x2*UP2x2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP1x3*UP2x3*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP1x3*UP2x3*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x3*UP2x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x3*UP2x1*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*UP1x1*UP2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*UP1x1*UP2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x3*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x3*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x3*UP2x1*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x3*UP2x1*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x2*UP2x1*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x2*UP2x1*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*UP1x1*UP2x2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*sw**2*UP1x1*UP2x2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAk*NMk*UP1x3*UP2x3*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP1x3*UP2x3*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_689 = Coupling(name = 'GC_689',
                  value = '(ee**2*complex(0,1)*UP2x1**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x2**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x2*UP2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x2*UP2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*UP2x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*sw**2*UP2x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*UP2x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*UP2x1**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*UP2x2**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*UP2x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*UP2x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*UP2x1*UP2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*UP2x1*UP2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x2*UP2x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x2*UP2x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*UP2x1*UP2x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*UP2x1*UP2x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*UP2x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*sw**2*UP2x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_690 = Coupling(name = 'GC_690',
                  value = '(3*ee**2*complex(0,1)*US1x1**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US1x2*US3x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x2*US1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x2*US1x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US1x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US1x2*US3x1*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US1x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x1*US1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x1*US1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US1x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x2*US1x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US1x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US1x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US1x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US1x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US1x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US1x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_691 = Coupling(name = 'GC_691',
                  value = '(3*ee**2*complex(0,1)*US1x1*US2x1*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x3*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x3*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x1*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x1*US3x2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x1*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US2x2*US3x2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x3*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x3*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x2*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x2*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x1*US3x1*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US2x1*US3x1*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x1*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US2x2*US3x1*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x3*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x3*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US2x1*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US2x1*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US2x1*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2*US2x2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x3*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x3*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x3*US2x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x3*US2x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*US2x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*US2x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x1*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x1*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x3*US2x2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x3*US2x2*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*US1x2*US2x3*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US2x3*US3x1)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*US1x3*US2x1*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x3*US2x1*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*US1x1*US2x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US2x3*US3x2)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*US1x2*US2x1*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US2x1*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMAl*NMl*US1x1*US2x2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US2x2*US3x3)/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMAk*NMk*US1x3*US2x3*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3*US2x3*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_692 = Coupling(name = 'GC_692',
                  value = '(3*ee**2*complex(0,1)*US2x1**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US2x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x2**2*US3x1*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x2**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3**2*US3x1*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x1*US2x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1*US2x2*US3x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US2x2*US3x2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US2x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US2x3**2*US3x2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x1*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x2*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x2*US2x3*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x1*US2x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1*US2x2*US3x1*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US2x2*US3x1*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US2x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US2x3**2*US3x1*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US2x2**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*US2x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x1*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x1*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x2*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x2*US2x3*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x2*US2x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x2*US2x3*US3x1*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x1*US2x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x1*US2x3*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x1*US2x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x1*US2x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US2x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US2x3**2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_693 = Coupling(name = 'GC_693',
                  value = '(3*ee**2*complex(0,1)*US1x1*US3x1**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x2*US3x1*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US3x1*US3x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x2*US3x1*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US3x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x3*US3x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x3*US3x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US3x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US3x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US3x1**2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x2*US3x1**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US3x1**2*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x1*US3x1*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*US3x1*US3x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x1*US3x1*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US1x2*US3x2**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US1x3*US3x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US1x3*US3x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US1x3*US3x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US1x3*US3x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US1x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US1x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x3*US3x1*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x3*US3x1*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x2*US3x1*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x2*US3x1*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US1x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US1x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US1x3*US3x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US1x3*US3x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_694 = Coupling(name = 'GC_694',
                  value = '(3*ee**2*complex(0,1)*US2x1*US3x1**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x2*US3x1*US3x2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x2*US3x1*US3x2*vd)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x2*US3x1*US3x2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1*US3x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x3*US3x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x3*US3x1*US3x3*vd)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x3*US3x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x3*US3x2*US3x3*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US2x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US2x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x2*US3x1**2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x2*US3x1**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x2*US3x1**2*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x1*US3x1*US3x2*vu)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1*US3x1*US3x2*vu)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x1*US3x1*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US2x2*US3x2**2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (2*complex(0,1)*NMk*NMl*US2x3*US3x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*sw**2*US2x3*US3x1*US3x3*vu)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMl**2*US2x3*US3x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMl**2*sw**2*US2x3*US3x2*US3x3*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*NMl*US2x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*NMl*sw**2*US2x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x3*US3x1*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x3*US3x1*US3x2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x2*US3x1*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x2*US3x1*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*US2x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*sw**2*US2x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US2x3*US3x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US2x3*US3x3**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_695 = Coupling(name = 'GC_695',
                  value = '(3*ee**2*complex(0,1)*US3x1**3*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (3*complex(0,1)*NMl**2*US3x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US3x1*US3x2**2*vd)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x1*US3x2**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US3x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x1*US3x3**2*vd)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US3x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US3x2*US3x3**2*vd)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US3x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) - (3*ee**2*complex(0,1)*US3x1**2*US3x2*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x1**2*US3x2*vu)/((-1 + sw)*(1 + sw)) + (3*ee**2*complex(0,1)*US3x2**3*vu)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (3*complex(0,1)*NMk*NMl*US3x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMk*NMl*sw**2*US3x1*US3x3**2*vu)/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMl**2*US3x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMl**2*sw**2*US3x2*US3x3**2*vu)/((-1 + sw)*(1 + sw)) - (3*complex(0,1)*NMAl*NMl*US3x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (3*complex(0,1)*NMAl*NMl*sw**2*US3x1*US3x2*US3x3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAk*NMk*US3x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAk*NMk*sw**2*US3x3**3*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_696 = Coupling(name = 'GC_696',
                  value = '-((ee*complex(0,1)*I87x11*VV1x1)/sw)',
                  order = {'QED':1})

GC_697 = Coupling(name = 'GC_697',
                  value = '-((ee*complex(0,1)*I87x22*VV1x1)/sw)',
                  order = {'QED':1})

GC_698 = Coupling(name = 'GC_698',
                  value = '-((ee*complex(0,1)*I87x33*VV1x1)/sw)',
                  order = {'QED':1})

GC_699 = Coupling(name = 'GC_699',
                  value = '-((ee*complex(0,1)*I89x11*VV1x1)/sw)',
                  order = {'QED':1})

GC_700 = Coupling(name = 'GC_700',
                  value = '-((ee*complex(0,1)*I89x22*VV1x1)/sw)',
                  order = {'QED':1})

GC_701 = Coupling(name = 'GC_701',
                  value = 'complex(0,1)*I8x33*VV1x2',
                  order = {'QED':1})

GC_702 = Coupling(name = 'GC_702',
                  value = 'complex(0,1)*I8x36*VV1x2',
                  order = {'QED':1})

GC_703 = Coupling(name = 'GC_703',
                  value = '-((ee*complex(0,1)*I89x33*VV1x1)/sw) + complex(0,1)*I90x33*VV1x2',
                  order = {'QED':1})

GC_704 = Coupling(name = 'GC_704',
                  value = '-((ee*complex(0,1)*I89x36*VV1x1)/sw) + complex(0,1)*I90x36*VV1x2',
                  order = {'QED':1})

GC_705 = Coupling(name = 'GC_705',
                  value = '(ee*UP1x1*UU1x2*VV1x1)/(sw*cmath.sqrt(2)) + (ee*UP1x2*UU1x1*VV1x2)/(sw*cmath.sqrt(2)) - (NMl*UP1x3*UU1x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_706 = Coupling(name = 'GC_706',
                  value = '(ee*UP2x1*UU1x2*VV1x1)/(sw*cmath.sqrt(2)) + (ee*UP2x2*UU1x1*VV1x2)/(sw*cmath.sqrt(2)) - (NMl*UP2x3*UU1x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_707 = Coupling(name = 'GC_707',
                  value = '-((ee*complex(0,1)*US1x1*UU1x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*UU1x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*UU1x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_708 = Coupling(name = 'GC_708',
                  value = '-((ee*complex(0,1)*US2x1*UU1x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*UU1x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*UU1x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_709 = Coupling(name = 'GC_709',
                  value = '-((ee*complex(0,1)*US3x1*UU1x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*UU1x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*UU1x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_710 = Coupling(name = 'GC_710',
                  value = '(ee*UP1x1*UU2x2*VV1x1)/(sw*cmath.sqrt(2)) + (ee*UP1x2*UU2x1*VV1x2)/(sw*cmath.sqrt(2)) - (NMl*UP1x3*UU2x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_711 = Coupling(name = 'GC_711',
                  value = '(ee*UP2x1*UU2x2*VV1x1)/(sw*cmath.sqrt(2)) + (ee*UP2x2*UU2x1*VV1x2)/(sw*cmath.sqrt(2)) - (NMl*UP2x3*UU2x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_712 = Coupling(name = 'GC_712',
                  value = '-((ee*complex(0,1)*US1x1*UU2x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*UU2x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*UU2x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_713 = Coupling(name = 'GC_713',
                  value = '-((ee*complex(0,1)*US2x1*UU2x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*UU2x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*UU2x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_714 = Coupling(name = 'GC_714',
                  value = '-((ee*complex(0,1)*US3x1*UU2x2*VV1x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*UU2x1*VV1x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*UU2x2*VV1x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_715 = Coupling(name = 'GC_715',
                  value = '-((ee*complex(0,1)*I87x11*VV2x1)/sw)',
                  order = {'QED':1})

GC_716 = Coupling(name = 'GC_716',
                  value = '-((ee*complex(0,1)*I87x22*VV2x1)/sw)',
                  order = {'QED':1})

GC_717 = Coupling(name = 'GC_717',
                  value = '-((ee*complex(0,1)*I87x33*VV2x1)/sw)',
                  order = {'QED':1})

GC_718 = Coupling(name = 'GC_718',
                  value = '-((ee*complex(0,1)*I89x11*VV2x1)/sw)',
                  order = {'QED':1})

GC_719 = Coupling(name = 'GC_719',
                  value = '-((ee*complex(0,1)*I89x22*VV2x1)/sw)',
                  order = {'QED':1})

GC_720 = Coupling(name = 'GC_720',
                  value = 'complex(0,1)*I8x33*VV2x2',
                  order = {'QED':1})

GC_721 = Coupling(name = 'GC_721',
                  value = 'complex(0,1)*I8x36*VV2x2',
                  order = {'QED':1})

GC_722 = Coupling(name = 'GC_722',
                  value = '-((ee*complex(0,1)*I89x33*VV2x1)/sw) + complex(0,1)*I90x33*VV2x2',
                  order = {'QED':1})

GC_723 = Coupling(name = 'GC_723',
                  value = '-((ee*complex(0,1)*I89x36*VV2x1)/sw) + complex(0,1)*I90x36*VV2x2',
                  order = {'QED':1})

GC_724 = Coupling(name = 'GC_724',
                  value = '(ee*UP1x1*UU1x2*VV2x1)/(sw*cmath.sqrt(2)) + (ee*UP1x2*UU1x1*VV2x2)/(sw*cmath.sqrt(2)) - (NMl*UP1x3*UU1x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_725 = Coupling(name = 'GC_725',
                  value = '(ee*UP2x1*UU1x2*VV2x1)/(sw*cmath.sqrt(2)) + (ee*UP2x2*UU1x1*VV2x2)/(sw*cmath.sqrt(2)) - (NMl*UP2x3*UU1x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_726 = Coupling(name = 'GC_726',
                  value = '-((ee*complex(0,1)*US1x1*UU1x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*UU1x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*UU1x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_727 = Coupling(name = 'GC_727',
                  value = '-((ee*complex(0,1)*US2x1*UU1x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*UU1x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*UU1x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_728 = Coupling(name = 'GC_728',
                  value = '-((ee*complex(0,1)*US3x1*UU1x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*UU1x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*UU1x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_729 = Coupling(name = 'GC_729',
                  value = '(ee*UP1x1*UU2x2*VV2x1)/(sw*cmath.sqrt(2)) + (ee*UP1x2*UU2x1*VV2x2)/(sw*cmath.sqrt(2)) - (NMl*UP1x3*UU2x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_730 = Coupling(name = 'GC_730',
                  value = '(ee*UP2x1*UU2x2*VV2x1)/(sw*cmath.sqrt(2)) + (ee*UP2x2*UU2x1*VV2x2)/(sw*cmath.sqrt(2)) - (NMl*UP2x3*UU2x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_731 = Coupling(name = 'GC_731',
                  value = '-((ee*complex(0,1)*US1x1*UU2x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*UU2x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*UU2x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_732 = Coupling(name = 'GC_732',
                  value = '-((ee*complex(0,1)*US2x1*UU2x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*UU2x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*UU2x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_733 = Coupling(name = 'GC_733',
                  value = '-((ee*complex(0,1)*US3x1*UU2x2*VV2x1)/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*UU2x1*VV2x2)/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*UU2x2*VV2x2)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_734 = Coupling(name = 'GC_734',
                  value = '(UP1x1*yd3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_735 = Coupling(name = 'GC_735',
                  value = '(UP2x1*yd3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_736 = Coupling(name = 'GC_736',
                  value = '-((complex(0,1)*US1x1*yd3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_737 = Coupling(name = 'GC_737',
                  value = '-((complex(0,1)*US2x1*yd3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_738 = Coupling(name = 'GC_738',
                  value = '-((complex(0,1)*US3x1*yd3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_739 = Coupling(name = 'GC_739',
                  value = '(UP1x1*ye3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_740 = Coupling(name = 'GC_740',
                  value = '(UP2x1*ye3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_741 = Coupling(name = 'GC_741',
                  value = '-((complex(0,1)*US1x1*ye3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_742 = Coupling(name = 'GC_742',
                  value = '-((complex(0,1)*US2x1*ye3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_743 = Coupling(name = 'GC_743',
                  value = '-((complex(0,1)*US3x1*ye3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_744 = Coupling(name = 'GC_744',
                  value = '(UP1x2*yu3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_745 = Coupling(name = 'GC_745',
                  value = '(UP2x2*yu3x3)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_746 = Coupling(name = 'GC_746',
                  value = '-((complex(0,1)*US1x2*yu3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_747 = Coupling(name = 'GC_747',
                  value = '-((complex(0,1)*US2x2*yu3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_748 = Coupling(name = 'GC_748',
                  value = '-((complex(0,1)*US3x2*yu3x3)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_749 = Coupling(name = 'GC_749',
                  value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_750 = Coupling(name = 'GC_750',
                  value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_751 = Coupling(name = 'GC_751',
                  value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_752 = Coupling(name = 'GC_752',
                  value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_753 = Coupling(name = 'GC_753',
                  value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_754 = Coupling(name = 'GC_754',
                  value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_755 = Coupling(name = 'GC_755',
                  value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_756 = Coupling(name = 'GC_756',
                  value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_757 = Coupling(name = 'GC_757',
                  value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_758 = Coupling(name = 'GC_758',
                  value = '(complex(0,1)*I11x33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_759 = Coupling(name = 'GC_759',
                  value = '(complex(0,1)*I11x36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_760 = Coupling(name = 'GC_760',
                  value = '(complex(0,1)*I22x33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_761 = Coupling(name = 'GC_761',
                  value = '(complex(0,1)*I22x36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_762 = Coupling(name = 'GC_762',
                  value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN1x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN1x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_763 = Coupling(name = 'GC_763',
                  value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN1x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN1x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_764 = Coupling(name = 'GC_764',
                  value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_765 = Coupling(name = 'GC_765',
                  value = '(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_766 = Coupling(name = 'GC_766',
                  value = '(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_767 = Coupling(name = 'GC_767',
                  value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_768 = Coupling(name = 'GC_768',
                  value = '(cw*ee*complex(0,1)*NN5x3*complexconjugate(NN1x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN5x4*complexconjugate(NN1x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_769 = Coupling(name = 'GC_769',
                  value = '(complex(0,1)*I50x33*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x33*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_770 = Coupling(name = 'GC_770',
                  value = '(complex(0,1)*I50x36*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x36*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN1x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN1x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_771 = Coupling(name = 'GC_771',
                  value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN1x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN1x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_772 = Coupling(name = 'GC_772',
                  value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN1x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN1x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_773 = Coupling(name = 'GC_773',
                  value = '-((cw*ee*UP1x1*complexconjugate(NN1x1)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw))) + (ee*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x1)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP1x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_774 = Coupling(name = 'GC_774',
                  value = '-((cw*ee*UP2x1*complexconjugate(NN1x1)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw))) + (ee*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x1)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP2x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_775 = Coupling(name = 'GC_775',
                  value = '-((cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x1)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x1)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_776 = Coupling(name = 'GC_776',
                  value = '-((cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x1)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x1)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_777 = Coupling(name = 'GC_777',
                  value = '-((cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x1)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x2)*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x1)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x2)*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x3)*complexconjugate(NN1x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x3)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x4)*complexconjugate(NN1x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN1x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_778 = Coupling(name = 'GC_778',
                  value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_779 = Coupling(name = 'GC_779',
                  value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_780 = Coupling(name = 'GC_780',
                  value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_781 = Coupling(name = 'GC_781',
                  value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_782 = Coupling(name = 'GC_782',
                  value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_783 = Coupling(name = 'GC_783',
                  value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_784 = Coupling(name = 'GC_784',
                  value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_785 = Coupling(name = 'GC_785',
                  value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_786 = Coupling(name = 'GC_786',
                  value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_787 = Coupling(name = 'GC_787',
                  value = '(complex(0,1)*I11x33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_788 = Coupling(name = 'GC_788',
                  value = '(complex(0,1)*I11x36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_789 = Coupling(name = 'GC_789',
                  value = '(complex(0,1)*I22x33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_790 = Coupling(name = 'GC_790',
                  value = '(complex(0,1)*I22x36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN2x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_791 = Coupling(name = 'GC_791',
                  value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN2x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN2x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_792 = Coupling(name = 'GC_792',
                  value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN2x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN2x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_793 = Coupling(name = 'GC_793',
                  value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_794 = Coupling(name = 'GC_794',
                  value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_795 = Coupling(name = 'GC_795',
                  value = '(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_796 = Coupling(name = 'GC_796',
                  value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_797 = Coupling(name = 'GC_797',
                  value = '(cw*ee*complex(0,1)*NN5x3*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN5x4*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_798 = Coupling(name = 'GC_798',
                  value = '(complex(0,1)*I50x33*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x33*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_799 = Coupling(name = 'GC_799',
                  value = '(complex(0,1)*I50x36*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x36*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN2x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_800 = Coupling(name = 'GC_800',
                  value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN2x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN2x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_801 = Coupling(name = 'GC_801',
                  value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN2x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN2x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_802 = Coupling(name = 'GC_802',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN1x1)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x1)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_803 = Coupling(name = 'GC_803',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN1x1)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x1)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_804 = Coupling(name = 'GC_804',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x1)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x1)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_805 = Coupling(name = 'GC_805',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x1)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x1)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_806 = Coupling(name = 'GC_806',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN2x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x3)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x4)*complexconjugate(NN2x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x1)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x2)*complexconjugate(NN2x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x1)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x2)*complexconjugate(NN2x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x4)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x5)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x3)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x5)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x3)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x4)*complexconjugate(NN2x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN1x5)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_807 = Coupling(name = 'GC_807',
                  value = '-((cw*ee*UP1x1*complexconjugate(NN2x1)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw))) + (ee*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x1)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP1x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_808 = Coupling(name = 'GC_808',
                  value = '-((cw*ee*UP2x1*complexconjugate(NN2x1)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw))) + (ee*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x1)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP2x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_809 = Coupling(name = 'GC_809',
                  value = '-((cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x1)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x1)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_810 = Coupling(name = 'GC_810',
                  value = '-((cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x1)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x1)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_811 = Coupling(name = 'GC_811',
                  value = '-((cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x1)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x2)*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x1)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x2)*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x3)*complexconjugate(NN2x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x3)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x4)*complexconjugate(NN2x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN2x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_812 = Coupling(name = 'GC_812',
                  value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_813 = Coupling(name = 'GC_813',
                  value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_814 = Coupling(name = 'GC_814',
                  value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_815 = Coupling(name = 'GC_815',
                  value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_816 = Coupling(name = 'GC_816',
                  value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_817 = Coupling(name = 'GC_817',
                  value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_818 = Coupling(name = 'GC_818',
                  value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_819 = Coupling(name = 'GC_819',
                  value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_820 = Coupling(name = 'GC_820',
                  value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_821 = Coupling(name = 'GC_821',
                  value = '(complex(0,1)*I11x33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_822 = Coupling(name = 'GC_822',
                  value = '(complex(0,1)*I11x36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_823 = Coupling(name = 'GC_823',
                  value = '(complex(0,1)*I22x33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_824 = Coupling(name = 'GC_824',
                  value = '(complex(0,1)*I22x36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN3x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_825 = Coupling(name = 'GC_825',
                  value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN3x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN3x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_826 = Coupling(name = 'GC_826',
                  value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN3x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN3x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_827 = Coupling(name = 'GC_827',
                  value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_828 = Coupling(name = 'GC_828',
                  value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_829 = Coupling(name = 'GC_829',
                  value = '-(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_830 = Coupling(name = 'GC_830',
                  value = '(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_831 = Coupling(name = 'GC_831',
                  value = '(cw*ee*complex(0,1)*NN5x3*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN5x4*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_832 = Coupling(name = 'GC_832',
                  value = '(complex(0,1)*I50x33*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x33*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_833 = Coupling(name = 'GC_833',
                  value = '(complex(0,1)*I50x36*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x36*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN3x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN3x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN3x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_834 = Coupling(name = 'GC_834',
                  value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN3x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN3x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_835 = Coupling(name = 'GC_835',
                  value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN3x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN3x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_836 = Coupling(name = 'GC_836',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN1x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_837 = Coupling(name = 'GC_837',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN1x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_838 = Coupling(name = 'GC_838',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_839 = Coupling(name = 'GC_839',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_840 = Coupling(name = 'GC_840',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN1x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_841 = Coupling(name = 'GC_841',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN2x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_842 = Coupling(name = 'GC_842',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN2x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_843 = Coupling(name = 'GC_843',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_844 = Coupling(name = 'GC_844',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_845 = Coupling(name = 'GC_845',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN3x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x3)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x4)*complexconjugate(NN3x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x1)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x2)*complexconjugate(NN3x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x1)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x2)*complexconjugate(NN3x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x4)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x5)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x3)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x5)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x3)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x4)*complexconjugate(NN3x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN2x5)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_846 = Coupling(name = 'GC_846',
                  value = '-((cw*ee*UP1x1*complexconjugate(NN3x1)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw))) + (ee*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN3x1)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP1x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_847 = Coupling(name = 'GC_847',
                  value = '-((cw*ee*UP2x1*complexconjugate(NN3x1)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw))) + (ee*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN3x1)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP2x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_848 = Coupling(name = 'GC_848',
                  value = '-((cw*ee*complex(0,1)*US1x1*complexconjugate(NN3x1)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US1x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN3x1)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_849 = Coupling(name = 'GC_849',
                  value = '-((cw*ee*complex(0,1)*US2x1*complexconjugate(NN3x1)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US2x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN3x1)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_850 = Coupling(name = 'GC_850',
                  value = '-((cw*ee*complex(0,1)*US3x1*complexconjugate(NN3x1)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US3x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN3x2)*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN3x1)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN3x2)*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN3x3)*complexconjugate(NN3x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN3x3)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN3x4)*complexconjugate(NN3x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN3x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_851 = Coupling(name = 'GC_851',
                  value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_852 = Coupling(name = 'GC_852',
                  value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_853 = Coupling(name = 'GC_853',
                  value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_854 = Coupling(name = 'GC_854',
                  value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_855 = Coupling(name = 'GC_855',
                  value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_856 = Coupling(name = 'GC_856',
                  value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_857 = Coupling(name = 'GC_857',
                  value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_858 = Coupling(name = 'GC_858',
                  value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_859 = Coupling(name = 'GC_859',
                  value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_860 = Coupling(name = 'GC_860',
                  value = '(complex(0,1)*I11x33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_861 = Coupling(name = 'GC_861',
                  value = '(complex(0,1)*I11x36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_862 = Coupling(name = 'GC_862',
                  value = '(complex(0,1)*I22x33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_863 = Coupling(name = 'GC_863',
                  value = '(complex(0,1)*I22x36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN4x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_864 = Coupling(name = 'GC_864',
                  value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN4x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN4x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_865 = Coupling(name = 'GC_865',
                  value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN4x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN4x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_866 = Coupling(name = 'GC_866',
                  value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_867 = Coupling(name = 'GC_867',
                  value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_868 = Coupling(name = 'GC_868',
                  value = '-(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_869 = Coupling(name = 'GC_869',
                  value = '-(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_870 = Coupling(name = 'GC_870',
                  value = '(cw*ee*complex(0,1)*NN5x3*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*NN5x4*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_871 = Coupling(name = 'GC_871',
                  value = '(complex(0,1)*I50x33*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x33*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_872 = Coupling(name = 'GC_872',
                  value = '(complex(0,1)*I50x36*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x36*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN4x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN4x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN4x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_873 = Coupling(name = 'GC_873',
                  value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN4x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN4x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_874 = Coupling(name = 'GC_874',
                  value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN4x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN4x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_875 = Coupling(name = 'GC_875',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN1x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_876 = Coupling(name = 'GC_876',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN1x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_877 = Coupling(name = 'GC_877',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_878 = Coupling(name = 'GC_878',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_879 = Coupling(name = 'GC_879',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN1x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_880 = Coupling(name = 'GC_880',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN2x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_881 = Coupling(name = 'GC_881',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN2x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_882 = Coupling(name = 'GC_882',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_883 = Coupling(name = 'GC_883',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_884 = Coupling(name = 'GC_884',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN2x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_885 = Coupling(name = 'GC_885',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN3x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN3x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_886 = Coupling(name = 'GC_886',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN3x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN3x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_887 = Coupling(name = 'GC_887',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN3x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN3x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN3x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN3x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_888 = Coupling(name = 'GC_888',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN3x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN3x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN3x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN3x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_889 = Coupling(name = 'GC_889',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN3x3)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN3x4)*complexconjugate(NN4x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN3x3)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN3x4)*complexconjugate(NN4x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN3x1)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN3x2)*complexconjugate(NN4x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN3x1)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN3x2)*complexconjugate(NN4x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN3x4)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN3x5)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN3x3)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN3x5)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN3x3)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN3x4)*complexconjugate(NN4x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN3x5)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_890 = Coupling(name = 'GC_890',
                  value = '-((cw*ee*UP1x1*complexconjugate(NN4x1)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw))) + (ee*UP1x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN4x1)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP1x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_891 = Coupling(name = 'GC_891',
                  value = '-((cw*ee*UP2x1*complexconjugate(NN4x1)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw))) + (ee*UP2x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN4x1)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP2x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_892 = Coupling(name = 'GC_892',
                  value = '-((cw*ee*complex(0,1)*US1x1*complexconjugate(NN4x1)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US1x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN4x1)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_893 = Coupling(name = 'GC_893',
                  value = '-((cw*ee*complex(0,1)*US2x1*complexconjugate(NN4x1)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US2x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN4x1)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_894 = Coupling(name = 'GC_894',
                  value = '-((cw*ee*complex(0,1)*US3x1*complexconjugate(NN4x1)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US3x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN4x2)*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN4x1)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN4x2)*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN4x3)*complexconjugate(NN4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN4x3)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN4x4)*complexconjugate(NN4x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN4x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_895 = Coupling(name = 'GC_895',
                  value = '(cw*ee*complex(0,1)*Rd1x1*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd1x1*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd1x1*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_896 = Coupling(name = 'GC_896',
                  value = '(cw*ee*complex(0,1)*Rd2x2*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd2x2*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd2x2*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_897 = Coupling(name = 'GC_897',
                  value = '-((cw*ee*complex(0,1)*Rl1x1*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl1x1*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl1x1*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_898 = Coupling(name = 'GC_898',
                  value = '-((cw*ee*complex(0,1)*Rl2x2*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*Rl2x2*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl2x2*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_899 = Coupling(name = 'GC_899',
                  value = '-((cw*ee*complex(0,1)*Rn1x1*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn1x1*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn1x1*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_900 = Coupling(name = 'GC_900',
                  value = '-((cw*ee*complex(0,1)*Rn2x2*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn2x2*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn2x2*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_901 = Coupling(name = 'GC_901',
                  value = '-((cw*ee*complex(0,1)*Rn3x3*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*Rn3x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rn3x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_902 = Coupling(name = 'GC_902',
                  value = '(cw*ee*complex(0,1)*Ru1x1*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru1x1*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru1x1*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_903 = Coupling(name = 'GC_903',
                  value = '(cw*ee*complex(0,1)*Ru2x2*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru2x2*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru2x2*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_904 = Coupling(name = 'GC_904',
                  value = '(complex(0,1)*I11x33*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x33*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd3x3*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd3x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd3x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_905 = Coupling(name = 'GC_905',
                  value = '(complex(0,1)*I11x36*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I11x36*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Rd6x3*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rd6x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rd6x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_906 = Coupling(name = 'GC_906',
                  value = '(complex(0,1)*I22x33*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x33*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl3x3*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl3x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl3x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_907 = Coupling(name = 'GC_907',
                  value = '(complex(0,1)*I22x36*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I22x36*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*Rl6x3*complexconjugate(NN5x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Rl6x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Rl6x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_908 = Coupling(name = 'GC_908',
                  value = '-((ee*complex(0,1)*UU1x1*complexconjugate(NN5x2))/sw) - (ee*complex(0,1)*UU1x2*complexconjugate(NN5x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_909 = Coupling(name = 'GC_909',
                  value = '-((ee*complex(0,1)*UU2x1*complexconjugate(NN5x2))/sw) - (ee*complex(0,1)*UU2x2*complexconjugate(NN5x3))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_910 = Coupling(name = 'GC_910',
                  value = '-(cw*ee*complex(0,1)*NN1x3*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN1x4*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_911 = Coupling(name = 'GC_911',
                  value = '-(cw*ee*complex(0,1)*NN2x3*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN2x4*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_912 = Coupling(name = 'GC_912',
                  value = '-(cw*ee*complex(0,1)*NN3x3*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN3x4*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_913 = Coupling(name = 'GC_913',
                  value = '-(cw*ee*complex(0,1)*NN4x3*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN4x4*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_914 = Coupling(name = 'GC_914',
                  value = '-(cw*ee*complex(0,1)*NN5x3*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*NN5x4*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw))',
                  order = {'QED':1})

GC_915 = Coupling(name = 'GC_915',
                  value = '(complex(0,1)*I50x33*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x33*sw**2*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru3x3*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru3x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru3x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_916 = Coupling(name = 'GC_916',
                  value = '(complex(0,1)*I50x36*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I50x36*sw**2*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*Ru6x3*complexconjugate(NN5x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*Ru6x3*complexconjugate(NN5x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*Ru6x3*sw*complexconjugate(NN5x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_917 = Coupling(name = 'GC_917',
                  value = '-((ee*complex(0,1)*VV1x1*complexconjugate(NN5x2))/sw) + (ee*complex(0,1)*VV1x2*complexconjugate(NN5x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_918 = Coupling(name = 'GC_918',
                  value = '-((ee*complex(0,1)*VV2x1*complexconjugate(NN5x2))/sw) + (ee*complex(0,1)*VV2x2*complexconjugate(NN5x4))/(sw*cmath.sqrt(2))',
                  order = {'QED':1})

GC_919 = Coupling(name = 'GC_919',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN1x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN1x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_920 = Coupling(name = 'GC_920',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN1x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN1x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_921 = Coupling(name = 'GC_921',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN1x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN1x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_922 = Coupling(name = 'GC_922',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN1x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN1x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_923 = Coupling(name = 'GC_923',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN1x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN1x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN1x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN1x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN1x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN1x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN1x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN1x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_924 = Coupling(name = 'GC_924',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN2x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN2x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_925 = Coupling(name = 'GC_925',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN2x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN2x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_926 = Coupling(name = 'GC_926',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN2x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN2x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_927 = Coupling(name = 'GC_927',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN2x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN2x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_928 = Coupling(name = 'GC_928',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN2x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN2x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN2x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN2x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN2x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN2x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN2x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN2x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_929 = Coupling(name = 'GC_929',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN3x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN3x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_930 = Coupling(name = 'GC_930',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN3x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN3x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_931 = Coupling(name = 'GC_931',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN3x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN3x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN3x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN3x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_932 = Coupling(name = 'GC_932',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN3x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN3x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN3x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN3x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_933 = Coupling(name = 'GC_933',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN3x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN3x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN3x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN3x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN3x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN3x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN3x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN3x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN3x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN3x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN3x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN3x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN3x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN3x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN3x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_934 = Coupling(name = 'GC_934',
                  value = '-(cw*ee*UP1x1*complexconjugate(NN4x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN4x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP1x1*complexconjugate(NN4x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP1x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN4x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP1x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP1x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP1x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_935 = Coupling(name = 'GC_935',
                  value = '-(cw*ee*UP2x1*complexconjugate(NN4x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN4x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*UP2x1*complexconjugate(NN4x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*UP2x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN4x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (NMl*UP2x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMl*sw**2*UP2x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (NMk*UP2x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_936 = Coupling(name = 'GC_936',
                  value = '-(cw*ee*complex(0,1)*US1x1*complexconjugate(NN4x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN4x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US1x1*complexconjugate(NN4x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US1x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN4x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_937 = Coupling(name = 'GC_937',
                  value = '-(cw*ee*complex(0,1)*US2x1*complexconjugate(NN4x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN4x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US2x1*complexconjugate(NN4x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US2x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN4x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_938 = Coupling(name = 'GC_938',
                  value = '-(cw*ee*complex(0,1)*US3x1*complexconjugate(NN4x3)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN4x4)*complexconjugate(NN5x1))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN4x3)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN4x4)*complexconjugate(NN5x2))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*US3x1*complexconjugate(NN4x1)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (ee*complex(0,1)*US3x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN4x2)*complexconjugate(NN5x3))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN4x1)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN4x2)*complexconjugate(NN5x4))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN4x4)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN4x5)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN4x3)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN4x5)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN4x3)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN4x4)*complexconjugate(NN5x5))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN4x5)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_939 = Coupling(name = 'GC_939',
                  value = '-((cw*ee*UP1x1*complexconjugate(NN5x1)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw))) + (ee*UP1x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP1x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP1x2*complexconjugate(NN5x1)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (ee*UP1x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP1x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) + (NMl*UP1x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP1x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP1x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP1x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP1x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_940 = Coupling(name = 'GC_940',
                  value = '-((cw*ee*UP2x1*complexconjugate(NN5x1)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw))) + (ee*UP2x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*sw*(1 + sw)) - (ee*sw*UP2x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*UP2x2*complexconjugate(NN5x1)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (ee*UP2x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*sw*(1 + sw)) + (ee*sw*UP2x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) + (NMl*UP2x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMl*UP2x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMl*sw**2*UP2x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (NMk*UP2x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (NMk*sw**2*UP2x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_941 = Coupling(name = 'GC_941',
                  value = '-((cw*ee*complex(0,1)*US1x1*complexconjugate(NN5x1)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US1x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US1x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US1x2*complexconjugate(NN5x1)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US1x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US1x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US1x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US1x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US1x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US1x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_942 = Coupling(name = 'GC_942',
                  value = '-((cw*ee*complex(0,1)*US2x1*complexconjugate(NN5x1)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US2x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US2x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US2x2*complexconjugate(NN5x1)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US2x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US2x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US2x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US2x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US2x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US2x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_943 = Coupling(name = 'GC_943',
                  value = '-((cw*ee*complex(0,1)*US3x1*complexconjugate(NN5x1)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw))) + (ee*complex(0,1)*US3x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*US3x1*complexconjugate(NN5x2)*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*US3x2*complexconjugate(NN5x1)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (ee*complex(0,1)*US3x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*sw*(1 + sw)) + (ee*complex(0,1)*sw*US3x2*complexconjugate(NN5x2)*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x3*complexconjugate(NN5x3)*complexconjugate(NN5x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x2*complexconjugate(NN5x3)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*US3x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl*sw**2*US3x1*complexconjugate(NN5x4)*complexconjugate(NN5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMk*US3x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMk*sw**2*US3x3*complexconjugate(NN5x5)**2*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_944 = Coupling(name = 'GC_944',
                  value = '-(complex(0,1)*G*complexconjugate(Rd1x1)*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_945 = Coupling(name = 'GC_945',
                  value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_946 = Coupling(name = 'GC_946',
                  value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_947 = Coupling(name = 'GC_947',
                  value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_948 = Coupling(name = 'GC_948',
                  value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_949 = Coupling(name = 'GC_949',
                  value = '(cw*ee*complex(0,1)*NN5x1*complexconjugate(Rd1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rd1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rd1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_950 = Coupling(name = 'GC_950',
                  value = '-(complex(0,1)*G*complexconjugate(Rd2x2)*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_951 = Coupling(name = 'GC_951',
                  value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_952 = Coupling(name = 'GC_952',
                  value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_953 = Coupling(name = 'GC_953',
                  value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_954 = Coupling(name = 'GC_954',
                  value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_955 = Coupling(name = 'GC_955',
                  value = '(cw*ee*complex(0,1)*NN5x1*complexconjugate(Rd2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rd2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rd2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_956 = Coupling(name = 'GC_956',
                  value = '-(complex(0,1)*G*complexconjugate(Rd3x3)*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_957 = Coupling(name = 'GC_957',
                  value = '(complex(0,1)*I3x33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_958 = Coupling(name = 'GC_958',
                  value = '(complex(0,1)*I3x33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_959 = Coupling(name = 'GC_959',
                  value = '(complex(0,1)*I3x33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_960 = Coupling(name = 'GC_960',
                  value = '(complex(0,1)*I3x33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_961 = Coupling(name = 'GC_961',
                  value = '(complex(0,1)*I3x33*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x33*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*complexconjugate(Rd3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rd3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rd3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_962 = Coupling(name = 'GC_962',
                  value = 'complex(0,1)*G*complexconjugate(Rd3x6)*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_963 = Coupling(name = 'GC_963',
                  value = '(complex(0,1)*I4x33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_964 = Coupling(name = 'GC_964',
                  value = '(complex(0,1)*I4x33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_965 = Coupling(name = 'GC_965',
                  value = '(complex(0,1)*I4x33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_966 = Coupling(name = 'GC_966',
                  value = '(complex(0,1)*I4x33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_967 = Coupling(name = 'GC_967',
                  value = '(complex(0,1)*I4x33*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x33*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rd3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_968 = Coupling(name = 'GC_968',
                  value = 'complex(0,1)*G*complexconjugate(Rd4x4)*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_969 = Coupling(name = 'GC_969',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_970 = Coupling(name = 'GC_970',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_971 = Coupling(name = 'GC_971',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_972 = Coupling(name = 'GC_972',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_973 = Coupling(name = 'GC_973',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rd4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_974 = Coupling(name = 'GC_974',
                  value = 'complex(0,1)*G*complexconjugate(Rd5x5)*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_975 = Coupling(name = 'GC_975',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_976 = Coupling(name = 'GC_976',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_977 = Coupling(name = 'GC_977',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_978 = Coupling(name = 'GC_978',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_979 = Coupling(name = 'GC_979',
                  value = '(cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rd5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_980 = Coupling(name = 'GC_980',
                  value = '-(complex(0,1)*G*complexconjugate(Rd6x3)*cmath.sqrt(2))',
                  order = {'QCD':1})

GC_981 = Coupling(name = 'GC_981',
                  value = '(complex(0,1)*I3x36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_982 = Coupling(name = 'GC_982',
                  value = '(complex(0,1)*I3x36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_983 = Coupling(name = 'GC_983',
                  value = '(complex(0,1)*I3x36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_984 = Coupling(name = 'GC_984',
                  value = '(complex(0,1)*I3x36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_985 = Coupling(name = 'GC_985',
                  value = '(complex(0,1)*I3x36*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I3x36*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*complexconjugate(Rd6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rd6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rd6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_986 = Coupling(name = 'GC_986',
                  value = 'complex(0,1)*G*complexconjugate(Rd6x6)*cmath.sqrt(2)',
                  order = {'QCD':1})

GC_987 = Coupling(name = 'GC_987',
                  value = '(complex(0,1)*I4x36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_988 = Coupling(name = 'GC_988',
                  value = '(complex(0,1)*I4x36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_989 = Coupling(name = 'GC_989',
                  value = '(complex(0,1)*I4x36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_990 = Coupling(name = 'GC_990',
                  value = '(complex(0,1)*I4x36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_991 = Coupling(name = 'GC_991',
                  value = '(complex(0,1)*I4x36*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I4x36*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rd6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                  order = {'QED':1})

GC_992 = Coupling(name = 'GC_992',
                  value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_993 = Coupling(name = 'GC_993',
                  value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_994 = Coupling(name = 'GC_994',
                  value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_995 = Coupling(name = 'GC_995',
                  value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_996 = Coupling(name = 'GC_996',
                  value = '-((cw*ee*complex(0,1)*NN5x1*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN5x2*complexconjugate(Rl1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rl1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_997 = Coupling(name = 'GC_997',
                  value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_998 = Coupling(name = 'GC_998',
                  value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_999 = Coupling(name = 'GC_999',
                  value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                  order = {'QED':1})

GC_1000 = Coupling(name = 'GC_1000',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1001 = Coupling(name = 'GC_1001',
                   value = '-((cw*ee*complex(0,1)*NN5x1*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) - (ee*complex(0,1)*NN5x2*complexconjugate(Rl2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rl2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1002 = Coupling(name = 'GC_1002',
                   value = '(complex(0,1)*I17x33*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x33*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1003 = Coupling(name = 'GC_1003',
                   value = '(complex(0,1)*I17x33*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x33*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1004 = Coupling(name = 'GC_1004',
                   value = '(complex(0,1)*I17x33*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x33*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1005 = Coupling(name = 'GC_1005',
                   value = '(complex(0,1)*I17x33*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x33*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1006 = Coupling(name = 'GC_1006',
                   value = '(complex(0,1)*I17x33*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x33*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN5x1*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rl3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rl3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1007 = Coupling(name = 'GC_1007',
                   value = '(complex(0,1)*I18x33*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x33*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1008 = Coupling(name = 'GC_1008',
                   value = '(complex(0,1)*I18x33*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x33*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1009 = Coupling(name = 'GC_1009',
                   value = '(complex(0,1)*I18x33*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x33*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1010 = Coupling(name = 'GC_1010',
                   value = '(complex(0,1)*I18x33*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x33*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1011 = Coupling(name = 'GC_1011',
                   value = '(complex(0,1)*I18x33*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x33*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rl3x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1012 = Coupling(name = 'GC_1012',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1013 = Coupling(name = 'GC_1013',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1014 = Coupling(name = 'GC_1014',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1015 = Coupling(name = 'GC_1015',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1016 = Coupling(name = 'GC_1016',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rl4x4)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1017 = Coupling(name = 'GC_1017',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1018 = Coupling(name = 'GC_1018',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1019 = Coupling(name = 'GC_1019',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1020 = Coupling(name = 'GC_1020',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1021 = Coupling(name = 'GC_1021',
                   value = '(cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rl5x5)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1022 = Coupling(name = 'GC_1022',
                   value = '(complex(0,1)*I17x36*NN1x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x36*NN1x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1023 = Coupling(name = 'GC_1023',
                   value = '(complex(0,1)*I17x36*NN2x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x36*NN2x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1024 = Coupling(name = 'GC_1024',
                   value = '(complex(0,1)*I17x36*NN3x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x36*NN3x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1025 = Coupling(name = 'GC_1025',
                   value = '(complex(0,1)*I17x36*NN4x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x36*NN4x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1026 = Coupling(name = 'GC_1026',
                   value = '(complex(0,1)*I17x36*NN5x3)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I17x36*NN5x3*sw**2)/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN5x1*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*complexconjugate(Rl6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rl6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1027 = Coupling(name = 'GC_1027',
                   value = '(complex(0,1)*I18x36*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x36*sw**2*complexconjugate(NN1x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1028 = Coupling(name = 'GC_1028',
                   value = '(complex(0,1)*I18x36*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x36*sw**2*complexconjugate(NN2x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1029 = Coupling(name = 'GC_1029',
                   value = '(complex(0,1)*I18x36*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x36*sw**2*complexconjugate(NN3x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1030 = Coupling(name = 'GC_1030',
                   value = '(complex(0,1)*I18x36*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x36*sw**2*complexconjugate(NN4x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1031 = Coupling(name = 'GC_1031',
                   value = '(complex(0,1)*I18x36*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I18x36*sw**2*complexconjugate(NN5x3))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Rl6x6)*cmath.sqrt(2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1032 = Coupling(name = 'GC_1032',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1033 = Coupling(name = 'GC_1033',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1034 = Coupling(name = 'GC_1034',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1035 = Coupling(name = 'GC_1035',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1036 = Coupling(name = 'GC_1036',
                   value = '-((cw*ee*complex(0,1)*NN5x1*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN5x2*complexconjugate(Rn1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rn1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1037 = Coupling(name = 'GC_1037',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1038 = Coupling(name = 'GC_1038',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1039 = Coupling(name = 'GC_1039',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1040 = Coupling(name = 'GC_1040',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1041 = Coupling(name = 'GC_1041',
                   value = '-((cw*ee*complex(0,1)*NN5x1*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN5x2*complexconjugate(Rn2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rn2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1042 = Coupling(name = 'GC_1042',
                   value = '-((cw*ee*complex(0,1)*NN1x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN1x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1043 = Coupling(name = 'GC_1043',
                   value = '-((cw*ee*complex(0,1)*NN2x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN2x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1044 = Coupling(name = 'GC_1044',
                   value = '-((cw*ee*complex(0,1)*NN3x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN3x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1045 = Coupling(name = 'GC_1045',
                   value = '-((cw*ee*complex(0,1)*NN4x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN4x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1046 = Coupling(name = 'GC_1046',
                   value = '-((cw*ee*complex(0,1)*NN5x1*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))) + (ee*complex(0,1)*NN5x2*complexconjugate(Rn3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Rn3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1047 = Coupling(name = 'GC_1047',
                   value = '-(complex(0,1)*G*complexconjugate(Ru1x1)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1048 = Coupling(name = 'GC_1048',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1049 = Coupling(name = 'GC_1049',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1050 = Coupling(name = 'GC_1050',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1051 = Coupling(name = 'GC_1051',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1052 = Coupling(name = 'GC_1052',
                   value = '(cw*ee*complex(0,1)*NN5x1*complexconjugate(Ru1x1))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*complexconjugate(Ru1x1))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Ru1x1))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1053 = Coupling(name = 'GC_1053',
                   value = '-(complex(0,1)*G*complexconjugate(Ru2x2)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1054 = Coupling(name = 'GC_1054',
                   value = '(cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1055 = Coupling(name = 'GC_1055',
                   value = '(cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1056 = Coupling(name = 'GC_1056',
                   value = '(cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1057 = Coupling(name = 'GC_1057',
                   value = '(cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1058 = Coupling(name = 'GC_1058',
                   value = '(cw*ee*complex(0,1)*NN5x1*complexconjugate(Ru2x2))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*complexconjugate(Ru2x2))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Ru2x2))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1059 = Coupling(name = 'GC_1059',
                   value = '-(complex(0,1)*G*complexconjugate(Ru3x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1060 = Coupling(name = 'GC_1060',
                   value = '(complex(0,1)*I37x33*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x33*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1061 = Coupling(name = 'GC_1061',
                   value = '(complex(0,1)*I37x33*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x33*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1062 = Coupling(name = 'GC_1062',
                   value = '(complex(0,1)*I37x33*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x33*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1063 = Coupling(name = 'GC_1063',
                   value = '(complex(0,1)*I37x33*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x33*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1064 = Coupling(name = 'GC_1064',
                   value = '(complex(0,1)*I37x33*NN5x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x33*NN5x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*complexconjugate(Ru3x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*complexconjugate(Ru3x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Ru3x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1065 = Coupling(name = 'GC_1065',
                   value = 'complex(0,1)*G*complexconjugate(Ru3x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1066 = Coupling(name = 'GC_1066',
                   value = '(complex(0,1)*I38x33*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x33*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1067 = Coupling(name = 'GC_1067',
                   value = '(complex(0,1)*I38x33*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x33*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1068 = Coupling(name = 'GC_1068',
                   value = '(complex(0,1)*I38x33*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x33*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1069 = Coupling(name = 'GC_1069',
                   value = '(complex(0,1)*I38x33*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x33*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1070 = Coupling(name = 'GC_1070',
                   value = '(complex(0,1)*I38x33*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x33*sw**2*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Ru3x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1071 = Coupling(name = 'GC_1071',
                   value = 'complex(0,1)*G*complexconjugate(Ru4x4)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1072 = Coupling(name = 'GC_1072',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1073 = Coupling(name = 'GC_1073',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1074 = Coupling(name = 'GC_1074',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1075 = Coupling(name = 'GC_1075',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1076 = Coupling(name = 'GC_1076',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Ru4x4)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1077 = Coupling(name = 'GC_1077',
                   value = 'complex(0,1)*G*complexconjugate(Ru5x5)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1078 = Coupling(name = 'GC_1078',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1079 = Coupling(name = 'GC_1079',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1080 = Coupling(name = 'GC_1080',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1081 = Coupling(name = 'GC_1081',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1082 = Coupling(name = 'GC_1082',
                   value = '(-2*cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Ru5x5)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1083 = Coupling(name = 'GC_1083',
                   value = '-(complex(0,1)*G*complexconjugate(Ru6x3)*cmath.sqrt(2))',
                   order = {'QCD':1})

GC_1084 = Coupling(name = 'GC_1084',
                   value = '(complex(0,1)*I37x36*NN1x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x36*NN1x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1085 = Coupling(name = 'GC_1085',
                   value = '(complex(0,1)*I37x36*NN2x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x36*NN2x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1086 = Coupling(name = 'GC_1086',
                   value = '(complex(0,1)*I37x36*NN3x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x36*NN3x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1087 = Coupling(name = 'GC_1087',
                   value = '(complex(0,1)*I37x36*NN4x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x36*NN4x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1088 = Coupling(name = 'GC_1088',
                   value = '(complex(0,1)*I37x36*NN5x4)/((-1 + sw)*(1 + sw)) - (complex(0,1)*I37x36*NN5x4*sw**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*complexconjugate(Ru6x3))/(3.*(-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*complexconjugate(Ru6x3))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*complexconjugate(Ru6x3))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1089 = Coupling(name = 'GC_1089',
                   value = 'complex(0,1)*G*complexconjugate(Ru6x6)*cmath.sqrt(2)',
                   order = {'QCD':1})

GC_1090 = Coupling(name = 'GC_1090',
                   value = '(complex(0,1)*I38x36*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x36*sw**2*complexconjugate(NN1x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1091 = Coupling(name = 'GC_1091',
                   value = '(complex(0,1)*I38x36*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x36*sw**2*complexconjugate(NN2x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1092 = Coupling(name = 'GC_1092',
                   value = '(complex(0,1)*I38x36*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x36*sw**2*complexconjugate(NN3x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1093 = Coupling(name = 'GC_1093',
                   value = '(complex(0,1)*I38x36*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x36*sw**2*complexconjugate(NN4x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1094 = Coupling(name = 'GC_1094',
                   value = '(complex(0,1)*I38x36*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (complex(0,1)*I38x36*sw**2*complexconjugate(NN5x4))/((-1 + sw)*(1 + sw)) - (2*cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(Ru6x6)*cmath.sqrt(2))/(3.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1095 = Coupling(name = 'GC_1095',
                   value = '-((ee*complex(0,1)*I23x11*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1096 = Coupling(name = 'GC_1096',
                   value = '-((ee*complex(0,1)*I23x22*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1097 = Coupling(name = 'GC_1097',
                   value = '-((ee*complex(0,1)*I7x11*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1098 = Coupling(name = 'GC_1098',
                   value = '-((ee*complex(0,1)*I7x22*complexconjugate(UU1x1))/sw)',
                   order = {'QED':1})

GC_1099 = Coupling(name = 'GC_1099',
                   value = 'complex(0,1)*I88x33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1100 = Coupling(name = 'GC_1100',
                   value = 'complex(0,1)*I91x33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1101 = Coupling(name = 'GC_1101',
                   value = 'complex(0,1)*I91x36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1102 = Coupling(name = 'GC_1102',
                   value = '-((ee*complex(0,1)*I23x33*complexconjugate(UU1x1))/sw) + complex(0,1)*I24x33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1103 = Coupling(name = 'GC_1103',
                   value = '-((ee*complex(0,1)*I23x36*complexconjugate(UU1x1))/sw) + complex(0,1)*I24x36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1104 = Coupling(name = 'GC_1104',
                   value = '-((ee*complex(0,1)*I7x33*complexconjugate(UU1x1))/sw) + complex(0,1)*I9x33*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1105 = Coupling(name = 'GC_1105',
                   value = '-((ee*complex(0,1)*I7x36*complexconjugate(UU1x1))/sw) + complex(0,1)*I9x36*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1106 = Coupling(name = 'GC_1106',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN1x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1107 = Coupling(name = 'GC_1107',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN2x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1108 = Coupling(name = 'GC_1108',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN3x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1109 = Coupling(name = 'GC_1109',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN4x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1110 = Coupling(name = 'GC_1110',
                   value = '-((ee*complex(0,1)*NN5x2*complexconjugate(UU1x1))/sw) - (ee*complex(0,1)*NN5x3*complexconjugate(UU1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1111 = Coupling(name = 'GC_1111',
                   value = 'ee*complex(0,1)*UU1x1*complexconjugate(UU1x1) + ee*complex(0,1)*UU1x2*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1112 = Coupling(name = 'GC_1112',
                   value = '-((cw*ee*complex(0,1)*UU1x1*complexconjugate(UU1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU1x1*complexconjugate(UU1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU1x2*complexconjugate(UU1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU1x2*complexconjugate(UU1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1113 = Coupling(name = 'GC_1113',
                   value = 'ee*complex(0,1)*UU2x1*complexconjugate(UU1x1) + ee*complex(0,1)*UU2x2*complexconjugate(UU1x2)',
                   order = {'QED':1})

GC_1114 = Coupling(name = 'GC_1114',
                   value = '-((cw*ee*complex(0,1)*UU2x1*complexconjugate(UU1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU2x1*complexconjugate(UU1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU2x2*complexconjugate(UU1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU2x2*complexconjugate(UU1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1115 = Coupling(name = 'GC_1115',
                   value = '-((ee*complex(0,1)*I23x11*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1116 = Coupling(name = 'GC_1116',
                   value = '-((ee*complex(0,1)*I23x22*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1117 = Coupling(name = 'GC_1117',
                   value = '-((ee*complex(0,1)*I7x11*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1118 = Coupling(name = 'GC_1118',
                   value = '-((ee*complex(0,1)*I7x22*complexconjugate(UU2x1))/sw)',
                   order = {'QED':1})

GC_1119 = Coupling(name = 'GC_1119',
                   value = 'complex(0,1)*I88x33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1120 = Coupling(name = 'GC_1120',
                   value = 'complex(0,1)*I91x33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1121 = Coupling(name = 'GC_1121',
                   value = 'complex(0,1)*I91x36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1122 = Coupling(name = 'GC_1122',
                   value = '-((ee*complex(0,1)*I23x33*complexconjugate(UU2x1))/sw) + complex(0,1)*I24x33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1123 = Coupling(name = 'GC_1123',
                   value = '-((ee*complex(0,1)*I23x36*complexconjugate(UU2x1))/sw) + complex(0,1)*I24x36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1124 = Coupling(name = 'GC_1124',
                   value = '-((ee*complex(0,1)*I7x33*complexconjugate(UU2x1))/sw) + complex(0,1)*I9x33*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1125 = Coupling(name = 'GC_1125',
                   value = '-((ee*complex(0,1)*I7x36*complexconjugate(UU2x1))/sw) + complex(0,1)*I9x36*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1126 = Coupling(name = 'GC_1126',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN1x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1127 = Coupling(name = 'GC_1127',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN2x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1128 = Coupling(name = 'GC_1128',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN3x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1129 = Coupling(name = 'GC_1129',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN4x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1130 = Coupling(name = 'GC_1130',
                   value = '-((ee*complex(0,1)*NN5x2*complexconjugate(UU2x1))/sw) - (ee*complex(0,1)*NN5x3*complexconjugate(UU2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1131 = Coupling(name = 'GC_1131',
                   value = 'ee*complex(0,1)*UU1x1*complexconjugate(UU2x1) + ee*complex(0,1)*UU1x2*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1132 = Coupling(name = 'GC_1132',
                   value = '-((cw*ee*complex(0,1)*UU1x1*complexconjugate(UU2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU1x1*complexconjugate(UU2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU1x2*complexconjugate(UU2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU1x2*complexconjugate(UU2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1133 = Coupling(name = 'GC_1133',
                   value = 'ee*complex(0,1)*UU2x1*complexconjugate(UU2x1) + ee*complex(0,1)*UU2x2*complexconjugate(UU2x2)',
                   order = {'QED':1})

GC_1134 = Coupling(name = 'GC_1134',
                   value = '-((cw*ee*complex(0,1)*UU2x1*complexconjugate(UU2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*UU2x1*complexconjugate(UU2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*UU2x2*complexconjugate(UU2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*UU2x2*complexconjugate(UU2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1135 = Coupling(name = 'GC_1135',
                   value = '-((ee*complex(0,1)*I31x11*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1136 = Coupling(name = 'GC_1136',
                   value = '-((ee*complex(0,1)*I31x22*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1137 = Coupling(name = 'GC_1137',
                   value = '-((ee*complex(0,1)*I31x33*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1138 = Coupling(name = 'GC_1138',
                   value = '-((ee*complex(0,1)*I51x11*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1139 = Coupling(name = 'GC_1139',
                   value = '-((ee*complex(0,1)*I51x22*complexconjugate(VV1x1))/sw)',
                   order = {'QED':1})

GC_1140 = Coupling(name = 'GC_1140',
                   value = 'complex(0,1)*I84x33*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1141 = Coupling(name = 'GC_1141',
                   value = 'complex(0,1)*I84x36*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1142 = Coupling(name = 'GC_1142',
                   value = '-((ee*complex(0,1)*I51x33*complexconjugate(VV1x1))/sw) + complex(0,1)*I53x33*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1143 = Coupling(name = 'GC_1143',
                   value = '-((ee*complex(0,1)*I51x36*complexconjugate(VV1x1))/sw) + complex(0,1)*I53x36*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1144 = Coupling(name = 'GC_1144',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN1x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1145 = Coupling(name = 'GC_1145',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN2x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1146 = Coupling(name = 'GC_1146',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN3x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1147 = Coupling(name = 'GC_1147',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN4x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1148 = Coupling(name = 'GC_1148',
                   value = '-((ee*complex(0,1)*NN5x2*complexconjugate(VV1x1))/sw) + (ee*complex(0,1)*NN5x4*complexconjugate(VV1x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1149 = Coupling(name = 'GC_1149',
                   value = 'ee*complex(0,1)*VV1x1*complexconjugate(VV1x1) + ee*complex(0,1)*VV1x2*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1150 = Coupling(name = 'GC_1150',
                   value = '-((cw*ee*complex(0,1)*VV1x1*complexconjugate(VV1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV1x1*complexconjugate(VV1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV1x2*complexconjugate(VV1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV1x2*complexconjugate(VV1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1151 = Coupling(name = 'GC_1151',
                   value = 'ee*complex(0,1)*VV2x1*complexconjugate(VV1x1) + ee*complex(0,1)*VV2x2*complexconjugate(VV1x2)',
                   order = {'QED':1})

GC_1152 = Coupling(name = 'GC_1152',
                   value = '-((cw*ee*complex(0,1)*VV2x1*complexconjugate(VV1x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV2x1*complexconjugate(VV1x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV2x2*complexconjugate(VV1x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV2x2*complexconjugate(VV1x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1153 = Coupling(name = 'GC_1153',
                   value = '-((ee*UP1x1*complexconjugate(UU1x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*UP1x2*complexconjugate(UU1x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(UU1x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1154 = Coupling(name = 'GC_1154',
                   value = '-((ee*UP2x1*complexconjugate(UU1x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*UP2x2*complexconjugate(UU1x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(UU1x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1155 = Coupling(name = 'GC_1155',
                   value = '-((ee*complex(0,1)*US1x1*complexconjugate(UU1x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*complexconjugate(UU1x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(UU1x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1156 = Coupling(name = 'GC_1156',
                   value = '-((ee*complex(0,1)*US2x1*complexconjugate(UU1x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*complexconjugate(UU1x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(UU1x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1157 = Coupling(name = 'GC_1157',
                   value = '-((ee*complex(0,1)*US3x1*complexconjugate(UU1x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*complexconjugate(UU1x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(UU1x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1158 = Coupling(name = 'GC_1158',
                   value = '-((ee*UP1x1*complexconjugate(UU2x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*UP1x2*complexconjugate(UU2x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(UU2x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1159 = Coupling(name = 'GC_1159',
                   value = '-((ee*UP2x1*complexconjugate(UU2x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*UP2x2*complexconjugate(UU2x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(UU2x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1160 = Coupling(name = 'GC_1160',
                   value = '-((ee*complex(0,1)*US1x1*complexconjugate(UU2x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*complexconjugate(UU2x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(UU2x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1161 = Coupling(name = 'GC_1161',
                   value = '-((ee*complex(0,1)*US2x1*complexconjugate(UU2x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*complexconjugate(UU2x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(UU2x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1162 = Coupling(name = 'GC_1162',
                   value = '-((ee*complex(0,1)*US3x1*complexconjugate(UU2x2)*complexconjugate(VV1x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*complexconjugate(UU2x1)*complexconjugate(VV1x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(UU2x2)*complexconjugate(VV1x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1163 = Coupling(name = 'GC_1163',
                   value = '-((ee*complex(0,1)*I31x11*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1164 = Coupling(name = 'GC_1164',
                   value = '-((ee*complex(0,1)*I31x22*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1165 = Coupling(name = 'GC_1165',
                   value = '-((ee*complex(0,1)*I31x33*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1166 = Coupling(name = 'GC_1166',
                   value = '-((ee*complex(0,1)*I51x11*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1167 = Coupling(name = 'GC_1167',
                   value = '-((ee*complex(0,1)*I51x22*complexconjugate(VV2x1))/sw)',
                   order = {'QED':1})

GC_1168 = Coupling(name = 'GC_1168',
                   value = 'complex(0,1)*I84x33*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1169 = Coupling(name = 'GC_1169',
                   value = 'complex(0,1)*I84x36*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1170 = Coupling(name = 'GC_1170',
                   value = '-((ee*complex(0,1)*I51x33*complexconjugate(VV2x1))/sw) + complex(0,1)*I53x33*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1171 = Coupling(name = 'GC_1171',
                   value = '-((ee*complex(0,1)*I51x36*complexconjugate(VV2x1))/sw) + complex(0,1)*I53x36*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1172 = Coupling(name = 'GC_1172',
                   value = '-((ee*complex(0,1)*NN1x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN1x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1173 = Coupling(name = 'GC_1173',
                   value = '-((ee*complex(0,1)*NN2x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN2x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1174 = Coupling(name = 'GC_1174',
                   value = '-((ee*complex(0,1)*NN3x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN3x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1175 = Coupling(name = 'GC_1175',
                   value = '-((ee*complex(0,1)*NN4x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN4x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1176 = Coupling(name = 'GC_1176',
                   value = '-((ee*complex(0,1)*NN5x2*complexconjugate(VV2x1))/sw) + (ee*complex(0,1)*NN5x4*complexconjugate(VV2x2))/(sw*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1177 = Coupling(name = 'GC_1177',
                   value = 'ee*complex(0,1)*VV1x1*complexconjugate(VV2x1) + ee*complex(0,1)*VV1x2*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1178 = Coupling(name = 'GC_1178',
                   value = '-((cw*ee*complex(0,1)*VV1x1*complexconjugate(VV2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV1x1*complexconjugate(VV2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV1x2*complexconjugate(VV2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV1x2*complexconjugate(VV2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1179 = Coupling(name = 'GC_1179',
                   value = 'ee*complex(0,1)*VV2x1*complexconjugate(VV2x1) + ee*complex(0,1)*VV2x2*complexconjugate(VV2x2)',
                   order = {'QED':1})

GC_1180 = Coupling(name = 'GC_1180',
                   value = '-((cw*ee*complex(0,1)*VV2x1*complexconjugate(VV2x1))/((-1 + sw)*sw*(1 + sw))) + (cw*ee*complex(0,1)*sw*VV2x1*complexconjugate(VV2x1))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*VV2x2*complexconjugate(VV2x2))/(2.*(-1 + sw)*sw*(1 + sw)) + (cw*ee*complex(0,1)*sw*VV2x2*complexconjugate(VV2x2))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1181 = Coupling(name = 'GC_1181',
                   value = '-((ee*UP1x1*complexconjugate(UU1x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*UP1x2*complexconjugate(UU1x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(UU1x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1182 = Coupling(name = 'GC_1182',
                   value = '-((ee*UP2x1*complexconjugate(UU1x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*UP2x2*complexconjugate(UU1x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(UU1x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1183 = Coupling(name = 'GC_1183',
                   value = '-((ee*complex(0,1)*US1x1*complexconjugate(UU1x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*complexconjugate(UU1x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(UU1x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1184 = Coupling(name = 'GC_1184',
                   value = '-((ee*complex(0,1)*US2x1*complexconjugate(UU1x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*complexconjugate(UU1x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(UU1x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1185 = Coupling(name = 'GC_1185',
                   value = '-((ee*complex(0,1)*US3x1*complexconjugate(UU1x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*complexconjugate(UU1x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(UU1x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1186 = Coupling(name = 'GC_1186',
                   value = '-((ee*UP1x1*complexconjugate(UU2x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*UP1x2*complexconjugate(UU2x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) + (NMl*UP1x3*complexconjugate(UU2x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1187 = Coupling(name = 'GC_1187',
                   value = '-((ee*UP2x1*complexconjugate(UU2x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*UP2x2*complexconjugate(UU2x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) + (NMl*UP2x3*complexconjugate(UU2x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1188 = Coupling(name = 'GC_1188',
                   value = '-((ee*complex(0,1)*US1x1*complexconjugate(UU2x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US1x2*complexconjugate(UU2x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US1x3*complexconjugate(UU2x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1189 = Coupling(name = 'GC_1189',
                   value = '-((ee*complex(0,1)*US2x1*complexconjugate(UU2x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US2x2*complexconjugate(UU2x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US2x3*complexconjugate(UU2x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1190 = Coupling(name = 'GC_1190',
                   value = '-((ee*complex(0,1)*US3x1*complexconjugate(UU2x2)*complexconjugate(VV2x1))/(sw*cmath.sqrt(2))) - (ee*complex(0,1)*US3x2*complexconjugate(UU2x1)*complexconjugate(VV2x2))/(sw*cmath.sqrt(2)) - (complex(0,1)*NMl*US3x3*complexconjugate(UU2x2)*complexconjugate(VV2x2))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1191 = Coupling(name = 'GC_1191',
                   value = '-((UP1x1*complexconjugate(yd3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1192 = Coupling(name = 'GC_1192',
                   value = '-((UP2x1*complexconjugate(yd3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1193 = Coupling(name = 'GC_1193',
                   value = '-((complex(0,1)*US1x1*complexconjugate(yd3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1194 = Coupling(name = 'GC_1194',
                   value = '-((complex(0,1)*US2x1*complexconjugate(yd3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1195 = Coupling(name = 'GC_1195',
                   value = '-((complex(0,1)*US3x1*complexconjugate(yd3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1196 = Coupling(name = 'GC_1196',
                   value = '-((UP1x1*complexconjugate(ye3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1197 = Coupling(name = 'GC_1197',
                   value = '-((UP2x1*complexconjugate(ye3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1198 = Coupling(name = 'GC_1198',
                   value = '-((complex(0,1)*US1x1*complexconjugate(ye3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1199 = Coupling(name = 'GC_1199',
                   value = '-((complex(0,1)*US2x1*complexconjugate(ye3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1200 = Coupling(name = 'GC_1200',
                   value = '-((complex(0,1)*US3x1*complexconjugate(ye3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1201 = Coupling(name = 'GC_1201',
                   value = '-((UP1x2*complexconjugate(yu3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1202 = Coupling(name = 'GC_1202',
                   value = '-((UP2x2*complexconjugate(yu3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1203 = Coupling(name = 'GC_1203',
                   value = '-((complex(0,1)*US1x2*complexconjugate(yu3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1204 = Coupling(name = 'GC_1204',
                   value = '-((complex(0,1)*US2x2*complexconjugate(yu3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1205 = Coupling(name = 'GC_1205',
                   value = '-((complex(0,1)*US3x2*complexconjugate(yu3x3))/cmath.sqrt(2))',
                   order = {'QED':1})

GC_1206 = Coupling(name = 'GC_1206',
                   value = 'complex(0,1)*I15x33*cmath.cos(beta)',
                   order = {'QED':1})

GC_1207 = Coupling(name = 'GC_1207',
                   value = 'complex(0,1)*I1x33*cmath.cos(beta)',
                   order = {'QED':1})

GC_1208 = Coupling(name = 'GC_1208',
                   value = '(complex(0,1)*I30x33*NMl*vs*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1209 = Coupling(name = 'GC_1209',
                   value = '(complex(0,1)*I30x36*NMl*vs*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1210 = Coupling(name = 'GC_1210',
                   value = '(complex(0,1)*I34x33*NMl*vs*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1211 = Coupling(name = 'GC_1211',
                   value = '(complex(0,1)*I34x36*NMl*vs*cmath.cos(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1212 = Coupling(name = 'GC_1212',
                   value = 'complex(0,1)*I14x33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1213 = Coupling(name = 'GC_1213',
                   value = 'complex(0,1)*I16x33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1214 = Coupling(name = 'GC_1214',
                   value = 'complex(0,1)*I2x33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1215 = Coupling(name = 'GC_1215',
                   value = 'complex(0,1)*I99x33*cmath.sin(beta)',
                   order = {'QED':1})

GC_1216 = Coupling(name = 'GC_1216',
                   value = '(ee*UP1x2*cmath.cos(beta))/(2.*sw) + (ee*UP1x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1217 = Coupling(name = 'GC_1217',
                   value = '-(ee**2*UP1x2*cmath.cos(beta))/(2.*sw) - (ee**2*UP1x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1218 = Coupling(name = 'GC_1218',
                   value = '(ee**2*UP1x2*cmath.cos(beta))/(2.*sw) + (ee**2*UP1x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1219 = Coupling(name = 'GC_1219',
                   value = '-(cw*ee**2*UP1x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*UP1x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1220 = Coupling(name = 'GC_1220',
                   value = '(cw*ee**2*UP1x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*UP1x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1221 = Coupling(name = 'GC_1221',
                   value = '(ee*UP2x2*cmath.cos(beta))/(2.*sw) + (ee*UP2x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1222 = Coupling(name = 'GC_1222',
                   value = '-(ee**2*UP2x2*cmath.cos(beta))/(2.*sw) - (ee**2*UP2x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1223 = Coupling(name = 'GC_1223',
                   value = '(ee**2*UP2x2*cmath.cos(beta))/(2.*sw) + (ee**2*UP2x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1224 = Coupling(name = 'GC_1224',
                   value = '-(cw*ee**2*UP2x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*UP2x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1225 = Coupling(name = 'GC_1225',
                   value = '(cw*ee**2*UP2x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) + (cw*ee**2*UP2x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1226 = Coupling(name = 'GC_1226',
                   value = '-(ee*complex(0,1)*US1x2*cmath.cos(beta))/(2.*sw) + (ee*complex(0,1)*US1x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1227 = Coupling(name = 'GC_1227',
                   value = '(ee**2*complex(0,1)*US1x2*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*US1x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1228 = Coupling(name = 'GC_1228',
                   value = '(cw*ee**2*complex(0,1)*US1x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*US1x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1229 = Coupling(name = 'GC_1229',
                   value = '-(ee*complex(0,1)*US2x2*cmath.cos(beta))/(2.*sw) + (ee*complex(0,1)*US2x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1230 = Coupling(name = 'GC_1230',
                   value = '(ee**2*complex(0,1)*US2x2*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*US2x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1231 = Coupling(name = 'GC_1231',
                   value = '(cw*ee**2*complex(0,1)*US2x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*US2x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1232 = Coupling(name = 'GC_1232',
                   value = '-(ee*complex(0,1)*US3x2*cmath.cos(beta))/(2.*sw) + (ee*complex(0,1)*US3x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1233 = Coupling(name = 'GC_1233',
                   value = '(ee**2*complex(0,1)*US3x2*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*US3x1*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':2})

GC_1234 = Coupling(name = 'GC_1234',
                   value = '(cw*ee**2*complex(0,1)*US3x2*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*US3x1*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1235 = Coupling(name = 'GC_1235',
                   value = '(complex(0,1)*NMl*NN1x5*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*sw**2*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1236 = Coupling(name = 'GC_1236',
                   value = '(complex(0,1)*NMl*NN2x5*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*sw**2*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1237 = Coupling(name = 'GC_1237',
                   value = '(complex(0,1)*NMl*NN3x5*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*sw**2*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1238 = Coupling(name = 'GC_1238',
                   value = '(complex(0,1)*NMl*NN4x5*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*sw**2*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1239 = Coupling(name = 'GC_1239',
                   value = '(complex(0,1)*NMl*NN5x5*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x5*sw**2*UU1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN5x3*UU1x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x3*sw*UU1x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN5x1*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*UU1x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*UU1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1240 = Coupling(name = 'GC_1240',
                   value = '(complex(0,1)*NMl*NN1x5*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*sw**2*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN1x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN1x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1241 = Coupling(name = 'GC_1241',
                   value = '(complex(0,1)*NMl*NN2x5*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*sw**2*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN2x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN2x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1242 = Coupling(name = 'GC_1242',
                   value = '(complex(0,1)*NMl*NN3x5*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*sw**2*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN3x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN3x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1243 = Coupling(name = 'GC_1243',
                   value = '(complex(0,1)*NMl*NN4x5*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*sw**2*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN4x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN4x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1244 = Coupling(name = 'GC_1244',
                   value = '(complex(0,1)*NMl*NN5x5*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x5*sw**2*UU2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*NN5x3*UU2x1*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x3*sw*UU2x1*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*NN5x1*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*UU2x2*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*sw*UU2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1245 = Coupling(name = 'GC_1245',
                   value = '-(ee**2*complex(0,1)*I27x11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I27x11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1246 = Coupling(name = 'GC_1246',
                   value = '-(ee**2*complex(0,1)*I27x22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I27x22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1247 = Coupling(name = 'GC_1247',
                   value = '-(ee**2*complex(0,1)*I27x33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I28x33*cmath.sin(beta) + (complex(0,1)*I29x33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I27x33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1248 = Coupling(name = 'GC_1248',
                   value = '-(ee**2*complex(0,1)*I27x36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I28x36*cmath.sin(beta) + (complex(0,1)*I29x36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I27x36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1249 = Coupling(name = 'GC_1249',
                   value = '-(ee**2*complex(0,1)*I33x11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I33x11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1250 = Coupling(name = 'GC_1250',
                   value = '-(ee**2*complex(0,1)*I33x22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I33x22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1251 = Coupling(name = 'GC_1251',
                   value = '-(ee**2*complex(0,1)*I33x33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I35x33*cmath.sin(beta) + (complex(0,1)*I36x33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I33x33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1252 = Coupling(name = 'GC_1252',
                   value = '-(ee**2*complex(0,1)*I33x36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I35x36*cmath.sin(beta) + (complex(0,1)*I36x36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I33x36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1253 = Coupling(name = 'GC_1253',
                   value = '-(ee**2*complex(0,1)*I41x11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I41x11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1254 = Coupling(name = 'GC_1254',
                   value = '-(ee**2*complex(0,1)*I41x22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I41x22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1255 = Coupling(name = 'GC_1255',
                   value = '-(ee**2*complex(0,1)*I54x11*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I54x11*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1256 = Coupling(name = 'GC_1256',
                   value = '-(ee**2*complex(0,1)*I54x22*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) - (ee**2*complex(0,1)*I54x22*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1257 = Coupling(name = 'GC_1257',
                   value = '(ee**2*complex(0,1)*vu*cmath.cos(beta))/(2.*sw) - (ee**2*complex(0,1)*vd*cmath.sin(beta))/(2.*sw)',
                   order = {'QED':1})

GC_1258 = Coupling(name = 'GC_1258',
                   value = '(cw*ee**2*complex(0,1)*vu*cmath.cos(beta))/(2.*(-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*vd*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1259 = Coupling(name = 'GC_1259',
                   value = '(complex(0,1)*I45x33*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I43x33*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1260 = Coupling(name = 'GC_1260',
                   value = '(complex(0,1)*I45x36*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I43x36*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1261 = Coupling(name = 'GC_1261',
                   value = '(complex(0,1)*I45x63*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I43x63*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1262 = Coupling(name = 'GC_1262',
                   value = '(complex(0,1)*I45x66*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I43x66*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1263 = Coupling(name = 'GC_1263',
                   value = '(complex(0,1)*I55x33*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I60x33*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1264 = Coupling(name = 'GC_1264',
                   value = '(complex(0,1)*I55x36*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I60x36*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1265 = Coupling(name = 'GC_1265',
                   value = '(complex(0,1)*I55x63*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I60x63*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1266 = Coupling(name = 'GC_1266',
                   value = '(complex(0,1)*I55x66*NMl*vs*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I60x66*NMl*vs*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':2})

GC_1267 = Coupling(name = 'GC_1267',
                   value = 'complex(0,1)*I42x33*cmath.cos(beta) + (complex(0,1)*I47x33*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I48x33*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I44x33*cmath.sin(beta) + (complex(0,1)*I46x33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I47x33*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1268 = Coupling(name = 'GC_1268',
                   value = 'complex(0,1)*I42x36*cmath.cos(beta) + (complex(0,1)*I47x36*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I48x36*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I44x36*cmath.sin(beta) + (complex(0,1)*I46x36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I47x36*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1269 = Coupling(name = 'GC_1269',
                   value = 'complex(0,1)*I42x63*cmath.cos(beta) + (complex(0,1)*I47x63*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I48x63*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x63*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I44x63*cmath.sin(beta) + (complex(0,1)*I46x63*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x63*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I47x63*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1270 = Coupling(name = 'GC_1270',
                   value = 'complex(0,1)*I42x66*cmath.cos(beta) + (complex(0,1)*I47x66*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I48x66*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x66*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I44x66*cmath.sin(beta) + (complex(0,1)*I46x66*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I41x66*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I47x66*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1271 = Coupling(name = 'GC_1271',
                   value = 'complex(0,1)*I57x33*cmath.cos(beta) + (complex(0,1)*I61x33*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I59x33*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x33*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I56x33*cmath.sin(beta) + (complex(0,1)*I58x33*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x33*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I61x33*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1272 = Coupling(name = 'GC_1272',
                   value = 'complex(0,1)*I57x36*cmath.cos(beta) + (complex(0,1)*I61x36*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I59x36*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x36*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I56x36*cmath.sin(beta) + (complex(0,1)*I58x36*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x36*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I61x36*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1273 = Coupling(name = 'GC_1273',
                   value = 'complex(0,1)*I57x63*cmath.cos(beta) + (complex(0,1)*I61x63*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I59x63*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x63*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I56x63*cmath.sin(beta) + (complex(0,1)*I58x63*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x63*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I61x63*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1274 = Coupling(name = 'GC_1274',
                   value = 'complex(0,1)*I57x66*cmath.cos(beta) + (complex(0,1)*I61x66*vd*cmath.cos(beta))/cmath.sqrt(2) + (complex(0,1)*I59x66*vu*cmath.cos(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x66*vu*cmath.cos(beta))/(2.*sw**2*cmath.sqrt(2)) + complex(0,1)*I56x66*cmath.sin(beta) + (complex(0,1)*I58x66*vd*cmath.sin(beta))/cmath.sqrt(2) - (ee**2*complex(0,1)*I54x66*vd*cmath.sin(beta))/(2.*sw**2*cmath.sqrt(2)) + (complex(0,1)*I61x66*vu*cmath.sin(beta))/cmath.sqrt(2)',
                   order = {'QED':1})

GC_1275 = Coupling(name = 'GC_1275',
                   value = '(ee*complex(0,1)*NN1x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*sw**2*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1276 = Coupling(name = 'GC_1276',
                   value = '(ee*complex(0,1)*NN2x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*sw**2*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1277 = Coupling(name = 'GC_1277',
                   value = '(ee*complex(0,1)*NN3x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*sw**2*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1278 = Coupling(name = 'GC_1278',
                   value = '(ee*complex(0,1)*NN4x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*sw**2*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1279 = Coupling(name = 'GC_1279',
                   value = '(ee*complex(0,1)*NN5x4*VV1x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x4*sw*VV1x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*VV1x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*VV1x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN5x5*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x5*sw**2*VV1x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1280 = Coupling(name = 'GC_1280',
                   value = '(ee*complex(0,1)*NN1x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN1x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN1x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN1x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN1x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN1x5*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN1x5*sw**2*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1281 = Coupling(name = 'GC_1281',
                   value = '(ee*complex(0,1)*NN2x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN2x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN2x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN2x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN2x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN2x5*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN2x5*sw**2*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1282 = Coupling(name = 'GC_1282',
                   value = '(ee*complex(0,1)*NN3x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN3x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN3x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN3x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN3x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN3x5*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN3x5*sw**2*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1283 = Coupling(name = 'GC_1283',
                   value = '(ee*complex(0,1)*NN4x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN4x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN4x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN4x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN4x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN4x5*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN4x5*sw**2*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1284 = Coupling(name = 'GC_1284',
                   value = '(ee*complex(0,1)*NN5x4*VV2x1*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*NN5x4*sw*VV2x1*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*NN5x1*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*NN5x2*VV2x2*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*NN5x2*sw*VV2x2*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*NN5x5*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*NN5x5*sw**2*VV2x2*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1285 = Coupling(name = 'GC_1285',
                   value = '(complex(0,1)*NMl*complexconjugate(NN1x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN1x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1286 = Coupling(name = 'GC_1286',
                   value = '(complex(0,1)*NMl*complexconjugate(NN2x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN2x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1287 = Coupling(name = 'GC_1287',
                   value = '(complex(0,1)*NMl*complexconjugate(NN3x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN3x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1288 = Coupling(name = 'GC_1288',
                   value = '(complex(0,1)*NMl*complexconjugate(NN4x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN4x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1289 = Coupling(name = 'GC_1289',
                   value = '(complex(0,1)*NMl*complexconjugate(NN5x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN5x5)*complexconjugate(UU1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN5x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN5x3)*complexconjugate(UU1x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN5x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN5x2)*complexconjugate(UU1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1290 = Coupling(name = 'GC_1290',
                   value = '(complex(0,1)*NMl*complexconjugate(NN1x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN1x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1291 = Coupling(name = 'GC_1291',
                   value = '(complex(0,1)*NMl*complexconjugate(NN2x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN2x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1292 = Coupling(name = 'GC_1292',
                   value = '(complex(0,1)*NMl*complexconjugate(NN3x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN3x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1293 = Coupling(name = 'GC_1293',
                   value = '(complex(0,1)*NMl*complexconjugate(NN4x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN4x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1294 = Coupling(name = 'GC_1294',
                   value = '(complex(0,1)*NMl*complexconjugate(NN5x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN5x5)*complexconjugate(UU2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (ee*complex(0,1)*complexconjugate(NN5x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN5x3)*complexconjugate(UU2x1)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*complexconjugate(NN5x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*sw*complexconjugate(NN5x2)*complexconjugate(UU2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2))',
                   order = {'QED':1})

GC_1295 = Coupling(name = 'GC_1295',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN1x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN1x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1296 = Coupling(name = 'GC_1296',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN2x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN2x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1297 = Coupling(name = 'GC_1297',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN3x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN3x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1298 = Coupling(name = 'GC_1298',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN4x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN4x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1299 = Coupling(name = 'GC_1299',
                   value = '(ee*complex(0,1)*complexconjugate(NN5x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN5x4)*complexconjugate(VV1x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN5x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN5x2)*complexconjugate(VV1x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN5x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN5x5)*complexconjugate(VV1x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1300 = Coupling(name = 'GC_1300',
                   value = '(ee*complex(0,1)*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN1x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN1x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN1x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN1x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN1x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1301 = Coupling(name = 'GC_1301',
                   value = '(ee*complex(0,1)*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN2x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN2x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN2x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN2x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN2x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1302 = Coupling(name = 'GC_1302',
                   value = '(ee*complex(0,1)*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN3x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN3x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN3x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN3x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN3x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1303 = Coupling(name = 'GC_1303',
                   value = '(ee*complex(0,1)*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN4x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN4x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN4x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN4x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN4x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1304 = Coupling(name = 'GC_1304',
                   value = '(ee*complex(0,1)*complexconjugate(NN5x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)) - (ee*complex(0,1)*sw*complexconjugate(NN5x4)*complexconjugate(VV2x1)*cmath.cos(beta))/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*complexconjugate(NN5x1)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (ee*complex(0,1)*complexconjugate(NN5x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*sw*(1 + sw)*cmath.sqrt(2)) - (ee*complex(0,1)*sw*complexconjugate(NN5x2)*complexconjugate(VV2x2)*cmath.cos(beta))/((-1 + sw)*(1 + sw)*cmath.sqrt(2)) + (complex(0,1)*NMl*complexconjugate(NN5x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl*sw**2*complexconjugate(NN5x5)*complexconjugate(VV2x2)*cmath.sin(beta))/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1305 = Coupling(name = 'GC_1305',
                   value = '-(ee*complex(0,1)*cmath.cos(beta)**2) - ee*complex(0,1)*cmath.sin(beta)**2',
                   order = {'QED':1})

GC_1306 = Coupling(name = 'GC_1306',
                   value = '2*ee**2*complex(0,1)*cmath.cos(beta)**2 + 2*ee**2*complex(0,1)*cmath.sin(beta)**2',
                   order = {'QED':2})

GC_1307 = Coupling(name = 'GC_1307',
                   value = '(ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*sw**2) + (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*sw**2)',
                   order = {'QED':2})

GC_1308 = Coupling(name = 'GC_1308',
                   value = '(cw*ee*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*sw*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (cw*ee*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw*(1 + sw)) - (cw*ee*complex(0,1)*sw*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':1})

GC_1309 = Coupling(name = 'GC_1309',
                   value = '-((cw*ee**2*complex(0,1)*cmath.cos(beta)**2)/((-1 + sw)*sw*(1 + sw))) + (2*cw*ee**2*complex(0,1)*sw*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (cw*ee**2*complex(0,1)*cmath.sin(beta)**2)/((-1 + sw)*sw*(1 + sw)) + (2*cw*ee**2*complex(0,1)*sw*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1310 = Coupling(name = 'GC_1310',
                   value = '(2*ee**2*complex(0,1)*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.cos(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*sw**2*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (2*ee**2*complex(0,1)*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*cmath.sin(beta)**2)/(2.*(-1 + sw)*sw**2*(1 + sw)) - (2*ee**2*complex(0,1)*sw**2*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1311 = Coupling(name = 'GC_1311',
                   value = '(complex(0,1)*NMl**2*US1x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*US1x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*US1x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US1x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US1x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1312 = Coupling(name = 'GC_1312',
                   value = '(complex(0,1)*NMl**2*US2x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*US2x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*US2x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US2x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US2x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1313 = Coupling(name = 'GC_1313',
                   value = '(complex(0,1)*NMl**2*US3x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US3x3*vs*cmath.cos(beta)**2)/((-1 + sw)*(1 + sw)) + (2*complex(0,1)*NMk*NMl*US3x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (2*complex(0,1)*NMk*NMl*sw**2*US3x3*vs*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMl**2*US3x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*sw**2*US3x3*vs*cmath.sin(beta)**2)/((-1 + sw)*(1 + sw))',
                   order = {'QED':2})

GC_1314 = Coupling(name = 'GC_1314',
                   value = '-(ee**2*complex(0,1)*US1x1*vd*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US1x1*vd*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*US1x2*vu*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US1x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US1x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US1x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US1x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US1x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US1x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US1x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*US1x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*US1x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US1x1*vd*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US1x2*vu*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US1x2*vu*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1315 = Coupling(name = 'GC_1315',
                   value = '-(ee**2*complex(0,1)*US2x1*vd*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US2x1*vd*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*US2x2*vu*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US2x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US2x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US2x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US2x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US2x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US2x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US2x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*US2x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*US2x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US2x1*vd*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US2x2*vu*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US2x2*vu*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

GC_1316 = Coupling(name = 'GC_1316',
                   value = '-(ee**2*complex(0,1)*US3x1*vd*cmath.cos(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US3x1*vd*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) + (ee**2*complex(0,1)*US3x2*vu*cmath.cos(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US3x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US3x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US3x2*vd*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US3x2*vd*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (ee**2*complex(0,1)*US3x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*(1 + sw)) - (complex(0,1)*NMl**2*US3x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US3x1*vu*cmath.cos(beta)*cmath.sin(beta))/(2.*(-1 + sw)*sw**2*(1 + sw)) + (complex(0,1)*NMl**2*sw**2*US3x1*vu*cmath.cos(beta)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (complex(0,1)*NMAl*NMl*US3x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) - (complex(0,1)*NMAl*NMl*sw**2*US3x3*cmath.cos(beta)*cmath.sqrt(2)*cmath.sin(beta))/((-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US3x1*vd*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw)) - (ee**2*complex(0,1)*US3x2*vu*cmath.sin(beta)**2)/(2.*(-1 + sw)*(1 + sw)) + (ee**2*complex(0,1)*US3x2*vu*cmath.sin(beta)**2)/(4.*(-1 + sw)*sw**2*(1 + sw))',
                   order = {'QED':1})

