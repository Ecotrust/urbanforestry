import math
from sympy.solvers import solve
from sympy import Symbol


class Equation(object):
    '''
    Creates an Equation class, which records the parameters to predict some 
    metric of a tree given a single variable measured from a tree. 
    
    Keyword arguments:
    ind_var -- the independent variable used in the equation, as a string
    pred_var -- the dependent variable predicted by the equation, as a string
    eq_form -- the form of equation (e.g., 'linear', 'quadratic'...) as a string
    param_a -- the first equation parameter (usually the y-intercept)
    param_b... e -- the second through fifth equation parameters
    '''
    def __init__(self, eq_form, param_a, param_b, param_c = None, 
                 param_d = None, param_e = None):
        self.eq_form = eq_form
        self.a = param_a
        self.b = param_b
        self.c = param_c
        self.d = param_d
        self.e = param_e
        
    def calc(self, measurement, inverse = False):
        '''
        Using the implied independent variable, calculates the dependent variable.
        Alternatively, if inverse = True, solves for the independent variable.
        
        Keyword arguments:
        measured_thing -- the value of the measured tree attribute
        inverse -- inverts the equation to solve for the independent variable
        '''
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        
        if not inverse: # we're solving for y using x
            x = float(measurement)
            if self.eq_form == 'linear' or self.eq_form == 'lin':
                value = a + b*x
            elif self.eq_form == 'quadratic' or self.eq_form == 'quad':
                value = a + b*x + c*x**2
            elif self.eq_form == 'cubic' or self.eq_form == 'cub':
                value = a + b*x + c*x**2 + d*x**3
            elif self.eq_form == 'quartic'  or self.eq_form == 'quart':
                value = a + b*x + c*x**2 + d*x**3 + e*x**4

            # note that math.log is natural log                
            elif self.eq_form == 'loglog1':
                value = math.exp(a + b*math.log(math.log(x+1)) + c/2.0)
            elif self.eq_form == 'loglog2':
                value = math.exp(a + b*math.log(math.log(x+1)) + math.sqrt(x)*(c/2.0))
            elif self.eq_form == 'loglog3':
                value = math.exp(a + b*math.log(math.log(x+1)) + x*(c/2.0))
            elif self.eq_form == 'loglog4':
                value = math.exp(a + b*math.log(math.log(x+1)) + x**2 * (c/2.0))
            elif self.eq_form == 'expow1':
                value = math.exp(a + b*x + c/2.0)
            elif self.eq_form == 'expow2':
                value = math.exp(a + b*x + math.sqrt(x)*(c/2.0))
            elif self.eq_form == 'expow3':
                value = math.exp(a + b*x + x*(c/2.0))
            elif self.eq_form == 'expow4':
                value = math.exp(a + b*x + (x**2)*(c/2.0))
        
        elif inverse: # we're solving for x using y
            y = float(measurement)
            x = Symbol('x')
            if self.eq_form == 'linear':
                value = solve (a + b*x - y, x)
            elif self.eq_form == 'quadratic':
                value = solve(a + b*x + c*x**2 - y, x)
            elif self.eq_form == 'cubic':
                value = solve(a + b*x + c*x**2 + d*x**3 - y, x)
            elif self.eq_form == 'quartic':
                value = solve(a + b*x + c*x**2 + d*x**3 + e*x**4 - y, x)

            # note that math.log is natural log
            elif self.eq_form == 'loglog1':                
                value = solve(math.exp(a + b*math.log(math.log(x+1)) + c/2.0) - y, x)
            elif self.eq_form == 'loglog2':                
                value = solve(math.exp(a + b*math.log(math.log(x+1)) + math.sqrt(x)*(c/2.0)) - y, x)
            elif self.eq_form == 'loglog3':                
                value = solve(math.exp(a + b*math.log(math.log(x+1)) + x*(c/2.0)) - y, x)
            elif self.eq_form == 'loglog4':                
                value = solve(math.exp(a + b*math.log(math.log(x+1)) + x**2 * (c/2.0)) - y, x)
            elif self.eq_form == 'expow1':
                value = solve(math.exp(a + b*x + c/2.0) - y, x)
            elif self.eq_form == 'expow2':
                value = solve(math.exp(a + b*x + math.sqrt(x)*(c/2.0)) - y, x)
            elif self.eq_form == 'expow3':
                value = solve(math.exp(a + b*x + x*(c/2.0)) - y, x)
            elif self.eq_form == 'expow4':
                value = solve(math.exp(a + b*x + (x**2)*(c/2.0)) - y, x)            
        
        if type(value) is list:
            value = value[0]
        return value
        

class Species(object):
    '''
    Creates a Species class, with option to record the species name. 
    
    Keyword arguments:
    sp_name -- optional name as a string, such as common name or latin name.
    '''
    def __init__(self, sp_name = None):
        self.sp_name = sp_name
        
    def add_eq(self, ind_var, pred_var, eq_form, param_a, param_b, 
               param_c = None, param_d = None, param_e = None):
        '''
        Adds an Equation object as an attribute of a Species object. The name 
        of the attribute is set as <ind_var>_to_<pred_var> (e.g., dbh_to_ht).
        
        Keyword arguments:
        ind_var -- the independent variable used in the equation, as a string
        pred_var -- the dependent variable predicted by the equation, as a string
        eq_form -- the form of equation (e.g., 'linear', 'quadratic'...) as a string
        param_a -- the first equation parameter (usually the y-intercept)
        param_b... e -- the second through fifth equation parameters
        '''
        # All equation forms have at least two parameters, a and b, only some 
        # equations have additional parameters
        setattr(self, ind_var + '_to_' + pred_var, 
                Equation(eq_form, param_a, param_b, param_c, param_d, param_e))

        
'''
# Example usage to create an equation for a species
       
# First, create the Species object
ACRU = Species('Acer rubrum')

# Then, add a an equation (in this case, dbh_to_ht)
ACRU.add_eq('dbh', 'ht', 'quad', 2.15761, 0.29592, -0.00158)

# If you wanted to calculate ht from a tree with known dbh of 6.5cm:
ACRU.dbh_to_ht.calc(6.5)        

# If you wanted to calculate dbh from a tree with known ht of 10m:
ACRU.dbh_to_ht.calc(10, inverse=True) 
'''        
        
# The Tree Equations for the Pacific Northwest Region
# ...
        
        
ACMA= Species('Acer macrophyllum')
ACMA.add_eq('dbh','age','quad',0.57316,0.51132,0.00191)
ACMA.add_eq('cdia','dbh','cub',0.2992,1.76512,0.55846,-0.01836)
ACMA.add_eq('dbh','crown dia','quad',0.39369,0.24616,-0.00081)
ACMA.add_eq('dbh','crown ht','loglogw1',0.28858,1.74158,0.03472)
ACMA.add_eq('age','dbh','cub',-0.21317,1.97429,-0.01051,0.00003)
ACMA.add_eq('dbh','leaf area','loglogw1',-1.54932,5.67086,0.32352)
ACMA.add_eq('dbh','tree ht','loglogw1',0.85132,1.48977,0.02171)

ACPL= Species('Acer platanoides')
ACPL.add_eq('dbh','age','cub',0.69436,0.34048,0.00733,-0.00004)
ACPL.add_eq('cdia','dbh','cub',0.75329,1.40261,0.35382,-0.00942)
ACPL.add_eq('dbh','crown dia','cub',0.0926,0.36724,-0.00266,0.00001)
ACPL.add_eq('dbh','crown ht','quad',1.83357,0.23724,-0.00089)
ACPL.add_eq('age','dbh','cub',-0.85211,2.4427,-0.02993,0.00022)
ACPL.add_eq('dbh','leaf area','loglogw1',-0.65552,5.15935,0.25353)
ACPL.add_eq('dbh','tree ht','quad',2.56416,0.3118,-0.00127)

ACRU= Species('Acer rubrum')
ACRU.add_eq('dbh','age','quad',0.19481,0.65916,-0.00084)
ACRU.add_eq('cdia','dbh','cub',0.64872,1.37806,0.27303,-0.00666)
ACRU.add_eq('dbh','crown dia','cub',0.16933,0.43289,-0.00526,0.00003)
ACRU.add_eq('dbh','crown ht','cub',0.57374,0.52078,-0.00953,0.00007)
ACRU.add_eq('age','dbh','lin',-0.28784,1.61264)
ACRU.add_eq('dbh','leaf area','loglogw1',-0.4824,4.97236,0.24582)
ACRU.add_eq('dbh','tree ht','cub',1.95274,0.61834,-0.01179,0.00008)

ACSA2= Species('Acer saccharum')
ACSA2.add_eq('dbh','age','quad',0.35741,0.33851,0.00463)
ACSA2.add_eq('cdia','dbh','cub',0.33495,3.17841,0.28395,-0.0112)
ACSA2.add_eq('dbh','crown dia','quad',-0.04297,0.26749,-0.00063)
ACSA2.add_eq('dbh','crown ht','quad',1.2475,0.32945,-0.00126)
ACSA2.add_eq('age','dbh','cub',1.47953,1.65885,0.00513,-0.00015)
ACSA2.add_eq('dbh','leaf area','loglogw3',-1.68905,5.95471,0.00136)
ACSA2.add_eq('dbh','tree ht','cub',1.10482,0.62678,-0.0077,0.00004)

BEPE= Species('Betula pendula')
BEPE.add_eq('dbh','age','quad',0.39331,0.52785,0.00437)
BEPE.add_eq('cdia','dbh','cub',0.82226,0.2215,0.8505,-0.04466)
BEPE.add_eq('dbh','crown dia','cub',0.35044,0.37618,-0.00631,0.00007)
BEPE.add_eq('dbh','crown ht','quad',1.14308,0.53014,-0.00455)
BEPE.add_eq('age','dbh','quad',-0.91466,2.11608,-0.02273)
BEPE.add_eq('dbh','leaf area','quad',-0.05315,0.15892,0.2975)
BEPE.add_eq('dbh','tree ht','quad',1.42669,0.70065,-0.00668)

CABEF= Species('Carpinus betulus Fastigiata')
CABEF.add_eq('dbh','age','cub',-0.97624,1.02991,-0.01942,0.00016)
CABEF.add_eq('cdia','dbh','loglogw1',2.31114,1.4967,0.03034)
CABEF.add_eq('dbh','crown dia','cub',0.124,0.11118,0.00861,-0.00015)
CABEF.add_eq('dbh','crown ht','loglogw1',0.21394,1.66187,0.02607)
CABEF.add_eq('age','dbh','cub',-0.99936,2.16074,-0.08853,0.00412)
CABEF.add_eq('dbh','leaf area','loglogw1',-3.09187,6.99008,0.17664)
CABEF.add_eq('dbh','tree ht','loglogw1',0.68427,1.47106,0.01592)

CADE2= Species('Calocedrus decurrens')
CADE2.add_eq('dbh','age','cub',-3.92586,1.19536,-0.01523,0.00007)
CADE2.add_eq('cdia','dbh','cub',-0.3338,-4.50349,14.53807,-2.22352)
CADE2.add_eq('dbh','crown dia','loglogw1',-0.81693,1.39977,0.05206)
CADE2.add_eq('dbh','crown ht','cub',1.04793,0.20865,0.00281,-0.00004)
CADE2.add_eq('age','dbh','lin',-0.5418,2.23052)
CADE2.add_eq('dbh','leaf area','loglogw1',-3.08102,5.42583,0.12706)
CADE2.add_eq('dbh','tree ht','quad',-0.43078,0.46333,-0.00309)

CRLA80= Species('Crataegus laevigata')
CRLA80.add_eq('dbh','age','loglogw1',0.66848,2.26341,0.21909)
CRLA80.add_eq('cdia','dbh','cub',1.02603,-0.06687,0.77348,-0.03728)
CRLA80.add_eq('dbh','crown dia','quad',0.17454,0.40208,-0.00444)
CRLA80.add_eq('dbh','crown ht','quad',0.93095,0.25848,-0.0016)
CRLA80.add_eq('age','dbh','cub',0.20544,0.81745,0.01594,-0.00034)
CRLA80.add_eq('dbh','leaf area','cub',1.26793,-1.41462,0.39699,-0.00419)
CRLA80.add_eq('dbh','tree ht','quad',1.88852,0.33423,-0.00252)

FASYAT= Species('Fagus sylvatica atropunicea')
FASYAT.add_eq('dbh','age','cub',1.00624,0.24853,0.0104,-0.0001)
FASYAT.add_eq('cdia','dbh','cub',0.8162,1.20719,0.26657,-0.00642)
FASYAT.add_eq('dbh','crown dia','cub',-0.02818,0.46803,-0.0055,0.00003)
FASYAT.add_eq('dbh','crown ht','cub',0.241,0.5817,-0.00797,0.00005)
FASYAT.add_eq('age','dbh','cub',-2.64454,3.28563,-0.09364,0.00158)
FASYAT.add_eq('dbh','leaf area','loglogw1',-1.00002,5.53664,0.09336)
FASYAT.add_eq('dbh','tree ht','cub',0.25817,0.79152,-0.0123,0.00007)

FRLA= Species('Fraxinus latifolia')
FRLA.add_eq('dbh','age','quad',0.54234,0.5756,0.00247)
FRLA.add_eq('cdia','dbh','loglogw1',1.76367,2.39953,0.06678)
FRLA.add_eq('dbh','crown dia','cub',0.18904,0.37563,-0.00399,0.00002)
FRLA.add_eq('dbh','crown ht','quad',0.53509,0.37254,-0.00163)
FRLA.add_eq('age','dbh','cub',-1.63283,2.29867,-0.02735,0.00016)
FRLA.add_eq('dbh','leaf area','loglogw2',-1.27832,5.54199,0.05757)
FRLA.add_eq('dbh','tree ht','quad',1.40753,0.49125,-0.0026)

ILOP= Species('Ilex opaca')
ILOP.add_eq('dbh','age','cub',-7.63703,3.39474,-0.1696,0.00288)
ILOP.add_eq('cdia','dbh','loglogw1',1.90024,2.00291,0.03392)
ILOP.add_eq('dbh','crown dia','quad',0.07202,0.28414,-0.00229)
ILOP.add_eq('dbh','crown ht','loglogw1',-0.77752,2.29877,0.05786)
ILOP.add_eq('age','dbh','cub',2.53924,0.07682,0.08877,-0.00172)
ILOP.add_eq('dbh','leaf area','loglogw1',-7.00874,8.2918,0.53597)
ILOP.add_eq('dbh','tree ht','loglogw1',-0.07365,1.89671,0.06102)

LIST= Species('Liquidambar styraciflua')
LIST.add_eq('dbh','age','cub',0.75623,0.27517,0.01399,-0.00016)
LIST.add_eq('cdia','dbh','lin',-1.81154,4.27467)
LIST.add_eq('dbh','crown dia','quad',0.41201,0.27728,-0.00101)
LIST.add_eq('dbh','crown ht','quad',0.89761,0.40549,-0.00264)
LIST.add_eq('age','dbh','lin',2.02364,1.64102)
LIST.add_eq('dbh','leaf area','cub',0.61572,-0.63818,0.47404,-0.00342)
LIST.add_eq('dbh','tree ht','quad',1.37134,0.5091,-0.00376)

MOAL= Species('Morus alba')
MOAL.add_eq('dbh','age','cub',-0.84292,1.65967,-0.0534,0.00062)
MOAL.add_eq('cdia','dbh','cub',0.69882,0.38155,0.30043,-0.00873)
MOAL.add_eq('dbh','crown dia','quad',0.403768894,0.452307942,-0.002497943)
MOAL.add_eq('dbh','crown ht','quad',0.69217,0.33321,-0.00245)
MOAL.add_eq('age','dbh','cub',-4.50256,4.49267,-0.14464,0.00188)
MOAL.add_eq('dbh','leaf area','loglogw1',0.21087,4.87716,0.11207)
MOAL.add_eq('dbh','tree ht','quad',1.50551,0.42182,-0.00304)

PHCA= Species('Phoenix canariensis')
PHCA.add_eq('dbh','crown dia','quad',0.807,1.028,-0.00301)
PHCA.add_eq('dbh','crown ht','quad',-0.205,0.897,-0.0212)
PHCA.add_eq('dbh','leaf area','cub',-12.237,14.0084,0.4132,-0.0292)
PHCA.add_eq('dbh','tree ht','quad',-0.365,0.453,-0.00235)

PHDA4= Species('Phoenix dactylifera')
PHDA4.add_eq('dbh','crown dia','cub',0.01576,1.002,-0.073,0.00205)
PHDA4.add_eq('dbh','crown ht','cub',0.00664,0.335,0.00438,-0.00027)
PHDA4.add_eq('dbh','leaf area','cub',-0.00432,2.477,-0.0315,0.00115)
PHDA4.add_eq('dbh','tree ht','cub',0,0.259,0.00662,-0.00009)

PICO5= Species('Pinus contorta var. bolanderi')
PICO5.add_eq('dbh','age','cub',1.04669,-0.1142,0.03967,-0.00051)
PICO5.add_eq('cdia','dbh','loglogw1',1.95606,2.16085,0.03133)
PICO5.add_eq('dbh','crown dia','quad',0.35277,0.23532,-0.00137)
PICO5.add_eq('dbh','crown ht','lin',1.54381,0.16978)
PICO5.add_eq('age','dbh','loglogw2',1.57707,1.6895,0.00427)
PICO5.add_eq('dbh','leaf area','loglogw1',-0.76822,4.59751,0.14668)
PICO5.add_eq('dbh','tree ht','lin',1.55686,0.23113)

POTR2= Species('Populus balsamifera subsp. trichocarpa')
POTR2.add_eq('dbh','age','cub',1.49435,-0.11181,0.1546,-0.00359)
POTR2.add_eq('cdia','dbh','lin',-0.50224,3.16188)
POTR2.add_eq('dbh','crown dia','loglogw1',0.12163,1.62583,0.07523)
POTR2.add_eq('dbh','crown ht','quad',1.31743,0.49369,-0.00184)
POTR2.add_eq('age','dbh','cub',-0.52236,1.19511,0.01779,-0.00011)
POTR2.add_eq('dbh','leaf area','loglogw1',-0.36109,4.93446,0.49534)
POTR2.add_eq('dbh','tree ht','quad',1.60968,0.53395,-0.00187)

PRCEKW= Species('Prunus cerasifera Thundercloud')
PRCEKW.add_eq('dbh','age','cub',0.90691,0.16681,0.02458,-0.00026)
PRCEKW.add_eq('cdia','dbh','loglogw1',1.61707,2.40368,0.04956)
PRCEKW.add_eq('dbh','crown dia','quad',0.17605,0.37473,-0.00359)
PRCEKW.add_eq('dbh','crown ht','quad',1.13295,0.25941,-0.00271)
PRCEKW.add_eq('age','dbh','quad',0.78364,0.88584,-0.00949)
PRCEKW.add_eq('dbh','leaf area','loglogw1',-0.8793,4.90311,0.23211)
PRCEKW.add_eq('dbh','tree ht','quad',1.75135,0.35285,-0.00369)

PRSE2= Species('Prunus serrulata')
PRSE2.add_eq('dbh','age','cub',0.55903,0.81962,-0.00508,0.00002)
PRSE2.add_eq('cdia','dbh','quad',-2.77143,7.51392,-0.17692)
PRSE2.add_eq('dbh','crown dia','quad',0.40377,0.45231,-0.0025)
PRSE2.add_eq('dbh','crown ht','lin',1.10845,0.09734)
PRSE2.add_eq('age','dbh','lin',-0.65197,1.81238)
PRSE2.add_eq('dbh','leaf area','loglogw1',-1.93065,5.12856,0.74048)
PRSE2.add_eq('dbh','tree ht','cub',1.54418,0.31365,-0.00561,0.00004)

PSME= Species('Pseudotsuga menziesii')
PSME.add_eq('dbh','age','quad',0.8747,0.357,0.00306)
PSME.add_eq('cdia','dbh','cub',1.53581,0.2103,0.70805,-0.02168)
PSME.add_eq('dbh','crown dia','loglogw2',-0.82603,2.35246,0.00584)
PSME.add_eq('dbh','crown ht','quad',0.2225,0.49037,-0.0024)
PSME.add_eq('age','dbh','quad',-1.19688,2.36185,-0.01267)
PSME.add_eq('dbh','leaf area','loglogw1',-3.15072,6.98931,0.57231)
PSME.add_eq('dbh','tree ht','quad',0.04179,0.5776,-0.00274)

PYAN= Species('Malus angustifolia')
PYAN.add_eq('dbh','age','cub',1.494352444,-0.111810012,0.154595537,-0.003591378)
PYAN.add_eq('cdia','dbh','lin',-0.145962016,3.083196595)
PYAN.add_eq('dbh','crown dia','quad',0.017297264,0.489633526,-0.009720787)
PYAN.add_eq('dbh','crown ht','cub',1.10172,-0.26591,0.05319,-0.00153)
PYAN.add_eq('age','dbh','quad',0.78363989,0.885844908,-0.009487773)
PYAN.add_eq('dbh','leaf area','loglogw1',-1.65725,5.18953,0.71529)
PYAN.add_eq('dbh','tree ht','cub',1.67077,0.07431,0.03374,-0.00119)

PYKA= Species('Pyrus kawakamii')
PYKA.add_eq('dbh','age','cub',-0.73974,0.56603,0.01225,-0.00018)
PYKA.add_eq('cdia','dbh','lin',-4.89161,4.3895)
PYKA.add_eq('dbh','crown dia','quad',0.99604,0.29207,-0.00136)
PYKA.add_eq('dbh','crown ht','lin',1.64501,0.11198)
PYKA.add_eq('age','dbh','lin',2.83176,1.30003)
PYKA.add_eq('dbh','leaf area','loglogw1',-0.56606,4.36954,0.3546)
PYKA.add_eq('dbh','tree ht','quad',1.81068,0.32025,-0.00236)

QUAG= Species('Quercus agrifolia')
QUAG.add_eq('dbh','age','quad',-0.79896,0.80279,-0.00169)
QUAG.add_eq('cdia','dbh','loglogw1',1.62756,2.54237,0.0562)
QUAG.add_eq('dbh','crown dia','quad',0.35867,0.27404,-0.00079)
QUAG.add_eq('dbh','crown ht','cub',0.10736,0.28514,-0.00312,0.00001)
QUAG.add_eq('age','dbh','lin',1.55732,1.48691)
QUAG.add_eq('dbh','leaf area','quad',-2.44397,2.13252,0.06884)
QUAG.add_eq('dbh','tree ht','cub',2.00172,0.32851,-0.00331,0.00001)

QURU= Species('Quercus rubra')
QURU.add_eq('dbh','age','quad',0.64315,0.34124,0.00265)
QURU.add_eq('cdia','dbh','cub',0.69338,1.08675,0.26902,-0.0058)
QURU.add_eq('dbh','crown dia','loglogw1',-0.33383,2.2058,0.0115)
QURU.add_eq('dbh','crown ht','lin',1.80318,0.21653)
QURU.add_eq('age','dbh','cub',0.12978,1.83846,0.02001,-0.00036)
QURU.add_eq('dbh','leaf area','cub',1.5132,-1.40346,0.45025,-0.00179)
QURU.add_eq('dbh','tree ht','cub',4.32011,0.02074,0.00596,-0.00004)

TIAM= Species('Tilia americana')
TIAM.add_eq('dbh','age','cub',0.25539,0.67674,-0.0057,0.00009)
TIAM.add_eq('cdia','dbh','quad',0.14685,3.4123,0.09475)
TIAM.add_eq('dbh','crown dia','quad',0.00315,0.30261,-0.00136)
TIAM.add_eq('dbh','crown ht','lin',1.79973,0.20395)
TIAM.add_eq('age','dbh','cub',-0.57558,1.6899,0.00442,-0.00016)
TIAM.add_eq('dbh','leaf area','loglogw1',-0.53808,5.14273,0.16444)
TIAM.add_eq('dbh','tree ht','quad',2.16179,0.30424,-0.0008)

TICO= Species('Tilia cordata')
TICO.add_eq('dbh','age','cub',0.77726,0.4449,0.00919,-0.00006)
TICO.add_eq('cdia','dbh','quad',0.91801,2.66134,0.13245)
TICO.add_eq('dbh','crown dia','quad',-0.01436,0.29841,-0.0012)
TICO.add_eq('dbh','crown ht','cub',0.85835,0.37393,-0.00368,0.00002)
TICO.add_eq('age','dbh','cub',-0.74132,1.99578,-0.02803,0.00027)
TICO.add_eq('dbh','leaf area','cub',0.66057,-0.60429,0.39734,-0.00279)
TICO.add_eq('dbh','tree ht','quad',2.18773,0.32224,-0.0011)

ULAM= Species('Ulmus americana')
ULAM.add_eq('dbh','age','cub',0.90064,0.34281,0.00785,-0.00004)
ULAM.add_eq('cdia','dbh','cub',0.11354,1.14531,0.31746,-0.0076)
ULAM.add_eq('dbh','crown dia','cub',0.33815,0.43976,-0.00425,0.00002)
ULAM.add_eq('dbh','crown ht','cub',0.78744,0.55014,-0.00515,0.00002)
ULAM.add_eq('age','dbh','quad',-0.70738,1.81652,-0.00501)
ULAM.add_eq('dbh','leaf area','loglogw1',-0.52583,5.28335,0.13193)
ULAM.add_eq('dbh','tree ht','cub',1.86129,0.61188,-0.00496,0.00002)

WARO= Species('Washingtonia robusta')
WARO.add_eq('dbh','crown dia','quad',0.09678,0.411,-0.00834)
WARO.add_eq('dbh','crown ht','quad',0.115,0.433,-0.0094)
WARO.add_eq('dbh','leaf area','cub',0.116,0.09026,0.22,-0.00547)
WARO.add_eq('dbh','tree ht','quart',-0.163,0.783,0.00079,-0.00015)
