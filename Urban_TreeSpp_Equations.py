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
ACRU.dbh_to_ht.calc(6.5, inverse= False)        

# If you wanted to calculate dbh from a tree with known ht of 10m:
ACRU.dbh_to_ht.calc(10, inverse=True) 
'''        
        
# The Tree Equations for the Pacific Northwest Region
# ...
