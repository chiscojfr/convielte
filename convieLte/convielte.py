import ply.lex as lex
import ply.yacc as yacc
import sys


#LEXER de los prefijos para conversiones metricas simples
tokens = [

    'FLOAT',
    'INT',
    'NAME',
    'FEMTO',
    'PICO',
    'NANO',
    'MICRO',
    'MILLI',
    'CENTI',
    'DECI',
    'DEKA',
    'HECTO',
    'KILO',
    'MEGA',
    'GIGA',
    'TERA',
    'EQUALS',
    'PLUS',
    'MINUS',
    'DIVIDE',
    'MULTIPLY',
    'LESSTHAN',
    'GREATERTHAN',
    'LPAR',
    'RPAR',
    'COMMA',
    'HELP'

]

t_PLUS = r'\+'
t_MINUS = r'\-'
t_MULTIPLY = r'\*'
t_DIVIDE = r'\/'
t_EQUALS= r'\='
t_ignore = r' '
t_LESSTHAN = r'\<'
t_GREATERTHAN = r'\>'
t_LPAR = r'\('
t_RPAR = r'\)'
t_COMMA = r'\,'


def t_HELP(t):
    r'help'
    t.type = 'HELP'
    return t

def t_NAME(t):
    r'[a-zA-Z_][a-zA-Z_0-9]*'
    t.type = 'NAME'
    return t


def t_FEMTO(t):
    r'\d+(?:\.\d+f | f)'
    t.value = float(t.value[:-1])*(10**-15)
    return t

def t_PICO(t):
    r'\d+(?:\.\d+p | p)'
    t.value = float(t.value[:-1])*(10**-12)
    return t

def t_NANO(t):
    r'\d+(?:\.\d+n | n)'
    t.value = float(t.value[:-1])*(10**-9)
    return t

def t_MICRO(t):
    r'\d+(?:\.\d+u | u)'
    t.value = float(t.value[:-1])*(10**-6)
    return t

def t_MILLI(t):
    r'\d+(?:\.\d+m | m)'
    t.value = float(t.value[:-1])*(10**-3)
    return t

def t_CENTI(t):
    r'\d+(?:\.\d+c | c)'
    t.value = float(t.value[:-1])*(10**-2)
    return t

def t_DEKA(t):
    r'\d+(?:\.\d+da | da)'
    t.value = float(t.value[:-2])*(10)
    return t

def t_DECI(t):
    r'\d+(?:\.\d+d | d)'
    t.value = float(t.value[:-1])*(10**-1)
    return t

def t_HECTO(t):
    r'\d+(?:\.\d+h | h)'
    t.value = float(t.value[:-1])*(10**2)
    return t

def t_KILO(t):
    r'\d+(?:\.\d+k | k)'
    t.value = float(t.value[:-1])*(10**3)
    return t

def t_MEGA(t):
    r'\d+(?:\.\d+M | M)'
    t.value = float(t.value[:-1])*(10**6)
    return t

def t_GIGA(t):
    r'\d+(?:\.\d+G | G)'
    t.value = float(t.value[:-1])*(10**9)
    return t

def t_TERA(t):

    r'\d+(?:\.\d+T | T)'
    t.value = float(t.value[:-1])*(10**12)
    return t

def t_FLOAT(t):
    r'\d+\.\d+'
    t.value = float(t.value)
    return t

def t_INT(t):
    r'\d+'
    t.value = int(t.value)
    return t

def t_error(t):
    print("Illegal characters!")
    t.lexer.skip(1)

lexer = lex.lex()



#PARSER de prefijos metricos simples

precedence = (

    ('left', 'PLUS', 'MINUS'),
    ('left', 'MULTIPLY', 'DIVIDE')

)


def p_conv(p):
    '''
    conv : expression
         | var_assign
         | empty
    '''
    print(run(p[1]))

def p_help(p):
    '''
    expression : HELP
    '''
    if p[1] == 'help':
        p[0]= "Welcome to ConvieLte, a programming language to simplify your conversion problems and help you calculate values." \
              "\n\nThe possible types of conversion are: \n" \
              "\nconvielte<prx>(NUMBER,PREFIX) --> Any number in base format conversion to" \
              " FEMTO(f), PICO (p), NANO (n), MICRO (u), MILLI (m), \n\t\t\t\t\t\t\t\tCENTI (c), DECI (d), DEKA (da), HECTO (h), KILO (k), MEGA (M), GIGA (G), TERA (T) " \
              "\nconvielte<L>(VALUE,TYPE) --> Any number in base L to oz, gal, m^3" \
              "\nconvielte<oz>(VALUE,TYPE) --> Any number in base oz to L, gal, m^3" \
              "\nconvielte<m>(VALUE,TYPE) --> Any number in base m to in, ft, yd" \
              "\nconvielte<in>(VALUE,TYPE) --> Any number in base in to m, cm, mm, ft, yd" \
              "\nconvielte<ft>(VALUE,TYPE) --> Any number in base ft to in, m, yd" \
              "\nconvielte<yd>(VALUE,TYPE) --> Any number in base yd to in, ft, m" \
              "\nconvielte<rad>(VALUE,TYPE) --> Any number in rad to deg" \
              "\nconvielte<deg>(VALUE,TYPE) --> Any number in deg to rad" \
              "\nconvielte<km>(VALUE,TYPE) --> Any number in km to mi" \
              "\nconvielte<mi>(VALUE,TYPE) --> Any number in mi to km" \
              "\nconvielte<F>(VALUE,TYPE) --> Any number in F to C, K" \
              "\nconvielte<C>(VALUE,TYPE) --> Any number in C to F, K" \
              "\nconvielte<K>(VALUE,TYPE) --> Any number in K to F, C" \
              "\nconvielte<g>(VALUE,TYPE) --> Any number in g to mg, kg, lb" \
              "\nconvielte<lb>(VALUE,TYPE) --> Any number in lb to g, kg, t" \
              "\nconvielte<t>(VALUE,TYPE) --> Any number in t to kg, lb" \
              "\nconvielte<s>(VALUE,TYPE) --> Any number in s to hr, min, day" \
              "\nconvielte<hr>(VALUE,TYPE) --> Any number in hr to s, min, day" \
              "\nconvielte<min>(VALUE,TYPE) --> Any number in min to s, hr, day" \
              "\nconvielte<day>(VALUE,TYPE) --> Any number in day to s, min, hr" \
              "\n\nThe available formulas are: \n" \
              "\nohmsLaw<V>(VALUE OF CURRENT,VALUE OF RESISTANCE) --> Calculates voltage using Ohm's Law (V = i * R)" \
              "\nohmsLaw<I>(VALUE OF VOLTAGE,VALUE OF RESISTANCE) --> Calculates current using Ohm's Law (i = V / R)" \
              "\nohmsLaw<V>(VALUE OF CURRENT,VALUE OF RESISTANCE) --> Calculates voltage using Ohm's Law (V = i * R)" \
              "\nohmsLaw<R>(VALUE OF VOLTAGE,VALUE OF CURRENT) --> Calculates resistance using Ohm's Law (R = V / i)" \
              "\nelectricPotencyVI<P>(VALUE OF VOLTAGE,VALUE OF CURRENT) --> Calculates power (P = V * i)" \
              "\nelectricPotencyVI<I>(VALUE OF POWER,VALUE OF VOLTAGE) --> Calculates current (i = P / V)" \
              "\nelectricPotencyVI<V>(VALUE OF POWER,VALUE OF CURRENT) --> Calculates voltage (V = P / i)" \
              "\nelectricPotencyRV<P>(VALUE OF VOLTAGE SQUARED,VALUE OF RESISTANCE) --> Calculates power (P = V^2 / R)" \
              "\nelectricPotencyRV<R>(VALUE OF VOLTAGE SQUARED,VALUE OF POWER) --> Calculates resistance (R = V^2 / P)" \
              "\nelectricPotencyRI<P>(VALUE OF CURRENT SQUARED,VALUE OF RESISTANCE) --> Calculates power (P = i^2 / R)" \
              "\nelectricPotencyRI<R>(VALUE OF POWER,VALUE OF CURRENT SQUARED) --> Calculates resistance (R = P / i^2)" \
              "\nelectricPotencyEt<P>(VALUE OF ENERGY,VALUE OF TIME) --> Calculates power (P = E / t)" \
              "\nelectricPotencyEt<E>(VALUE OF POWER,VALUE OF TIME) --> Calculates energy (E = P * t)" \
              "\nelectricPotencyEt<t>(VALUE OF ENERGY,VALUE OF POWER) --> Calculates time (t = E / P)" \
              "\nforceEquation<force>(VALUE OF ACCELERATION,VALUE OF MASS) --> Calculates force (F = a * m)" \
              "\nforceEquation<acceleration>(VALUE OF FORCE,VALUE OF MASS) --> Calculates acceleration (a = F / m)" \
              "\nforceEquation<mass>(VALUE OF FORCE,VALUE OF ACCELERATION) --> Calculates mass (m = F / a)" \
              "\npositionEquation<position>(VALUE OF INITIAL POSITION,VALUE OF VELOCITY*TIME) --> Calculates horizontal position (Ph = Pi + V*t)" \
              "\npositionEquation<time>(VALUE OF HORIZONTAL POSITION,VALUE OF INITIAL POSITION,VALUE OF VELOCITY) --> Calculates time (t = Ph-Pi / V)" \
              "\npositionEquation<velocity>(VALUE OF HORIZONTAL POSITION,VALUE OF INITIAL POSITION,VALUE OF TIME) --> Calculates velocity (V = Ph-Pi / t)" \
              "\nfreefallEquation<height>(VALUE OF VELOCITY,VALUE OF TIME,VALUE OF TIME SQUARED) --> Calculates height (h = V*t - 9.8*(t^2) / 2)" \
              "\nfreefallEquation<velocity>(VALUE OF HEIGHT,VALUE OF TIME SQUARED,VALUE OF TIME) --> Calculates height (V = h + 9.8*(t^2) / 2) / t)" \
              "\nworkEquation<W>(VALUE OF FORCE,VALUE OF DISPLACEMENT) --> Calculates work (W = F * x)" \
              "\nworkEquation<F>(VALUE OF WORK,VALUE OF DISPLACEMENT) --> Calculates force (F = W / x)" \
              "\nworkEquation<displacement>(VALUE OF WORK,VALUE OF FORCE) --> Calculates displacement (x = W / F)" \
              "\nkineticEquation<E>(VALUE OF MASS, VALUE OF VELOCITY SQUARED) --> Calculates the kinetic energy (Ec = (m * V^2) / 2)" \
              "\nkineticEquation<mass>(VALUE OF KINETIC ENERGY,VALUE OF VELOCITY SQUARED) --> Calculates the mass (m = 2 * (Ec / V^2))" \
              "\npotentialEquation<E>(VALUE OF MASS,VALUE OF HEIGHT) --> Calculates the potential energy (Ep = m * h * 9.8)" \
              "\npotentialEquation<mass>(VALUE OF POTENTIAL ENERGY,VALUE OF HEIGHT) --> Calculates the mass (m = Ep / (9.8 * h))" \
              "\npotentialEquation<height>(VALUE OF POTENTIAL ENERGY,VALUE OF MASS) --> Calculates the height (h = Ep / (9.8 * m))"

names = {}

def p_var_assign(p):
    '''
    var_assign : NAME EQUALS expression
    '''
    p[0]= ('=', p[1], p[3])
    names[p[1]] = p[3]


def p_expression(p):
    '''
    expression : expression MULTIPLY expression
                    | expression DIVIDE expression
                    | expression PLUS expression
                    | expression MINUS expression
    '''

def p_expression_conv(p):
    '''
    expression : FLOAT
                | INT
                | FEMTO
                | PICO
                | NANO
                | MICRO
                | MILLI
                | CENTI
                | DECI
                | DEKA
                | HECTO
                | KILO
                | MEGA
                | GIGA
                | TERA
    '''
    p[0]=p[1]


def p_expression_var(p):
    '''
    expression : NAME
    '''
    p[0]=p[1]
    try:
        p[0] = names[p[1]]
    except LookupError:
        p[0] = p[1]


def p_expression_functions(p):
    '''
    expression : expression LESSTHAN expression GREATERTHAN LPAR expression COMMA expression RPAR
               | expression LESSTHAN expression GREATERTHAN LPAR expression COMMA expression COMMA expression RPAR
    '''

    if p[1] == 'convielte':

         # Liquid Measure Converter
        if p[3] == 'L':
            if p[8] == 'oz':
                p[0] = str(float(p[6]) * 33.814) + ' oz'  # L to oz conversion
            elif p[8] == 'gal':
                p[0] = str(float(p[6]) * 0.264) + ' gal'  # L to gallons conversion
            elif p[8] == 'm3':
                p[0] = str(float(p[6]) * 0.001) + ' m^3'  # L to m^3 conversion

        if p[3] == 'oz':
            if p[8] == 'L':
                p[0] = str(float(p[6]) * 0.0295735) + ' L'  # oz to L conversion
            elif p[8] == 'gal':
                p[0] = str(float(p[6]) * 0.0078125) + ' gal'  # oz to gallons conversion
            elif p[8] == 'ml':
                p[0] = str(float(p[6]) * 29.5735) + ' ml'  # oz to ml conversion

        # Measure Converter
        if p[3] == 'm':
            if p[8] == 'in':
                p[0] = str(float(p[6]) * 39.37) + ' in'  # m to in conversion
            if p[8] == 'ft':
                p[0] = str(float(p[6]) * 3.28) + ' ft'  # m to ft conversion
            if p[8] == 'yd':
                p[0] = str(float(p[6]) * 1.094) + ' yd'  # m to yd conversion

        if p[3] == 'in':
            if p[8] == 'm':
                p[0] = str(float(p[6]) * 0.0254) + ' m'  # in to m conversion
            if p[8] == 'cm':
                p[0] = str(float(p[6]) * 2.54) + ' cm'  # in to cm conversion
            if p[8] == 'mm':
                p[0] = str(float(p[6]) * 25.4) + ' mm'  # in to cm conversion
            if p[8] == 'ft':
                p[0] = str(float(p[6]) * 0.083) + ' ft'  # in to ft conversion
            if p[8] == 'yd':
                p[0] = str(float(p[6]) * 0.028) + ' yd'  # in to yd conversion

        if p[3] == 'ft':
            if p[8] == 'in':
                p[0] = str(float(p[6]) * 12) + ' in'  # ft to in conversion
            if p[8] == 'm':
                p[0] = str(float(p[6]) * 0.3048) + ' m'  # ft to m conversion
            if p[8] == 'yd':
                p[0] = str(float(p[6]) * 0.33) + ' yd'  # ft to yd conversion

        if p[3] == 'yd':
            if p[8] == 'in':
                p[0] = str(float(p[6]) * 36) + ' in'  # yd to in conversion
            if p[8] == 'ft':
                p[0] = str(float(p[6]) * 3) + ' ft'  # yd to ft conversion
            if p[8] == 'm':
                p[0] = str(float(p[6]) * 0.9144) + ' m'  # yd to m conversion

                
        #Plane Angle Converter
        if p[3] == 'rad':
            if p[8] == 'deg':
                p[0] = str(float(p[6]) * 57.296) + ' deg'  # rad to deg conversion
        if p[3] == 'deg':
            if p[8] == 'rad':
                p[0] = str(float(p[6]) * 0.017) + ' rad'  # rad to deg conversion
       
        # Distance Converter
        if p[3] == 'km':
            if p[8] == 'mi':
                p[0] = str(float(p[6]) * 0.6214) + ' mi'  # km to mi conversion

        if p[3] == 'mi':
            if p[8] == 'km':
                p[0] = str(float(p[6]) * 1.60934) + ' km'  # mi to km conversion

        # Temperature Converter
        if p[3] == 'F':
            if p[8] == 'C':
                p[0] = str(float(p[6] - 32) * 5 / 9) + ' C'  # F to C conversion
            elif p[8] == 'K':
                p[0] = str(float(p[6] + 459.67) * 5 / 9) + 'K' # F to K conversion

        if p[3] == 'C':
            if p[8] == 'F':
                p[0] = str(float(p[6] * 9 / 5) + 32) + ' F'  # C to F conversion
            elif p[8] == 'K':
                p[0] = str(float(p[6] + 273.15)) + 'K' # C to K conversion

        if p[3] == 'K':
            if p[8] == 'F':
                p[0] = str(float(p[6] * (9 / 5)) - 459.67) + 'F' # K to F conversion
            elif p[8] == 'C':
                p[0] = str(float(p[6] - 273.15)) + 'C' # K to C conversion

        # Weight Converter
        if p[3] == 'g':
            if p[8] == 'mg':
                p[0] = str(float(p[6] * 1000)) + ' mg' # g to mg conversion
            elif p[8] == 'kg':
                p[0] = str(float(p[6] / 1000)) + ' kg' # g to kg conversion
            elif p[8] == 'lb':
                p[0] = str(float(p[6] / 453.59237)) + ' lb' # g to lb conversion

        if p[3] == "lb":
            if p[8] == 'g':
                p[0] = str(float(p[6] * 453.59237)) + ' g' # lb to g conversion
            elif p[8] == 'kg':
                p[0] = str(float(p[6] * 0.45359237)) + ' kg' # lb to kg conversion
            elif p[8] == 't':
                p[0] = str(float(p[6] * 0.00045359237)) + ' t' # lb to t conversion

        if p[3] == "t":
            if p[8] == 'kg':
                p[0] = str(float(p[6] * 1000)) + ' kg' # t to kg conversion
            elif p[8] == "lb":
                p[0] = str(float(p[6] / 0.00045359237)) + ' lb' # t to lb conversion
                
                
        # Time Converter
        if p[3] == 's':
            if p[8] == 'hr':
                p[0] = str(float(p[6] * 0.00028)) + ' hr'  # s to hr conversion
            elif p[8] == 'min':
                p[0] = str(float(p[6] * 60)) + ' min'  # s to min conversion
            elif p[8] == 'day':
                p[0] = str(float(p[6] * 0.042 )) + ' day'  # s to day conversion

        if p[3] == 'hr':
            if p[8] == 's':
                p[0] = str(float(p[6] * 3600)) + ' s'  # hr to s conversion
            elif p[8] == 'min':
                p[0] = str(float(p[6] * 0.017)) + ' min'  # hr to min conversion
            elif p[8] == 'day':
                p[0] = str(float(p[6] * 0.0000116)) + ' day'  # hr to day conversion

        if p[3] == 'min':
            if p[8] == 'hr':
                p[0] = str(float(p[6] * 0.017)) + ' hr'  # min to hr conversion
            elif p[8] == 's':
                p[0] = str(float(p[6] * 60)) + ' s'  # min to s conversion
            elif p[8] == 'day':
                p[0] = str(float(p[6] * 0.00069)) + ' day'  # min to day conversion

        if p[3] == 'day':
            if p[8] == 'hr':
                p[0] = str(float(p[6] * 24)) + ' hr'  # day to hr conversion
            elif p[8] == 's':
                p[0] = str(float(p[6] * 86400)) + ' s'  # day to s conversion
            elif p[8] == 'min':
                p[0] = str(float(p[6] * 1440)) + ' min'  # day to min conversion

        # Prefix Converter
        if p[3] == 'prx':
            if p[8] == 'E':
                p[0] = str(float(p[6]) / 1000000000000000000) + ' E'  # to Exa conversion
            if p[8] == 'P':
                p[0] = str(float(p[6]) / 1000000000000000) + ' P'  # to Peta conversion
            if p[8] == 'T':
                p[0] = str(float(p[6]) / 1000000000000) + ' T'  # to Tera conversion
            if p[8] == 'G':
                p[0] = str(float(p[6]) / 1000000000) + ' G'  # to Giga conversion
            if p[8] == 'M':
                p[0] = str(float(p[6]) / 1000000) + ' M'  # to Mega conversion
            if p[8] == 'k':
                p[0] = str(float(p[6]) / 1000) + ' k'  # to Kilo conversion
            if p[8] == 'h':
                p[0] = str(float(p[6]) / 100) + ' h'  # to hecto conversion
            if p[8] == 'da':
                p[0] = str(float(p[6]) / 10) + ' da'  # to deca conversion

            if p[8] == 'd':
                p[0] = str(float(p[6]) / 0.1) + ' d'  # to deci conversion
            if p[8] == 'c':
                p[0] = str(float(p[6]) / 0.01) + ' c'  # to centi conversion
            if p[8] == 'm':
                p[0] = str(float(p[6]) / 0.001) + ' m'  # to mili conversion
            if p[8] == 'u':
                p[0] = str(float(p[6]) / 0.000001) + ' u'  # to micro conversion
            if p[8] == 'n':
                p[0] = str(float(p[6]) / 0.000000001) + ' n'  # to nano conversion
            if p[8] == 'p':
                p[0] = str(float(p[6]) / 0.000000000001) + ' p'  # to pico conversion
            if p[8] == 'f':
                p[0] = str(float(p[6]) / 0.000000000000001) + ' f'  # to femto conversion
            if p[8] == 'a':
                p[0] = str(float(p[6]) / 0.000000000000000001) + ' a'  # to atto conversion

    # Electric Equations
    elif p[1] == 'ohmsLaw':
        if p[3] == 'V':
            p[0] = str(float(p[6]) * float(p[8])) + ' V'  # V = i * R
        if p[3] == 'I':
            p[0] = str(float(p[6]) / float(p[8])) + ' A'  # i = V / R
        if p[3] == 'R':
            p[0] = str(float(p[6]) / float(p[8])) + ' Ohms'  # R = V / i

    elif p[1] == 'electricPotencyVI':
        if p[3] == 'P':
            p[0] = str(float(p[6]) * float(p[8])) + ' W' # P = V*I
        if p[3] == 'I':
            p[0] = str(float(p[6]) / float(p[8])) + ' A'  # i = P/V
        if p[3] == 'V':
            p[0] = str(float(p[6]) / float(p[8])) + ' V'  # i = P/I

    elif p[1] == 'electricPotencyRV':
        if p[3] == 'P':
            p[0] = str(float(p[6]) * float(p[6]) / float(p[8])) + ' W' # P = (V^2)/R
        if p[3] == 'R':
            p[0] = str(float(p[6]) * float(p[6]) / float(p[8])) + ' Ohms' # R = (V^2)/P

    elif p[1] == 'electricPotencyRI':
        if p[3] == 'P':
            p[0] = str(float(p[6]) * float(p[6]) * float(p[8])) + ' W'  # P = (I^2)*R
        if p[3] == 'R':
            p[0] = str(float(p[8])/ (float(p[6]) * float(p[6]))) + ' Ohms' # R = P/(I^2)

    elif p[1] == 'electricPotencyEt':
        if p[3] == 'P':
            p[0] = str(float(p[6]) / float(p[8])) + ' W'  # P = E/t
        if p[3] == 'E':
            p[0] = str(float(p[6]) * float(p[8])) + ' J'  # E = P*t
        if p[3] == 't':
            p[0] = str(float(p[6]) / float(p[8])) + ' s'  # t = E/P

    # Physics Equations
    elif p[1] == 'forceEquation':
        if p[3] == 'force':
            p[0] = str(float(p[6]) * float(p[8])) + 'N'  # Force = Acceleration * Mass
        if p[3] == 'acceleration':
            p[0] = str(float(p[6]) / float(p[8])) + 'm / s^2'  # Acceleration = Force / Mass
        if p[3] == 'mass':
            p[0] = str(float(p[6]) / float(p[8])) + 'kg'  # Mass = Force / Acceleration

    elif p[1] == 'positionEquation': # X = Xo + Vo*t
        if p[3] == 'position':
            p[0] = str(float(p[6]) + float(p[8]) * float(p[10])) + ' m' # posicion horizontal = posicion inicial + Velocidad*tiempo
        if p[3] == 'time':
            p[0] = str((float(p[6]) - float(p[8])) / float(p[10])) + ' s' # tiempo = (posicion horizontal - posicion inicial)/velocidad
        if p[3] == 'velocity':
            p[0] = str((float(p[6]) - float(p[8])) / float(p[10])) + ' m/s' # velocidad = posicion horizontal - posicion inicial)/tiempo

    elif p[1] == 'freefallEquation': # Y = Vo*t(1) - 1/2*g. t^2
        if p[3] == 'height':
            p[0] = str(float(p[6]) * float(p[8]) - 9.8*float(p[8])*float(p[8])/2) + ' m' # Altura = velocidad*t - 9.8*(t^2)/2
        if p[3] == 'velocity':
            p[0] = str((float(p[6]) + 9.8 * float(p[8])*float(p[8]) / 2)/float(p[8])) + ' m/s' # velocidad = (altura + 9.8*(t^2)/2)/t

    elif p[1] == 'workEquation':
        if p[3] == 'W':
            p[0] = str(float(p[6]) * float(p[8])) + ' J'  # W = F*X
        if p[3] == 'F':
            p[0] = str(float(p[6]) / float(p[8])) + ' N' # F = W/X
        if p[3] == 'displacement':
            p[0] = str(float(p[6]) / float(p[8])) + ' m' # X = W/F

    elif p[1] == 'kineticEquation': # Ec = 1/2 m V^2
        if p[3] == 'E':
            p[0] = str(float(p[6]) * float(p[8]) / 2) + ' J' # Ec = (m*V^2)/2
        if p[3] == 'mass':
            p[0] = str(2 * float(p[6]) / (float(p[8]) * float(p[8]))) + ' kg' # m = 2*Ec/(V^2)

    elif p[1] == 'potentialEquation': # Ep = m*g*h
        if p[3] == 'E':
            p[0] = str(float(p[6]) * float(p[8]) * 9.8) + ' J' # Ep = m*g*h
        if p[3] == 'mass':
            p[0] = str(float(p[6]) / (float(p[8]) * 9.8)) + ' kg' # m = Ep/g*h
        if p[3] == 'height':
            p[0] = str(float(p[6]) / (float(p[8]) * 9.8)) + ' m' # h = Ep/g*m


def p_empty(p):
    '''
    empty :
    '''
    p[0] = None

def p_error(p):
    print("Syntax Error!")

parser = yacc.yacc()
env = {}

def run(p):
    global env
    if type(p) == tuple:
        if p[0]=='+':
            return run(p[1]) + run(p[2])
        elif p[0]=='-':
            return run(p[1]) - run(p[2])
        elif p[0]=='/':
            return run(p[1]) / run(p[2])
        elif p[0]=='*':
            return run(p[1]) * run(p[2])
        elif p[0]=='=':
            env[p[1]] = run(p[2])
            return env
        elif p[0]=='var':
            if p[1] not in env:
                return "Undefined variable found!"
            return env[p[1]]
    else:
        return p

while True:
    try:
        s = input('>> ')
    except EOFError:
        break
    parser.parse(s)