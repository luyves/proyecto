### Tipos duales ###
global pylab, pl
import pylab as pl

class Dual(object):
    def __init__(self,fun,der=0):
        self.fun = fun
        self.der = der
    def __str__(self):
        return "Dual(%s,%s)" % (str(self.fun),str(self.der))    
    def __real__(self, a):
        return Dual(a,0)
    def __add__(self, b):
        if type(self) < type(b):
            b = Dual(b)
        return Dual(self.fun+b.fun,self.der+b.der)
    def __radd__(self, b):
        return self.__add__(b)
    def __sub__(self, b):
        if type(self) < type(b):
            b = Dual(b)
        return Dual(self.fun-b.fun,self.der-b.der)
    def __rsub__(self, b):
        return self.__sub__(b)
    def __mul__(self, b):
        if type(self) < type(b):
            return Dual(b*self.fun,b*self.der)
        return Dual(self.fun*b.fun, (self.fun*b.der) + (b.fun*self.der))
    def __rmul__(self, b):
        return self.__mul__(b)
    def __div__(self, b):
        if type(self) < type(b) and b != 0:
            return Dual(self.fun/float(b), self.der/float(b))
        return Dual(self.fun/b.fun, (self.der - (self.fun/float(b.fun))*b.der)/float(b.fun))
    def __rdiv__(self,b):
        return self.__div__(b)
    def __pow__(self,b):
        return Dual(self.fun**b, b*(self.fun**(b-1))*self.der)

### Funciones para tipos duales
def xdual(x0):
    return Dual(x0,1)
def fact(x):
    if x == 0 or x == 1:
        prod = 1
    else:
        prod = x*fact(x-1)
    return prod


### Series de Taylor ###

class Taylor(object):
    def __init__(self, pol):
        self.pol = pol
    def __str__(self):
        return "Taylor(%s)" % (str(self.pol))
    def __add__(self, b):
        if type(b) != (list and Taylor):
            b = Taylor([b])
        nuevoPol = [i+j for i,j in zip(promocionTaylor(self,b).pol,promocionTaylor(b,self).pol)]
        return Taylor(nuevoPol)
    def __radd__(self, b):
        return self.__add__(b)
    def __sub__(self, b):
        if type(b) != (list and Taylor):
            b = Taylor([b])
        nuevoPol = [i-j for i,j in zip(promocionTaylor(self,b).pol,promocionTaylor(b,self).pol)]
        return Taylor(nuevoPol)
    def __rsub__(self, b):
        return self.__sub__(b)
    def __mul__(self, b):
        if type(b) != (list and Taylor):
            b = Taylor([b])
        n = gradoMaxTaylor(self) + gradoMaxTaylor(b) - 1
        r = Taylor([0]*n)
        A = promocionTaylor(self,r)
        B = promocionTaylor(b,r)
        for k in range(n):
            suma = 0
            for j in range(k+1):
                suma += A.pol[j]*B.pol[k-j]
            r.pol[k] = suma
        return r
    def __rmul__(self, b):
        return self.__mul__(b)
    def __div__(self, b):
        if type(b) != (list and Taylor):
            b = Taylor([b])
        A = promocionTaylor(self,b)
        B = promocionTaylor(b,self)
        n = gradoMaxTaylor(A)
        r = Taylor([0]*n)
        s = 0 # indice desde donde empezamos
        while B.pol[s] == 0: # checamos si el primer termino no es nulo
            s += 1;
        r.pol[0] = float(A.pol[s])/B.pol[s];
        for k in range(s+1,n):
            suma = 0
            for j in range(k):
                suma += r.pol[j]*B.pol[k-j]
            r.pol[k-s] = float(A.pol[k]-suma)/B.pol[s]    
        return r
    def __rdiv__(self,b):
        return self.__div__(b)
    def __pow__(self,n):
        res = Taylor([1.])
        for k in range(n):
            res = res*self
        return res

### Funciones para polinomios de Taylor

def gradoMaxTaylor(a,b=Taylor([0])):
    return max(len(a.pol), len(b.pol))
def promocionTaylor(a,b):
    n = gradoMaxTaylor(a,b)-gradoMaxTaylor(a)
    nuevoPol = a.pol[:]
    nuevoPol.extend([0]*n)
    return Taylor(nuevoPol)
def autoPromocion(a,n):
    return promocionTaylor(a,Taylor([0]*n))

### Metodos funcionales compartidos ###

def sqrt(self):
    if type(self) == Dual:
        return Dual(pl.sqrt(self.fun), (0.5)*(self.der/pl.sqrt(self.fun)))
    elif type(self) == Taylor:
        return self**(0.5)
    else:
        return pl.sqrt(self)
def exp(self, nmax=5):
    if type(self) == Dual:
        return Dual(pl.exp(self.fun), self.der*pl.exp(self.fun))
    elif type(self) == Taylor:
        self = autoPromocion(self,max(gradoMaxTaylor(self),nmax))
        n = gradoMaxTaylor(self)
        exp_t = Taylor([0]*n)
        exp_t.pol[0] = exp(self.pol[0])
        for k in range(1,n):
            suma = 0
            for j in range(k):
                suma += (k-j)*self.pol[k-j]*exp_t.pol[j]
            exp_t.pol[k] = suma*(1./k)
        return exp_t
    else:
        return pl.exp(self)
def log(self, nmax=5):
    if type(self) == Dual:
        return Dual(pl.log(self.fun), self.der/self.fun)
    elif type(self) == Taylor:
        self = autoPromocion(self,max(gradoMaxTaylor(self),nmax))
        n = gradoMaxTaylor(self)
        L = Taylor([0]*n)
        s = 0 # indice no nulo desde donde empezamos
        while self.pol[s] == 0:
            s += 1
        L.pol[s] = log(self.pol[s])
        for k in range(s+1,n):
            suma = 0
            for j in range(s+1,k):
                suma += (j)*L.pol[j]*self.pol[k-j]
            L.pol[k] = (1./self.pol[s])*(self.pol[k]-suma/(k-s))
        return L
    else:
        return pl.log(self)
def sin(self,nmax=5):
    if type(self) == Dual:
        return Dual(pl.sin(self.fun), self.der*pl.cos(self.fun))
    elif type(self) == Taylor:
        self = autoPromocion(self,max(gradoMaxTaylor(self),nmax))
        n = gradoMaxTaylor(self)
        S = Taylor([0]*n)
        for k in range(n):
            S += ((-1.)**k)*(self**(2*k+1))/fact(2*k+1)
        return S
    else:
        return pl.sin(self)
def cos(self,nmax=5):
    if type(self) == Dual:
        return Dual(pl.cos(self.fun), -self.der*pl.sin(self.fun)) 
    elif type(self) == Taylor:
        self = autoPromocion(self,max(gradoMaxTaylor(self),nmax))
        n = gradoMaxTaylor(self)
        C = Taylor([0]*n)
        for k in range(n):
            C += ((-1.)**k)*(self**(2*k))/fact(2*k)
        return C
    else:
        return pl.cos(self)
def tan(self):
    if type(self) == Dual:
        return Dual(pl.tan(self.fun), self.der/(pl.cos(self.fun)**2))
    elif type(self) == Taylor:
        return sin(self)/cos(self)
    else:
        return pl.tan(self)
def cot(self):
    if type(self) == Dual:
        return Dual(1./pl.tan(self.fun), -1.*self.der/(pl.sin(self.fun)**2))
    elif type(self) == Taylor:
        return cos(self)/sin(self)
    else:
        return 1./pl.tan(self)
def sec(self):
    if type(self) == Dual:
        return Dual(1./pl.cos(self.fun), self.der*pl.tan(self.fun)/pl.cos(self.fun))
    elif type(self) == Taylor:
        return 1./cos(self)
    else:
        return 1./pl.cos(self)
def csc(self):
    if type(self) == Dual:
        return Dual(1./pl.sin(self.fun), -1.*self.der/(pl.tan(self.fun)*pl.sin(self.fun)))
    elif type(self) == Taylor:
        return 1./sin(self)
    else:
        return 1./pl.sin(self)

### Integrador para Duffing ###

def paso(a, eps= 1e-8, p=21):
    p = min(p, gradoMaxTaylor(a))-1
    while p != 0 and a.pol[p] == 0:
        p = p-1
        print p
    h = 0.5*(eps/abs(a.pol[p]))**(1./p)
    return h
def evHorner(a, t0):
    n = gradoMaxTaylor(a)
    x1 = a.pol[-1]
    for k in range(1,n):
        x1 = a.pol[n-k-1] + x1*t0
    return x1
def funcDuffing(x,y,t,omega=1.):
    return 0.3*cos(omega*t)+x-x**3-0.2*y
def coefsDuffing(x0,y0,t0,F,omega=1,n=21):
    x = [x0] + [0.]*(n-1)
    y = [y0] + [0.]*(n-1)
    for k in range(n-1):
        bx = Taylor(x)
        by = Taylor(y)
        f = F(bx,by,t0,omega)
        y[k+1] = f.pol[k]/(k+1)
        x[k+1] = y[k]/(k+1)
    return Taylor(x), Taylor(y)
def intDuffing(x0,y0,t0,tf,f,omega=1,eps=1e-20):
    t = [t0]
    x = [x0]
    y = [y0]
    while t0<tf and x0<1000:
        X,Y = coefsDuffing(x0,y0,t0,f,omega,n=21)
        hx,hy = paso(X,eps),paso(Y,eps)
        h = min(hx,hy)
        t0 += h
        x0,y0 = evHorner(X,h),evHorner(Y,h)
        t.append(t0)
        x.append(x0)
        y.append(y0)
    return t,x,y
