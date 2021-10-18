import numpy as np 

class PointInfinity(object):
    """ Wrapper to indicate the point at infinity serving as the identity element. """
    def __init__(self):
        pass

class ArithmeticAdapter(object):
    """ Adapter class to implement arithmetic over the required Galois field F_p. 
        NOTE: ArithmeticAdapater is an interal class hidden in ECPoint from the user. The class attributes p and n will be automatically set by the EllipticCurve() constructor if any char is given at all hence cannot be run with p and n equal None.
    """
    p = None # Charateristic 
    n = None # Dimension
    def __init__(self, x: int):
        self.coor = x % ECPoint.p

    def __add__(self, y):
        return ArithmeticAdapter(
            (self.coor + y.coor) % ECPoint.p
            )
    
    def __mul__(self, y):
        return ArithmeticAdapterr(
            (self.coor * y.coor) % ECPoint.p
            )
    
    def __sub__(self, y):
        return ArithmeticAdapter(
            (self.coor - y.coor) % self.__class__.p
            )
    
    def __div__(self, y):
        # Since a finite field F_p^n forms a multiplicative group
        # we know that the inverse of a is a^(p^n-2).
        size = np.power(self.__class__.p, self.__n__)
        y_inv = np.power(y.coor, size-2)
        return ArithmeticAdapter(
            (self.coor * y_inv) % self.__class__.p
        )
    
    def __eq__(self, y):
        if self.coor == y.coor:
            return True 
        return False

class ECPoint:
    """ Wrapper to encapsulate point coordinates and arithmetic over Galois Field."""
    p = 0 
    n = 1
    def __init__(self, x,y):
        if self.__class__.p:
            self.x = ArithmeticAdapter(x)
            self.y = ArithmeticAdapter(y)
        else:
            self.x = x 
            self.y = y 
        
    def __eq__(self, Q):
        if self.x == Q.x and self.y == Q.y:
            return True 
        return False
    
    def __str__(self):
        if self.p:
            print(f"Point: ({self.x.coor},{self.y.coor})")


class EllipticCurve:
    """ Implements EC structure.
        EC is defined by equation: y^2 = x^3 + Ax + B 
    """
    def __init__(self,A,B, p=0,n=1):
        if p:
            ECPoint.p = p 
            ECPoint.n = n 
        self.A = A 
        self.B = B 
    
    def addition(self,P,Q):
        # P + O = O + P = P 
        if isinstance(Q, PointInfinity):
            return ECPoint(P.x, P.y)
        
        if isinstance(P, PointInfinity):
            return ECPoint(Q.x,Q.y)

        if P == Q:
            lambd = (3*np.power(self.x,2) + self.A) / 2* self.y

        # P + (-P) = O 
        if P.x == Q.x and P.y == -Q.y:
            return PointInfinity()
        
        lambd = (Q.y - P.y) / (Q.x - P.x)
        x_r = np.power(lambd, 2) - P.x - Q.x 
        y_r = lambd*(P.x - x_r) - P.y
        return ECPoint(x_r,-y_r)
    
    def multiply(P, n):
        # Recursive multiplication algorithm requiring log(n)
        # iterations of doublning and addition to compute multiplication.
        if n == 0:
            return 0 
        elif n == 1:
            return P 
        elif n % 2 == 1:
            return self.addition(P, self.multiply(P, n-1))
        else:
            return self.multiply(self.addition(P,P), n/2)


""" 
Required interface:
curve = EllipticCurve(A,B)
curve1 = EllipticCurve(A1,B1,p)
P = ECPoint(X,Y)
Q = ECPoint(X,Y)
curve.addition(P,Q)
curve.addition(P,Q) # Working over F_p
"""


