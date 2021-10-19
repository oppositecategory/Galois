import numpy as np
import matplotlib.pyplot as plt 

class PointInfinity(object):
    """ Wrapper to indicate the point at infinity serving as the identity element. """
    def __init__(self):
        pass

class ECPoint:
    """ Wrapper to encapsulate point coordinates."""
    def __init__(self, x,y):
        self.x = x 
        self.y = y 

    def __eq__(self, Q):
        if self.x == Q.x and self.y == Q.y:
            return True 
        return False
    
    def __str__(self):
        return f"Point: ({self.x},{self.y})"
        

class EllipticCurve:
    """ Implements EC structure.
        EC is defined by equation: y^2 = x^3 + Ax + B 
        args:
            - p = The field characteristic the Elliptic Curve is defined over. i.e E(F_p)
                  If p=None/0 then we work over the real numbers.
            - n = The dimension of the field. 
            - A,B = Parameters of the Elliptic Curve.
    """
    def __init__(self,A,B, p=0,n=1):
        self.p = p 
        self.n = n 
        if self.p:
            self.size = pow(self.p,self.n)
            self.A = A % p 
            self.B = B % p 
        else:
            self.size = None
            self.A = A 
            self.B = B 
    
    def _add_mod_p(self, P: ECPoint, Q : ECPoint):
        E = lambda x,y: pow(y,2) - pow(x,3) - self.A*x - self.B
        if E(P.x,P.y) !=0 or E(Q.x,Q.y) != 0:
            raise RuntimeError("Points must lie on the curve.")

        if isinstance(Q, PointInfinity):
            return ECPoint(P.x, P.y)
        
        if isinstance(P, PointInfinity):
            return ECPoint(Q.x,Q.y)

        if P == Q:
            # To invert we use the simple observation that a^m=1
            # where m=p^n in finite fields. 
            numerator = (3*pow(P.x,2) + self.A) % self.p 
            inv = pow(2*P.y,self.size-2,self.p)
            lambd = (numerator * inv) % self.p 
        else:
            inv = pow(Q.x-P.x,self.size-2, self.p)
            lambd = ((Q.y - P.y) * inv) % self.p

        # P + (-P) = O 
        if P.x == Q.x and P.y == -Q.y:
            return PointInfinity()
        
        x_r = (np.power(lambd,2) - P.x - Q.x) % self.p
        y_r = (lambd*(P.x - x_r) - P.y) % self.p 
        return ECPoint(x_r,y_r)

    
    def _add(self,P:ECPoint,Q:ECPoint):
        E = lambda x,y: pow(y,2) - pow(x,3) - self.A*x - self.B
        if E(P.x,P.y) !=0 or E(Q.x,Q.y) != 0:
            raise RuntimeError("Points must lie on the curve.")
        # P + O = O + P = P 
        if isinstance(Q, PointInfinity):
            return ECPoint(P.x, P.y)
        
        if isinstance(P, PointInfinity):
            return ECPoint(Q.x,Q.y)

        if P == Q:
            lambd = (3*(P.x*P.x) + self.A) / 2 * P.y
        else:
            lambd = (Q.y - P.y) / (Q.x - P.x)

        # P + (-P) = O 
        if P.x == Q.x and P.y == -Q.y:
            return PointInfinity()
        
        x_r = (lambd*lambd) - P.x - Q.x 
        y_r = lambd*(P.x - x_r) - P.y
        return ECPoint(x_r,-y_r) 
    
    def addition(self, P : ECPoint, Q : ECPoint):
        if self.p: 
            return self._add_mod_p(P,Q)
        else:
            return self._add(P, Q)
        
    
    def multiply(self, P, n):
        """Recursive multiplication algorithm requiring log(n)
           iterations of doublning and addition to compute multiplication.   
        """
        if n == 0:
            return 0 
        elif n == 1:
            return P 
        elif n % 2 == 1:
            return self.addition(P, self.multiply(P, n-1))
        else:
            return self.multiply(self.addition(P,P), n/2)

    def _plot_elliptic_curve_real(self):
        # Source: https://stackoverflow.com/questions/19756043/python-matplotlib-elliptic-curves
        y, x = np.ogrid[-5:5:100j, -5:5:100j]
        E = pow(y,2) - self.A * x - self.B
        plt.contour(x.ravel(),y.ravel(), E, [0])
        plt.grid()
        plt.show()
    
    def _plot_elliptic_curve_finite_field(self):
        x = y = np.array([i for i in range(self.size)])
        E = (pow(y,2) - self.A * x - self.B) % self.p
        points_modp = np.where(E==0)[0]
        plt.scatter(x[points_modp],y[points_modp])
        plt.grid()
        plt.show()

    def plot_elliptic_curve(self):
        if self.p:
            self._plot_elliptic_curve_finite_field()
        else:
            self._plot_elliptic_curve_real()

curve = EllipticCurve(A=-4,B=20)



""" 
Required interface:
curve = EllipticCurve(A,B)
curve1 = EllipticCurve(A1,B1,p)
P = ECPoint(X,Y)
Q = ECPoint(X,Y)
curve1.addition(P,Q) # Working over F_p
curve.addition(P,Q) 
"""


