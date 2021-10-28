import numpy as np 
from abc import ABC, abstractmethod

from ec import * 


class M_221(EllipticCurve):
    """ Elliptic curve M-221.
        y^2 = x^3 + 117050x^2+x
    """
    def __init__(self):
        p = np.power(2,221)-3
        super().__init__(self,A=117050,B=0,p=p)
    
    def generate_public_key(self, d):
        #return self.multiply(self.G,d)
        pass 
    

class Curve25519(EllipticCurve):
    """ Elliptic Curve 25519.
        y^2 = x^3+486662x^2+x
    """
    def __init__(self):
        p = np.power(2,255) - 19
        super().__init__(self,A=486662,B=0,p=p)
        self.G = ECPoint(
            9, 14781619447589544791020593568409986887264606134616475288964881837755586237401
        )
        self.n = np.pow(2,252) + 2774231777737235353585193779088364849
    
    def generate_public_key(self, d):
        return EllipticCurve.multiply(self.G, d)


class ECDH:
    """ Implements a party in the Elliptic Curve Diffie Hellman Protocol.
    """ 
    def __init__(self, 
                curve : EllipticCurve,
                private_key : ECPoint):
        # Internal parameters of the party.
        self.curve = curve 
        self.private_key = private_key
        self.public_key = self.curve.multiply(self.curve.G, 
                                              self.private_key)

        # Public key of the other party in the communication.
        self.party_pub_key = None 

        self.shared_secret = None 
    
    def set_public_key(self, key):
        # Assoicate the client with the given key for further communication.
        self.party_pub_key = key 

    def _compute_shared_secret(self):
        assert self.party_pub_key 
        self.shared_secret = self.curve.generate_public_key(self.private_key)

    def validate_communication(self, secret):
        if not self.shared_secret:
            self._compute_shared_secret()
        
        if secret == self.shared_secret:
            return True 
        return False 
        
        

        
     


