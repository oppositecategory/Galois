import numpy as np 
from abc import ABC, abstractmethod

from ec import * 

from Crypto.Cipher import AES


class M_221(EllipticCurve):
    """ Elliptic curve M-221.
        y^2 = x^3 + 117050x^2+x
    """
    def __init__(self):
        p = np.power(2,221)-3
        super().__init__(self,A=117050,B=0,p=p)

class Curve25519:
    """ Elliptic Curve 25519.
        y^2 = x^3+486662x^2+x
    """
    def __init__(self):
        self.c = EllipticCurve(A=486662,B=0,p=np.power(2,255) - 19)

        # Generator 
        self.G = ECPoint(
            9, 
            14781619447589544791020593568409986887264606134616475288964881837755586237401
        )

        # Order of E(F_p)
        self.n = np.power(2,252) + 277423177773723535358519377908836484

    def test(self):
        return self.c.addition(self.G, self.G)
    
    def generate_public_key(self, priv):
        return self.c.multiply(self.G, priv)


class SymmetricECC:
    """ Elliptic Curve Diffie Hellman Protocol. 
        Class implements a ECDH scheme where one use a base elliptic curve for generating a shared key.
        The shared key is used later for symmetric encryption. 
    """ 
    def __init__(self, 
                curve : EllipticCurve,
                private_key : int):

        # Internal parameters of the party.
        self.curve = curve 
        self.private_key = private_key

        # Public key is d * Q, where Q is the pre-agreed generator, d is the private key.    
        self.public_key = self.curve.generate_public_key(self.private_key)

        # Public key of the other party in the communication.
        self.party_pub_key = None 

        self.shared_key = None 
    
    def set_public_key(self, key):
        # Assoicate the client with the given key for further communication.
        self.party_pub_key = key 

    def get_public_key(self):
        return self.public_key

    def _establish_shared_key(self):
        # Establishes a shared secret key. 
        self.shared_key = self.curve.c.multiply(self.party_pub_key, self.private_key).x

    def encrypt_message(self, messge):
        # Message encryption using symmetric key.
        assert self.party_pub_key , "Missing public key of the other end in the communication."
        if not self.shared_key:
            self._establish_shared_key()

        # AES symmetric key encryption using the shared secret key. 
        cipher = AES.new(self.shared_key, AES.MODE_EAX)
        ciphertext, tag = cipher.encrypt_and_digest(message)

        return ciphertext



alice = SymmetricECC(Curve25519(), 3)
bob = SymmetricECC(Curve25519(), 5)



#alice.set_public_key(bob.get_public_key())
#alice.encrypt_message('hi')













    







        
     


