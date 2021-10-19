from galois import ec 


curve = ec.EllipticCurve(A=2,B=3,p=97,n=1) # Over F_97 
curve1 = ec.EllipticCurve(2, 3) # Over the reals
P = ec.ECPoint(2, 3)
