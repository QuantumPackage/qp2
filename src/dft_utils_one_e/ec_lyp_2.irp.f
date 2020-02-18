double precision function ec_lyp2(RhoA,RhoB,GA,GB,GAB)
 include 'constants.include.F'
 implicit none
 double precision, intent(in) :: RhoA,RhoB,GA,GB,GAB
 double precision :: Tol,caa,cab,cac,cad,cae,RA,RB,comega,cdelta,cLaa,cLbb,cLab,E
 ec_lyp2 = 0.d0
 Tol=1D-14
 E=2.718281828459045D0
 caa=0.04918D0
 cab=0.132D0
 cac=0.2533D0
 cad=0.349D0
 cae=(2D0**(11D0/3D0))*((3D0/10D0)*((3D0*(Pi**2D0))**(2D0/3D0)))


  RA = MAX(RhoA,0D0)
  RB = MAX(RhoB,0D0)
  IF ((RA.gt.Tol).OR.(RB.gt.Tol)) THEN
   IF ((RA.gt.Tol).AND.(RB.gt.Tol)) THEN
    comega = 1D0/(E**(cac/(RA+RB)**(1D0/3D0))*(RA+RB)**(10D0/3D0)*(cad+(RA+RB)**(1D0/3D0)))
    cdelta = (cac+cad+(cac*cad)/(RA+RB)**(1D0/3D0))/(cad+(RA+RB)**(1D0/3D0))
    cLaa = (cab*comega*RB*(RA-3D0*cdelta*RA-9D0*RB-((-11D0+cdelta)*RA**2D0)/(RA+RB)))/9D0
    cLbb = (cab*comega*RA*(-9D0*RA+(RB*(RA-3D0*cdelta*RA-4D0*(-3D0+cdelta)*RB))/(RA+RB)))/9D0
    cLab = cab*comega*(((47D0-7D0*cdelta)*RA*RB)/9D0-(4D0*(RA+RB)**2D0)/3D0)
    ec_lyp2 = -(caa*(cLaa*GA+cLab*GAB+cLbb*GB+cab*cae*comega*RA*RB*(RA**(8D0/3D0)+RB**(8D0/3D0))+(4D0*RA*RB)/(RA+RB+cad*(RA+RB)**(2D0/3D0))))
  endif
 endif
end
