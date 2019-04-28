#!/usr/bin/python

usage = """
*** Notes and Comments ***

1. The calculations are based mainly on two papers:

   "Analytical method for calculation of stresses and material 
   birefringence in polarization-maintaining optical fiber"
   by P. L. Chu and R. A. Sammut, J. Lightwave Technology, vol. LT-2,
   No. 5, Oct. 1984, pp.650-662
   
   "Design and fabrication of low-loss and low-crosstalk polarization-
   maintaining optical fibers," by Y. Sasaki, et al., 
   J. Lightwave Technology, vol. LT-4, No. 8, Aug. 1986, pp.1097-1102.


   Additional information on stress calculations and birefringence 
   can be found in:

   "Theory of elasticity" by S. P. Timoshenko & J. N. Goodier (1970)
   McGraw Hill, NY

   "Advanced strength and applied elasticity"
   by A. Ugural & S. Fenster (2003), Chap. 3, p.115, Prentice Hall
   
   "Stress-induced index profile distortion in optical waveguides"
   by G. W. Schere, Appl. Opt., Vol. 19, No. 12, June 1980,
   pp.2000-2006
   
   "Optical and mechanical effects of frozen-in stresses and strains
   in optical fibers"
   by A. D. Yablon, IEEE J. of QE, Vol. 10, No. 2, 2004, pp.300-311
   
2. The stress results stay the same as long as the geometry ratios are
   the same; In most of the plots, the radius is normalized with
   respect to the radius of the cladding.

3. The "Calculate best config (r fixed)" button is used to calculate
   the optimum stress rod diameter for a panda fiber that maximizes
   the birefringence, given a fixed core to cladding ratio and a fixed
   distance r from the center of the fiber to the edge of the SAP.
   (For a fixed stress rod diameter, the birefringence increases as r
   decreases).

4. Thermal expansion coefficient estimates are based on the values
   given in:
   
   "Single-mode fiber optics: principles and applications" 
   by Luc B. Jeunhomme, CRC Press (2nd edition, 1989)

5. Note that the change in refractive index has the opposite sign of
   the stress; i.e., a negative stress (compression) induces an
   *increase* in the refractive index.

6. Interpretation of stress: positive stress corresponds to tension
   ("pulling apart"), while negative stress corresponds to
   compression. The stress rods regions have a higher temperature
   expansion coefficient, thus the volume shrinks as temperature drops
   back to ambient. However, the surrounding area restricts (prevents)
   it from doing so, thus it forces that region to *expand*; hence the
   stress in the SAP region is *positive*, even though (in its natural
   state) it would shrink.  The region between the SAPs is also under
   positive stress (tension) along the x-axis, but is slightly
   compressed along the y-axis.  However the change in refractive
   index along each axis has a more complicated relationship with the
   stress: it depends also on the stress perpendicular to that
   particular axis.  The end result is that in the region near the
   core, it is under tension along the x-axis, but because of stress
   along the y-axis (which is negative), the overall delta n actually
   *increases*, and thus the x-axis is in fact the *slow* axis.
   
7. The ellipses in the stress visualization plot represents the
   compression and/or tension of the volume elements, with the
   directions of the major and minor axes of the ellipses representing
   the directions of the eigenvectors that diagonalize the stress
   tensor. The lengths of the major and minor axes are exaggerated to
   show the effect, and should not be taken literally.  The actual
   amount of compress and tension in glass is very small.
   
8. The boundary condition plot shows how well the boundary conditions
   are met after truncating the Airy stress function to N terms (the
   exact solution is an infinite series).  In practice N=20 gives a
   very good match for the panda fiber. The parameter N can be set in
   the main control panel.
   
9. For the bow-tie fiber, both the Airy stress function and the
   contribution from the SAPs are in the form of infinite series.  The
   boundary conditions as calculated will be perfectly matched, simply
   because we truncate both series with the same number of terms.
   However, the series associated with the SAPs is slow converging
   (the terms decay as 1/n), so requires many more terms to get a
   satisfactory result (the default is set to N=150). On top of this,
   the Gibb's phenomenon also kicks in around the edges of the SAP, so
   there are some "hot spots" that should be interpreted with care.
   """

# ----------------------------------------------------------------------------------------------------
# tooltips TBD
