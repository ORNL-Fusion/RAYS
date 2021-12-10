# RAYS
Plasma ray tracing code.  A more accurate description would be skeleton for a plasma ray
tracing code.  This is not a production code.  The primary objective is to provide a starting
point for discussion of what a modern ray tracing code, built on a green-field site might
look like.  A secondary objective is to provide a very simple, stripped down framework in
which some basic questions can be looked at such as, is some particular ODE solver better for
rays tracing applications than another?, is it feasible or perhaps advantageous to compute
the needed derivatives of the dispersion relation numerically?, do spline fits to equilibrium
data, which have limited continuous derivatives, introduce errors?

With respect to the first objective we should consider what are the disiderata for such a 
code?  What would we be trying to accomplish?  I suggest:

* Ease of use
* Flexibility - easy support for multiple geometries, dispersion models, equilibrium models, etc
* Extensibility - easy to add new features without disturbing existing ones
* Verifiability - ease of comparison with simple situations for which rigorous answers are available
* Maintainability - straightforward, human readable coding with different modules as parallel in structure as possible
* Speed and accuracy

What is here is a rework of the MORAYS fortran 90 code written in the 1990s by Caiyi Wang
which was based on the ancient RAYS code.  I have stripped out all graphics and 
preprocessing, and for now all reference to warm plasma effects.  It uses fortran 90 modules
extensively and focuses on defining generic interfaces which are supported by modules which
may be very specific.  There are many parts in the closet left from the previous editions,
which I believe can be easily adapted and reattached to this architecture.  But for now
the focus is on simplicity.  Other capabilities can be added as there is interest, and as we
see whether this architecture really is a good approach to build on.


