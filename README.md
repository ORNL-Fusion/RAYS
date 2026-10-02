# RAYS
Although referred to as a ray tracing ‘code’ or ‘mini-App’, RAYS might more
properly be considered as a framework for constructing geometrical optics codes
or applications. It consists of a set of statically linked libraries that
provide ray tracing services, several codes that employ the libraries to carry
out ray tracing modeling, plus a collection of graphics applications for
displaying the results.  All of the services of the static libraries are user
accessible and can easily be used directly to provide geometrical optics
capabilities inside other codes, such as transport models or stability codes.
The work has been supported by the RF SciDAC project, by direct subcontracts
with ORNL, and by the AToM SciDAC project.

Ray tracing has been an important tool for radio-frequency (RF) applications in
Fusion, and is likely to continue to be for a long time.  Ray tracing can give
good quantitative answers for short wavelength modes and often for longer
wavelengths, where geometrical optics approximations are not well satisfied.
Ray tracing is easy to implement for complicated geometries where full-wave
treatment is essentially impossible.  We do have sophisticated full-wave codes
(AORSA, TORIC, TORLH) but these are orders of magnitude more computationally
expensive than geometrical optics, when geometrical optics is applicable.  The
intent in developing RAYS was to essentially start on a green-field-site and
take advantage of recent advances in programming languages, programming
standards, and computer capabilities so as to provide a framework for
development into the future.  Specific design goals were:


* Ease of use – support for multiple geometries, dispersion models etc in a single
implementation, simple I/O structure
* Extensibility – ease of adding new features without disturbing existing ones,
simple hierarchical organization
* Maintainability – minimal use of external libraries, simple build system
* Verifiability - ease of comparison with analytically solvable cases using
single code build
* Speed and accuracy – consideration of parallelism from beginning of design
phase, numerical precision control

The vision is to go beyond the ‘hero code’ model that can only be modified or
developed by the original author, but that is so transparently laid out that
anyone who understands the physics and has modest programming experience, can be
a developer.


The design process was informed by lessons learned from development of the
Integrated Plasma Simulator (IPS) carried out under the SWIM SciDAC project.
That is, by careful consideration of the identification of separable component
functions, and development of robust interfaces between components which need
never (at least rarely) modified.  The largest of these separation decisions was
to adopt a Kernel/Front-End structure.  The function of the kernel is the actual
integration of the ray equations, which are a Hamiltonian system.  This task is
generic, the structure of which is independent of geometry, frequency regime,
wave dispersion model, or differential equation solver method.  On the other
hand, the applications of the ray tracing are not
at all generic.  The selection of rays to be traced depends on the geometry, the
launcher characteristics, and on the particular physics being studied.
Similarly, the kinds of analysis and graphical display of results depend
entirely on the geometry and particular needs of the study.  The solution
adopted by RAYS is to develop one kernel and many front ends.  A front-end
consists of an application code, a post-processing system, and appropriate
graphical display capabilities.  Users can tailor these to fit their needs.  The
kernel is compiled into a statically linked library, RAYS_lib, as is the core of
the post-processing system, post_process_lib.  To minimize reliance on external
dependencies, additional mathematical and other support services are provided by
other linkable static libraries.  At present there are only two external
dependencies – OpenMP and netCDF.  RAYS is open source and is accessible at the
ORNL Fusion github repository – https://github.com/ORNL-Fusion/RAYS.


![RAYS-and-kernel](readme_images/RAYS and kernel-fig)
