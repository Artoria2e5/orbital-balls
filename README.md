orbital-balls
=============

orbital-balls is probably a solar system sandbox that runs on the Arduboy. It may or may not be a Game Jam 3 entry.

orbital-balls uses patched conics, with major simplifications.

The Renderer
------------

The renderer requests points for mean anomaly at an increment of pi/4, approximating the orbit by an 8-gon (with some bounds checking). Each segment is enhanced by a quadratic interpolator, allowing for some 24-gon display.

SOI handling
------------

Patched conics simulation requires SOI handling -- when does the body move to another center mass?

The "exit to parent" SOI handling is exact. It involves working the position evaluation back in reverse.

The "enter sibling SOI" handling is based on two quadratic interpolations given above. That provides a quartic expression for square of distance, which we can minimize with single-variable calculus and a cubic equation solver. Instead of using fixed increments, however, the interpolation windows are of uniform size.

To make the results consistent, the state vectors used for "enter sibling" will always come from interpolation. In both cases, some addition/subtraction on the state vectors are required.

Numerical precision
-------------------

Oh, hell yeah, we are using single-precision floats here! Every position is recorded relative to the parent body, so perhaps we are fine?

License
-------

Things made by me, that is everything except Vector.h and poly34.\*, are licensed under CC0 1.0.
