# Summer Projects Presentation

- brief introduction to exterior calculus
  - a k-form is a thing that can be integrated over a k-manifold
  - in 3D 1-forms and 2-forms have three components at each point, and are often
    treated as vector fields.

- how to discretize?
  - we want to keep the integrable structure
  - sampling the field at some points is most obvious, but it's unclear how to
    integrate
  - instead, put components on respective pieces of a cell complex. now
    integration is clear.
  - this can also be thought of as sampling by integrating over the piece,
    instead of just sampling a point.
  - this discretization preserves stokes theorem (give an example of the
    exterior derivative)
  - this can be thought of an approximation of the continuum where we take the
    limit as the grid gets fine. or, this can be thought of as an interesting
    discrete structure in its own right.
  - other EC operations like hodge duality and the laplacian can be formulated
    in DEC, but are less natural.

- fourier methods
  - start by introducing the point charge potential on a grid problem.
    - trick: a function on Z^n is just a sum of delta functions, so we can use
      the normal fourier transform.
    - in fact, a discrete function is just something times the comb function,
      and a periodic function is something convolved with the comb, explaining
      why FT sends periodic to discrete and discrete to periodic.
    - we get this expression for the FT of the potential, which we can expand
      into a series, or, more usefully, integrate over one period to recover the
      values.

  - next example: div free flow in 2d with point source
    - can just take your favorite continuous streamfunction and sample and take
      exterior derivative. eg, theta
    - alternatively, you can impose a gauge condition, and then we get this
      system of eqns whose solution looks very similar to the first problem.

  - a more useful example: stokes
    - like before we get a small system of eqns in freq space, which we can
      invert
    - I got a computer to do this. Here's a graph of the speed compared with a
      standard linear solver method.
    - Since there are N^3 cells and factoring an mxm matrix generally takes
      O(m^3) time, the graph agrees with the theoretical O(N^9) growth for
      setting up QR, O(N^6) for solving with QR, and O(N^3 log N) for both parts
      of fourier.

- and also: here's an electromegnetism simulation. It looks cool!
