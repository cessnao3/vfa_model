# Very Flexible Aircraft Model

Matlab model of a simplified Very Flexible Aircraft (VFA) model by Travis E. Gibson, Anuradha Annaswamy, and Eugene Lavretsky.

* "Closed-loop Reference Model adaptive control : with application to very flexible aircraft," http://hdl.handle.net/1721.1/87974
* "Modeling for Control of Very Flexible Aircraft," https://doi.org/10.2514/6.2011-6202

## File List

| Name                | Description                                                      |
|---------------------|------------------------------------------------------------------|
| atmosphere4.m       | Standard atmospheric model for the 1976 NASA Standard Atmosphere |
| find_steady_state.m | Finds steady-state parameters for the VFA model using fmincon    |
| odefunc.m           | Helper function to call *vfa_deriv.m*                            |
| vfa_deriv.m         | Nonlinear VFA model                                              |

## References

* Gibson, T. E. (2014). Closed-loop Reference Model adaptive control : with application to very flexible aircraft. Massachusetts Institute of Technology, Department of Mechanical Engineering. Cambridge, Massachusetts: Massachusetts Institute of Technology. Retrieved from http://hdl.handle.net/1721.1/87974
* Gibson, T., Annaswamy, A., & Lavretsky, E. (2011). Modeling for Control of Very Flexible Aircraft. AIAA Guidance, Navigation, and Control Conference. Portland, Oregon. doi:10.2514/6.2011-6202
