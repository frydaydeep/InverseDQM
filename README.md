### A realization of InverseDQM

<!--
**frydaydeep/frydaydeep** is a āØ _special_ āØ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- š­ Iām currently working on ...
- š± Iām currently learning ...
- šÆ Iām looking to collaborate on ...
- š¤ Iām looking for help with ...
- š¬ Ask me about ...
- š« How to reach me: ...
- š Pronouns: ...
- ā” Fun fact: ...
-->

A numerical method is proposed to calculate the integral from a series of discrete functional value,
derived directly from differential quadrature method (DQM).

Imagine different order of a differential equation as different steps of a ladder, and your foot on the 
step where you can get hold of a series of discrete function value. DQM tells you how to climb down, and 
inverseDQM proposed here tells you how to climb up as long as you have the right amount of boundary condtions.

You can find a example "IDQMforEularBeam.jl", and a pdf explainning the theory.
*try to play with the DQM function in the code first and see how it works before you try IDQM. 

A more general function version is on its way where you input your discrete data and boundary conditions you have
and the order you'd like to have for results. It returns the discrete data in your desired order.
