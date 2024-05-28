I created this matlab function in order to acoustically model a plate reverb I built on, but have made some changes so that it can model most plates. 

It returns four sets of modes, the quasi-longitudinal modes (which are nearly inaudible, but still present,) the transverse shear modes (which are also most likely not very present) and the bending modes, (which tend to be the most audible modes.) These returned variables are calculated using a fairly simple analytical method. Finally, through modelling the plate via the finite difference method, and taking the eigenvalues of the stiffness matrix, we can get a final set of modes (although this seems to be a less reliable method.)

Longitudinal waves move in the same direction of the sound, the transverse shear waves move perpendicular to the sound, and the bending waves "combine" the two. I suggest reading the first chapter Acoustics and Psychoacoustics by David M. Howard and Jamie A. S. Angus if you want a solid understanding of how these waves behave.

I calculated the modes analytically by calculating the speed of sound in the medium, and used this value to calculate the transverse and longitudinal waves in a theoretical bar. I then extrapolated the bar example into two dimensions. For the bending waves, I used the function in the textbook mentioned above.

I also played around with generating a model from the equations used in classical plate theory. I believe these visualizations to be roughly correct, but the frequencies calculated are most likely inaccurate, as I mentioned before.

**I experimentally found that the calculated bending modes were the most accurate** , and I hope that this can be useful in anyone elses future plate modelling adventures!

The testing file is just for my own testing of the file
