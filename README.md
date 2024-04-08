I created this matlab function in order to acoustically model a plate reverb I worked on, but have made some changes so that it can model any simply supported plate. 

It returns three sets of modes, the quasi-longitudinal modes (which are nearly inaudible, but still present,) the transverse shear modes (which are also most likely not very present) and the bending modes, (which tend to be the most audible modes.) These returned variables are calculated using a fairly simple analytical method.

Longitudinal waves move in the same direction of the sound, the transverse shear waves move perpendicular to the sound, and the bending waves "combine" the two. I suggest reading the first chapter Acoustics and Psychoacoustics by David M. Howard and Jamie A. S. Angus if you want a solid understanding of what these waves are.

I calculated the modes analytically by calculating the speed of sound in the medium, and used this value to calculate the transverse and longitudinal waves in a theoretical bar. I then extrapolated the bar example into two dimensions. For the bending waves, I used the function in the textbook mentioned above.

I also played around with generating a model from the equations used in classical plate theory. I believe these visualizations to be roughly correct, but the frequencies calculated are most likely inaccurate.

I experimentally found that the bending modes were the most accurate, and I hope that this can be used!

