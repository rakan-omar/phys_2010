# phys_2010

okay so there is an explosion of a slightly randomised force, somewhere between quarters and three-quarters of the way
through projectile motion

sim3d and sim3din2d are three dimensional (the latter shows 2 dimensional graphs, so it's easier to read), both are working but the results are disastrous.

the bomb is assumed to be a circle. and the fragments (with masses randomised to be equal to the sum M) are clean cut sectors of the circle. then, their centre of mass is along theradius that starts from the middle of the sector of each fragment (by symmetry). then we can calculate the angle from the horizontal to each fragment's centre of mass, which gives us it's direction of motion

the speed of each is given by conservation of momentum (x and y components seperately):

speed = (impulse + (initial_speed*fragment_mass)) / fragment_mass

the current results show that the centre of mass is always close or approximately equal to the original trajectory (even for
large forces) but not necessarily exactly the same.
BTW, these results also hold if we move the explosion to before the appex

I heavily recommend playing around with the inital speed, initial bomb mass, and explosion force (range).

if you can't see the graph because of the labels, comment out the last "legend" line
