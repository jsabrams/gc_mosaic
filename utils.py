import numpy as np

def spacing_fn(x, y, en = 1.633922, et = 1.6666, es = 1.126081, ei = 1.488036, s0 = 1.0/120):
	"""
	Specifies the average ganglion cell spacing from its neighbors when the ganglion cell is at (x,y)
	in degrees relative to the fovea
	Inputs:	x: an n x m array of x positions relative to the fovea (usually n x 1)
			y: an n x m array of y positions relative to the fovea (usually n x 1)
			en: nasal visual field half spacing (default half spacings are from Drasdo et al., 2007)
			et: temporal visual field half spacing
			es: superior visual field half spacing
			ei: inferior visual field half spacing
			s0: the spacing in the fovea (default is 1/120 of a degree)
	Output: Q: an n x m array of spacings in degrees (usually n x 1) 
	"""
	
	#Spacing for the four quadrants
	Q1 = s0 * (np.sqrt(((x**2) / (et**2)) + ((y**2) / (es**2))) + 1)
	Q2 = s0 * (np.sqrt(((x**2) / (en**2)) + ((y**2) / (es**2))) + 1)
	Q3 = s0 * (np.sqrt(((x**2) / (en**2)) + ((y**2) / (ei**2))) + 1)
	Q4 = s0 * (np.sqrt(((x**2) / (et**2)) + ((y**2) / (ei**2))) + 1)

	#Piece together the quadrants
	Q = np.zeros(x.shape);
	Q[np.logical_and(x >= 0, y >= 0)] = Q1[np.logical_and(x >= 0, y >= 0)]
	Q[np.logical_and(x < 0, y >= 0)] = Q2[np.logical_and(x < 0, y >= 0)]
	Q[np.logical_and(x < 0, y < 0)] = Q3[np.logical_and(x < 0, y < 0)]
	Q[np.logical_and(x >= 0, y < 0)] = Q4[np.logical_and(x >= 0, y < 0)]

	return Q

def cart2pol(x, y):
	"""
	Convert cartesian coordinates into polar coordinates
	Inputs: 	x: n x 1 array of x positions
				y: n x 1 array of y positions
	Outputs:	theta: n x 1 array of angles
				rho: n x 1 array of distances
	"""
    theta = np.arctan2(y, x)
    rho = np.sqrt(x**2 + y**2)
    
    return (theta, rho)

def pol2cart(theta, rho):
	"""
	Convert polar coordinates into cartesian coordinates
	Inputs:		theta: n x 1 array of angles
				rho: n x 1 array of distances 	
	Outputs:	x: n x 1 array of x positions
				y: n x 1 array of y positions
	"""
	x = rho * np.cos(theta)
	y = rho * np.sin(theta)
	
	return (x, y)

def intersect_circle(x1, y1, x2, y2, en = 1.633922, et = 1.6666, es = 1.126081, ei = 1.488036, s0 = 1.0/120):
	"""
	Solve for the intersection of two circles whose centers are (x1,y1) and (x2,y2), defaults are needed for
	the spacing function
	Inputs:	x1: center x coordinate of the first circle
			y1: center y coordinate of the first circle
			x2: center x coordinate of the second circle
			y2: center y coordinate of the second circle
			en: nasal visual field half spacing (default half spacings are from Drasdo et al., 2007)
			et: temporal visual field half spacing
			es: superior visual field half spacing
			ei: inferior visual field half spacing
			s0: the spacing in the fovea (default is 1/120 of a degree)
	Output:	xI_1: x coordinate for intersection one
			yI_1: y coordinate for intersection one
			xI_2: x coordinate for intersection two
			yI_2: y coordinate for intersection two
	"""
	sp_1 = spacing_fn(x1, y1, en = 1.633922, et = 1.6666, es = 1.126081, ei = 1.488036, s0 = 1.0/120)
	sp_2 = spacing_fn(x2, y2, en = 1.633922, et = 1.6666, es = 1.126081, ei = 1.488036, s0 = 1.0/120)

	x_diff = x2 - x1
	y_diff = y2 - y1

	cent_dist = np.sqrt(x_diff**2 + y_diff**2)
	c1_to_ln = (cent_dist**2 + sp1**2 - sp2**2) / (2 * cent_dist) #Dist from center 1 to the line joining the points of intersection

	xI_1 = x1 + x_diff * c1_to_ln / cent_dist + (y_diff / cent_dist) * np.sqrt(sp1**2 - sp2**2)
	yI_1 = y1 + y_diff * c1_to_ln / cent_dist - (x_diff / cent_dist) * np.sqrt(sp1**2 - sp2**2)

	xI_2 = x1 + x_diff * c1_to_ln / cent_dist - (y_diff / cent_dist) * np.sqrt(sp1**2 - sp2**2)
	yI_2 = y1 + y_diff * c1_to_ln / cent_dist + (x_diff / cent_dist) * np.sqrt(sp1**2 - sp2**2)

	return ()

def circle_spacer(cx, cy):
	"""
	Takes a vector of coordinates (cx, cy) and output a modified circle in which the points are more equally spaced
	Inputs:		cx: n x 1 array of x coordinates
				cy: n x 1 array of y coordinates
	Outputs:	new_x: n x 1 array of new x coordinates
				new_y: n x 1 array of new y coordinates
	"""
	theta, rho = cart2pol(cx, cy)

	#Get the angle between the first and last coordinates
	if theta[0] < -np.pi / 2.0 and theta[-1] > np.pi / 2.0:
		theta_0_e = np.mod(theta[0], 2.0 * np.pi) - theta[-1]
	else:
		theta_0_e = theta[0] - theta[-1]

	#Get the angle between the 1st and 2nd coordinates
	if theta[1] < -np.pi / 2.0 and theta[0] > np.pi / 2.0:
		theta_1_2 = np.mod(theta[1], 2.0 * np.pi) - theta[0]
	else:
		theta_1_2 = theta[1] - theta[0]

	#Get the increment that we subtract from each interval
	incr = theta_1_2 - theta_0_e

	#Subtract the proper increments
	I = np.linspace(1, len(cx), num=len(cx))
	theta = theta - (I - 1) * incr / len(cx)

	new_x, new_y = pol2cart(theta, rho)

	return (new_x, new_y)
