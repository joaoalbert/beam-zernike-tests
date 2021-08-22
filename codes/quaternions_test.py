import numpy as np

class quaternion(object):

	def __init__(self, x, y, z, w):
	
		self.quaternion = np.array(x,y,z,w)
		
	def multiply(q1, q2):
	
		x = q1.quaternion[0] * q2.quaternion[0]
		
		return quaternion()
