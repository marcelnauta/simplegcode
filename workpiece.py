import copy
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np

import constants as c

class WorkpiecePreview(object):
    def __init__(self, size, origin = [0,0,0]):
        self.size = size
        self.origin = origin
        
    def plot_three_axis(self, points):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Draw lines individually to vary color
        num_points = np.size(points, 0)
        colors = pl.cm.gist_earth(np.linspace(0,0.9,num_points))
        for n in range(1, num_points):
            ax.plot(points[n-1:n+1,c.AXIS_X], points[n-1:n+1,c.AXIS_Y], points[n-1:n+1,c.AXIS_Z], '-', color=colors[n])
        self._draw_bounding_box(ax)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        plt.show()
    
    def _draw_bounding_box(self, axes):
        for axis_1, axis_2, axis_3 in [[c.AXIS_X, c.AXIS_Y, c.AXIS_Z],
                                       [c.AXIS_Y, c.AXIS_Z, c.AXIS_X],
                                       [c.AXIS_Z, c.AXIS_X, c.AXIS_Y]]:
            for is_max_2 in [0, 1]:                            
                for is_max_3 in [0, 1]:
                    line = [np.array([o,o]) for o in self.origin]
                    line[axis_1][1] += self.size[axis_1]
                    for i in range(2):
                        line[axis_2][i] += is_max_2 * self.size[axis_2]
                        line[axis_3][i] += is_max_3 * self.size[axis_3]
                    axes.plot3D(line[c.AXIS_X], line[c.AXIS_Y], -line[c.AXIS_Z], color="b")
        
        