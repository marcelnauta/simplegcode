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

    def estimate_time(self, points, speed):
        distance = np.sum(np.sqrt((points[:-1,c.AXIS_X] - points[1:,c.AXIS_X])**2
                                 + (points[:-1,c.AXIS_Y] - points[1:,c.AXIS_Y])**2
                                 + (points[:-1,c.AXIS_Z] - points[1:,c.AXIS_Z])**2))
        return distance / speed

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

    def plot_birds_eye(self, points):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')

        # Draw lines individually to vary color
        num_points = np.size(points, 0)
        max_depth = np.max(points[:,c.AXIS_Z])

        if 0: # color code depth (too slow)
            colors = pl.cm.gist_earth(0.9 * points[:,c.AXIS_Z] / max_depth)
            for n in range(1, num_points):
                ax.plot(points[n-1:n+1,c.AXIS_X], points[n-1:n+1,c.AXIS_Y], '-', color=colors[n])
        else:
            ax.plot(points[:,c.AXIS_X], points[:,c.AXIS_Y], '-')
        self._draw_bounding_box_2d(ax)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        plt.show()

    def _draw_bounding_box_2d(self, axes):
        min_x = self.origin[c.AXIS_X]
        min_y = self.origin[c.AXIS_Y]
        max_x = self.origin[c.AXIS_X] + self.size[c.AXIS_X]
        max_y = self.origin[c.AXIS_Y] + self.size[c.AXIS_Y]
        axes.plot([min_x, max_x, max_x, min_x, min_x],
                  [min_y, min_y, max_y, max_y, min_y],
                  color="b")

