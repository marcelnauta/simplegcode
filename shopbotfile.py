import copy
import os
import math
import numpy as np

import constants as c
from helpers import divide_into_equal_passes, center_equally_spaced_points_in_range
from helpers import divide_with_clearout_depths
from routerbit import RouterBit
import templates
import toolpaths
import workpiece

class ShopBotFile(object):
    def __init__(self, filename = 'tmp.sbp'):
        self.outfile = open(filename, 'w')
        self._write(templates.get_header())
        self.movement_speeds = c.MovementSpeeds()
        self.ramp_speeds = c.RampSpeeds()
        self.set_speed()
        self.set_ramps()

    def _write(self, out_str):
        self.outfile.write(out_str)
        print(out_str)

    def close(self):
        self._write(templates.get_footer())
        self.outfile.close()

    def set_speed(self, **kwargs):
        self.movement_speeds.set(**kwargs)
        self._write(self.movement_speeds.get_command())

    def set_ramps(self, **kwargs):
        self.ramp_speeds.set(**kwargs)
        self._write(self.ramp_speeds.get_command())

    def move_to(self, point, command_prefix = 'M'):
        if len(point) == 4: # At Maker Labs, the CNC is 4 axis, but the 4th axis acts as the 5th
            point = (point[0], point[1], point[2], 0.0, point[3])
        point_as_str = ['{0:.6f}'.format(x) for x in point]
        command_tag = command_prefix + '{0:d}'.format(len(point))
        self._write(','.join([command_tag] + point_as_str) + '\n')

    def jog_to(self, point):
        self.move_to(point, command_prefix = 'J')

    def add_points(self, point_list):
        [self.move_to(point) for point in point_list]

    def pause(self, time):
        self._write('PAUSE {}\n'.format(time))


if __name__ == '__main__':
    #make_dovetail_domino()
    dovetail_bit = RouterBit(0.5)
    endmill_bit = RouterBit(0.5, clearout_depth_inches = 2.01)

    move_to_path = toolpaths.MoveToPath([12.0, 2.0, 0.0, math.pi])
    shift_2d = toolpaths.Shift2D(1.0, 1.0)
    move_to_path.add_transform(shift_2d)

    hole_path = toolpaths.HolePath(bit = endmill_bit, diameter = 0.75, depth = 1.55,
                                   direction = c.DIRECTION_CONVENTIONAL)


    multiple_holes_path = toolpaths.MultipleHolesPath(hole_path)
    x_offsets = center_equally_spaced_points_in_range(0, 42, 2.5, 2.0)
    y_offsets = [0.0] #center_equally_spaced_points_in_range(0, 42, 2.5, 1.0)
    for x_idx, offset_x in enumerate(x_offsets):
        direction = 1 if x_idx % 2 == 0 else -1
        for offset_y in y_offsets[::direction]: #y_offsets[::direction]:
            multiple_holes_path.add_location(offset_x = offset_x, offset_y = offset_y)

    shop_bot_file = ShopBotFile('custom_size_hole.sbp')
    shop_bot_file.set_ramps(small_circle_diameter = 0.2, xy_move_ramp_speed = 0.8)

    shop_bot_file.add_points(multiple_holes_path.get_points())
    shop_bot_file.close()

    points = np.array(multiple_holes_path.get_points())

    workpiece_preview = workpiece.WorkpiecePreview([2, 2, 2], origin=[-1, -1, 0])
    workpiece_preview.plot_three_axis(points)

    print(points)