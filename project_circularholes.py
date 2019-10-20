import argparse
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
from shopbotfile import ShopBotFile
import workpiece

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Drill an array of constant sized holes.')

    parser.add_argument('--bit-diameter', default=0.5, type=float,
                        help='Diameter of the endmill bit being used.')
    parser.add_argument('--bit-clearout-depth', default=2.0+1/32.0, type=float,
                        help='Maximum depth before a clearout operation is performed.')

    parser.add_argument('--hole-diameter', default=1.0, type=float,
                        help='Diameter of the holes being drilled.')

    parser.add_argument('--hole-depth', default=0.25, type=float,
                        help='Depth of the holes being drilled.')

    parser.add_argument('--x-spacing', default=2.5, type=float,
                        help='Distance between holes in the x-direction')
    parser.add_argument('--x-length', default=0.0, type=float,
                        help='Total length of the workpiece in the x-direction.')
    parser.add_argument('--x-min-edge', default=0.0, type=float,
                        help='Minimum distance of the holes from the edge of the workpiece in the x-direction.')

    parser.add_argument('--y-spacing', default=2.5, type=float,
                        help='Distance between holes in the y-direction')
    parser.add_argument('--y-length', default=0.0, type=float,
                        help='Total length of the workpiece in the y-direction.')
    parser.add_argument('--y-min-edge', default=0.0, type=float,
                        help='Minimum distance of the holes from the edge of the workpiece in the y-direction.')

    args = parser.parse_args()

    endmill_bit = RouterBit(args.bit_diameter, clearout_depth_inches = args.bit_clearout_depth)

    hole_path = toolpaths.HolePath(bit = endmill_bit,
                                   diameter = args.hole_diameter,
                                   depth = args.hole_depth,
                                   direction = c.DIRECTION_CONVENTIONAL)

    multiple_holes_path = toolpaths.MultipleHolesPath(hole_path)
    x_offsets = center_equally_spaced_points_in_range(0, args.x_length, args.x_spacing, args.x_min_edge)
    y_offsets = center_equally_spaced_points_in_range(0, args.y_length, args.y_spacing, args.y_min_edge)
    print('X locations: {}'.format(x_offsets))
    print('Y locations: {}'.format(y_offsets))
    for x_idx, offset_x in enumerate(x_offsets):
        direction = 1 if x_idx % 2 == 0 else -1
        for offset_y in y_offsets[::direction]:
            multiple_holes_path.add_location(offset_x = offset_x, offset_y = offset_y)

    shop_bot_file = ShopBotFile('custom_size_hole.sbp')
    shop_bot_file.set_ramps(small_circle_diameter = 0.2, xy_move_ramp_speed = 0.8)

    shop_bot_file.add_points(multiple_holes_path.get_points())
    shop_bot_file.close()

    points = np.array(multiple_holes_path.get_points())

    width_x = args.x_length
    origin_x = 0.0
    if width_x == 0.0:
        width_x = 2.0 * args.hole_diameter
        origin_x = -0.5 * width_x

    width_y = args.y_length
    origin_y = 0.0
    if width_y == 0.0:
        width_y = 2.0 * args.hole_diameter
        origin_y = -0.5 * width_y
    workpiece_preview = workpiece.WorkpiecePreview([width_x, width_y, args.hole_depth],
                                                   origin=[origin_x, origin_y, 0])
    workpiece_preview.plot_birds_eye(points)
    workpiece_preview.plot_three_axis(points)